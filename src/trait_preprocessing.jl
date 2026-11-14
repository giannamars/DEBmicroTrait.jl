struct TraitPreprocessor
    names::Vector{Symbol}
    id_column::Union{Nothing,Int}
    binary_columns::Vector{Int}
    count_columns::Vector{Int}
    continuous_columns::Vector{Int}
    count_log_min::Vector{Float64}
    count_log_max::Vector{Float64}
    continuous_mean::Vector{Float64}
    continuous_std::Vector{Float64}
end

function _read_trait_csv(path::AbstractString; delimiter = ',', header = true)
    raw_table = DelimitedFiles.readdlm(path, delimiter, Any, '\n'; header = header)
    if header
        data, raw_header = raw_table
        header_names = strip.(vec(raw_header))
        if !isempty(header_names)
            header_names[1] = replace(header_names[1], '\ufeff' => "")
        end
        names = Symbol.(header_names)
        if !isempty(names) && isempty(String(names[1]))
            data = data[:, 2:end]
            names = names[2:end]
        end
        return data, names
    end

    data = raw_table
    names = [Symbol("column_$(index)") for index in 1:size(data, 2)]
    return data, names
end

function _normalize_trait_columns(columns, names::Vector{Symbol})
    columns === nothing && return Int[]
    items = columns isa Tuple ? collect(columns) : collect(columns isa AbstractVector ? columns : (columns,))
    indices = Int[]
    for item in items
        if item isa Integer
            1 <= item <= length(names) || throw(ArgumentError("column index $(item) is out of bounds"))
            push!(indices, Int(item))
        elseif item isa Symbol
            index = findfirst(==(item), names)
            index === nothing && throw(ArgumentError("unknown column name $(item)"))
            push!(indices, index)
        elseif item isa AbstractString
            symbol = Symbol(item)
            index = findfirst(==(symbol), names)
            index === nothing && throw(ArgumentError("unknown column name $(item)"))
            push!(indices, index)
        else
            throw(ArgumentError("column identifiers must be integers, symbols, or strings"))
        end
    end
    return indices
end

function _normalize_single_column(column, names::Vector{Symbol})
    column === nothing && return nothing
    indices = _normalize_trait_columns((column,), names)
    length(indices) == 1 || throw(ArgumentError("expected exactly one id column"))
    return first(indices)
end

function _normalize_schema_list(values, names::Vector{Symbol})
    values === nothing && return Int[]
    return _normalize_trait_columns(values, names)
end

function _read_trait_schema(schema_path::AbstractString, names::Vector{Symbol})
    schema = TOML.parsefile(schema_path)
    id_column = haskey(schema, "id_column") ? _normalize_single_column(schema["id_column"], names) : nothing
    binary_columns = haskey(schema, "binary_columns") ? _normalize_schema_list(schema["binary_columns"], names) : Int[]
    count_columns = haskey(schema, "count_columns") ? _normalize_schema_list(schema["count_columns"], names) : Int[]
    continuous_columns = haskey(schema, "continuous_columns") ? _normalize_schema_list(schema["continuous_columns"], names) : Int[]

    return (; id_column, binary_columns, count_columns, continuous_columns)
end

function _detect_delimiter(path::AbstractString)
    open(path, "r") do io
        line = readline(io)
        occursin('\t', line) && !occursin(',', line) && return '\t'
        occursin(';', line) && !occursin(',', line) && return ';'
        return ','
    end
end

function _load_relabund_map(path::AbstractString; id_column = :id, value_column = Symbol("rel.abund"), delimiter = ',', header = true)
    detected_delimiter = delimiter == ',' ? _detect_delimiter(path) : delimiter
    data, names = _read_trait_csv(path; delimiter = detected_delimiter, header = header)
    id_index = _normalize_single_column(id_column, names)
    value_index = _normalize_single_column(value_column, names)

    relabund_map = Dict{String,Float64}()
    for row_index in 1:size(data, 1)
        id_value = data[row_index, id_index]
        id_value === missing && continue
        abundance = _parse_trait_cell(data[row_index, value_index])
        abundance === missing && continue
        relabund_map[string(id_value)] = Float64(abundance)
    end

    return relabund_map
end

function align_relabund(ids::AbstractVector, path::AbstractString; id_column = :id, value_column = Symbol("rel.abund"), delimiter = ',', header = true, missing_value = 0.0)
    relabund_map = _load_relabund_map(path; id_column = id_column, value_column = value_column, delimiter = delimiter, header = header)
    values = Vector{Float64}(undef, length(ids))
    for index in eachindex(ids)
        identifier = string(ids[index])
        values[index] = get(relabund_map, identifier, Float64(missing_value))
    end
    return values
end

function _parse_trait_cell(cell)
    if cell === missing
        return missing
    elseif cell isa Missing
        return missing
    elseif cell isa Bool
        return cell ? 1.0 : 0.0
    elseif cell isa Number
        return Float64(cell)
    elseif cell isa AbstractString
        text = strip(cell)
        isempty(text) && return missing
        lowered = lowercase(text)
        lowered in ("na", "nan", "missing", "null") && return missing
        lowered == "true" && return 1.0
        lowered == "false" && return 0.0
        parsed = tryparse(Float64, text)
        parsed === nothing && throw(ArgumentError("could not parse trait value `$(text)` as a number"))
        return parsed
    else
        return Float64(cell)
    end
end

function _fit_count_statistics(values::AbstractVector{<:Real})
    transformed = log1p.(max.(Float64.(values), 0.0))
    minimum_value = isempty(transformed) ? 0.0 : minimum(transformed)
    maximum_value = isempty(transformed) ? 0.0 : maximum(transformed)
    return minimum_value, maximum_value
end

function _fit_continuous_statistics(values::AbstractVector{<:Real})
    finite_values = filter(isfinite, Float64.(values))
    isempty(finite_values) && return 0.0, 0.0
    length(finite_values) == 1 && return first(finite_values), 0.0

    mean_value = sum(finite_values) / length(finite_values)
    squared_error = sum((value - mean_value)^2 for value in finite_values)
    standard_deviation = sqrt(squared_error / (length(finite_values) - 1))
    return mean_value, standard_deviation
end

"""
    fit_trait_preprocessor(path; id_column, binary_columns, count_columns, continuous_columns, schema_file)

Read a trait CSV once, fit per-column scaling statistics, and return a `TraitPreprocessor`
that can be reused to transform the same schema again. If `schema_file` is provided, it
controls the column groups. Otherwise, all three trait column groups must be provided
explicitly.
"""
function fit_trait_preprocessor(path::AbstractString;
    id_column = nothing,
    binary_columns = nothing,
    count_columns = nothing,
    continuous_columns = nothing,
    schema_file = nothing,
    delimiter = ',',
    header = true,
)
    data, names = _read_trait_csv(path; delimiter = delimiter, header = header)
    nrows, ncols = size(data)

    if schema_file !== nothing
        schema = _read_trait_schema(schema_file, names)
        id_index = schema.id_column
        binary_indices = schema.binary_columns
        count_indices = schema.count_columns
        continuous_indices = schema.continuous_columns
    else
        id_index = _normalize_single_column(id_column, names)
        if binary_columns !== nothing && count_columns !== nothing && continuous_columns !== nothing
            binary_indices = _normalize_trait_columns(binary_columns, names)
            count_indices = _normalize_trait_columns(count_columns, names)
            continuous_indices = _normalize_trait_columns(continuous_columns, names)
        else
            throw(ArgumentError("provide schema_file or all three explicit trait column groups"))
        end
    end

    all_indices = vcat(binary_indices, count_indices, continuous_indices)
    length(unique(all_indices)) == length(all_indices) || throw(ArgumentError("trait columns must not overlap"))
    id_index === nothing || (id_index in all_indices && throw(ArgumentError("id_column must not overlap with trait columns")))
    expected_columns = id_index === nothing ? collect(1:ncols) : setdiff(collect(1:ncols), [id_index])
    sort(all_indices) == expected_columns || throw(ArgumentError("binary, count, and continuous columns must cover every non-id CSV column exactly once"))

    count_log_min = Float64[]
    count_log_max = Float64[]
    for column_index in count_indices
        column_values = Float64[]
        for row_index in 1:nrows
            cell = _parse_trait_cell(data[row_index, column_index])
            value = cell === missing ? 0.0 : Float64(cell)
            isfinite(value) || (value = 0.0)
            value < 0 && throw(ArgumentError("count columns must not contain negative values"))
            push!(column_values, value)
        end
        minimum_value, maximum_value = _fit_count_statistics(column_values)
        push!(count_log_min, minimum_value)
        push!(count_log_max, maximum_value)
    end

    continuous_mean = Float64[]
    continuous_std = Float64[]
    for column_index in continuous_indices
        column_values = Float64[]
        for row_index in 1:nrows
            cell = _parse_trait_cell(data[row_index, column_index])
            cell === missing && continue
            value = Float64(cell)
            isfinite(value) && push!(column_values, value)
        end
        mean_value, std_value = _fit_continuous_statistics(column_values)
        push!(continuous_mean, mean_value)
        push!(continuous_std, std_value)
    end

    return TraitPreprocessor(collect(names), id_index, binary_indices, count_indices, continuous_indices, count_log_min, count_log_max, continuous_mean, continuous_std)
end

function _transform_binary_column(data::AbstractVector)
    values = Vector{Float64}(undef, length(data))
    for index in eachindex(data)
        cell = _parse_trait_cell(data[index])
        value = cell === missing ? 0.0 : Float64(cell)
        isfinite(value) || (value = 0.0)
        (value == 0.0 || value == 1.0) || throw(ArgumentError("binary columns must contain only 0/1 values"))
        values[index] = value
    end
    return values
end

function _transform_count_column(data::AbstractVector, minimum_value::Real, maximum_value::Real)
    values = Vector{Float64}(undef, length(data))
    for index in eachindex(data)
        cell = _parse_trait_cell(data[index])
        value = cell === missing ? 0.0 : Float64(cell)
        isfinite(value) || (value = 0.0)
        value < 0 && throw(ArgumentError("count columns must not contain negative values"))
        values[index] = log1p(value)
    end

    scale = maximum_value - minimum_value
    if scale == 0.0
        fill!(values, 0.0)
    else
        @inbounds for index in eachindex(values)
            values[index] = (values[index] - minimum_value) / scale
        end
    end
    return values
end

function _transform_continuous_column(data::AbstractVector, mean_value::Real, std_value::Real)
    values = Vector{Float64}(undef, length(data))
    if std_value == 0.0
        fill!(values, 0.0)
        return values
    end

    @inbounds for index in eachindex(data)
        cell = _parse_trait_cell(data[index])
        if cell === missing
            values[index] = 0.0
        else
            standardized = (Float64(cell) - mean_value) / std_value
            values[index] = isfinite(standardized) ? standardized : 0.0
        end
    end
    return values
end

function transform_trait_csv(path::AbstractString, preprocessor::TraitPreprocessor; delimiter = ',', header = true)
    data, names = _read_trait_csv(path; delimiter = delimiter, header = header)
    names == preprocessor.names || throw(ArgumentError("CSV column names do not match the fitted preprocessor"))

    nrows, ncols = size(data)
    output_columns = preprocessor.id_column === nothing ? ncols : ncols - 1
    transformed = Matrix{Float64}(undef, nrows, output_columns)
    ids = preprocessor.id_column === nothing ? nothing : Vector{String}(undef, nrows)

    if preprocessor.id_column !== nothing
        id_column = preprocessor.id_column
        for row_index in 1:nrows
            value = data[row_index, id_column]
            ids[row_index] = value === missing ? "" : string(value)
        end
    end

    output_index = 1
    count_lookup = Dict(column_index => position for (position, column_index) in pairs(preprocessor.count_columns))
    continuous_lookup = Dict(column_index => position for (position, column_index) in pairs(preprocessor.continuous_columns))

    for column_index in 1:ncols
        column_index == preprocessor.id_column && continue
        if column_index in preprocessor.binary_columns
            transformed[:, output_index] = _transform_binary_column(view(data, :, column_index))
        elseif haskey(count_lookup, column_index)
            position = count_lookup[column_index]
            transformed[:, output_index] = _transform_count_column(view(data, :, column_index), preprocessor.count_log_min[position], preprocessor.count_log_max[position])
        elseif haskey(continuous_lookup, column_index)
            position = continuous_lookup[column_index]
            transformed[:, output_index] = _transform_continuous_column(view(data, :, column_index), preprocessor.continuous_mean[position], preprocessor.continuous_std[position])
        else
            throw(ArgumentError("unexpected column index $(column_index)"))
        end
        output_index += 1
    end

    return (; data = transformed, ids = ids)
end

"""
    load_and_preprocess_traits(path; id_column, binary_columns, count_columns, continuous_columns, schema_file)

Read a CSV of trait measurements and return a named tuple with the transformed matrix and
the fitted preprocessing metadata. If `id_column` is provided, it is returned separately in
the `ids` field and excluded from the feature matrix. If the three trait column groups are
omitted, `schema_file` must be provided.
"""
function load_and_preprocess_traits(path::AbstractString;
    id_column = nothing,
    binary_columns = nothing,
    count_columns = nothing,
    continuous_columns = nothing,
    schema_file = nothing,
    relabund_file = nothing,
    relabund_id_column = :id,
    relabund_value_column = Symbol("rel.abund"),
    relabund_missing_value = 0.0,
    delimiter = ',',
    header = true,
)
    preprocessor = fit_trait_preprocessor(path; id_column = id_column, binary_columns = binary_columns, count_columns = count_columns, continuous_columns = continuous_columns, schema_file = schema_file, delimiter = delimiter, header = header)
    transformed = transform_trait_csv(path, preprocessor; delimiter = delimiter, header = header)
    relabund = relabund_file === nothing ? nothing : align_relabund(transformed.ids, relabund_file; id_column = relabund_id_column, value_column = relabund_value_column, delimiter = delimiter, header = header, missing_value = relabund_missing_value)
    return (; data = transformed.data, ids = transformed.ids, relabund, preprocessor)
end

function save_trait_schema(path::AbstractString, preprocessor::TraitPreprocessor)
    open(path, "w") do io
        if preprocessor.id_column !== nothing
            println(io, "id_column = \"$(preprocessor.names[preprocessor.id_column])\"")
        end
        println(io, "binary_columns = [")
        for column_index in preprocessor.binary_columns
            println(io, "  \"$(preprocessor.names[column_index])\",")
        end
        println(io, "]")
        println(io, "count_columns = [")
        for column_index in preprocessor.count_columns
            println(io, "  \"$(preprocessor.names[column_index])\",")
        end
        println(io, "]")
        println(io, "continuous_columns = [")
        for column_index in preprocessor.continuous_columns
            println(io, "  \"$(preprocessor.names[column_index])\",")
        end
        println(io, "]")
    end
    return path
end

function trait_preprocessing_summary(result)
    hasproperty(result, :data) || throw(ArgumentError("result must have a data field"))
    hasproperty(result, :preprocessor) || throw(ArgumentError("result must have a preprocessor field"))

    data = result.data
    preprocessor = result.preprocessor
    ids = hasproperty(result, :ids) ? result.ids : nothing
    relabund = hasproperty(result, :relabund) ? result.relabund : nothing

    return (; 
        n_rows = size(data, 1),
        n_features = size(data, 2),
        n_ids = ids === nothing ? 0 : length(ids),
        has_relabund = relabund !== nothing,
        has_id = preprocessor.id_column !== nothing,
        n_binary = length(preprocessor.binary_columns),
        n_count = length(preprocessor.count_columns),
        n_continuous = length(preprocessor.continuous_columns),
        matrix_size = size(data),
    )
end

function _csv_escape(value)
    text = string(value)
    occursin(',', text) || occursin('"', text) || occursin('\n', text) || occursin('\r', text) || return text
    return '"' * replace(text, '"' => "\"\"") * '"'
end

function save_preprocessed_traits(path::AbstractString, result)
    hasproperty(result, :data) || throw(ArgumentError("result must have a data field"))
    hasproperty(result, :preprocessor) || throw(ArgumentError("result must have a preprocessor field"))

    data = result.data
    preprocessor = result.preprocessor
    ids = hasproperty(result, :ids) ? result.ids : nothing

    feature_names = String[]
    for column_index in 1:length(preprocessor.names)
        column_index == preprocessor.id_column && continue
        push!(feature_names, String(preprocessor.names[column_index]))
    end

    open(path, "w") do io
        if ids === nothing
            println(io, join(_csv_escape.(feature_names), ','))
            for row_index in 1:size(data, 1)
                println(io, join(_csv_escape.(data[row_index, :]), ','))
            end
        else
            println(io, join((_csv_escape("id"), _csv_escape.(feature_names)...), ','))
            for row_index in 1:size(data, 1)
                row_values = Any[ids[row_index]; data[row_index, :]]
                println(io, join(_csv_escape.(row_values), ','))
            end
        end
    end

    return path
end
