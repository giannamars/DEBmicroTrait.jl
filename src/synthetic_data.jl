function hysomi_parameter_constraint(theta::AbstractVector)
    length(theta) >= 11 || throw(ArgumentError("parameter vector must contain at least 11 values"))
    mu = theta[1]
    r_B = theta[5]
    d_B = theta[6]
    m_max = theta[11]
    return m_max < mu && mu < d_B && d_B < r_B
end

function hysomi_latin_hypercube(ranges::AbstractMatrix, n::Integer; rng = nothing)
    n > 0 || throw(ArgumentError("n must be positive"))
    size(ranges, 2) == 2 || throw(ArgumentError("ranges must have two columns"))

    dimensions = size(ranges, 1)
    samples = Matrix{Float64}(undef, n, dimensions)

    for dimension in 1:dimensions
        lower = ranges[dimension, 1]
        upper = ranges[dimension, 2]
        draws = rng === nothing ? rand(n) : rand(rng, n)
        permutation = sortperm(rng === nothing ? rand(n) : rand(rng, n))
        strata = (collect(0:n-1) .+ draws) ./ n
        strata = strata[permutation]
        samples[:, dimension] = lower .+ strata .* (upper - lower)
    end

    return samples
end

function hysomi_constrained_sampling(
    ranges::AbstractMatrix,
    n::Integer;
    rng = nothing,
    oversample_factor::Integer = 7,
    constraint = hysomi_parameter_constraint,
    max_attempts::Integer = 10,
)
    n > 0 || throw(ArgumentError("n must be positive"))
    oversample_factor > 0 || throw(ArgumentError("oversample_factor must be positive"))
    max_attempts > 0 || throw(ArgumentError("max_attempts must be positive"))

    samples = Matrix{Float64}(undef, n, size(ranges, 1))
    accepted = 0

    for _ in 1:max_attempts
        pool = hysomi_latin_hypercube(ranges, max(n, oversample_factor * n); rng = rng)
        for row in eachrow(pool)
            constraint(row) || continue
            accepted += 1
            @views samples[accepted, :] .= row
            accepted == n && return samples
        end
    end

    throw(ArgumentError("could not draw $n constrained samples after $max_attempts attempts"))
end

hysomi_traits_from_parameters(trait_map::TraitEncoder, parameters) = trait_vector(trait_map, parameters)
hysomi_traits_from_parameters(trait_map::Function, parameters) = collect(trait_map(parameters))

function hysomi_generate_synthetic_sample(trait_map, teacher_encoder::TraitEncoder, parameters, soilpar, u0::AbstractVector, tspan::Tuple{<:Real, <:Real}; kwargs...)
    traits = hysomi_traits_from_parameters(trait_map, parameters)
    sample = hysomi_generate_sample(teacher_encoder, traits, soilpar, u0, tspan; kwargs...)
    return (; sample..., parameters = collect(parameters))
end

function hysomi_process_filter(target::AbstractMatrix;
    min_biomass = 0.05,
    max_biomass = 3.125,
    min_carbon = 1.5,
    max_carbon = 65.0,
    min_total_biomass = 0.5,
    min_peak_index = 15,
    max_peak_index = 48,
    min_width = 5,
    max_width = 48,
    threshold = 0.5,
)
    size(target, 1) == length(HYSOMI_STATE_NAMES) || throw(ArgumentError("target must have six state rows"))
    all(isfinite, target) || return false

    B_a = view(target, 1, :)
    B_i = view(target, 2, :)
    C_l_L = view(target, 3, :)
    C_H = view(target, 4, :)
    C_s_L = view(target, 5, :)

    Btot = B_a .+ B_i
    Ctot = C_l_L .+ C_H .+ C_s_L
    B_a_max = maximum(B_a)
    peak_index = argmax(B_a) - 1
    width = count(value -> value > threshold, B_a ./ max(B_a_max, eps(eltype(B_a))))

    return all(
        (
            all(value -> min_biomass < value < max_biomass, Btot),
            all(value -> min_carbon < value < max_carbon, Ctot),
            maximum(Btot) > min_total_biomass,
            B_i[end] > 0.5 * B_a_max,
            peak_index < max_peak_index,
            peak_index > min_peak_index,
            width < max_width,
            width > min_width,
        ),
    )
end

function hysomi_process_filter(dataset::AbstractVector; kwargs...)
    mask = [hysomi_process_filter(sample.target; kwargs...) for sample in dataset]
    filtered = [dataset[index] for index in eachindex(dataset) if mask[index]]
    return (; filtered, mask)
end

function hysomi_generate_synthetic_dataset(
    trait_map,
    teacher_encoder::TraitEncoder,
    parameter_ranges::AbstractMatrix,
    n::Integer;
    soilpar,
    u0::AbstractVector,
    tspan::Tuple{<:Real, <:Real},
    rng = nothing,
    oversample_factor::Integer = 7,
    constraint = hysomi_parameter_constraint,
    kwargs...,
)
    parameters = hysomi_constrained_sampling(
        parameter_ranges,
        n;
        rng = rng,
        oversample_factor = oversample_factor,
        constraint = constraint,
    )

    return [
        hysomi_generate_synthetic_sample(trait_map, teacher_encoder, view(parameters, index, :), soilpar, u0, tspan; kwargs...)
        for index in 1:size(parameters, 1)
    ]
end

function hysomi_generate_python_style_synthetic_dataset(
    trait_map,
    teacher_encoder::TraitEncoder,
    parameter_ranges::AbstractMatrix,
    n::Integer;
    soilpar,
    u0::AbstractVector,
    tspan::Tuple{<:Real, <:Real},
    rng = nothing,
    oversample_factor::Integer = 7,
    constraint = hysomi_parameter_constraint,
    kwargs...,
)
    n > 0 || throw(ArgumentError("n must be positive"))

    raw_parameters = hysomi_constrained_sampling(
        parameter_ranges,
        max(n, oversample_factor * n);
        rng = rng,
        oversample_factor = oversample_factor,
        constraint = constraint,
    )

    raw_dataset = [
        hysomi_generate_synthetic_sample(trait_map, teacher_encoder, view(raw_parameters, index, :), soilpar, u0, tspan; kwargs...)
        for index in 1:size(raw_parameters, 1)
    ]

    mask = [hysomi_process_filter(sample.target) for sample in raw_dataset]
    filtered_parameters = [raw_dataset[index].parameters for index in eachindex(raw_dataset) if mask[index]]
    filtered_dataset = [raw_dataset[index] for index in eachindex(raw_dataset) if mask[index]]

    return (; raw_parameters, raw_dataset, mask, filtered_parameters, filtered_dataset)
end

function _hysomi_write_serialized(path::AbstractString, value)
    open(path, "w") do io
        Serialization.serialize(io, value)
    end
    return path
end

function hysomi_save_synthetic_dataset_bundle(out_dir::AbstractString, bundle; config = nothing, overwrite::Bool = true)
    mkpath(out_dir)

    paths = Dict{String,String}()
    artifact_paths = (
        bundle = joinpath(out_dir, "synthetic_data_bundle.jls"),
        raw_parameters = joinpath(out_dir, "raw_parameters.jls"),
        raw_dataset = joinpath(out_dir, "raw_dataset.jls"),
        mask = joinpath(out_dir, "mask.jls"),
        filtered_parameters = joinpath(out_dir, "filtered_parameters.jls"),
        filtered_dataset = joinpath(out_dir, "filtered_dataset.jls"),
    )

    for (name, path) in pairs(artifact_paths)
        if !overwrite && isfile(path)
            throw(ArgumentError("refusing to overwrite existing file $path"))
        end
        value = name == :bundle ? bundle : getproperty(bundle, name)
        _hysomi_write_serialized(path, value)
        paths[String(name)] = path
    end

    if config !== nothing
        config_path = joinpath(out_dir, "synthetic_data.toml")
        if !overwrite && isfile(config_path)
            throw(ArgumentError("refusing to overwrite existing file $config_path"))
        end
        open(config_path, "w") do io
            TOML.print(io, _hysomi_config_value(config))
        end
        paths["config"] = config_path
    end

    sorted_paths = [Symbol(key) => value for (key, value) in sort(collect(paths); by = first)]
    summary = (
        out_dir = out_dir,
        raw_count = length(bundle.raw_dataset),
        filtered_count = length(bundle.filtered_dataset),
        paths = (; sorted_paths...),
    )

    _hysomi_write_serialized(joinpath(out_dir, "synthetic_data_summary.jls"), summary)
    return summary
end

function hysomi_load_synthetic_dataset_bundle(out_dir::AbstractString)
    bundle_path = joinpath(out_dir, "synthetic_data_bundle.jls")
    summary_path = joinpath(out_dir, "synthetic_data_summary.jls")
    config_path = joinpath(out_dir, "synthetic_data.toml")

    isfile(bundle_path) || throw(ArgumentError("missing synthetic dataset bundle at $bundle_path"))
    isfile(summary_path) || throw(ArgumentError("missing synthetic dataset summary at $summary_path"))

    bundle = open(bundle_path, "r") do io
        Serialization.deserialize(io)
    end
    summary = open(summary_path, "r") do io
        Serialization.deserialize(io)
    end

    config = isfile(config_path) ? TOML.parsefile(config_path) : nothing
    return (; bundle, summary, config)
end
