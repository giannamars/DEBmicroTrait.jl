const HYSOMI_TRAIT_INPUT_DIM = 38

hysomi_trait_encoder(weights, biases; hidden_activation = relu, output_activation = sigmoid, dropout_probability = 0.0, batchnorm_statistics = nothing) =
    TraitEncoder(weights, biases, hidden_activation;
        output_activation = output_activation,
        parameter_names = HYSOMI_BIOPAR_NAMES,
        parameter_ranges = HYSOMI_BIOPAR_RANGES,
        dropout_probability = dropout_probability,
        batchnorm_statistics = batchnorm_statistics,
    )

function hysomi_parameter_constraint_loss(parameters::HysomiParameters)
    biopar = parameters.biopar
    return constraint_penalty(biopar[11], biopar[1]) + constraint_penalty(biopar[1], biopar[6]) + constraint_penalty(biopar[6], biopar[5])
end

function hysomi_parameter_constraint_loss(parameters::NamedTuple)
    return hysomi_parameter_constraint_loss(HysomiParameters(parameters))
end

function hysomi_parameter_constraint_loss(parameters::AbstractVector)
    length(parameters) >= 6 || throw(ArgumentError("parameter vector must contain at least the 13 biopar values"))
    return constraint_penalty(parameters[11], parameters[1]) + constraint_penalty(parameters[1], parameters[6]) + constraint_penalty(parameters[6], parameters[5])
end

function softargmax1d(values; beta = 1e4)
    if values isa AbstractVector
        shifted = values .- maximum(values)
        probabilities = exp.(beta .* shifted)
        probabilities = probabilities ./ sum(probabilities)
        positions = collect(0:(length(values) - 1))
        return sum(probabilities .* positions)
    elseif ndims(values) == 2
        shifted = values .- maximum(values; dims = 1)
        probabilities = exp.(beta .* shifted)
        probabilities = probabilities ./ sum(probabilities; dims = 1)
        positions = collect(0:(size(values, 1) - 1))
        return vec(sum(probabilities .* positions, dims = 1))
    else
        throw(ArgumentError("softargmax1d expects a vector or a matrix"))
    end
end

peak_width(values; threshold = 0.5, sharpness = 1e3) = sum(1.0 ./ (1.0 .+ exp.(-sharpness .* (values .- threshold))))

function constraint_penalty(lhs, rhs)
    value_type = lhs isa Number ? typeof(lhs) : eltype(lhs)
    penalty = (lhs .- rhs) ./ (lhs .+ eps(value_type))
    if penalty isa Number
        return max(penalty, zero(penalty))
    end
    penalty = max.(penalty, zero(eltype(penalty)))
    return sum(penalty) / length(penalty)
end

function hysomi_process_constraint_loss(prediction; sharpness = 1e8)
    ndims(prediction) == 2 || throw(ArgumentError("prediction must be a matrix with states along the first dimension"))
    size(prediction, 1) == length(HYSOMI_STATE_NAMES) || throw(ArgumentError("prediction must have six state rows"))

    B_a = view(prediction, 1, :)
    B_i = view(prediction, 2, :)
    C_l_L = view(prediction, 3, :)
    C_H = view(prediction, 4, :)
    C_s_L = view(prediction, 5, :)

    Btot = B_a .+ B_i
    Ctot = C_l_L .+ C_H .+ C_s_L
    B_a_max = maximum(B_a)
    B_a_ratio = B_a ./ max(B_a_max, eps(eltype(B_a)))
    B_a_width = peak_width(B_a_ratio; threshold = 0.5, sharpness = sharpness)

    penalties = (
        constraint_penalty(Btot, 3.180),
        constraint_penalty(0.05, Btot),
        constraint_penalty(Ctot, 65.0),
        constraint_penalty(1.5, Ctot),
        constraint_penalty(0.5, maximum(Btot)),
        constraint_penalty(0.5 * B_a_max, B_i[end]),
        constraint_penalty(softargmax1d(B_a; beta = sharpness), 48.0),
        constraint_penalty(15.0, softargmax1d(B_a; beta = sharpness)),
        constraint_penalty(B_a_width, 48.0),
        constraint_penalty(5.0, B_a_width),
        constraint_penalty(0.9 * C_H[1], C_H[end]),
    )

    return sum(penalties)
end

function hysomi_trajectory_matrix(solution)
    if hasproperty(solution, :u) && !isempty(solution.u)
        return hcat(solution.u...)
    elseif solution isa AbstractVector
        return reshape(collect(solution), :, 1)
    elseif solution isa AbstractMatrix
        return solution
    else
        throw(ArgumentError("solution must be an ODE solution or a state-by-time matrix"))
    end
end

function hysomi_trait_parameters_from_state(template::TraitEncoder, state, traits)
    state_vector = collect(state)
    activation = collect(traits)
    offset = 1
    layer_count = length(template.weights)

    for index in eachindex(template.weights)
        weight = template.weights[index]
        count = length(weight)
        weight_matrix = reshape(view(state_vector, offset:offset + count - 1), size(weight))
        offset += count

        bias = template.biases[index]
        count = length(bias)
        bias_vector = view(state_vector, offset:offset + count - 1)
        offset += count

        activation = weight_matrix * activation .+ bias_vector
        if index < layer_count
            activation = template.hidden_activation.(activation)
            if template.batchnorm_statistics !== nothing
                activation = _apply_batchnorm(activation, template.batchnorm_statistics[index])
            end
            activation = _apply_dropout(activation, template.dropout_probability; training = false)
        else
            activation = template.output_activation.(activation)
        end
    end

    offset - 1 == length(state_vector) || throw(ArgumentError("state does not match the template encoder shape"))

    if template.parameter_ranges !== nothing
        activation = _apply_parameter_ranges(activation, template.parameter_ranges)
    end

    return NamedTuple{template.parameter_names}(Tuple(activation))
end

function hysomi_loss_components_from_parameters(biopar, soilpar, u0::AbstractVector, tspan::Tuple{<:Real,<:Real}, target;
    data_weight = 1.0,
    parameter_weight = 0.0,
    process_weight = 0.0,
    alg = nothing,
    sensealg = nothing,
    saveat = nothing,
    state_indices = eachindex(HYSOMI_STATE_NAMES),
    kwargs...,
)
    model = build_hysomi_model(hysomi_state_layout(), build_hysomi_parameters(biopar, soilpar))
    if sensealg === nothing
        solution = simulate(model, u0, tspan; saveat = saveat, alg = alg, kwargs...)
    else
        solution = simulate(model, u0, tspan; saveat = saveat, alg = alg, sensealg = sensealg, kwargs...)
    end
    prediction = hysomi_trajectory_matrix(solution)

    data_loss = hysomi_data_loss(prediction, target; state_indices = state_indices)
    parameter_loss = hysomi_parameter_constraint_loss(model.parameters)
    process_loss = hysomi_process_constraint_loss(prediction)
    total_loss = data_weight * data_loss + parameter_weight * parameter_loss + process_weight * process_loss

    return (; total_loss, data_loss, parameter_loss, process_loss, prediction, model)
end

function hysomi_data_loss(prediction, target; state_indices = eachindex(HYSOMI_STATE_NAMES))
    pred_matrix = hysomi_trajectory_matrix(prediction)
    target_matrix = hysomi_trajectory_matrix(target)

    if size(pred_matrix, 2) != size(target_matrix, 2)
        return sum(abs2, pred_matrix[state_indices, end] .- vec(target_matrix[state_indices, end]))
    end

    size(pred_matrix) == size(target_matrix) || throw(ArgumentError("prediction and target must have matching shapes, or target must be a terminal state vector"))
    return sum(abs2, pred_matrix[state_indices, :] .- target_matrix[state_indices, :])
end

function hysomi_loss_components(encoder::TraitEncoder, traits::AbstractVector, soilpar, u0::AbstractVector, tspan::Tuple{<:Real,<:Real}, target;
    data_weight = 1.0,
    parameter_weight = 0.0,
    process_weight = 0.0,
    alg = nothing,
    sensealg = nothing,
    saveat = nothing,
    state_indices = eachindex(HYSOMI_STATE_NAMES),
    kwargs...,
)
    biopar = trait_parameters(encoder, traits)
    return hysomi_loss_components_from_parameters(biopar, soilpar, u0, tspan, target;
        data_weight = data_weight,
        parameter_weight = parameter_weight,
        process_weight = process_weight,
        alg = alg,
        sensealg = sensealg,
        saveat = saveat,
        state_indices = state_indices,
        kwargs...,
    )
end

function hysomi_batch_loss(encoder::TraitEncoder, batch; kwargs...)
    traits, target, soilpar, u0, tspan = batch
    return hysomi_loss_components(encoder, traits, soilpar, u0, tspan, target; kwargs...)
end

function hysomi_encoder_state(encoder::TraitEncoder)
    return vcat((vec(weight) for weight in encoder.weights)..., (vec(bias) for bias in encoder.biases)...)
end

function hysomi_rebuild_encoder(template::TraitEncoder, state)
    state_vector = collect(state)
    weights = Vector{AbstractMatrix{eltype(state_vector)}}(undef, length(template.weights))
    biases = Vector{AbstractVector{eltype(state_vector)}}(undef, length(template.biases))

    offset = 1
    for index in eachindex(template.weights)
        weight = template.weights[index]
        count = length(weight)
        weights[index] = reshape(view(state_vector, offset:offset + count - 1), size(weight))
        offset += count
    end

    for index in eachindex(template.biases)
        bias = template.biases[index]
        count = length(bias)
        biases[index] = collect(view(state_vector, offset:offset + count - 1))
        offset += count
    end

    offset - 1 == length(state_vector) || throw(ArgumentError("state does not match the template encoder shape"))

    return TraitEncoder(
        Tuple(weights),
        Tuple(biases),
        template.hidden_activation;
        output_activation = template.output_activation,
        parameter_names = template.parameter_names,
        parameter_ranges = template.parameter_ranges,
        dropout_probability = template.dropout_probability,
        batchnorm_statistics = template.batchnorm_statistics,
    )
end

function hysomi_generate_sample(teacher_encoder::TraitEncoder, traits::AbstractVector, soilpar, u0::AbstractVector, tspan::Tuple{<:Real, <:Real};
    alg = nothing,
    sensealg = nothing,
    saveat = nothing,
    kwargs...,
)
    biopar = trait_parameters(teacher_encoder, traits)
    model = build_hysomi_model(hysomi_state_layout(), build_hysomi_parameters(biopar, soilpar))
    if sensealg === nothing
        solution = simulate(model, u0, tspan; saveat = saveat, alg = alg, kwargs...)
    else
        solution = simulate(model, u0, tspan; saveat = saveat, alg = alg, sensealg = sensealg, kwargs...)
    end

    return (; traits = collect(traits), target = hysomi_trajectory_matrix(solution), soilpar, u0 = collect(u0), tspan)
end

function hysomi_generate_dataset(teacher_encoder::TraitEncoder, traits_batch::AbstractVector; soilpar, u0::AbstractVector, tspan::Tuple{<:Real, <:Real}, kwargs...)
    return [hysomi_generate_sample(teacher_encoder, traits, soilpar, u0, tspan; kwargs...) for traits in traits_batch]
end

function hysomi_batch_loss(state, batch::AbstractVector, template::TraitEncoder; kwargs...)
    total = zero(Float64)
    for sample in batch
        biopar = hysomi_trait_parameters_from_state(template, state, sample.traits)
        components = hysomi_loss_components_from_parameters(biopar, sample.soilpar, sample.u0, sample.tspan, sample.target; kwargs...)
        total += components.total_loss
    end
    return total / length(batch)
end

function hysomi_batch_indices(n::Integer, batch_size::Integer; shuffle::Bool = true, rng = nothing)
    n > 0 || throw(ArgumentError("n must be positive"))
    batch_size > 0 || throw(ArgumentError("batch_size must be positive"))

    indices = collect(1:n)
    if shuffle
        order = rng === nothing ? sortperm(rand(n)) : sortperm(rand(rng, n))
        indices = indices[order]
    end

    batches = Vector{Vector{Int}}()
    for first_index in 1:batch_size:n
        last_index = min(first_index + batch_size - 1, n)
        push!(batches, indices[first_index:last_index])
    end
    return batches
end

function hysomi_split_dataset(dataset::AbstractVector; validation_fraction::Real = 0.2, shuffle::Bool = true, rng = nothing)
    0 <= validation_fraction < 1 || throw(ArgumentError("validation_fraction must be in [0, 1)"))
    n = length(dataset)
    n > 0 || throw(ArgumentError("dataset must not be empty"))

    indices = collect(1:n)
    if shuffle
        order = rng === nothing ? sortperm(rand(n)) : sortperm(rand(rng, n))
        indices = indices[order]
    end

    validation_count = validation_fraction == 0 ? 0 : max(1, round(Int, n * validation_fraction))
    validation_count = min(validation_count, max(n - 1, 0))
    split_index = n - validation_count
    split_index = max(1, split_index)
    train_indices = indices[1:split_index]
    validation_indices = indices[split_index+1:end]

    return dataset[train_indices], dataset[validation_indices]
end

function hysomi_train_epoch(state, dataset::AbstractVector, template::TraitEncoder;
    batch_size::Integer = 1,
    step_size = 1e-3,
    backend = Zygote.withgradient,
    shuffle::Bool = true,
    rng = nothing,
    kwargs...,
)
    batch_size > 0 || throw(ArgumentError("batch_size must be positive"))
    batches = hysomi_batch_indices(length(dataset), batch_size; shuffle = shuffle, rng = rng)

    current_state = collect(state)
    batch_losses = Float64[]
    for batch_indices in batches
        batch = dataset[batch_indices]
        result = hysomi_training_step(current_state, batch, template; step_size = step_size, backend = backend, kwargs...)
        push!(batch_losses, result.value)
        current_state = result.updated_state
    end

    return (; state = current_state, batch_losses, loss = isempty(batch_losses) ? NaN : sum(batch_losses) / length(batch_losses))
end

function hysomi_train(state, dataset::AbstractVector, template::TraitEncoder;
    epochs::Integer = 1,
    batch_size::Integer = 1,
    step_size = 1e-3,
    backend = Zygote.withgradient,
    validation_set = nothing,
    validation_fraction::Real = 0.0,
    shuffle::Bool = true,
    rng = nothing,
    kwargs...,
)
    epochs > 0 || throw(ArgumentError("epochs must be positive"))

    train_dataset = dataset
    validation_data = validation_set
    if validation_data === nothing && validation_fraction > 0
        train_dataset, validation_data = hysomi_split_dataset(dataset; validation_fraction = validation_fraction, shuffle = shuffle, rng = rng)
    end

    train_losses = Vector{Float64}(undef, epochs)
    validation_losses = validation_data === nothing ? Float64[] : Vector{Float64}(undef, epochs)
    current_state = collect(state)

    for epoch in 1:epochs
        epoch_result = hysomi_train_epoch(current_state, train_dataset, template;
            batch_size = batch_size,
            step_size = step_size,
            backend = backend,
            shuffle = shuffle,
            rng = rng,
            kwargs...,
        )
        current_state = epoch_result.state
        train_losses[epoch] = epoch_result.loss
        if validation_data !== nothing
            validation_losses[epoch] = hysomi_batch_loss(current_state, validation_data, template; kwargs...)
        end
    end

    return (; state = current_state, train_losses, validation_losses)
end

function _hysomi_config_value(value)
    if value isa NamedTuple
        return Dict(string(key) => _hysomi_config_value(getproperty(value, key)) for key in keys(value))
    elseif value isa AbstractDict
        return Dict(string(key) => _hysomi_config_value(entry) for (key, entry) in value)
    elseif value isa AbstractArray
        return [_hysomi_config_value(entry) for entry in value]
    elseif value isa Symbol
        return String(value)
    else
        return value
    end
end

function _hysomi_dict(value)
    if value isa NamedTuple
        return Dict{String,Any}(string(key) => _hysomi_dict(getproperty(value, key)) for key in keys(value))
    elseif value isa AbstractDict
        return Dict{String,Any}(string(key) => _hysomi_dict(entry) for (key, entry) in value)
    elseif value isa AbstractArray
        return [_hysomi_dict(entry) for entry in value]
    elseif value isa Symbol
        return String(value)
    else
        return value
    end
end

function _hysomi_config_string(value)
    if value isa NamedTuple
        parts = [string(key, "=", _hysomi_config_string(getproperty(value, key))) for key in keys(value)]
        return "(" * join(parts, ",") * ")"
    elseif value isa AbstractDict
        entries = sort!(collect(value); by = first)
        parts = [string(key, "=", _hysomi_config_string(entry)) for (key, entry) in entries]
        return "{" * join(parts, ",") * "}"
    elseif value isa AbstractArray
        return "[" * join(_hysomi_config_string.(value), ",") * "]"
    else
        return sprint(show, value)
    end
end

function hysomi_training_hash(config; length_id::Integer = 12)
    length_id > 0 || throw(ArgumentError("length_id must be positive"))
    digest = bytes2hex(SHA.sha256(codeunits(_hysomi_config_string(config))))
    return digest[1:min(length_id, lastindex(digest))]
end

function hysomi_hyperpar_to_toml(NN_hyperpar, loss_hyperpar, training_hyperpar, training_outputs_library_path; length_id::Integer = 12)
    config = (; NN_hyperpars = NN_hyperpar, loss_hyperpar = loss_hyperpar, training_hyperpar = training_hyperpar)
    hash_id = hysomi_training_hash(config; length_id = length_id)
    out_dir = joinpath(String(training_outputs_library_path), "train_" * hash_id)
    mkpath(out_dir)
    path_to_file = joinpath(out_dir, "hyper_dict.toml")
    open(path_to_file, "w") do io
        TOML.print(io, _hysomi_config_value(config))
    end
    return (; config, hash_id, out_dir, path_to_file)
end

function hysomi_load_training_config(config_path::AbstractString)
    return TOML.parsefile(config_path)
end

function hysomi_save_training_result(out_dir::AbstractString, result)
    mkpath(out_dir)
    open(joinpath(out_dir, "training_result.jls"), "w") do io
        Serialization.serialize(io, result)
    end
    return joinpath(out_dir, "training_result.jls")
end

function hysomi_load_serialized_dataset(dataset_path::AbstractString)
    open(dataset_path, "r") do io
        return Serialization.deserialize(io)
    end
end

function _hysomi_training_hyperparameters(training_hyperpar)
    epochs = Int(get(training_hyperpar, "n_epochs", get(training_hyperpar, "epochs", 1)))
    batch_size = Int(get(training_hyperpar, "batch_size", 1))
    step_size = get(training_hyperpar, "lr", get(training_hyperpar, "step_size", 1e-3))
    validation_fraction = get(training_hyperpar, "validation_fraction", 0.0)
    shuffle = get(training_hyperpar, "shuffle", true)
    seed = haskey(training_hyperpar, "seed") ? training_hyperpar["seed"] : nothing
    return (; epochs, batch_size, step_size, validation_fraction, shuffle, seed)
end

function hysomi_train_from_config(config_path::AbstractString, template::TraitEncoder;
    dataset_loader = hysomi_load_serialized_dataset,
    initial_state = hysomi_encoder_state(template),
    out_dir = nothing,
    save_result::Bool = false,
    backend = Zygote.withgradient,
    kwargs...,
)
    config = hysomi_load_training_config(config_path)
    training_hyperpar = get(config, "training_hyperpar", config)
    haskey(training_hyperpar, "dataset_path") || throw(ArgumentError("training config must contain dataset_path"))

    if out_dir === nothing && haskey(training_hyperpar, "output_dir")
        out_dir = training_hyperpar["output_dir"]
    end

    dataset = dataset_loader(String(training_hyperpar["dataset_path"]))
    hyper = _hysomi_training_hyperparameters(training_hyperpar)
    state = haskey(training_hyperpar, "initial_state") ? collect(training_hyperpar["initial_state"]) : collect(initial_state)

    result = hysomi_train(state, dataset, template;
        epochs = hyper.epochs,
        batch_size = hyper.batch_size,
        step_size = hyper.step_size,
        backend = backend,
        validation_fraction = hyper.validation_fraction,
        shuffle = hyper.shuffle,
        kwargs...,
    )

    if save_result && out_dir !== nothing
        hysomi_save_training_result(String(out_dir), result)
    end

    return (; config, result, out_dir)
end

function hysomi_train_from_config(config::AbstractDict, template::TraitEncoder; kwargs...)
    temp_path = tempname()
    open(temp_path, "w") do io
        TOML.print(io, _hysomi_config_value(config))
    end
    try
        return hysomi_train_from_config(temp_path, template; kwargs...)
    finally
        isfile(temp_path) && rm(temp_path; force = true)
    end
end

function hysomi_run_experiment(
    teacher_encoder::TraitEncoder,
    template::TraitEncoder,
    train_traits::AbstractVector,
    test_traits::AbstractVector = train_traits;
    NN_hyperpar = (model = "hysomi_trait_encoder", input_dim = HYSOMI_TRAIT_INPUT_DIM),
    loss_hyperpar = (data_sel = collect(eachindex(HYSOMI_STATE_NAMES)), k = 1e3),
    training_hyperpar = NamedTuple(),
    training_outputs_library_path = pwd(),
    soilpar,
    u0::AbstractVector,
    tspan::Tuple{<:Real,<:Real},
    save_artifacts::Bool = true,
    backend = nothing,
    kwargs...,
)
    backend === nothing && throw(ArgumentError("pass a backend, for example Zygote.withgradient"))
    bundle = hysomi_hyperpar_to_toml(NN_hyperpar, loss_hyperpar, training_hyperpar, training_outputs_library_path)

    train_dataset = hysomi_generate_dataset(teacher_encoder, train_traits; soilpar = soilpar, u0 = u0, tspan = tspan, kwargs...)
    test_dataset = hysomi_generate_dataset(teacher_encoder, test_traits; soilpar = soilpar, u0 = u0, tspan = tspan, kwargs...)

    train_dataset_path = joinpath(bundle.out_dir, "train_dataset.jls")
    test_dataset_path = joinpath(bundle.out_dir, "test_dataset.jls")
    open(train_dataset_path, "w") do io
        Serialization.serialize(io, train_dataset)
    end
    open(test_dataset_path, "w") do io
        Serialization.serialize(io, test_dataset)
    end

    training_dict = Dict{String,Any}(_hysomi_dict(training_hyperpar))
    training_dict["dataset_path"] = train_dataset_path
    training_dict["output_dir"] = bundle.out_dir

    config = Dict{String,Any}(
        "NN_hyperpars" => _hysomi_dict(NN_hyperpar),
        "loss_hyperpar" => _hysomi_dict(loss_hyperpar),
        "training_hyperpar" => training_dict,
    )
    config_path = bundle.path_to_file
    open(config_path, "w") do io
        TOML.print(io, config)
    end

    train_result = hysomi_train_from_config(
        config_path,
        template;
        dataset_loader = hysomi_load_serialized_dataset,
        initial_state = hysomi_encoder_state(template),
        out_dir = bundle.out_dir,
        save_result = save_artifacts,
        backend = backend,
        kwargs...,
    )

    evaluation_loss = hysomi_batch_loss(train_result.result.state, test_dataset, template; kwargs...)
    trained_encoder = hysomi_rebuild_encoder(template, train_result.result.state)

    if save_artifacts
        open(joinpath(bundle.out_dir, "evaluation_result.jls"), "w") do io
            Serialization.serialize(io, (; evaluation_loss))
        end
    end

    return (; bundle, train_dataset, test_dataset, train_result, evaluation_loss, trained_encoder)
end

function hysomi_update_state(state, gradient, step_size)
    return state .- step_size .* gradient
end

function hysomi_training_step(state, batch::AbstractVector, template::TraitEncoder; step_size = 1e-3, backend = Zygote.withgradient, kwargs...)
    value, gradient = backend(current_state -> hysomi_batch_loss(current_state, batch, template; kwargs...), state)
    updated_state = hysomi_update_state(state, gradient[1], step_size)
    updated_value = hysomi_batch_loss(updated_state, batch, template; kwargs...)
    return (; value, gradient = gradient[1], updated_state, updated_value)
end

function hysomi_training_epoch(state, dataset::AbstractVector, template::TraitEncoder; batch_size = 1, step_size = 1e-3, backend = Zygote.withgradient, kwargs...)
    epoch_result = hysomi_train_epoch(state, dataset, template; batch_size = batch_size, step_size = step_size, backend = backend, kwargs...)
    return (; state = epoch_result.state, losses = epoch_result.batch_losses)
end
