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

softargmax1d(values; beta = 1e4) = begin
    shifted = values .- maximum(values)
    probabilities = exp.(beta .* shifted)
    probabilities = probabilities ./ sum(probabilities)
    positions = collect(0:length(values) - 1)
    return sum(probabilities .* positions)
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
    batch_size > 0 || throw(ArgumentError("batch_size must be positive"))
    losses = Float64[]
    current_state = state
    for first_index in 1:batch_size:length(dataset)
        last_index = min(first_index + batch_size - 1, length(dataset))
        batch = dataset[first_index:last_index]
        result = hysomi_training_step(current_state, batch, template; step_size = step_size, backend = backend, kwargs...)
        push!(losses, result.value)
        current_state = result.updated_state
    end
    return (; state = current_state, losses)
end
