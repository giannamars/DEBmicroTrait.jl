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
    return constraint_penalty(biopar.m_max, biopar.mu) + constraint_penalty(biopar.mu, biopar.r_deact) + constraint_penalty(biopar.r_deact, biopar.r_react)
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
    probabilities ./= sum(probabilities)
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

function hysomi_data_loss(prediction, target; state_indices = eachindex(HYSOMI_STATE_NAMES))
    pred_matrix = hysomi_trajectory_matrix(prediction)
    target_matrix = hysomi_trajectory_matrix(target)

    if size(target_matrix, 2) == 1 && size(pred_matrix, 2) != 1
        return sum(abs2, pred_matrix[state_indices, end] .- vec(target_matrix[state_indices, 1]))
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
)
    biopar = trait_parameters(encoder, traits)
    model = build_hysomi_model(hysomi_state_layout(), build_hysomi_parameters(biopar, soilpar))
    solution = simulate(model, u0, tspan; saveat = saveat, alg = alg, sensealg = sensealg)
    prediction = hysomi_trajectory_matrix(solution)

    data_loss = hysomi_data_loss(prediction, target; state_indices = state_indices)
    parameter_loss = hysomi_parameter_constraint_loss(model.parameters)
    process_loss = hysomi_process_constraint_loss(prediction)
    total_loss = data_weight * data_loss + parameter_weight * parameter_loss + process_weight * process_loss

    return (; total_loss, data_loss, parameter_loss, process_loss, prediction, model)
end

function hysomi_batch_loss(encoder::TraitEncoder, batch; kwargs...)
    traits, target, soilpar, u0, tspan = batch
    return hysomi_loss_components(encoder, traits, soilpar, u0, tspan, target; kwargs...)
end
