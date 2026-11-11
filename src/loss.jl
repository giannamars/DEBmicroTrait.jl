function loss(prediction, target, model; mass_penalty = 0.0, reference_mass = nothing)
    data_term = sum(abs2, prediction .- target)

    if mass_penalty == 0.0
        return data_term
    end

    reference = reference_mass === nothing ? mass_balance(target, model) : reference_mass
    balance_term = abs2(mass_balance(prediction, model) - reference)
    return data_term + mass_penalty * balance_term
end

_mean(values) = sum(values) / length(values)
_mean_abs(values) = _mean(abs.(values))
_mean_abs2(values) = _mean(abs2.(values))

softargmax2d(values; beta = 1e7) = _softargmax2d(values; beta = beta)

function _softargmax2d(values; beta = 1e7)
    ndims(values) >= 2 || throw(ArgumentError("values must have at least two dimensions"))
    rows = size(values, ndims(values) - 1)
    cols = size(values, ndims(values))
    leading = prod(size(values)[1:end-2]; init = 1)
    reshaped = reshape(values, leading, rows * cols)
    probabilities = _softmax_lastdim(beta .* reshaped)

    row_indices = repeat(collect(0:(rows - 1)), inner = cols)
    col_indices = repeat(collect(0:(cols - 1)), outer = rows)

    result_row = ((rows - 1) .* probabilities .* row_indices') |> sum_dims_last
    result_col = ((cols - 1) .* probabilities .* col_indices') |> sum_dims_last

    return hcat(result_row, result_col)
end

function _softmax_lastdim(values)
    max_values = maximum(values; dims = ndims(values))
    shifted = values .- max_values
    weights = exp.(shifted)
    return weights ./ sum(weights; dims = ndims(values))
end

function sum_dims_last(values)
    return vec(sum(values; dims = ndims(values)))
end

function peak_width_v2(values, threshold; k = 1e3)
    if values isa AbstractVector
        return sum(1.0 ./ (1.0 .+ exp.(-k .* (values .- threshold))))
    elseif ndims(values) == 2
        return vec(sum(1.0 ./ (1.0 .+ exp.(-k .* (values .- threshold))); dims = 1))
    elseif ndims(values) == 3
        return vec(sum(1.0 ./ (1.0 .+ exp.(-k .* (values .- threshold))); dims = 2))
    else
        throw(ArgumentError("values must be a vector, matrix, or 3D array"))
    end
end

struct AutomaticWeightedLossL2{T}
    params::Vector{T}
end

AutomaticWeightedLossL2(num::Integer = 15) = AutomaticWeightedLossL2(ones(Float64, num))

function (awl::AutomaticWeightedLossL2)(losses::AbstractVector)
    length(losses) == length(awl.params) || throw(ArgumentError("loss vector must match the number of learned weights"))
    total = zero(promote_type(eltype(losses), eltype(awl.params)))
    for index in eachindex(losses)
        weight = awl.params[index]
        total += 0.5 / (weight^2) * losses[index] + log1p(weight^2)
    end
    return total
end

struct NoConstraintModule{T}
    params::Vector{T}
end

no_constraint_module(num::Integer = 2) = NoConstraintModule(ones(Float64, num))

function (no_constraint::NoConstraintModule)(losses::AbstractVector)
    isempty(losses) && throw(ArgumentError("losses must not be empty"))
    return losses[1]
end

struct LossConstraintAutoW3{I,K}
    data_sel::I
    k::K
end

loss_constraint_autoW3(data_sel, k) = LossConstraintAutoW3(collect(data_sel), k)

function _ensure_state_time_batch(values)
    if ndims(values) == 2
        return reshape(values, size(values, 1), size(values, 2), 1)
    elseif ndims(values) == 3
        return values
    else
        throw(ArgumentError("state tensors must have two or three dimensions"))
    end
end

function _ensure_parameter_batch(values)
    if values isa AbstractVector
        return reshape(values, length(values), 1)
    elseif ndims(values) == 2
        return values
    else
        throw(ArgumentError("parameter tensors must have one or two dimensions"))
    end
end

function _parameter_matrix(values, nparameters)
    if size(values, 1) == nparameters
        return values
    elseif size(values, 2) == nparameters
        return permutedims(values)
    else
        throw(ArgumentError("parameter array must contain $nparameters parameters along one axis"))
    end
end

function _parameter_transform_list(par_trans_list)
    if par_trans_list isa Tuple && length(par_trans_list) == 2
        return par_trans_list[1], par_trans_list[2]
    elseif par_trans_list isa AbstractVector && length(par_trans_list) == 2
        return par_trans_list[1], par_trans_list[2]
    else
        throw(ArgumentError("par_trans_list must contain a scale and an offset"))
    end
end

function _unnormalize_parameters(values, par_trans_list)
    scale, offset = _parameter_transform_list(par_trans_list)
    if values isa AbstractVector
        return (values .- offset) ./ scale
    else
        scale_matrix = reshape(scale, :, 1)
        offset_matrix = reshape(offset, :, 1)
        return (values .- offset_matrix) ./ scale_matrix
    end
end

function _constraint(a, b)
    penalty = (a .- b) ./ a
    if penalty isa Number
        return max(penalty, zero(penalty))
    else
        return _mean(max.(penalty, zero(eltype(penalty))))
    end
end

function (loss_module::LossConstraintAutoW3)(pred, data, par, true_par, par_trans_list)
    prediction = _ensure_state_time_batch(pred)
    target = _ensure_state_time_batch(data)
    parameters = _parameter_matrix(_ensure_parameter_batch(par), 13)
    true_parameters = _parameter_matrix(_ensure_parameter_batch(true_par), 13)

    length(loss_module.data_sel) > 0 || throw(ArgumentError("data_sel must not be empty"))
    maximum(loss_module.data_sel) <= size(prediction, 1) || throw(ArgumentError("data_sel is out of bounds"))

    selected_prediction = prediction[loss_module.data_sel, :, :]
    selected_target = target[loss_module.data_sel, :, :]

    relu_loss = x -> max(x, zero(x))
    data_loss = _mean_abs2(selected_prediction .- selected_target)

    loss_pred_l1 = _mean_abs(selected_prediction .- selected_target)
    loss_pred_unmasked = _mean_abs2(prediction[1:min(5, size(prediction, 1)), :, :] .- target[1:min(5, size(target, 1)), :, :])
    loss_pred_l1_unmasked = _mean_abs(prediction[1:min(5, size(prediction, 1)), :, :] .- target[1:min(5, size(target, 1)), :, :])

    par_unranged = _unnormalize_parameters(parameters, par_trans_list)
    true_par_unranged = _unnormalize_parameters(true_parameters, par_trans_list)
    hidden_par_loss = _mean_abs(par_unranged .- true_par_unranged)
    theta_l1_tensor = vec(sum(abs.(par_unranged .- true_par_unranged); dims = 2) ./ size(par_unranged, 2))

    if size(prediction, 1) < 6 || size(target, 1) < 6
        throw(ArgumentError("prediction and target must have at least six state rows"))
    end

    B_a = view(prediction, 1, :, :)
    B_i = view(prediction, 2, :, :)
    C_l_L = view(prediction, 3, :, :)
    C_H = view(prediction, 4, :, :)
    C_s_L = view(prediction, 5, :, :)

    Btot = B_a .+ B_i
    Ctot = C_l_L .+ C_H .+ C_s_L
    B_a_max = vec(maximum(B_a; dims = 1))
    B_a_ratio = B_a ./ max.(B_a_max', eps(eltype(B_a_max)))
    B_a_width = peak_width_v2(B_a_ratio, 0.5; k = loss_module.k)

    par_constraint_tensor = [
        _constraint(parameters[11, :], parameters[1, :]),
        _constraint(parameters[1, :], parameters[6, :]),
        _constraint(parameters[6, :], parameters[5, :]),
    ]

    process_constraint_tensor = [
        _constraint(Btot, 3.180),
        _constraint(0.05, Btot),
        _constraint(Ctot, 65.0),
        _constraint(1.5, Ctot),
        _constraint(0.5, vec(maximum(Btot; dims = 1))),
        _constraint(0.5 .* B_a_max, B_i[end, :]),
        _constraint(softargmax1d(B_a; beta = loss_module.k), 48.0),
        _constraint(15.0, softargmax1d(B_a; beta = loss_module.k)),
        _constraint(B_a_width, 48.0),
        _constraint(5.0, B_a_width),
        _constraint(0.9 .* C_H[1, :], C_H[end, :]),
    ]

    constraint_tensor = vcat(par_constraint_tensor, process_constraint_tensor)
    loss_tensor = vcat(data_loss, constraint_tensor)
    ana_loss = vcat(
        loss_pred_l1,
        hidden_par_loss,
        loss_pred_unmasked,
        loss_pred_l1_unmasked,
        _mean_abs2(prediction[6, :, :] .- target[6, :, :]),
        _mean_abs2(prediction[1, :, :] .- target[1, :, :]),
        _mean_abs2(prediction[2, :, :] .- target[2, :, :]),
        _mean_abs2(prediction[3, :, :] .- target[3, :, :]),
        _mean_abs2(prediction[4, :, :] .- target[4, :, :]),
        _mean_abs2(prediction[5, :, :] .- target[5, :, :]),
    )

    return (; loss_tensor, loss_pred = data_loss, constraint_tensor, ana_loss, theta_L1_tensor = theta_l1_tensor)
end

