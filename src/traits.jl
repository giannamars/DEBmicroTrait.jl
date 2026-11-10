struct TraitEncoder{W,B,H,O,N,R}
    weights::W
    biases::B
    hidden_activation::H
    output_activation::O
    parameter_names::N
    parameter_ranges::R
end

function _normalize_parameter_names(parameter_names, output_dim)
    names = parameter_names === nothing ? ntuple(i -> Symbol(:parameter_, i), output_dim) : Tuple(parameter_names)
    output_dim == length(names) || throw(ArgumentError("parameter_names must match the number of outputs"))
    return names
end

function _normalize_parameter_ranges(parameter_ranges, output_dim)
    parameter_ranges === nothing && return nothing
    if parameter_ranges isa AbstractMatrix
        size(parameter_ranges, 2) == 2 || throw(ArgumentError("parameter_ranges must have two columns"))
        size(parameter_ranges, 1) == output_dim || throw(ArgumentError("parameter_ranges must match the number of outputs"))
        return (lower = collect(parameter_ranges[:, 1]), upper = collect(parameter_ranges[:, 2]))
    elseif parameter_ranges isa NamedTuple
        hasproperty(parameter_ranges, :lower) && hasproperty(parameter_ranges, :upper) ||
            throw(ArgumentError("parameter_ranges must define lower and upper bounds"))
        length(parameter_ranges.lower) == output_dim || throw(ArgumentError("parameter_ranges lower bounds must match the number of outputs"))
        length(parameter_ranges.upper) == output_dim || throw(ArgumentError("parameter_ranges upper bounds must match the number of outputs"))
        return parameter_ranges
    else
        throw(ArgumentError("parameter_ranges must be nothing, a two-column matrix, or a NamedTuple with lower and upper"))
    end
end

function TraitEncoder(weights::AbstractMatrix, bias::AbstractVector, activation = identity; parameter_names = nothing)
    size(weights, 1) == length(bias) || throw(ArgumentError("weights and bias must agree in their output dimension"))
    names = _normalize_parameter_names(parameter_names, size(weights, 1))
    TraitEncoder((weights,), (bias,), identity, activation, names, nothing)
end

function TraitEncoder(weights, biases, hidden_activation = relu;
    output_activation = sigmoid,
    parameter_names = nothing,
    parameter_ranges = nothing,
)
    weight_layers = Tuple(weights)
    bias_layers = Tuple(biases)
    length(weight_layers) == length(bias_layers) || throw(ArgumentError("weights and biases must have the same number of layers"))
    isempty(weight_layers) && throw(ArgumentError("weights must contain at least one layer"))
    for index in 1:length(weight_layers)
        weight_layer = weight_layers[index]
        bias_layer = bias_layers[index]
        weight_layer isa AbstractMatrix || throw(ArgumentError("each weight layer must be a matrix"))
        bias_layer isa AbstractVector || throw(ArgumentError("each bias layer must be a vector"))
        size(weight_layer, 1) == length(bias_layer) || throw(ArgumentError("each weight matrix must match its bias vector"))
        if index < length(weight_layers)
            size(weight_layer, 2) == size(weight_layers[index + 1], 1) || throw(ArgumentError("adjacent layers must agree in their dimensions"))
        end
    end

    names = _normalize_parameter_names(parameter_names, size(weight_layers[end], 1))
    ranges = _normalize_parameter_ranges(parameter_ranges, size(weight_layers[end], 1))
    TraitEncoder(weight_layers, bias_layers, hidden_activation, output_activation, names, ranges)
end

unitize(values, lower, upper) = (values .- lower) ./ (upper .- lower)

devectorize(values, lower, upper) = lower .+ values .* (upper .- lower)

relu(x) = ifelse(x > zero(x), x, zero(x))
sigmoid(x) = inv(one(x) + exp(-x))

function _apply_parameter_ranges(values, parameter_ranges)
    parameter_ranges === nothing && return values
    return devectorize(values, parameter_ranges.lower, parameter_ranges.upper)
end

function trait_vector(encoder::TraitEncoder, traits::AbstractVector)
    activation = traits
    last_layer = length(encoder.weights)

    for layer in 1:(last_layer - 1)
        activation = encoder.hidden_activation.(encoder.weights[layer] * activation .+ encoder.biases[layer])
    end

    encoded = encoder.output_activation.(encoder.weights[last_layer] * activation .+ encoder.biases[last_layer])
    return _apply_parameter_ranges(encoded, encoder.parameter_ranges)
end

encode_traits(encoder::TraitEncoder, traits::AbstractVector) = trait_vector(encoder, traits)

function trait_parameters(encoder::TraitEncoder, traits::AbstractVector)
    encoded = trait_vector(encoder, traits)
    return NamedTuple{encoder.parameter_names}(Tuple(encoded))
end

softplus(x) = log1p(exp(x))
