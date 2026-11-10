struct TraitEncoder{W,B,F,N}
    weights::W
    bias::B
    activation::F
    parameter_names::N
end

function TraitEncoder(weights::AbstractMatrix, bias::AbstractVector, activation = identity; parameter_names = nothing)
    size(weights, 1) == length(bias) || throw(ArgumentError("weights and bias must agree in their output dimension"))
    names = parameter_names === nothing ? ntuple(i -> Symbol(:parameter_, i), size(weights, 1)) : Tuple(parameter_names)
    size(weights, 1) == length(names) || throw(ArgumentError("parameter_names must match the number of outputs"))
    TraitEncoder(weights, bias, activation, names)
end

unitize(values, lower, upper) = (values .- lower) ./ (upper .- lower)

devectorize(values, lower, upper) = lower .+ values .* (upper .- lower)

trait_vector(encoder::TraitEncoder, traits::AbstractVector) = encoder.activation.(encoder.weights * traits .+ encoder.bias)
encode_traits(encoder::TraitEncoder, traits::AbstractVector) = trait_vector(encoder, traits)

function trait_parameters(encoder::TraitEncoder, traits::AbstractVector)
    encoded = trait_vector(encoder, traits)
    return NamedTuple{encoder.parameter_names}(Tuple(encoded))
end

softplus(x) = log1p(exp(x))
