struct TraitEncoder{W<:AbstractMatrix,B<:AbstractVector,F}
    weights::W
    bias::B
    activation::F
end

unitize(values, lower, upper) = (values .- lower) ./ (upper .- lower)

devectorize(values, lower, upper) = lower .+ values .* (upper .- lower)

encode_traits(encoder::TraitEncoder, traits::AbstractVector) = encoder.activation.(encoder.weights * traits .+ encoder.bias)

softplus(x) = log1p(exp(x))

