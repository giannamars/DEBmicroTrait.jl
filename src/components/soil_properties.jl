function aqueous_diffusivity(molecular_weight::Array{Float64,1})
    # Perry and Hilton (1973)
    D_S = @. 4.36e-9*molecular_weight^(-0.386)
end
