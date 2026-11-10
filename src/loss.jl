function loss(prediction, target, model; mass_penalty = 0.0, reference_mass = nothing)
    data_term = sum(abs2, prediction .- target)

    if mass_penalty == 0.0
        return data_term
    end

    reference = reference_mass === nothing ? mass_balance(target, model) : reference_mass
    balance_term = abs2(mass_balance(prediction, model) - reference)
    return data_term + mass_penalty * balance_term
end

