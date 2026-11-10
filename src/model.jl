function rhs!(du, u, model::SoilModel, t)
    layout = model.state_layout
    parameters = model.parameters

    substrate = view(u, state_range(layout, :substrate))
    biomass = view(u, state_range(layout, :biomass))
    mineral = view(u, state_range(layout, :mineral))

    k_uptake = parameter(parameters, :k_uptake)
    yield = parameter(parameters, :yield)
    k_decay = parameter(parameters, :k_decay)

    uptake = k_uptake * substrate[1] * biomass[1]
    decay = k_decay * biomass[1]

    du[state_range(layout, :substrate)] .= -uptake
    du[state_range(layout, :biomass)] .= yield * uptake - decay
    du[state_range(layout, :mineral)] .= (1 - yield) * uptake + decay

    return du
end

observables(u, model::SoilModel) = copy(u)

mass_balance(u, model::SoilModel) = sum(u)

