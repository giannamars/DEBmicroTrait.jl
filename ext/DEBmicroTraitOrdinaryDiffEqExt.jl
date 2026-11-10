module DEBmicroTraitOrdinaryDiffEqExt

using DEBmicroTrait
using OrdinaryDiffEq
using SciMLBase
using SciMLSensitivity

function DEBmicroTrait.simulate(model::SoilModel, u0::AbstractVector, tspan::Tuple{<:Real,<:Real};
    saveat = nothing,
    alg = Tsit5(),
    sensealg = InterpolatingAdjoint(),
    solve_function = SciMLBase.solve,
)
    problem = ODEProblem(
        DEBmicroTrait.sciml_rhs!,
        DEBmicroTrait.initial_state(model.state_layout, u0),
        tspan,
        DEBmicroTrait.parameter_vector(model),
    )
    if saveat === nothing
        return solve_function(problem, alg; sensealg = sensealg)
    end
    return solve_function(problem, alg; saveat = saveat, sensealg = sensealg)
end

end
