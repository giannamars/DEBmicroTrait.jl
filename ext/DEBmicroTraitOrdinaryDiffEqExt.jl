module DEBmicroTraitOrdinaryDiffEqExt

using DEBmicroTrait
using OrdinaryDiffEq
using SciMLBase
using SciMLSensitivity

function DEBmicroTrait.build_problem(model::SoilModel, u0::AbstractVector, tspan::Tuple{<:Real,<:Real}; problem_constructor = ODEProblem)
    problem_constructor(DEBmicroTrait.rhs!, DEBmicroTrait.initial_state(model.state_layout, u0), tspan, model)
end

function DEBmicroTrait.simulate(model::SoilModel, u0::AbstractVector, tspan::Tuple{<:Real,<:Real};
    saveat = nothing,
    alg = Tsit5(),
    sensealg = InterpolatingAdjoint(),
    problem_constructor = ODEProblem,
    solve_function = SciMLBase.solve,
)
    problem = DEBmicroTrait.build_problem(model, u0, tspan; problem_constructor = problem_constructor)
    if saveat === nothing
        return solve_function(problem, alg; sensealg = sensealg)
    end
    return solve_function(problem, alg; saveat = saveat, sensealg = sensealg)
end

end
