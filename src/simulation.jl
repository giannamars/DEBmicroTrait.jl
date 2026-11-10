initial_state(layout::StateLayout, state::AbstractVector) = length(state) == layout.dim ? collect(state) : throw(ArgumentError("state length must match layout dimension"))

struct SciMLBackend{P,S,A,SA}
    problem_constructor::P
    solve_function::S
    algorithm::A
    sense_algorithm::SA
end

function build_problem(model::SoilModel, u0::AbstractVector, tspan::Tuple{<:Real,<:Real}; problem_constructor = nothing)
    problem_constructor === nothing && throw(ArgumentError("pass a SciML ODEProblem constructor to build_problem"))
    return problem_constructor(rhs!, initial_state(model.state_layout, u0), tspan, model)
end

function simulate(model::SoilModel, u0::AbstractVector, tspan::Tuple{<:Real,<:Real};
    problem_constructor = nothing,
    solve_function = nothing,
    algorithm = nothing,
    sense_algorithm = nothing,
    saveat = nothing,
)
    solve_function === nothing && throw(ArgumentError("pass a SciML solve function to simulate"))
    problem = build_problem(model, u0, tspan; problem_constructor = problem_constructor)

    if saveat === nothing
        return solve_function(problem, algorithm; sensealg = sense_algorithm)
    end

    return solve_function(problem, algorithm; saveat = saveat, sensealg = sense_algorithm)
end

trajectory_loss(model::SoilModel, u0::AbstractVector, tspan::Tuple{<:Real,<:Real}, target; kwargs...) = begin
    solution = simulate(model, u0, tspan; kwargs...)
    prediction = solution.u[end]
    return loss(prediction, target, model)
end

