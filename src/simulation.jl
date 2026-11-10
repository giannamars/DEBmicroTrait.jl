initial_state(layout::StateLayout, state::AbstractVector) = length(state) == layout.dim ? collect(state) : throw(ArgumentError("state length must match layout dimension"))

function parameter_vector(model::SoilModel)
    return [
        parameter(model.parameters, :k_uptake),
        parameter(model.parameters, :yield),
        parameter(model.parameters, :k_decay),
    ]
end

function sciml_rhs!(du, u, p, t)
    k_uptake, yield, k_decay = p
    uptake = k_uptake * u[1] * u[2]
    decay = k_decay * u[2]

    du[1] = -uptake
    du[2] = yield * uptake - decay
    du[3] = (1 - yield) * uptake + decay

    return du
end

struct SciMLBackend{P,S,A,SA}
    problem_constructor::P
    solve_function::S
    algorithm::A
    sense_algorithm::SA
end

simulate(args...; kwargs...) = throw(ArgumentError("pass a SciML solve function to simulate"))

trajectory_loss(model, u0::AbstractVector, tspan::Tuple{<:Real,<:Real}, target; kwargs...) = begin
    solution = simulate(model, u0, tspan; kwargs...)
    prediction = solution.u[end]
    return loss(prediction, target, model)
end

function gradient_smoke_test(model, u0::AbstractVector, tspan::Tuple{<:Real,<:Real}, target;
    gradient_function = nothing,
    loss_function = trajectory_loss,
    kwargs...,
)
    gradient_function === nothing && throw(ArgumentError("pass a gradient_function to run the smoke test"))
    return gradient_function(loss_function, model, u0, tspan, target; kwargs...)
end

