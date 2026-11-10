struct TinyTrainingProblem{L,U,T,Tg,K,S,Alg,SA}
    layout::L
    u0::U
    tspan::T
    target::Tg
    k_uptake::K
    k_decay::K
    step_size::S
    alg::Alg
    sensealg::SA
end

struct TinyGradientBackend{F}
    gradient_function::F
end

function TinyTrainingProblem(layout::StateLayout, u0::AbstractVector, tspan::Tuple{<:Real,<:Real}, target::AbstractVector;
    k_uptake = 0.5,
    k_decay = 0.1,
    step_size = 0.1,
    alg = nothing,
    sensealg = nothing,
)
    TinyTrainingProblem(layout, collect(u0), tspan, collect(target), k_uptake, k_decay, step_size, alg, sensealg)
end

function tiny_loss(problem::TinyTrainingProblem, yield)
    model = build_model(problem.layout, build_parameters((k_uptake = problem.k_uptake, yield = yield, k_decay = problem.k_decay)))
    simulation = simulate(model, problem.u0, problem.tspan; saveat = 1.0, alg = problem.alg, sensealg = problem.sensealg)
    return sum(abs2, simulation.u[end] .- problem.target)
end

function tiny_training_step(problem::TinyTrainingProblem, yield; backend::TinyGradientBackend)
    value, gradient = backend.gradient_function(y -> tiny_loss(problem, y), yield)
    updated_yield = yield - problem.step_size * gradient[1]
    updated_value = tiny_loss(problem, updated_yield)
    return (; value, gradient = gradient[1], updated_yield, updated_value)
end

function run_tiny_training_example(; initial_yield = 0.4, backend = nothing, alg = nothing, sensealg = nothing)
    layout = build_state_layout([:substrate, :biomass, :mineral], [1, 1, 1])
    problem = TinyTrainingProblem(layout, [1.0, 0.2, 0.0], (0.0, 10.0), [0.75, 0.25, 0.25];
        k_uptake = 0.5,
        k_decay = 0.1,
        step_size = 0.1,
        alg = alg,
        sensealg = sensealg,
    )
    backend === nothing && throw(ArgumentError("pass a TinyGradientBackend to run the tiny training example"))
    result = tiny_training_step(problem, initial_yield; backend = backend)
    return (; initial_yield, result...)
end

function training_step!(model::SoilModel, batch; mass_penalty = 0.0)
    prediction, target = batch
    return loss(prediction, target, model; mass_penalty = mass_penalty)
end
