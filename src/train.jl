struct TinyTrainingParameters{K,Y,D}
    k_uptake::K
    yield::Y
    k_decay::D
end

TinyTrainingParameters(values::NamedTuple) = TinyTrainingParameters(values.k_uptake, values.yield, values.k_decay)

struct TinyTrainingProblem{L,U,T,Tg,P,S,Alg,SA}
    layout::L
    u0::U
    tspan::T
    target::Tg
    initial_parameters::P
    step_size::S
    alg::Alg
    sensealg::SA
end

struct TinyGradientBackend{F}
    gradient_function::F
end

function TinyTrainingProblem(layout::StateLayout, u0::AbstractVector, tspan::Tuple{<:Real,<:Real}, target::AbstractVector;
    initial_parameters = TinyTrainingParameters(0.5, 0.4, 0.1),
    step_size = 0.1,
    alg = nothing,
    sensealg = nothing,
)
    normalized_parameters = initial_parameters isa TinyTrainingParameters ? initial_parameters : tiny_parameters(initial_parameters)
    TinyTrainingProblem(layout, collect(u0), tspan, collect(target), normalized_parameters, step_size, alg, sensealg)
end

parameter_vector(parameters::TinyTrainingParameters) = [parameters.k_uptake, parameters.yield, parameters.k_decay]

tiny_parameters(parameters::AbstractVector) = TinyTrainingParameters(parameters[1], parameters[2], parameters[3])
tiny_parameters(parameters::NamedTuple) = TinyTrainingParameters(parameters)
tiny_parameters(parameters::TinyTrainingParameters) = parameters

function tiny_model(problem::TinyTrainingProblem, parameters::TinyTrainingParameters)
    build_model(problem.layout, build_parameters((k_uptake = parameters.k_uptake, yield = parameters.yield, k_decay = parameters.k_decay)))
end

function tiny_loss(problem::TinyTrainingProblem, parameters::TinyTrainingParameters)
    model = tiny_model(problem, parameters)
    simulation = simulate(model, problem.u0, problem.tspan; saveat = 1.0, alg = problem.alg, sensealg = problem.sensealg)
    return sum(abs2, simulation.u[end] .- problem.target)
end

function tiny_training_step(problem::TinyTrainingProblem; backend::TinyGradientBackend, parameters::TinyTrainingParameters = problem.initial_parameters)
    value, gradient = backend.gradient_function(p -> tiny_loss(problem, tiny_parameters(p)), parameter_vector(parameters))
    updated_parameter_vector = parameter_vector(parameters) .- problem.step_size .* gradient[1]
    updated_parameters = tiny_parameters(updated_parameter_vector)
    updated_value = tiny_loss(problem, updated_parameters)
    return (; value, gradient = gradient[1], updated_parameters, updated_value)
end

function run_tiny_training_example(; initial_parameters = TinyTrainingParameters(0.5, 0.4, 0.1), backend = nothing, alg = nothing, sensealg = nothing)
    layout = build_state_layout([:substrate, :biomass, :mineral], [1, 1, 1])
    problem = TinyTrainingProblem(layout, [1.0, 0.2, 0.0], (0.0, 10.0), [0.75, 0.25, 0.25];
        initial_parameters = initial_parameters,
        step_size = 0.1,
        alg = alg,
        sensealg = sensealg,
    )
    backend === nothing && throw(ArgumentError("pass a TinyGradientBackend to run the tiny training example"))
    result = tiny_training_step(problem; backend = backend)
    return (; initial_parameters, result...)
end

function training_step!(model::SoilModel, batch; mass_penalty = 0.0)
    prediction, target = batch
    return loss(prediction, target, model; mass_penalty = mass_penalty)
end
