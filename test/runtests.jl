using Test
using DEBmicroTrait

if Base.find_package("OrdinaryDiffEq") !== nothing
    @eval using OrdinaryDiffEq
end

if Base.find_package("SciMLSensitivity") !== nothing
    @eval using SciMLSensitivity
end

if Base.find_package("Zygote") !== nothing
    @eval using Zygote
end

@testset "layout" begin
    layout = build_state_layout([:substrate, :biomass, :mineral], [1, 1, 1])
    @test layout.dim == 3
    @test state_range(layout, :biomass) == 2:2
end

@testset "trait encoding" begin
    encoder = TraitEncoder([2.0 0.0; 0.0 2.0], [1.0, -1.0], identity)
    traits = [3.0, 4.0]
    @test encode_traits(encoder, traits) == [7.0, 7.0]
    @test unitize([5.0], 0.0, 10.0) == [0.5]
    @test devectorize([0.5], 0.0, 10.0) == [5.0]
end

@testset "mass balance" begin
    layout = build_state_layout([:substrate, :biomass, :mineral], [1, 1, 1])
    model = build_model(layout, build_parameters((k_uptake = 0.5, yield = 0.4, k_decay = 0.1)))
    u = [1.0, 0.2, 0.0]
    du = zeros(3)

    rhs!(du, u, model, 0.0)

    @test du[1] < 0
    @test isapprox(sum(du), 0.0; atol = 1e-12)
    @test isapprox(mass_balance(u, model), 1.2)
end

@testset "loss" begin
    layout = build_state_layout([:substrate, :biomass, :mineral], [1, 1, 1])
    model = build_model(layout, build_parameters((k_uptake = 0.5, yield = 0.4, k_decay = 0.1)))
    prediction = [1.0, 2.0, 3.0]
    target = [1.0, 1.0, 1.0]

    @test loss(prediction, target, model) == 5.0
    @test training_step!(model, (prediction, target)) == 5.0
end

@testset "simulate" begin
    layout = build_state_layout([:substrate, :biomass, :mineral], [1, 1, 1])
    model = build_model(layout, build_parameters((k_uptake = 0.0, yield = 0.4, k_decay = 0.0)))
    u0 = [1.0, 0.2, 0.0]

    if Base.find_package("OrdinaryDiffEq") !== nothing && Base.find_package("SciMLSensitivity") !== nothing
        solution = simulate(model, u0, (0.0, 1.0); saveat = 1.0)
        @test length(solution.t) == 2
        @test solution.u[end] == u0
    else
        @test_throws ArgumentError simulate(model, u0, (0.0, 1.0))
    end
end

@testset "backend container" begin
    backend = SciMLBackend(nothing, nothing, nothing, nothing)
    @test backend.problem_constructor === nothing
end

@testset "gradient smoke scaffold" begin
    layout = build_state_layout([:substrate, :biomass, :mineral], [1, 1, 1])
    model = build_model(layout, build_parameters((k_uptake = 0.0, yield = 0.4, k_decay = 0.0)))
    u0 = [1.0, 0.2, 0.0]
    target = [1.0, 0.2, 0.0]

    smoke_gradient(loss_function, model, u0, tspan, target; kwargs...) = (loss_function, model, u0, tspan, target, kwargs)

    result = gradient_smoke_test(model, u0, (0.0, 1.0), target; gradient_function = smoke_gradient)

    @test result[1] === trajectory_loss
    @test result[2] === model
    @test result[3] == u0
end

@testset "real sciml example" begin
    @test isfile(joinpath(@__DIR__, "..", "examples", "sciml_tiny_example.jl"))
end

@testset "real differentiability smoke" begin
    required_packages = ("OrdinaryDiffEq", "SciMLSensitivity", "Zygote")
    if all(pkg -> Base.find_package(pkg) !== nothing, required_packages)
        include(joinpath(@__DIR__, "..", "examples", "sciml_gradient_smoke.jl"))
        value, gradient = run_gradient_smoke_example()

        @test isfinite(value)
        @test length(gradient) == 1
        @test isfinite(gradient[1])
    else
        @test true
    end
end

@testset "tiny training example" begin
    @test isfile(joinpath(@__DIR__, "..", "examples", "train_tiny_example.jl"))
end

@testset "tiny training api" begin
    if all(pkg -> Base.find_package(pkg) !== nothing, ("OrdinaryDiffEq", "SciMLSensitivity", "Zygote"))
        include(joinpath(@__DIR__, "..", "examples", "train_tiny_example.jl"))
        result = run_tiny_training_example(
            initial_parameters = TINY_TRAINING_INITIAL_PARAMETERS,
            backend = TINY_TRAINING_BACKEND,
            alg = TINY_TRAINING_ALG,
            sensealg = TINY_TRAINING_SENSEALG,
        )

        @test isfinite(result.value)
        @test length(result.gradient) == 3
        @test result.updated_parameters isa TinyTrainingParameters
        @test all(isfinite, result.gradient)
        @test all(isfinite, parameter_vector(result.updated_parameters))
        @test isfinite(result.updated_value)
    else
        @test true
    end
end
