using Test
using DEBmicroTrait
using LinearAlgebra

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

@testset "trait mlp encoding" begin
    encoder = TraitEncoder(
        (
            [0.8 0.1 -0.2; -0.3 0.9 0.4; 0.2 0.5 0.7; 0.6 -0.4 0.3],
            [0.7 -0.1 0.2 0.5; -0.2 0.6 0.1 -0.3; 0.3 0.2 -0.5 0.4],
        ),
        (
            [0.1, 0.0, -0.1, 0.05],
            [-0.2, 0.05, -0.3],
        ),
        relu;
        output_activation = sigmoid,
        parameter_names = (:k_uptake, :yield, :k_decay),
        parameter_ranges = [0.05 0.8; 0.1 0.9; 0.01 0.35],
    )
    traits = [0.5, 0.4, 0.1]
    encoded = trait_parameters(encoder, traits)

    @test keys(encoded) == (:k_uptake, :yield, :k_decay)
    @test all(>=(0.01), values(encoded))
    @test all(<=(0.9), values(encoded))
    @test encoded == trait_parameters(encoder, traits)
end

@testset "trait regularization hooks" begin
    encoder = TraitEncoder(
        (
            [0.8 0.1 -0.2; -0.3 0.9 0.4; 0.2 0.5 0.7; 0.6 -0.4 0.3],
            [0.7 -0.1 0.2 0.5; -0.2 0.6 0.1 -0.3; 0.3 0.2 -0.5 0.4],
        ),
        (
            [0.1, 0.0, -0.1, 0.05],
            [-0.2, 0.05, -0.3],
        ),
        relu;
        output_activation = sigmoid,
        parameter_names = (:k_uptake, :yield, :k_decay),
        parameter_ranges = [0.05 0.8; 0.1 0.9; 0.01 0.35],
        dropout_probability = 0.25,
        batchnorm_statistics = (
            (mean = zeros(4), variance = ones(4)),
        ),
    )
    traits = [0.5, 0.4, 0.1]

    encoded = trait_vector(encoder, traits)
    @test length(encoded) == 3
    @test all(isfinite, encoded)
    @test trait_vector(encoder, traits; training = false) == encoded

    trained = trait_vector(encoder, traits; training = true)
    @test length(trained) == 3
    @test all(isfinite, trained)
end

@testset "trait parameters" begin
    encoder = TraitEncoder(Matrix{Float64}(I, 3, 3), zeros(3), identity; parameter_names = (:k_uptake, :yield, :k_decay))
    traits = [0.5, 0.4, 0.1]
    named_parameters = trait_parameters(encoder, traits)

    @test named_parameters.k_uptake == 0.5
    @test named_parameters.yield == 0.4
    @test named_parameters.k_decay == 0.1
    @test TinyTrainingParameters(named_parameters) == TinyTrainingParameters(0.5, 0.4, 0.1)
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
            initial_traits = TINY_TRAINING_TRAITS,
            trait_encoder = TINY_TRAINING_ENCODER,
            backend = TINY_TRAINING_BACKEND,
            alg = TINY_TRAINING_ALG,
            sensealg = TINY_TRAINING_SENSEALG,
        )
        expected_parameters = TinyTrainingParameters(trait_parameters(TINY_TRAINING_ENCODER, TINY_TRAINING_TRAITS))

        @test isfinite(result.value)
        @test length(result.gradient) == 3
        @test result.updated_parameters isa TinyTrainingParameters
        @test result.initial_parameters == expected_parameters
        @test all(isfinite, result.gradient)
        @test all(isfinite, parameter_vector(result.updated_parameters))
        @test isfinite(result.updated_value)
    else
        @test true
    end
end

@testset "hysomi core" begin
    layout = hysomi_state_layout()
    @test layout.dim == 6
    @test state_range(layout, :CO2) == 6:6

    biopar = (
        mu = 0.1,
        K_L = 10.0,
        Y_mu = 0.4,
        S_tr = 0.2,
        r_react = 0.5,
        r_deact = 0.3,
        K_C_H = 0.7,
        vmax = 0.2,
        fm = 0.4,
        Y_m = 0.5,
        m_max = 0.6,
        beta = 0.2,
        tau = 0.3,
    )
    soilpar = (
        k_a = 0.1,
        k_d = 0.01,
        C_S_max = 1.0,
        I_C_L = 0.0,
        I_C_H = 0.0,
    )
    model = build_hysomi_model(layout, build_hysomi_parameters(biopar, soilpar))
    u = [0.1, 0.2, 1.0, 10.0, 0.5, 0.0]
    du = zeros(6)

    hysomi_rhs!(du, u, parameter_vector(model), 0.0)

    @test length(du) == 6
    @test all(isfinite, du)
    @test isapprox(sum(du), 0.0; atol = 1e-12)
    @test isapprox(hysomi_mass_balance(u, model), sum(u))
end

@testset "hysomi trait scaffold" begin
    encoder = TraitEncoder(
        (
            zeros(20, HYSOMI_TRAIT_INPUT_DIM),
            zeros(length(HYSOMI_BIOPAR_NAMES), 20),
        ),
        (
            zeros(20),
            zeros(length(HYSOMI_BIOPAR_NAMES)),
        ),
        relu;
        output_activation = sigmoid,
        parameter_names = HYSOMI_BIOPAR_NAMES,
        parameter_ranges = HYSOMI_BIOPAR_RANGES,
    )
    traits = zeros(HYSOMI_TRAIT_INPUT_DIM)
    biopar = trait_parameters(encoder, traits)

    @test keys(biopar) == HYSOMI_BIOPAR_NAMES
    @test length(values(biopar)) == length(HYSOMI_BIOPAR_NAMES)
    @test all(isfinite, values(biopar))
end

@testset "hysomi loss scaffold" begin
    required_packages = ("OrdinaryDiffEq", "SciMLSensitivity", "Zygote")
    if all(pkg -> Base.find_package(pkg) !== nothing, required_packages)
        encoder = TraitEncoder(
            (
                zeros(20, HYSOMI_TRAIT_INPUT_DIM),
                zeros(length(HYSOMI_BIOPAR_NAMES), 20),
            ),
            (
                zeros(20),
                zeros(length(HYSOMI_BIOPAR_NAMES)),
            ),
            relu;
            output_activation = sigmoid,
            parameter_names = HYSOMI_BIOPAR_NAMES,
            parameter_ranges = HYSOMI_BIOPAR_RANGES,
        )
        traits = zeros(HYSOMI_TRAIT_INPUT_DIM)
        soilpar = (k_a = 0.1, k_d = 0.01, C_S_max = 1.0, I_C_L = 0.0, I_C_H = 0.0)
        u0 = [0.1, 0.2, 1.0, 10.0, 0.5, 0.0]
        model = build_hysomi_model(hysomi_state_layout(), build_hysomi_parameters(trait_parameters(encoder, traits), soilpar))
        solution = simulate(model, u0, (0.0, 1.0); saveat = 1.0)
        target = hysomi_trajectory_matrix(solution)

        components = hysomi_loss_components(encoder, traits, soilpar, u0, (0.0, 1.0), target; saveat = 1.0)

        @test isfinite(components.total_loss)
        @test isapprox(components.data_loss, 0.0; atol = 1e-8)
        @test components.model isa HysomiModel
    else
        @test true
    end
end
