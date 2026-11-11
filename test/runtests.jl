using Test
using DEBmicroTrait
using LinearAlgebra

make_hysomi_smoke_encoder(weight_offset, bias_offset) = TraitEncoder(
    (
        reshape(collect(range(weight_offset, step = 0.01, length = HYSOMI_TRAIT_INPUT_DIM)), 1, HYSOMI_TRAIT_INPUT_DIM),
        reshape(collect(range(weight_offset + 0.5, step = 0.01, length = length(HYSOMI_BIOPAR_NAMES))), length(HYSOMI_BIOPAR_NAMES), 1),
    ),
    (
        collect(range(bias_offset, step = 0.01, length = 1)),
        collect(range(bias_offset + 0.05, step = 0.01, length = length(HYSOMI_BIOPAR_NAMES))),
    ),
    relu;
    output_activation = sigmoid,
    parameter_names = HYSOMI_BIOPAR_NAMES,
    parameter_ranges = HYSOMI_BIOPAR_RANGES,
)

make_hysomi_smoke_traits(prefix) = vcat(prefix, zeros(HYSOMI_TRAIT_INPUT_DIM - length(prefix)))

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

@testset "hybrid loss helpers" begin
    @test isapprox(softargmax1d([0.0, 1.0, 0.0]; beta = 1e6), 1.0; atol = 1e-6)
    @test size(softargmax2d([0.0 0.0; 0.0 1.0]; beta = 1e6)) == (1, 2)
    @test peak_width_v2([0.1, 0.6, 0.9], 0.5; k = 1e6) ≈ 2.0

    awl = AutomaticWeightedLossL2(3)
    @test isfinite(awl([1.0, 2.0, 3.0]))
    @test no_constraint_module(2)([4.0, 5.0]) == 4.0

    pred = zeros(6, 60, 2)
    data = zeros(6, 60, 2)
    for batch in 1:2
        center = batch == 1 ? 20 : 35
        for time in 1:60
            pred[1, time, batch] = exp(-((time - center) / 6)^2)
            pred[2, time, batch] = 1.5
            pred[3, time, batch] = 0.3
            pred[4, time, batch] = 1.0
            pred[5, time, batch] = 0.1
            pred[6, time, batch] = 0.2
        end
    end
    data .= pred

    par = hcat([fill(value, 2) for value in (1.0, 0.5, 0.4, 0.2, 3.0, 2.0, 0.7, 0.2, 0.4, 0.5, 0.3, 0.2, 0.3)]...)
    true_par = copy(par)
    par_trans_list = (ones(13), zeros(13))
    loss_module = loss_constraint_autoW3(1:6, 1e3)
    result = loss_module(pred, data, par, true_par, par_trans_list)

    @test length(result.loss_tensor) == 15
    @test length(result.constraint_tensor) == 14
    @test length(result.ana_loss) == 10
    @test length(result.theta_L1_tensor) == 13
    @test all(isfinite, result.loss_tensor)
    @test isfinite(result.loss_pred)
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

@testset "hysomi training helpers" begin
    template = make_hysomi_smoke_encoder(0.05, -0.02)
    state = hysomi_encoder_state(template)
    rebuilt = hysomi_rebuild_encoder(template, state)
    traits = make_hysomi_smoke_traits([0.5, 0.4, 0.1])

    @test trait_parameters(rebuilt, traits) == trait_parameters(template, traits)

    if all(pkg -> Base.find_package(pkg) !== nothing, ("OrdinaryDiffEq", "SciMLSensitivity", "Zygote"))
        teacher = make_hysomi_smoke_encoder(0.15, 0.03)
        soilpar = (
            k_a = 0.1,
            k_d = 0.01,
            C_S_max = 1.0,
            I_C_L = 0.0,
            I_C_H = 0.0,
        )
        dataset = hysomi_generate_dataset(
            teacher,
            [traits, make_hysomi_smoke_traits([0.2, 0.3, 0.7])];
            soilpar = soilpar,
            u0 = [0.1, 0.2, 1.0, 10.0, 0.5, 0.0],
            tspan = (0.0, 1.0),
            alg = Tsit5(),
            saveat = 1.0,
        )

        @test length(dataset) == 2
        @test size(dataset[1].target, 1) == length(HYSOMI_STATE_NAMES)
        @test size(dataset[1].target, 2) == 2

        step = hysomi_training_step(
            state,
            dataset,
            template;
            step_size = 1e-10,
            backend = Zygote.withgradient,
            alg = Tsit5(),
            sensealg = InterpolatingAdjoint(),
            saveat = 1.0,
        )

        epoch = hysomi_training_epoch(
            state,
            dataset,
            template;
            batch_size = 1,
            step_size = 1e-10,
            backend = Zygote.withgradient,
            alg = Tsit5(),
            sensealg = InterpolatingAdjoint(),
            saveat = 1.0,
        )

        @test isfinite(step.value)
        @test isfinite(step.updated_value)
        @test length(step.gradient) == length(state)
        @test all(isfinite, step.gradient)
        @test all(isfinite, step.updated_state)
        @test length(epoch.losses) == 2
        @test all(isfinite, epoch.losses)
    else
        @test true
    end
end

@testset "hysomi synthetic data" begin
    constraint_sample = [0.5, 0.2, 0.3, 0.4, 0.8, 0.7, 0.2, 0.1, 0.6, 0.9, 0.1, 0.0, 0.0]
    invalid_sample = [0.5, 0.2, 0.3, 0.4, 0.8, 0.9, 0.2, 0.1, 0.6, 0.9, 0.1, 0.0, 0.0]

    @test hysomi_parameter_constraint(constraint_sample)
    @test !hysomi_parameter_constraint(invalid_sample)

    valid_target = permutedims(hcat(
        vcat(fill(0.1, 20), fill(0.4, 10), fill(0.1, 30)),
        fill(0.3, 60),
        fill(0.8, 60),
        fill(0.6, 60),
        fill(0.2, 60),
        fill(0.0, 60),
    ))
    @test hysomi_process_filter(valid_target)

    invalid_target = copy(valid_target)
    invalid_target[3, 1] = 0.0
    @test !hysomi_process_filter(invalid_target)

    sampled = hysomi_constrained_sampling(HYSOMI_BIOPAR_RANGES, 3)
    @test size(sampled) == (3, length(HYSOMI_BIOPAR_NAMES))
    @test all(row -> hysomi_parameter_constraint(collect(row)), eachrow(sampled))

    if all(pkg -> Base.find_package(pkg) !== nothing, ("OrdinaryDiffEq", "SciMLSensitivity", "Zygote"))
        teacher = make_hysomi_smoke_encoder(0.15, 0.03)
        trait_map = _ -> make_hysomi_smoke_traits([0.5, 0.4, 0.1])
        synthetic = hysomi_generate_synthetic_dataset(
            trait_map,
            teacher,
            HYSOMI_BIOPAR_RANGES,
            2;
            soilpar = (
                k_a = 0.1,
                k_d = 0.01,
                C_S_max = 1.0,
                I_C_L = 0.0,
                I_C_H = 0.0,
            ),
            u0 = [0.1, 0.2, 1.0, 10.0, 0.5, 0.0],
            tspan = (0.0, 1.0),
            alg = Tsit5(),
            saveat = 1.0,
        )

        @test length(synthetic) == 2
        @test all(sample -> hasproperty(sample, :parameters), synthetic)
        @test all(sample -> length(sample.parameters) == length(HYSOMI_BIOPAR_NAMES), synthetic)
        @test all(sample -> size(sample.target, 1) == length(HYSOMI_STATE_NAMES), synthetic)
        @test all(sample -> all(isfinite, sample.target), synthetic)
    else
        @test true
    end
end

@testset "hysomi synthetic data example" begin
    if Base.find_package("OrdinaryDiffEq") !== nothing
        include(joinpath(@__DIR__, "..", "examples", "hysomi_synthetic_data_example.jl"))
        result = run_hysomi_synthetic_data_example(2)
        @test size(result.raw_parameters, 2) == length(HYSOMI_BIOPAR_NAMES)
        @test length(result.raw_dataset) == size(result.raw_parameters, 1)
        @test length(result.mask) == length(result.raw_dataset)
        @test length(result.filtered_dataset) == count(identity, result.mask)
        @test all(sample -> length(sample.traits) == HYSOMI_TRAIT_INPUT_DIM, result.raw_dataset)
        @test all(sample -> size(sample.target, 1) == length(HYSOMI_STATE_NAMES), result.raw_dataset)
    else
        @test true
    end
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
