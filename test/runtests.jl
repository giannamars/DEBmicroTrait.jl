using Test
using DEBmicroTrait

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

    @test_throws ArgumentError simulate(model, u0, (0.0, 1.0))
    @test_throws ArgumentError build_problem(model, u0, (0.0, 1.0))
end

@testset "backend container" begin
    backend = SciMLBackend(nothing, nothing, nothing, nothing)
    @test backend.problem_constructor === nothing
end
