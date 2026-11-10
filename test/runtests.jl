using Test
using DEBmicroTrait

@testset "core" begin
    layout = build_state_layout([:P, :D, :E], [2, 3, 1])
    @test layout.dim == 6
    params = build_parameters((k = 1.0,))
    model = build_model(layout, params)
    @test model.state_layout.dim == 6
    @test model.parameters.values.k == 1.0
end

@testset "model" begin
    layout = build_state_layout([:x], [2])
    model = build_model(layout, build_parameters((k = 1.0,)))
    u = [1.0, 2.0]
    du = similar(u)
    rhs!(du, u, model, 0.0)
    @test du == [-1.0, -2.0]
    @test observables(u, model) == u
end

@testset "traits" begin
    @test encode_traits([1.0, 2.0], [3.0]) == [4.0, 5.0]
end

@testset "loss" begin
    layout = build_state_layout([:x], [1])
    model = build_model(layout, build_parameters((k = 1.0,)))
    @test loss([2.0, 3.0], [1.0, 1.0], model) == 5.0
end

@testset "train" begin
    layout = build_state_layout([:x], [1])
    model = build_model(layout, build_parameters((k = 1.0,)))
    returned_model, batch = training_step!(model, (1, 2, 3))
    @test returned_model === model
    @test batch == (1, 2, 3)
end
