using DEBmicroTrait
using OrdinaryDiffEq
using SciMLSensitivity

layout = build_state_layout([:substrate, :biomass, :mineral], [1, 1, 1])
model = build_model(layout, build_parameters((k_uptake = 0.5, yield = 0.4, k_decay = 0.1)))
u0 = [1.0, 0.2, 0.0]

solution = simulate(
    model,
    u0,
    (0.0, 10.0);
    saveat = 1.0,
    alg = Tsit5(),
    sensealg = InterpolatingAdjoint(),
)

println(solution.u[end])
