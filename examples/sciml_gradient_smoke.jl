using DEBmicroTrait
using OrdinaryDiffEq
using SciMLSensitivity
using Zygote

function run_gradient_smoke_example()
    layout = build_state_layout([:substrate, :biomass, :mineral], [1, 1, 1])
    u0 = [1.0, 0.2, 0.0]
    tspan = (0.0, 10.0)
    target = [0.8, 0.25, 0.15]

    loss_from_yield(yield) = begin
        model = build_model(layout, build_parameters((k_uptake = 0.5, yield = yield, k_decay = 0.1)))
        solution = simulate(model, u0, tspan; saveat = 1.0, alg = Tsit5(), sensealg = InterpolatingAdjoint())
        return sum(abs2, solution.u[end] .- target)
    end

    value, back = Zygote.withgradient(loss_from_yield, 0.4)
    return value, back[1]
end

if abspath(PROGRAM_FILE) == @__FILE__
    value, gradient = run_gradient_smoke_example()
    println(value)
    println(gradient)
end
