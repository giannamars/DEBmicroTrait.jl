# DEBmicroTrait

Minimal Julia scaffold for a soil-pool process model with an optional SciML backend.

## Tiny SciML example

Once `OrdinaryDiffEq.jl`, `SciMLBase.jl`, and `SciMLSensitivity.jl` are available in your Julia environment, you can call `simulate` like this:

```julia
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
```

If the SciML packages are not installed, the core package still loads, and the backend-specific methods stay inactive.

See `examples/sciml_tiny_example.jl` for the same setup as a runnable script.

For a real differentiability smoke test, see `examples/sciml_gradient_smoke.jl`.

For a tiny training prototype, see `examples/train_tiny_example.jl`.
