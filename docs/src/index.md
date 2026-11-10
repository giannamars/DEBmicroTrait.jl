```@meta
CurrentModule = DEBmicroTrait
```

# DEBmicroTrait

Documentation for [DEBmicroTrait](https://github.com/giannamars/DEBmicroTrait.jl).

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

The same example is available as `examples/sciml_tiny_example.jl`.

For a real gradient smoke test, see `examples/sciml_gradient_smoke.jl`.

```@index
```

```@autodocs
Modules = [DEBmicroTrait]
```
