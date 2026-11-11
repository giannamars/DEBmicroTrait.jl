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

## Hysomi workflow

The Julia Hysomi scaffold mirrors the Python workflow in three small steps:

1. Generate synthetic data with `examples/hysomi_synthetic_data_example.jl`.
2. Run training with `examples/hysomi_start_training.jl`.
3. Inspect a saved run with `examples/hysomi_post_training_analysis.jl`.

The end-to-end helper is `hysomi_run_experiment`, which ties generation, training, and evaluation together.

```julia
using DEBmicroTrait
using OrdinaryDiffEq
using SciMLSensitivity
using Zygote

result = hysomi_run_experiment(
    HYSOMI_SMOKE_TEACHER,
    HYSOMI_SMOKE_TEMPLATE,
    HYSOMI_SMOKE_TRAITS,
    reverse(HYSOMI_SMOKE_TRAITS);
    NN_hyperpar = (model = "hysomi_trait_encoder", input_dim = HYSOMI_TRAIT_INPUT_DIM),
    loss_hyperpar = (data_sel = collect(eachindex(HYSOMI_STATE_NAMES)), k = 1e3, awl = "AutomaticWeightedLossL2"),
    training_hyperpar = (epochs = 2, batch_size = 1, lr = 1e-10, validation_fraction = 0.0, shuffle = false),
    training_outputs_library_path = "./outputs",
    soilpar = HYSOMI_SMOKE_SOILPAR,
    u0 = HYSOMI_SMOKE_U0,
    tspan = HYSOMI_SMOKE_TSPAN,
    backend = Zygote.withgradient,
    alg = HYSOMI_SMOKE_ALG,
    sensealg = HYSOMI_SMOKE_SENSEALG,
    saveat = 1.0,
)
```

The corresponding post-training analysis helper is `hysomi_post_training_analysis(output_library_path, hash_id)`.

## Project layout

- `src/hysomi.jl` — Hysomi state, parameters, and ODE right-hand side.
- `src/hysomi_training.jl` — training helpers, config-driven entrypoints, and analysis utilities.
- `src/synthetic_data.jl` — synthetic dataset generation and process filtering.
- `examples/hysomi_synthetic_data_example.jl` — standalone synthetic-data demo.
- `examples/hysomi_start_training.jl` — start-training style launcher with on-disk data folders.
- `examples/hysomi_post_training_analysis.jl` — post-training summary helper.
- `examples/hysomi_experiment.jl` — compact end-to-end wrapper example.
