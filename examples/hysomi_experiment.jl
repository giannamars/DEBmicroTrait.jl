using DEBmicroTrait
using OrdinaryDiffEq
using SciMLSensitivity
using Zygote

include(joinpath(@__DIR__, "hysomi_training_smoke.jl"))

function run_hysomi_experiment_example()
    return DEBmicroTrait.hysomi_run_experiment(
        HYSOMI_SMOKE_TEACHER,
        HYSOMI_SMOKE_TEMPLATE,
        HYSOMI_SMOKE_TRAITS,
        reverse(HYSOMI_SMOKE_TRAITS);
        NN_hyperpar = (model = "hysomi_trait_encoder", input_dim = HYSOMI_TRAIT_INPUT_DIM),
        loss_hyperpar = (data_sel = collect(eachindex(HYSOMI_STATE_NAMES)), k = 1e3),
        training_hyperpar = (epochs = 2, batch_size = 1, lr = 1e-10, validation_fraction = 0.0, shuffle = false),
        training_outputs_library_path = mktempdir(),
        soilpar = HYSOMI_SMOKE_SOILPAR,
        u0 = HYSOMI_SMOKE_U0,
        tspan = HYSOMI_SMOKE_TSPAN,
        save_artifacts = false,
        backend = Zygote.withgradient,
        alg = HYSOMI_SMOKE_ALG,
        sensealg = HYSOMI_SMOKE_SENSEALG,
        saveat = 1.0,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    result = run_hysomi_experiment_example()
    println(result.evaluation_loss)
end
