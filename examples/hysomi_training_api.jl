using DEBmicroTrait
using OrdinaryDiffEq
using SciMLSensitivity
using Zygote

include(joinpath(@__DIR__, "hysomi_training_smoke.jl"))

function run_hysomi_training_api_example()
    dataset = hysomi_generate_dataset(
        HYSOMI_SMOKE_TEACHER,
        HYSOMI_SMOKE_TRAITS;
        soilpar = HYSOMI_SMOKE_SOILPAR,
        u0 = HYSOMI_SMOKE_U0,
        tspan = HYSOMI_SMOKE_TSPAN,
        alg = HYSOMI_SMOKE_ALG,
        adaptive = false,
        dt = 0.01,
        saveat = 1.0,
    )
    state = hysomi_encoder_state(HYSOMI_SMOKE_TEMPLATE)
    return hysomi_train(
        state,
        dataset,
        HYSOMI_SMOKE_TEMPLATE;
        epochs = 2,
        batch_size = 1,
        validation_fraction = 0.5,
        step_size = 1e-10,
        backend = Zygote.withgradient,
        alg = HYSOMI_SMOKE_ALG,
        sensealg = HYSOMI_SMOKE_SENSEALG,
        saveat = 1.0,
        shuffle = false,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    result = run_hysomi_training_api_example()
    println(result.train_losses)
    println(result.validation_losses)
end
