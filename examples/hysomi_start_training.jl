using DEBmicroTrait
using OrdinaryDiffEq
using SciMLSensitivity
using Zygote
using Serialization

include(joinpath(@__DIR__, "hysomi_training_smoke.jl"))

const HYSOMI_START_DATA_PATH = joinpath(@__DIR__, "data", "synthetic_data")
const HYSOMI_START_OUTPUT_PATH = joinpath(@__DIR__, "data", "outputs")

function run_hysomi_start_training_example(; data_path = HYSOMI_START_DATA_PATH, output_path = HYSOMI_START_OUTPUT_PATH, save_artifacts = true)
    mkpath(data_path)
    mkpath(output_path)

    result = DEBmicroTrait.hysomi_run_experiment(
        HYSOMI_SMOKE_TEACHER,
        HYSOMI_SMOKE_TEMPLATE,
        HYSOMI_SMOKE_TRAITS,
        reverse(HYSOMI_SMOKE_TRAITS);
        NN_hyperpar = (model = "hysomi_trait_encoder", input_dim = HYSOMI_TRAIT_INPUT_DIM),
        loss_hyperpar = (data_sel = collect(eachindex(HYSOMI_STATE_NAMES)), k = 1e3),
        training_hyperpar = (epochs = 2, batch_size = 1, lr = 1e-10, validation_fraction = 0.0, shuffle = false),
        training_outputs_library_path = output_path,
        soilpar = HYSOMI_SMOKE_SOILPAR,
        u0 = HYSOMI_SMOKE_U0,
        tspan = HYSOMI_SMOKE_TSPAN,
        save_artifacts = save_artifacts,
        backend = Zygote.withgradient,
        alg = HYSOMI_SMOKE_ALG,
        sensealg = HYSOMI_SMOKE_SENSEALG,
        saveat = 1.0,
    )

    if save_artifacts
        open(joinpath(data_path, "train_dataset.jls"), "w") do io
            Serialization.serialize(io, result.train_dataset)
        end
        open(joinpath(data_path, "test_dataset.jls"), "w") do io
            Serialization.serialize(io, result.test_dataset)
        end
    end

    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    result = run_hysomi_start_training_example()
    println(result.bundle.out_dir)
    println(result.evaluation_loss)
end
