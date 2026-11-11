using DEBmicroTrait
using OrdinaryDiffEq

const HYSOMI_SYNTHETIC_OUTPUT_PATH = joinpath(@__DIR__, "data", "synthetic_data")

function hysomi_synthetic_trait_map(parameters)
    normalized = unitize(collect(parameters), HYSOMI_BIOPAR_RANGES[:, 1], HYSOMI_BIOPAR_RANGES[:, 2])
    traits = Vector{Float64}(undef, HYSOMI_TRAIT_INPUT_DIM)
    for index in eachindex(traits)
        traits[index] = normalized[mod1(index, length(normalized))]
    end
    return traits
end

function make_hysomi_synthetic_teacher()
    return TraitEncoder(
        (
            reshape(range(0.05, step = 0.01, length = HYSOMI_TRAIT_INPUT_DIM), 1, HYSOMI_TRAIT_INPUT_DIM),
            reshape(range(0.25, step = 0.01, length = length(HYSOMI_BIOPAR_NAMES)), length(HYSOMI_BIOPAR_NAMES), 1),
        ),
        (
            collect(range(-0.03, step = 0.005, length = 1)),
            collect(range(0.02, step = 0.01, length = length(HYSOMI_BIOPAR_NAMES))),
        ),
        relu;
        output_activation = sigmoid,
        parameter_names = HYSOMI_BIOPAR_NAMES,
        parameter_ranges = HYSOMI_BIOPAR_RANGES,
    )
end

const HYSOMI_SYNTHETIC_SOILPAR = (
    k_a = 0.1,
    k_d = 0.01,
    C_S_max = 0.5,
    I_C_L = 0.0,
    I_C_H = 0.0,
)

const HYSOMI_SYNTHETIC_U0 = [0.0000125, 0.2499875, 1.0, 10.0, 0.0, 0.0]
const HYSOMI_SYNTHETIC_TSPAN = (0.0, 7.0)
const HYSOMI_SYNTHETIC_ALG = Tsit5()

function run_hysomi_synthetic_data_example(n::Integer = 3; output_path::AbstractString = HYSOMI_SYNTHETIC_OUTPUT_PATH, save_artifacts::Bool = false)
    teacher = make_hysomi_synthetic_teacher()
    result = hysomi_generate_python_style_synthetic_dataset(
        hysomi_synthetic_trait_map,
        teacher,
        HYSOMI_BIOPAR_RANGES,
        n;
        soilpar = HYSOMI_SYNTHETIC_SOILPAR,
        u0 = HYSOMI_SYNTHETIC_U0,
        tspan = HYSOMI_SYNTHETIC_TSPAN,
        alg = HYSOMI_SYNTHETIC_ALG,
        saveat = 1.0,
    )

    if save_artifacts
        DEBmicroTrait.hysomi_save_synthetic_dataset_bundle(
            output_path,
            result;
            config = (
                n = n,
                trait_input_dim = HYSOMI_TRAIT_INPUT_DIM,
                soilpar = HYSOMI_SYNTHETIC_SOILPAR,
                u0 = HYSOMI_SYNTHETIC_U0,
                tspan = HYSOMI_SYNTHETIC_TSPAN,
            ),
        )
    end

    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    result = run_hysomi_synthetic_data_example(save_artifacts = true)
    println(size(result.raw_parameters))
    println(length(result.raw_dataset))
    println(count(identity, result.mask))
    println(length(result.filtered_dataset))
    println(result.raw_dataset[1].traits[1:5])
    println(result.raw_dataset[1].target[:, end])
    println(HYSOMI_SYNTHETIC_OUTPUT_PATH)
end
