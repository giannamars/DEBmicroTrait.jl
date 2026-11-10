using DEBmicroTrait
using OrdinaryDiffEq
using SciMLSensitivity
using Zygote

make_hysomi_smoke_encoder(weight_offset, bias_offset) = TraitEncoder(
    (
        reshape(range(weight_offset, step = 0.01, length = HYSOMI_TRAIT_INPUT_DIM), 1, HYSOMI_TRAIT_INPUT_DIM),
        reshape(range(weight_offset + 0.5, step = 0.01, length = length(HYSOMI_BIOPAR_NAMES)), length(HYSOMI_BIOPAR_NAMES), 1),
    ),
    (
        collect(range(bias_offset, step = 0.01, length = 1)),
        collect(range(bias_offset + 0.05, step = 0.01, length = length(HYSOMI_BIOPAR_NAMES))),
    ),
    relu;
    output_activation = sigmoid,
    parameter_names = HYSOMI_BIOPAR_NAMES,
    parameter_ranges = HYSOMI_BIOPAR_RANGES,
)

make_hysomi_smoke_traits(prefix) = vcat(prefix, zeros(HYSOMI_TRAIT_INPUT_DIM - length(prefix)))

const HYSOMI_SMOKE_TEMPLATE = make_hysomi_smoke_encoder(0.05, -0.02)
const HYSOMI_SMOKE_TEACHER = make_hysomi_smoke_encoder(0.15, 0.03)
const HYSOMI_SMOKE_TRAITS = [make_hysomi_smoke_traits([0.5, 0.4, 0.1]), make_hysomi_smoke_traits([0.2, 0.3, 0.7])]
const HYSOMI_SMOKE_SOILPAR = (
    k_a = 0.1,
    k_d = 0.01,
    C_S_max = 1.0,
    I_C_L = 0.0,
    I_C_H = 0.0,
)
const HYSOMI_SMOKE_U0 = [0.1, 0.2, 1.0, 10.0, 0.5, 0.0]
const HYSOMI_SMOKE_TSPAN = (0.0, 1.0)
const HYSOMI_SMOKE_ALG = Tsit5()
const HYSOMI_SMOKE_SENSEALG = InterpolatingAdjoint()

function run_hysomi_training_smoke_example()
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
    step = hysomi_training_step(
        state,
        dataset,
        HYSOMI_SMOKE_TEMPLATE;
        step_size = 1e-10,
        backend = Zygote.withgradient,
        alg = HYSOMI_SMOKE_ALG,
        sensealg = HYSOMI_SMOKE_SENSEALG,
        saveat = 1.0,
    )
    return (; dataset, step)
end

if abspath(PROGRAM_FILE) == @__FILE__
    result = run_hysomi_training_smoke_example()
    println(result.step.value)
    println(result.step.updated_value)
end
