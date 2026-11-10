using DEBmicroTrait
using OrdinaryDiffEq
using SciMLSensitivity
using Zygote

const TINY_TRAINING_BACKEND = TinyGradientBackend(Zygote.withgradient)
const TINY_TRAINING_ALG = Tsit5()
const TINY_TRAINING_SENSEALG = InterpolatingAdjoint()
const TINY_TRAINING_ENCODER = TraitEncoder(
    (
        [0.8 0.1 -0.2; -0.3 0.9 0.4; 0.2 0.5 0.7; 0.6 -0.4 0.3],
        [0.7 -0.1 0.2 0.5; -0.2 0.6 0.1 -0.3; 0.3 0.2 -0.5 0.4],
    ),
    (
        [0.1, 0.0, -0.1, 0.05],
        [-0.2, 0.05, -0.3],
    ),
    relu;
    output_activation = sigmoid,
    parameter_names = (:k_uptake, :yield, :k_decay),
    parameter_ranges = [0.05 0.8; 0.1 0.9; 0.01 0.35],
)
const TINY_TRAINING_TRAITS = [0.5, 0.4, 0.1]

result = run_tiny_training_example(
    initial_traits = TINY_TRAINING_TRAITS,
    trait_encoder = TINY_TRAINING_ENCODER,
    backend = TINY_TRAINING_BACKEND,
    alg = TINY_TRAINING_ALG,
    sensealg = TINY_TRAINING_SENSEALG,
)
println(result)
