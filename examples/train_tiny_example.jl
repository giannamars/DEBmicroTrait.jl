using DEBmicroTrait
using OrdinaryDiffEq
using SciMLSensitivity
using Zygote

const TINY_TRAINING_BACKEND = TinyGradientBackend(Zygote.withgradient)
const TINY_TRAINING_ALG = Tsit5()
const TINY_TRAINING_SENSEALG = InterpolatingAdjoint()
const TINY_TRAINING_INITIAL_PARAMETERS = TinyTrainingParameters(0.5, 0.4, 0.1)

result = run_tiny_training_example(
    initial_parameters = TINY_TRAINING_INITIAL_PARAMETERS,
    backend = TINY_TRAINING_BACKEND,
    alg = TINY_TRAINING_ALG,
    sensealg = TINY_TRAINING_SENSEALG,
)
println(result)
