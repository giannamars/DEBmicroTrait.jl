using DEBmicroTrait
using OrdinaryDiffEq
using SciMLSensitivity
using Zygote

const TINY_TRAINING_BACKEND = TinyGradientBackend(Zygote.withgradient)
const TINY_TRAINING_ALG = Tsit5()
const TINY_TRAINING_SENSEALG = InterpolatingAdjoint()

result = run_tiny_training_example(
    backend = TINY_TRAINING_BACKEND,
    alg = TINY_TRAINING_ALG,
    sensealg = TINY_TRAINING_SENSEALG,
)
println(result)
