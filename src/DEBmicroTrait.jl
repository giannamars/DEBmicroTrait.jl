module DEBmicroTrait

include("core.jl")
include("traits.jl")
include("model.jl")
include("simulation.jl")
include("loss.jl")
include("train.jl")

export StateLayout, Parameters, SoilModel, ModelConfig
export build_state_layout, build_parameters, build_model
export state_range, parameter, rhs!, observables, mass_balance
export TraitEncoder, unitize, devectorize, relu, sigmoid, trait_vector, encode_traits, trait_parameters, softplus
export initial_state, SciMLBackend, simulate, trajectory_loss
export gradient_smoke_test
export loss, training_step!, TinyTrainingProblem, TinyTrainingParameters, TinyGradientBackend, parameter_vector, tiny_loss, tiny_training_step, run_tiny_training_example

end
