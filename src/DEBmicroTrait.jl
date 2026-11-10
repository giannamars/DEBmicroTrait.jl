module DEBmicroTrait

include("core.jl")
include("traits.jl")
include("hysomi.jl")
include("hysomi_training.jl")
include("model.jl")
include("simulation.jl")
include("loss.jl")
include("train.jl")

export StateLayout, Parameters, SoilModel, ModelConfig
export build_state_layout, build_parameters, build_model
export HysomiParameters, HysomiModel, HYSOMI_BIOPAR_NAMES, HYSOMI_SOILPAR_NAMES, HYSOMI_STATE_NAMES
export HYSOMI_BIOPAR_RANGES, HYSOMI_SOILPAR_RANGES, build_hysomi_parameters, build_hysomi_model, hysomi_state_layout
export hysomi_biopar_vector, hysomi_soilpar_vector, hysomi_parameter_vector, hysomi_rhs!, hysomi_observables, hysomi_mass_balance
export HYSOMI_TRAIT_INPUT_DIM, hysomi_trait_encoder, hysomi_parameter_constraint_loss, softargmax1d, peak_width, constraint_penalty, hysomi_process_constraint_loss, hysomi_trajectory_matrix, hysomi_data_loss, hysomi_loss_components, hysomi_batch_loss
export state_range, parameter, rhs!, observables, mass_balance
export TraitEncoder, unitize, devectorize, relu, sigmoid, trait_vector, encode_traits, trait_parameters, softplus
export initial_state, SciMLBackend, simulate, trajectory_loss
export gradient_smoke_test
export loss, training_step!, TinyTrainingProblem, TinyTrainingParameters, TinyGradientBackend, parameter_vector, tiny_loss, tiny_training_step, run_tiny_training_example

end
