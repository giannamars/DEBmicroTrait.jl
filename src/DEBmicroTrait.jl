module DEBmicroTrait

import SHA
import Serialization
import TOML

include("core.jl")
include("traits.jl")
include("hysomi.jl")
include("hysomi_training.jl")
include("synthetic_data.jl")
include("model.jl")
include("simulation.jl")
include("loss.jl")
include("train.jl")

export StateLayout, Parameters, SoilModel, ModelConfig
export build_state_layout, build_parameters, build_model
export HysomiParameters, HysomiModel, HYSOMI_BIOPAR_NAMES, HYSOMI_SOILPAR_NAMES, HYSOMI_STATE_NAMES
export HYSOMI_BIOPAR_RANGES, HYSOMI_SOILPAR_RANGES, build_hysomi_parameters, build_hysomi_model, hysomi_state_layout
export hysomi_biopar_vector, hysomi_soilpar_vector, hysomi_parameter_vector, hysomi_rhs!, hysomi_observables, hysomi_mass_balance
export HYSOMI_TRAIT_INPUT_DIM, hysomi_trait_encoder, hysomi_parameter_constraint_loss, softargmax1d, softargmax2d, peak_width, peak_width_v2, constraint_penalty, hysomi_process_constraint_loss, hysomi_trajectory_matrix, hysomi_data_loss, hysomi_loss_components, hysomi_batch_loss, hysomi_encoder_state, hysomi_rebuild_encoder, hysomi_generate_sample, hysomi_generate_dataset, hysomi_parameter_constraint, hysomi_latin_hypercube, hysomi_constrained_sampling, hysomi_generate_synthetic_sample, hysomi_process_filter, hysomi_generate_synthetic_dataset, hysomi_generate_python_style_synthetic_dataset, hysomi_save_synthetic_dataset_bundle, hysomi_load_synthetic_dataset_bundle, AutomaticWeightedLossL2, no_constraint_module, loss_constraint_autoW3, hysomi_batch_indices, hysomi_split_dataset, hysomi_train_epoch, hysomi_train, hysomi_training_hash, hysomi_hyperpar_to_toml, hysomi_load_training_config, hysomi_save_training_result, hysomi_load_serialized_dataset, hysomi_train_from_config, hysomi_run_experiment, hysomi_update_state, hysomi_training_step, hysomi_training_epoch
export state_range, parameter, rhs!, observables, mass_balance
export TraitEncoder, unitize, devectorize, relu, sigmoid, trait_vector, encode_traits, trait_parameters, softplus
export initial_state, SciMLBackend, simulate, trajectory_loss
export gradient_smoke_test
export loss, training_step!, TinyTrainingProblem, TinyTrainingParameters, TinyGradientBackend, parameter_vector, tiny_loss, tiny_training_step, run_tiny_training_example

end
