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
export TraitEncoder, unitize, devectorize, encode_traits, softplus
export initial_state, SciMLBackend, build_problem, simulate, trajectory_loss
export loss, training_step!

end
