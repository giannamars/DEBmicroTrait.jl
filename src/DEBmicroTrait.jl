module DEBmicroTrait

include("core.jl")
include("traits.jl")
include("model.jl")
include("loss.jl")
include("train.jl")

export StateLayout, Parameters, SoilModel, ModelConfig
export build_state_layout, build_parameters, build_model
export state_range, parameter, rhs!, observables, mass_balance
export TraitEncoder, unitize, devectorize, encode_traits, softplus
export loss, training_step!

end
