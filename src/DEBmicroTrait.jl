module DEBmicroTrait

include("core.jl")
include("model.jl")
include("traits.jl")
include("loss.jl")
include("train.jl")

export StateLayout, Parameters, ModelConfig
export build_state_layout, build_parameters, build_model
export rhs!, observables, loss
export encode_traits, training_step!

end
