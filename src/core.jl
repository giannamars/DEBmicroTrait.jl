struct StateLayout
    names::Vector{Symbol}
    sizes::Vector{Int}
    dim::Int
end

struct Parameters
    values::NamedTuple
end

struct ModelConfig
    state_layout::StateLayout
    parameters::Parameters
end

build_state_layout(names::Vector{Symbol}, sizes::Vector{Int}) = StateLayout(names, sizes, sum(sizes))
build_parameters(values::NamedTuple) = Parameters(values)
build_model(state_layout::StateLayout, parameters::Parameters) = ModelConfig(state_layout, parameters)
