struct StateLayout
    names::Vector{Symbol}
    sizes::Vector{Int}
    ranges::Vector{UnitRange{Int}}
    dim::Int
end

struct Parameters{T}
    values::T
end

struct SoilModel
    state_layout::StateLayout
    parameters::Parameters
end

const ModelConfig = SoilModel

function build_state_layout(names::AbstractVector{Symbol}, sizes::AbstractVector{<:Integer})
    length(names) == length(sizes) || throw(ArgumentError("names and sizes must match"))

    ranges = UnitRange{Int}[]
    first_index = 1
    for size in sizes
        last_index = first_index + Int(size) - 1
        push!(ranges, first_index:last_index)
        first_index = last_index + 1
    end

    StateLayout(collect(names), Int.(collect(sizes)), ranges, first_index - 1)
end

build_parameters(values::NamedTuple) = Parameters(values)

build_model(state_layout::StateLayout, parameters::Parameters) = SoilModel(state_layout, parameters)

function state_range(layout::StateLayout, name::Symbol)
    index = findfirst(==(name), layout.names)
    index === nothing && throw(KeyError(name))
    return layout.ranges[index]
end

parameter(parameters::Parameters, name::Symbol) = getproperty(parameters.values, name)

