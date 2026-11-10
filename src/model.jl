function rhs!(du, u, p::ModelConfig, t)
    @inbounds for i in eachindex(u)
        du[i] = -u[i]
    end
    return du
end

observables(u, p::ModelConfig) = copy(u)
