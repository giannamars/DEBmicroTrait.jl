function normalized_substrate_flux!(F_c::Vector{T}, S::Vector{T}, K::Matrix{T}) where {T<:Real}
    (I,J) = size(K)
    @assert I == size(K,1)
    @assert J == size(K,2)
    @inbounds for l in 1:I, j in 1:J
            val = S[l]/K[l,j]
            F_c[j] = F_c[j] + val[1]
    end
    return F_c
end


function conjugate_substrate_flux!(F_r::Vector{T}, E::Vector{T}, K::Matrix{T}, N_SB::Matrix{T}) where {T<:Real}
    (I,J) = size(K)
    @assert J == size(E,1)
    @assert I == size(F_r,1)
    @inbounds for i in 1:I, l in 1:J
        val = N_SB[i,l]*E[l]/K[i,l]
        F_r[i] = val[1]
    end
    return F_r
end


function ECA_kinetics!(ECA::Matrix{T}, S::Vector{T}, E::Vector{T}, K::Matrix{T}, k2p::Vector{T}, N_SB::Matrix{T}) where {T<:Real}
    (I,J) = size(K)
    @assert I == size(S,1)
    @assert J == size(E,1)
    @assert (I,J) == size(K)
    F_c = zeros(J)
    F_r = zeros(I)
    F_r = conjugate_substrate_flux!(F_r, E, K, N_SB)
    F_c = normalized_substrate_flux!(F_c, S, K)
    @inbounds for i in 1:I, j in 1:J
            val = k2p[i]*N_SB[i,j]*E[j]*(S[i]/K[i,j])/(1 + F_r[i] + F_c[j])
            ECA[i,j] = val[1]
    end

    any(x->x==true, isnan.(ECA)) ? throw(DomainError("NaN in DEBmicroTrait.ECA_kinetics!")) : return ECA
end

function ECA_kinetics!(ECA::Matrix{T}, S::Vector{T}, E::Vector{T}, K::Matrix{T}, k2p::Matrix{T}, N_SB::Matrix{T}) where {T<:Real}
    (I,J) = size(K)
    @assert I == size(S,1)
    @assert J == size(E,1)
    @assert (I,J) == size(K)
    F_c = zeros(J)
    F_r = zeros(I)
    F_r = conjugate_substrate_flux!(F_r, E, K, N_SB)
    F_c = normalized_substrate_flux!(F_c, S, K)
    @inbounds for i in 1:I, j in 1:J
            val = k2p[i,j]*N_SB[i,j]*E[j]*(S[i]/K[i,j])/(1 + F_r[i] + F_c[j])
            ECA[i,j] = val[1]
    end

    any(x->x==true, isnan.(ECA)) ? throw(DomainError("NaN in DEBmicroTrait.ECA_kinetics!")) : return ECA
end

function conjugate_substrate_flux!(F_r::Vector{T}, E::Vector{T}, K::Matrix{T}) where {T<:Real}
    (I,J) = size(K)
    @assert J == size(E,1)
    @assert I == size(F_r,1)
    @inbounds for i in 1:I, l in 1:J
        val = E[l]/K[i,l]
        F_r[i] = F_r[i] + val[1]
    end
    return F_r
end

function normalized_substrate_flux!(F_c::Vector{T}, S::Vector{T}, K::Matrix{T}, α_kin::Matrix{T}) where {T<:Real}
    (I,J) = size(K)
    @assert I == size(K,1)
    @assert J == size(K,2)
    @inbounds for j in 1:J, l in 1:I
            val = α_kin[l,j]*S[l]/K[l,j]
            F_c[j] = F_c[j] + val[1]
    end
    return F_c
end

function ECA_kinetics!(ECA::Matrix{T}, S::Vector{T}, E::Vector{T}, K::Matrix{T}, k2p::Vector{T}, α_kin::Matrix{T}, t::String) where {T<:Real}
    (I,J) = size(ECA)
    @assert I == size(S,1)
    @assert J == size(E,1)
    F_c = zeros(J)
    F_r = zeros(I)
    F_r = conjugate_substrate_flux!(F_r, E, K)
    F_c = normalized_substrate_flux!(F_c, S, K, α_kin)
    @inbounds for i in 1:I, j in 1:J
            val = k2p[j]*E[j]*(α_kin[i,j]*S[i]/K[i,j])/(1 + F_r[i] + F_c[j])
            ECA[i,j] = val[1]
    end
    any(x->x==true, isnan.(ECA)) ? throw(DomainError("NaN in DEBmicroTrait.ECA_kinetics!")) : return ECA
end
