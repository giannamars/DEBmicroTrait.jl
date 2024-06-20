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


function SUPECA_kinetics!(SUPECA::Array{T}, A::Vector{T}, B::Vector{T}, E::Vector{T}, K::Array{T}, k2p::Array{T}, N_SB::Array{T}) where {T<:Real}
    (I,J,L) = size(K)
    @assert I == size(A,1)
    @assert J == size(B,1)
    @assert L == size(E,1)
    F_c_A = zeros(L)
    F_c_B = zeros(L)
    F_r_A = zeros(I)
    F_r_B = zeros(J)
    G_A   = zeros(I,L)
    G_B   = zeros(J,L)
    G_AB  = zeros(I,J,L)
    F_r_A = conjugate_substrate_flux!(F_r_A, E, K[:,1,:], N_SB[:,1,:])
    F_r_B = conjugate_substrate_flux!(F_r_B, E, K[1,:,:], N_SB[1,:,:])
    F_c_A = normalized_substrate_flux!(F_c_A, A, K[:,1,:])
    F_c_B = normalized_substrate_flux!(F_c_B, B, K[1,:,:])
    F_c_AB = @. F_c_A + F_c_B

    @inbounds for i in 1:I, j in 1:J, k in 1:L
            G_A[i,k] = F_c_A[k] + F_r_A[i]
            G_B[j,k] = F_c_B[k] + F_r_B[j]
            G_AB[i,j,k] = G_A[i,k] + G_B[j,k]
            val = k2p[i,j,k]*N_SB[i,j,k]*E[k]*A[i]/K[i,1,k]*B[j]/K[1,j,k]/(G_A[i,k]*G_B[j,k]/G_AB[i,j,k]*F_c_AB[k] + F_c_AB[k] - (F_c_A[k]*G_B[j,k] + G_A[i,k]*F_c_B[k] - G_A[i,k]*G_B[j,k])/G_AB[i,j,k])
            SUPECA[i,j,k] = val
    end

    any(x->x==true, isnan.(SUPECA)) ? throw(DomainError("NaN in DEBmicroTrait.SUPECA_kinetics!")) : return SUPECA
end