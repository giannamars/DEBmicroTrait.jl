abstract type AbstractDepolymerization end

struct Depolymerization{X0,VE,AKIN,KEP,PD} <: AbstractDepolymerization
    Χ_0::X0     # Number of C atoms per product molecule
    V_E::VE     # Specific hydrolysis rate
    α_kin::AKIN # Enzyme binding sites per polymer molecule
    K_EP::KEP   # Depolymerization half-saturation constant
    f_PD::PD    # Polymer-monomer association
end

for fn in fieldnames(Depolymerization)
    @eval $fn(p::Depolymerization) = p.$fn
end

depolymerization!(J_P::Vector{Float64}, p::Depolymerization, P::Vector{Float64}, X::Vector{Float64}) = begin
    n_substrates = size(P,1)
    n_consumers  = size(X,1)
    ECA  = ECA_kinetics!(zeros(n_substrates,n_consumers), P, X, K_EP(p), V_E(p), α_kin(p), "depoly")
    J_P = vcat(sum(ECA[1:n_substrates, :], dims=2)...)
    J_PD = sum(f_PD(p).*J_P, dims=1)[:]
    return J_P, J_PD
end
