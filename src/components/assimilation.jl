abstract type AbstractAssimilation end
abstract type AbstractAssimilationC <: AbstractAssimilation end

struct AssimilationC{NSB,KD,YDE,NC,NX} <: AbstractAssimilationC
    N_SB::NSB   # Mol transporters per mol of biomass carbon [mol/mol]
    K_D::KD     # Reference half-saturation constant [mol/m3]
    y_DE::YDE   # Yield on assimilation [mol/mol]
    N_C::NC     # Number of C atoms per monomer [-]
    N_X::NX     # Monomer id for return-on-investment [-]
end

for fn in fieldnames(AssimilationC)
    @eval $fn(p::AssimilationC) = p.$fn
end

assimilation!(J_DE::Vector{Float64}, p::AssimilationC, D::Vector{Float64}, V::Vector{Float64}) = begin
    n_substrates = size(D,1)
    n_consumers  = size(V,1)
    k2p = 180.0*60*60*ones(n_substrates, n_consumers)
    ECA  = ECA_kinetics!(zeros(n_substrates,n_consumers), D, V, K_D(p), k2p, N_SB(p))
    J_DE = vcat(sum((1.0 .- y_DE(p)).*N_C(p).*ECA[:, 1:n_consumers], dims=1)...)
end

uptake!(J_D::Vector{Float64}, p::AbstractAssimilation, D::Vector{Float64}, V::Vector{Float64}) = begin
    n_substrates = size(D,1)
    n_consumers  = size(V,1)
    k2p = 180.0*60*60*ones(n_substrates)
    ECA  = ECA_kinetics!(zeros(n_substrates,n_consumers), D, V, K_D(p), k2p, N_SB(p))
    J_D = vcat(sum(ECA[:, 1:n_consumers], dims=2)...).*N_C(p)
end

assimilation_production!(J_DE_CO2::Vector{Float64}, p::AbstractAssimilation, D::Vector{Float64}, V::Vector{Float64}) = begin
    n_substrates = size(D,1)
    n_consumers  = size(V,1)
    k2p = 180.0*60*60*ones(n_substrates)
    ECA  = ECA_kinetics!(zeros(n_substrates,n_consumers), D, V, K_D(p), k2p, N_SB(p))
    J_DE_CO2 = vcat(sum(y_DE(p).*N_C(p).*ECA[:, 1:n_consumers], dims=1)...)
end
