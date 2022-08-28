abstract type AbstractTurnover end

@columns struct Turnover{gV0,gV1,gX,gDads,gXads,FED,FVD,FVP,FV,FXD,FXP,FX} <: AbstractTurnover
    γ_V_0::gV0     | 1/h       | "Max. reserve and structure turnover rate"
    γ_V_1::gV1     | mol/m^3   | "Half-saturation constant for reserve and structure turnover rate"
    γ_X::gX        | 1/h       | "Enzyme turnover rate"
    γ_D_ads::gDads | 1/h       | "Adsorbed monomer turnover rate"
    γ_X_ads::gXads | 1/h       | "Adsorbed enzyme turnover rate"
    f_ED::FED      | _         | "Fraction of decayed reserve recycling to monomers"
    f_VD::FVD      | _         | "Fraction of decayed biomass recycling to monomer fractions"
    f_VP::FVP      | _         | "Fraction of decayed biomass recycling to polymer fractions"
    f_V::FV        | _         | "Fraction of decayed biomass recycling to monomers"
    f_XD::FXD      | _         | "Fraction of decayed enzyme recycling to monomer fractions"
    f_XP::FXP      | _         | "Fraction of decayed enzyme recycling to polymer fractions"
    f_X::FX        | _         | "Fraction of decayed enzyme recycling to polymers"
end

for fn in fieldnames(Turnover)
    @eval $fn(p::Turnover) = p.$fn
end

biomass_turnover!(J_B::Vector{Float64}, p::AbstractTurnover, B::Vector{Float64}) = begin
    J_B = γ_V_0(p).*B./(γ_V_1(p) .+ B)
end


enzyme_decay!(J_X::Vector{Float64}, p::AbstractTurnover, X::Vector{Float64}) = begin
    J_X = γ_X(p).*X
end

reserve_recycling!(J_ED::Vector{Float64}, p::AbstractTurnover, E::Vector{Float64}) = begin
    J_ED = f_ED(p).*sum(γ_V_0(p).*E./(γ_V_1(p) .+ E), dims=1)
end


biomass_recycling!(J_VD::Vector{Float64}, p::AbstractTurnover, V::Vector{Float64}) = begin
    J_VD = f_VD(p).*sum(f_V(p).*γ_V_0(p).*V./(γ_V_1(p) .+ V), dims=1)
    J_VP = f_VP(p).*sum((1 .- f_V(p)).*sum(γ_V_0(p).*V./(γ_V_1(p) .+ V), dims=1), dims=1)
    return J_VD, J_VP
end

enzyme_recycling!(J_XD::Vector{Float64}, p::AbstractTurnover, X::Vector{Float64}) = begin
    J_XD = (1 .- f_X(p)).*f_XD(p).*sum(γ_X(p).*X, dims=1)
    J_XP = f_X(p).*f_XP(p).*sum(γ_X(p).*X, dims=1)
    return J_XD, J_XP
end
