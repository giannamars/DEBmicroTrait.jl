abstract type AbstractMetabolism end
abstract type AbstractMetabolismC <: AbstractMetabolism end

@columns struct MetabolismC{KE,YEV,KM,YEM,AX,YEX,MGT} <: AbstractMetabolismC
    k_E::KE     | 1/h       | "Reserve export rate"
    y_EV::YEV   | mol/mol   | "Yield of structure on reserve"
    k_M::KM     | 1/h       | "Specific maintenance rate"
    y_EM::YEM   | mol/mol   | "Maintenance yield on reserve vs structure"
    α_X::AX     | _         | "Fraction of growth flux allocated to enzyme production"
    y_EX::YEX   | mol/mol   | "Yield on enzyme production"
    min_gt::MGT | 1/h       | "Minimum generation time"
end

for fn in fieldnames(MetabolismC)
    @eval $fn(p::MetabolismC) = p.$fn
end

growth!(r0, p::MetabolismC, E::Vector{Float64}, V::Vector{Float64}) = begin
    y_VM = y_EM(p)./y_EV(p)

    function f(r, j)
      m_E   = max(1e-8, E[j]/V[j])
      j_EC  = m_E*(k_E(p)[j] - r)
      j_EM  = k_M(p)[j]*y_EM(p)[j]
      jEM   = max(1e-8, min(j_EC, j_EM))
      jVM   = (j_EM - jEM)*y_VM[j]/y_EM(p)[j]
      j_EG  = (j_EC - jEM)/(1.0 + α_X(p)[j])
      j_G   = j_EG/y_EV(p)[j]
      res   = r - (j_G - jVM)
    end

    f1 = [r->f(r,j) for j in 1:size(r0,1)]
    r = Roots.find_zero.(f1, r0)
end

growth_production!(r, p::MetabolismC, E::Vector{Float64}, V::Vector{Float64}) = begin
      y_VM  = y_EM(p)./y_EV(p)
      m_E   = E./V
      j_EC  = m_E.*(k_E(p) .- r)
      j_EM  = k_M(p).*y_EM(p)
      jEM   = min(j_EC, j_EM)
      jVM   = (j_EM - jEM).*y_VM./y_EM(p)
      j_EG  = (j_EC - jEM)./(1.0 .+ α_X(p))
      x     = α_X(p).*j_EG./y_EX(p)
      rG_CO2  = (1 .- 1 ./y_EV(p)).*j_EG.*V
      rM_CO2  = @. (jEM + jVM)*V
      rX_CO2  = (1 .- 1 ./y_EX(p)).*α_X(p).*j_EG.*V
      return x, rG_CO2, rM_CO2, rX_CO2
end

enzyme_production!(x, p::MetabolismC, V::Vector{Float64}) = begin
    J_EX = sum(x.*V, dims=1)
end
