function aqueous_diffusivity(molecular_weight::Array{Float64,1})
    # Perry and Hilton (1973)
    D_S = @. 4.36e-9*molecular_weight^(-0.386)
end

function clapp_hornberger_par(pct_sand::Float64, pct_clay::Float64)
    chb = 2.91 + 0.195*pct_clay
    sat = 0.489 - 0.00126*pct_sand
    psisat = -10.0^(1.88-0.0131*pct_sand)
    ksat = 0.0070556*10.0^(-0.884+0.0153*pct_sand)
    return chb, sat, psisat, ksat
end

function moldrup_tau(por, epsi, theta, kappa)
    taug = epsi*(epsi/por)^(3.0/kappa)
    tauw = theta*((theta+1e-20)/por)^(kappa/3.0-1.0)
    return taug, tauw
end

function cosby_psi(s_sat, psisat, sat, chb)
    psi=max(psisat*(s_sat+1.e-10)^(-chb)*1.e-3*1.e4, -1e8*1e-3*1e4)
    dpsidvsm = -chb*psi/(s_sat*sat+1e-20)
    return psi, dpsidvsm
end

function cosby_hk(s_sat, ksat, chb)
    hk = ksat*(s_sat)^(2.0*chb+3.0)
end

function cosby_Dwpsi(s_sat, psisat, ksat, sat, chb)
    psi, dpsidvsm = cosby_psi(s_sat, psisat, sat, chb)
    hk = cosby_hk(s_sat, ksat, chb)
    Dwpsi = hk*1e-3*dpsidvsm*1e-4
end

function soil_affinity_properties(pct_sand::Float64, pct_clay::Float64, s_sat::Float64)
    chb, sat, psisat, ksat = clapp_hornberger_par(pct_sand, pct_clay)
    theta = s_sat*sat # phi_w
    epsi = sat-theta
    psi, dpsidvsm = cosby_psi(s_sat, psisat, sat, chb)
    film = max(exp(-13.65-0.857*log(-psi*1e-6)), 1e-7)
    τ_g, τ_w = moldrup_tau(sat, epsi, theta, chb)
    return theta, τ_g, τ_w, film
end

function effective_diffusivity(molecular_weight::Array{Float64,1}, pct_sand::Float64, pct_clay::Float64, s_sat::Float64)
    D_S = @. 4.36e-9*molecular_weight^(-0.386)
    ϕ_w, τ_g, τ_w, δ  = soil_affinity_properties(pct_sand, pct_clay, s_sat)
    D_eff = D_S*τ_w*ϕ_w
    return D_eff.*1e-9
end