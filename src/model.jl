using Distributions
########################################################################################################################
# init_batch_model
function init_batch_model(id_isolate, id_monomer, id_polymer, assimilation, enzymes, maintenance, protein_synthesis, turnover)
    #
    n_polymers = 1
    n_monomers = 1
    n_microbes = 1
    n_enzymes  = 1
    n_minerals = 0
    p_set      = Setup(n_polymers, n_monomers, n_microbes, n_enzymes, n_minerals)
    #
    k_E        = protein_synthesis["kE"][id_isolate]
    y_EV       = protein_synthesis["yEV"][id_isolate]
    k_M        = maintenance["kM"][id_isolate]
    y_EM       = assimilation["yEM"][id_isolate]
    α_X        = enzymes["alpha"][id_isolate].*100
    y_EX       = y_EV
    f_αX       = enzymes["f_alphaX"][id_polymer, id_isolate]
    min_gt     = protein_synthesis["mingt"][id_isolate]
    p_met      = MetabolismC([k_E], [y_EV], [k_M], [y_EM], [α_X], [y_EX], [f_αX], min_gt)
    #
    N_SB       = assimilation["NSB"][id_monomer,id_isolate]
    K_D        = assimilation["KD"][id_monomer,id_isolate]
    y_DE       = assimilation["yDE"][id_monomer]
    N_C        = assimilation["NC"][id_monomer]
    p_ass      = AssimilationC(N_SB*ones(1,1),K_D*ones(1,1),[y_DE],[N_C])
    #
    γ_V0       = turnover["gV0"][id_isolate]
    γ_V1       = turnover["gV1"]
    γ_X        = 1/(7*24)*ones(n_enzymes)
    #
    f_ED      = rand(Dirichlet(n_monomers,1))
    f_VD      = rand(Dirichlet(n_monomers,1))
    f_VP      = rand(Dirichlet(n_polymers,1))
    f_V       = ones(n_microbes)
    f_XD      = rand(Dirichlet(n_monomers,1))
    f_XP      = rand(Dirichlet(n_polymers,1))
    f_X       = ones(n_enzymes)
    p_turn     = Turnover([γ_V0],[γ_V1],γ_X,f_ED,f_VD,f_VP,f_V,f_XD,f_XP,f_X)
    #
    Χ_0        = [6, 2, 6, 5, 6, 6, 6][id_polymer] # C atoms per monomer [cellulose, chitin, mannan, xylan, xyloglucan, glucan, protein]
    ρ_P        = [1.55, 1.45, 1.5, 1.4, 1.45, 1.5, 1.35][id_polymer]
    R_E        = [4.1e-9, 4e-9, 3.9e-9, 3.7e-9, 4.0e-9, 4.1e-9, 5e-9][id_polymer] # enzyme radius
    V_E        = [3.0, 5.0, 4.0, 3.5, 3.0, 4.0, 5.0][id_polymer] # specific hydrolysis rate
    MW_E       = [55e3, 60e3, 50e3, 45e3, 50e3, 55e3, 60e3][id_polymer] # enzyme molecular weight
    D_E        = DEBmicroTrait.aqueous_diffusivity([MW_E])
    #
    R_P        = 20e-3*ones(n_polymers)
    #
    α_kin_P    = DEBmicroTrait.binding_sites([ρ_P], [R_E])
    K_EP_0     = DEBmicroTrait.polymer_affinity([R_E], R_P, [V_E], D_E)
    #
    f_PD       = zeros(n_polymers, n_monomers)
    for i in 1:n_polymers
        f_PD[i,:] = rand(Dirichlet(n_monomers,1))
    end
    p_depoly = Depolymerization([Χ_0], [V_E]*60*60*24, α_kin_P, K_EP_0, f_PD)
    #
    p = Params(p_set,p_met,p_ass,p_depoly,p_turn)
end
########################################################################################################################

# init_mixed_medium
function init_mixed_medium(id_isolate, assimilation, enzymes, maintenance, protein_synthesis, turnover)
    #
    n_polymers = 0
    n_monomers = 83
    n_microbes = 1
    n_enzymes  = 1
    n_minerals = 0
    p_set      = Setup(n_polymers, n_monomers, n_microbes, n_enzymes, n_minerals)
    #
    k_E        = protein_synthesis["kE"][id_isolate]
    y_EV       = protein_synthesis["yEV"][id_isolate]
    k_M        = maintenance["kM"][id_isolate]
    y_EM       = assimilation["yEM"][id_isolate]
    α_X        = enzymes["alpha"][id_isolate]
    y_EX       = y_EV
    f_αX       = ones(n_enzymes)
    min_gt     = protein_synthesis["mingt"][id_isolate]
    p_met      = MetabolismC([k_E], [y_EV], [k_M], [y_EM], [α_X], [y_EX], min_gt)
    #
    N_SB       = reshape(assimilation["NSB"][:,id_isolate], n_monomers, n_microbes)
    K_D        = reshape(assimilation["KD"][:,id_isolate], n_monomers, n_microbes)
    y_DE       = reshape(assimilation["yDE"][:,id_isolate], n_monomers, n_microbes)
    N_C        = assimilation["NC"].*4.5
    N_X        = ones(n_monomers)
    p_ass      = AssimilationC(N_SB,K_D,y_DE,N_C,N_X)
    #
    γ_V0       = turnover["gV0"][id_isolate]
    γ_V1       = turnover["gV1"]
    γ_X        = 1/(7*24)*ones(n_enzymes)
    #
    γ_D_ads    = zeros(n_monomers)
    γ_X_ads    = zeros(n_enzymes)
    f_ED       = zeros(n_monomers)
    f_V        = zeros(n_microbes)    # no structural biomass recycling to monomers
    f_VD       = ones(n_microbes)
    f_VP       = ones(n_microbes)
    f_X        = ones(n_enzymes)    # no enzymes recycling to monomers
    f_XD       = ones(n_enzymes)
    f_XP       = ones(n_enzymes)
    p_turn     = Turnover([γ_V0],[γ_V1],γ_X,γ_D_ads,γ_X_ads,f_ED,f_VD,f_VP,f_V,f_XD,f_XP,f_X)
    #
    p = Params(p_set,p_met,p_ass,nothing,p_turn)
end

# init_mixed_medium
function init_mixed_medium_r2(id_isolate, assimilation, enzymes, maintenance, protein_synthesis, turnover)
    #
    n_polymers = 0
    n_monomers = 43
    n_microbes = 1
    n_enzymes  = 1
    n_minerals = 0
    p_set      = Setup(n_polymers, n_monomers, n_microbes, n_enzymes, n_minerals)
    #
    k_E        = protein_synthesis["kE"][id_isolate]
    y_EV       = protein_synthesis["yEV"][id_isolate]
    k_M        = maintenance["kM"][id_isolate]
    y_EM       = assimilation["yEM"][id_isolate]
    α_X        = enzymes["alpha"][id_isolate]
    y_EX       = y_EV
    f_αX       = ones(n_enzymes)
    min_gt     = protein_synthesis["mingt"][id_isolate]
    p_met      = MetabolismC([k_E], [y_EV], [k_M], [y_EM], [α_X], [y_EX], min_gt)
    #
    N_SB       = reshape(assimilation["NSB"][:,id_isolate], n_monomers, n_microbes)
    K_D        = reshape(assimilation["KD"][:,id_isolate], n_monomers, n_microbes)
    y_DE       = reshape(assimilation["yDE"][:,id_isolate], n_monomers, n_microbes)
    N_C        = assimilation["NC"].*4.5
    N_X        = ones(n_monomers)
    p_ass      = AssimilationC(N_SB,K_D,y_DE,N_C,N_X)
    #
    γ_V0       = turnover["gV0"][id_isolate]
    γ_V1       = turnover["gV1"]
    γ_X        = 1/(7*24)*ones(n_enzymes)
    #
    γ_D_ads    = zeros(n_monomers)
    γ_X_ads    = zeros(n_enzymes)
    f_ED       = zeros(n_monomers)
    f_V        = zeros(n_microbes)    # no structural biomass recycling to monomers
    f_VD       = ones(n_microbes)
    f_VP       = ones(n_microbes)
    f_X        = ones(n_enzymes)    # no enzymes recycling to monomers
    f_XD       = ones(n_enzymes)
    f_XP       = ones(n_enzymes)
    p_turn     = Turnover([γ_V0],[γ_V1],γ_X,γ_D_ads,γ_X_ads,f_ED,f_VD,f_VP,f_V,f_XD,f_XP,f_X)
    #
    p = Params(p_set,p_met,p_ass,nothing,p_turn)
end

########################################################################################################################
# batch model
function batch_model!(du, u, p, t)
    P, D, E, V, X, CO2 = DEBmicroTrait.split_state_poly(u, p)
    # setup
    n_polymers                = p.setup_pars.n_polymers
    n_monomers                = p.setup_pars.n_monomers
    n_microbes                = p.setup_pars.n_microbes
    n_enzymes                 = p.setup_pars.n_enzymes
    # metabolism
    r                         = growth!(0.0*ones(n_microbes), p.metabolism_pars, E, V)
    x, rG_CO2, rM_CO2, rX_CO2 = growth_production!(r, p.metabolism_pars, E, V)
    J_EX                      = enzyme_production!(x, p.metabolism_pars, V)
    # assimilation
    J_DE         = assimilation!(zeros(n_microbes), p.assimilation_pars, D, V)
    J_DE_CO2     = assimilation_production!(zeros(n_microbes), p.assimilation_pars, D, V)
    J_D          = uptake!(zeros(n_monomers), p.assimilation_pars, D, V)
    # turnover
    J_ED         = reserve_recycling!(zeros(n_monomers), p.turnover_pars, E)
    J_X          = enzyme_decay!(zeros(n_enzymes), p.turnover_pars, X)
    J_XD, J_XP   = enzyme_recycling!(zeros(n_monomers), p.turnover_pars, X)
    J_V          = biomass_turnover!(zeros(n_microbes), p.turnover_pars, V)
    J_VD, J_VP   = biomass_recycling!(zeros(n_monomers), p.turnover_pars, V)
    J_E          = biomass_turnover!(zeros(n_microbes), p.turnover_pars, E)
    # depolymerization
    J_P, J_PD   = DEBmicroTrait.depolymerization!(zeros(n_polymers), p.depolymerization_pars, P, X)
    # system
    @. du[1:n_polymers] = - J_P + J_VP + J_XP
    @. du[1+n_polymers:n_polymers+n_monomers] = J_PD - J_D + J_ED + J_VD + J_XD
    @. du[1+n_polymers+n_monomers:n_polymers+n_monomers+n_microbes] =  J_DE - (p.metabolism_pars.k_E - r)*E - J_E
    @. du[1+n_polymers+n_monomers+n_microbes:n_polymers+n_monomers+2*n_microbes] = r*V - J_V
    @. du[1+n_polymers+n_monomers+2*n_microbes:n_polymers+n_monomers+2*n_microbes+n_enzymes] = J_EX - J_X
    @. du[1+n_polymers+n_monomers+2*n_microbes+n_enzymes:n_polymers+n_monomers+2*n_microbes+n_enzymes+n_microbes] = rG_CO2 + rX_CO2 + rM_CO2 + J_DE_CO2
    return du
end
########################################################################################################################
