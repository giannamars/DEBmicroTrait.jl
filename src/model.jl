########################################################################################################################
# init_batch_model
function init_batch_model(id_isolate, id_monomer, assimilation, enzymes, maintenance, protein_synthesis, turnover)
    #
    n_polymers = 0
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
    α_X        = enzymes["alpha"][id_isolate]
    y_EX       = y_EV
    f_αX       = ones(n_enzymes)
    min_gt     = protein_synthesis["mingt"][id_isolate]
    p_met      = MetabolismC([k_E], [y_EV], [k_M], [y_EM], [α_X], [y_EX], min_gt)
    #
    N_SB       = assimilation["NSB"][id_monomer,id_isolate]
    K_D        = assimilation["KD"][id_monomer,id_isolate]
    y_DE       = assimilation["yDE"][id_monomer]
    N_C        = assimilation["NC"][id_monomer]
    N_X        = 1
    p_ass      = AssimilationC(N_SB*ones(1,1),K_D*ones(1,1),[y_DE],[N_C],N_X)
    #
    γ_V0       = turnover["gV0"][id_isolate]
    γ_V1       = turnover["gV1"]
    γ_X        = 1/(7*24)*ones(n_enzymes)
    #
    γ_D_ads    = zeros(n_monomers)
    γ_X_ads    = zeros(n_enzymes)
    f_ED       = ones(n_monomers)
    f_V        = ones(n_microbes)    # all structural biomass recycling to monomers
    f_VD       = ones(n_microbes)
    f_VP       = ones(n_microbes)
    f_X        = zeros(n_enzymes)    # all enzymes recycling to monomers
    f_XD       = ones(n_enzymes)
    f_XP       = ones(n_enzymes)
    p_turn     = Turnover([γ_V0],[γ_V1],γ_X,γ_D_ads,γ_X_ads,f_ED,f_VD,f_VP,f_V,f_XD,f_XP,f_X)
    #
    p = Params(p_set,p_met,p_ass,nothing,p_turn)
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
    D, E, V, X, CO2 = DEBmicroTrait.split_state_batch(u, p)
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
    # system
    @. du[1+n_polymers:n_polymers+n_monomers] = - J_D + J_ED + J_VD + J_XD
    @. du[1+n_polymers+n_monomers:n_polymers+n_monomers+n_microbes] =  J_DE - (p.metabolism_pars.k_E - r)*E - J_E
    @. du[1+n_polymers+n_monomers+n_microbes:n_polymers+n_monomers+2*n_microbes] = r*V - J_V
    @. du[1+n_polymers+n_monomers+2*n_microbes:n_polymers+n_monomers+2*n_microbes+n_enzymes] = J_EX - J_X
    @. du[1+n_polymers+n_monomers+2*n_microbes+n_enzymes:n_polymers+n_monomers+2*n_microbes+n_enzymes+n_microbes] = rG_CO2 + rX_CO2 + rM_CO2 + J_DE_CO2
    return du
end
########################################################################################################################
