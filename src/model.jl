########################################################################################################################
# init_batch_model
function init_batch_model(assimilation, enzymes, maintenance, protein_synthesis, turnover, avena_ex, avena_root)
    #
    n_polymers = 0
    n_monomers = 83
    n_microbes = 39
    n_enzymes  = 1
    n_minerals = 0
    p_set      = Setup(n_polymers, n_monomers, n_microbes, n_enzymes, n_minerals)
    #
    k_E        = protein_synthesis["kE"];
    y_EV       = protein_synthesis["yEV"];
    k_M        = maintenance["kM"];
    y_EM       = assimilation["yEM"];
    α_X        = enzymes["alpha"];
    y_EX       = y_EV;
    f_αX       = rand(Dirichlet(n_enzymes, 1));
    min_gt     = protein_synthesis["mingt"];
    p_met      = MetabolismC(k_E, y_EV, k_M, y_EM, α_X, y_EX, f_αX, min_gt);
    #
    N_SB       = assimilation["NSB"]
    K_D        = assimilation["KD"]
    y_DE       = assimilation["yDE"]
    N_C        = assimilation["NC"]
    N_X        = 1
    p_ass      = AssimilationC(N_SB,K_D,y_DE,N_C,N_X)
    #
    γ_V0       = turnover["gV0"]
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
    p_turn     = Turnover(γ_V0,γ_V1,γ_X,γ_D_ads,γ_X_ads,f_ED,f_VD,f_VP,f_V,f_XD,f_XP,f_X)
    #
    γ_root                  = 5e-5;
    p_forcing = Forcing(avena_root["root_C"], γ_root, avena_ex["exudation_rate"]);
    #
    p = ParamsF(p_set,p_met,p_ass,nothing,p_turn,p_forcing)
    return p
end
########################################################################################################################

# init_mixed_medium
function init_mixed_medium(id_isolate, n_monomers, assimilation, enzymes, maintenance, protein_synthesis, turnover)
    #
    n_polymers = 0
    #n_monomers = 83
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
    p_met      = MetabolismC([k_E], [y_EV], [k_M], [y_EM], [α_X], [y_EX], f_αX, min_gt)
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
# init_rhizosphere_model
function init_rhizosphere_model(assimilation, enzymes, maintenance, protein_synthesis, turnover, avena_ex, avena_root)
    #
    n_polymers = 1;
    n_monomers = 6;
    n_microbes = 39;
    n_enzymes  = 1;
    n_minerals = 1;
    p_set      = Setup(n_polymers, n_monomers, n_microbes, n_enzymes, n_minerals);
    #
    k_E        = protein_synthesis["kE"];
    y_EV       = protein_synthesis["yEV"];
    k_M        = maintenance["kM"];
    y_EM       = assimilation["yEM"];
    α_X        = enzymes["alpha"];
    y_EX       = y_EV;
    f_αX       = rand(Dirichlet(n_enzymes, 1));
    min_gt     = protein_synthesis["mingt"];
    p_met      = MetabolismC(k_E, y_EV, k_M, y_EM, α_X, y_EX, f_αX, min_gt);
    #
    N_SB_D     = assimilation["NSB"];
    N_SB_M     = ones(n_monomers,n_minerals);
    N_SB       = hcat(N_SB_D, N_SB_M);
    K_D_0      = assimilation["KD"];
    K_M_0      = 2.4*ones(n_monomers,n_minerals);
    K_D        = hcat(K_D_0, K_M_0);
    y_DE       = assimilation["yDE"];
    N_C        = assimilation["NC"];
    N_X        = 1.0*ones(n_monomers);
    #M          = 50e-6*ones(n_minerals);
    M          = zeros(n_minerals);
    p_ass      = AssimilationCM(N_SB, K_D, y_DE, N_C, N_X, M);
    #
    γ_V0       = turnover["gV0"];
    #γ_V0       = turnover["gV0"]/100;
    γ_V1       = turnover["gV1"];
    γ_X        = 1/(7*24)*ones(n_enzymes);
    #γ_X        = 1/(7*24*10000)*ones(n_enzymes);
    γ_D_ads    = 1e-6*ones(n_monomers);
    γ_X_ads    = 1e-6*ones(n_enzymes);
    f_ED       = rand(Dirichlet(n_monomers, 1));
    #f_ED       = zeros(n_monomers);
    #f_V        = rand(n_microbes);    # fraction of structural biomass recycling to monomers
    #f_V        = rand(n_microbes);    # fraction of structural biomass recycling to monomers
    f_V        = zeros(n_microbes);    # fraction of structural biomass recycling to monomers
    f_VD       = rand(Dirichlet(n_monomers, 1));
    #f_VD       = zeros(n_monomers);
    f_VP       = rand(Dirichlet(n_polymers, 1));
    #f_VP       = zeros(n_polymers);
    #f_X        = rand(n_enzymes);   # fraction of enzymes recycling to monomers
    f_X        = zeros(n_enzymes);   # fraction of enzymes recycling to monomers
    f_XD       = rand(Dirichlet(n_monomers, 1));
    #f_XD       = zeros(n_monomers);
    f_XP       = rand(Dirichlet(n_polymers, 1));
    #f_XP       = zeros(n_polymers);
    p_turn     = Turnover(γ_V0,γ_V1,γ_X,γ_D_ads,γ_X_ads,f_ED,f_VD,f_VP,f_V,f_XD,f_XP,f_X);
    #
    Χ_0        = 10.0*ones(n_polymers);
    V_E        = 5.0*ones(n_enzymes);
    α_kin_P    = 3.79e-7*ones(n_polymers,n_enzymes);
    #α_kin_P    = 3.79e-3*ones(n_polymers,n_enzymes);
    K_EP_P     = 3.05*ones(n_polymers,n_enzymes);
    α_kin_M    = ones(n_minerals, n_enzymes);
    K_MX       = ones(n_minerals, n_enzymes);
    #α_kin_M    = 1e-6*ones(n_minerals, n_enzymes);
    #K_MX       = 1e-6*ones(n_minerals, n_enzymes);
    α_kin      = vcat(α_kin_P, α_kin_M);
    K_EP       = vcat(K_EP_P,K_MX);
    f_PD       = zeros(n_polymers, n_monomers);
    for i in 1:n_polymers
        f_PD[i,:] = rand(Dirichlet(n_monomers,1))
    end
    p_depoly   = DepolymerizationM(Χ_0, V_E, α_kin, K_EP, M, f_PD);
    #
    γ_root                  = 5e-5;
    p_forcing = Forcing(avena_root["root_C"], γ_root, avena_ex["exudation_rate"]./10.0);
    #
    p = ParamsF(p_set, p_met, p_ass, p_depoly, p_turn, p_forcing);
    return p
end
########################################################################################################################

########################################################################################################################
# init_ReSOM_model
function init_ReSOM_model(n_polymers, n_monomers, n_microbes, n_enzymes, n_minerals)
    #
    p_set      = Setup(n_polymers, n_monomers, n_microbes, n_enzymes, n_minerals)
    #
    k_E        = 1.34*ones(n_microbes)
    y_EV       = 1.51*ones(n_microbes)
    k_M        = 8.5e-5*ones(n_microbes)
    y_EM       = 1.0*ones(n_microbes)
    α_X        = 0.0027*ones(n_microbes)
    y_EX       = y_EV
    f_αX       = rand(Dirichlet(n_enzymes, 1))
    min_gt     = 4.36*ones(n_microbes)
    p_met      = MetabolismC(k_E, y_EV, k_M, y_EM, α_X, y_EX, f_αX, min_gt)
    #
    N_SB_D     = 3.95e-8*ones(n_monomers,n_microbes)
    N_SB_M     = ones(n_monomers,n_minerals)
    N_SB       = hcat(N_SB_D, N_SB_M)
    K_D_0      = 0.0061*ones(n_monomers,n_microbes)
    K_M_0      = 2.4*ones(n_monomers,n_minerals)
    K_D        = hcat(K_D_0, K_M_0)
    y_DE       = 0.2*ones(n_monomers, n_microbes)
    N_C        = 5.0*ones(n_monomers)
    N_X        = 1.0*ones(n_monomers)
    M          = 50e-6*ones(n_minerals)
    p_ass      = AssimilationCM(N_SB,K_D,y_DE,N_C,N_X,M)
    #
    γ_V0       = 0.26*ones(n_microbes)
    γ_V1       = 98.8*ones(n_microbes)
    γ_X        = 1/(7*24)*ones(n_enzymes)
    γ_D_ads    = rand(n_monomers)
    γ_X_ads    = rand(n_enzymes)
    f_ED       = rand(Dirichlet(n_monomers, 1))
    f_V        = rand(n_microbes)    # fraction of structural biomass recycling to monomers
    f_VD       = rand(Dirichlet(n_monomers, 1))
    f_VP       = rand(Dirichlet(n_polymers, 1))
    f_X        = rand(n_enzymes)    # fraction of enzymes recycling to monomers
    f_XD       = rand(Dirichlet(n_monomers, 1))
    f_XP       = rand(Dirichlet(n_polymers, 1))
    p_turn     = Turnover(γ_V0,γ_V1,γ_X,γ_D_ads,γ_X_ads,f_ED,f_VD,f_VP,f_V,f_XD,f_XP,f_X)
    #
    Χ_0        = 10.0*ones(n_polymers)
    V_E        = 5.0*ones(n_enzymes)
    α_kin_P    = 3.79e-7*ones(n_polymers,n_enzymes)
    K_EP_P     = 3.05*ones(n_polymers,n_enzymes)
    α_kin_M    = ones(n_minerals, n_enzymes)
    K_MX       = ones(n_minerals, n_enzymes)
    α_kin      = vcat(α_kin_P, α_kin_M)
    K_EP       = vcat(K_EP_P,K_MX)
    f_PD       = zeros(n_polymers, n_monomers)
    for i in 1:n_polymers
        f_PD[i,:] = rand(Dirichlet(n_monomers,1))
    end
    p_depoly   = DepolymerizationM(Χ_0, V_E, α_kin, K_EP, M, f_PD)
    #
    p = Params(p_set,p_met,p_ass,p_depoly,p_turn)
end
########################################################################################################################


########################################################################################################################
# rhizosphere model
function rhizosphere_model!(du, u, p, t)
    P, D, E, V, X, D_ads, X_ads, CO2 = split_state_ReSOM(u, p)
    # setup
    n_polymers                = p.setup_pars.n_polymers
    n_monomers                = p.setup_pars.n_monomers
    n_microbes                = p.setup_pars.n_microbes
    n_enzymes                 = p.setup_pars.n_enzymes
    n_minerals                = p.setup_pars.n_minerals
    # metabolism
    r                         = growth!(0.0*ones(n_microbes), p.metabolism_pars, E, V)
    x, rG_CO2, rM_CO2, rX_CO2 = growth_production!(r, p.metabolism_pars, E, V)
    J_EX                      = enzyme_production!(x, p.metabolism_pars, V)
    # assimilation
    J_DE         = assimilation!(zeros(n_monomers,n_microbes+n_minerals), p.assimilation_pars, D, V)
    J_DE_CO2     = assimilation_production!(zeros(n_monomers,n_microbes+n_minerals), p.assimilation_pars, D, V)
    J_D, J_DM    = uptake!(zeros(n_monomers,n_microbes+n_minerals), p.assimilation_pars, D, V)
    # turnover
    J_ED         = reserve_recycling!(zeros(n_monomers), p.turnover_pars, E)
    J_X          = enzyme_decay!(zeros(n_enzymes), p.turnover_pars, X)
    J_XD, J_XP   = enzyme_recycling!(zeros(n_monomers), p.turnover_pars, X)
    J_V          = biomass_turnover!(zeros(n_microbes), p.turnover_pars, V)
    J_VD, J_VP   = biomass_recycling!(zeros(n_monomers), p.turnover_pars, V)
    J_E          = biomass_turnover!(zeros(n_microbes), p.turnover_pars, E)
    J_D_ads      = adsorbed_monomer_decay!(zeros(n_monomers), p.turnover_pars, D_ads)
    J_X_ads      = adsorbed_enzyme_decay!(zeros(n_enzymes), p.turnover_pars, X_ads)
     # depolymerization
    J_P, J_PD, J_XM    = depolymerization!(zeros(n_polymers), p.depolymerization_pars, P, X)
    # system
    @. du[1:n_polymers] = - J_P + J_VP + J_XP
    @. du[1+n_polymers:n_polymers+n_monomers] = - J_D + J_PD - J_DM + J_ED + J_VD + J_XD + J_D_ads
    @. du[1+n_polymers+n_monomers:n_polymers+n_monomers+n_microbes] =  J_DE - (p.metabolism_pars.k_E - r)*E - J_E
    @. du[1+n_polymers+n_monomers+n_microbes:n_polymers+n_monomers+2*n_microbes] = r*V - J_V
    @. du[1+n_polymers+n_monomers+2*n_microbes:n_polymers+n_monomers+2*n_microbes+n_enzymes] = J_EX - J_X - J_XM + J_X_ads
    @. du[1+n_polymers+n_monomers+2*n_microbes+n_enzymes:n_polymers+n_monomers+2*n_microbes+n_enzymes+n_monomers] = J_DM - J_D_ads
    @. du[1+n_polymers+n_monomers+2*n_microbes+n_enzymes+n_monomers:n_polymers+n_monomers+2*n_microbes+2*n_enzymes+n_monomers] = J_XM - J_X_ads
    @. du[1+n_polymers+n_monomers+2*n_microbes+2*n_enzymes+n_monomers:n_polymers+n_monomers+2*n_microbes+2*n_enzymes+n_monomers+n_microbes] = rG_CO2 + rM_CO2 + rX_CO2 + J_DE_CO2
    return du
end
########################################################################################################################


# ReSOM model
########################################################################################################################

function ReSOM_model!(du, u, p, t)

    # setup
    P, D, E, V, X, D_ads, X_ads, CO2 = split_state_ReSOM(u, p)
    n_polymers                = p.setup_pars.n_polymers
    n_monomers                = p.setup_pars.n_monomers
    n_microbes                = p.setup_pars.n_microbes
    n_enzymes                 = p.setup_pars.n_enzymes
    n_minerals                = p.setup_pars.n_minerals

    # metabolism
    r                         = growth!(0.0*ones(n_microbes), p.metabolism_pars, E, V)
    x, rG_CO2, rM_CO2, rX_CO2 = growth_production!(r, p.metabolism_pars, E, V)
    J_EX                      = enzyme_production!(x, p.metabolism_pars, V)

    # assimilation
    J_DE         = assimilation!(zeros(n_monomers,n_microbes+n_minerals), p.assimilation_pars, D, V)
    J_DE_CO2     = assimilation_production!(zeros(n_monomers,n_microbes+n_minerals), p.assimilation_pars, D, V)
    J_D, J_DM    = uptake!(zeros(n_monomers,n_microbes+n_minerals), p.assimilation_pars, D, V)

    # turnover
    J_ED         = reserve_recycling!(zeros(n_monomers), p.turnover_pars, E)
    J_X          = enzyme_decay!(zeros(n_enzymes), p.turnover_pars, X)
    J_XD, J_XP   = enzyme_recycling!(zeros(n_monomers), p.turnover_pars, X)
    J_V          = biomass_turnover!(zeros(n_microbes), p.turnover_pars, V)
    J_VD, J_VP   = biomass_recycling!(zeros(n_monomers), p.turnover_pars, V)
    J_E          = biomass_turnover!(zeros(n_microbes), p.turnover_pars, E)
    J_D_ads      = adsorbed_monomer_decay!(zeros(n_monomers), p.turnover_pars, D_ads)
    J_X_ads      = adsorbed_enzyme_decay!(zeros(n_enzymes), p.turnover_pars, X_ads)

    # depolymerization
    J_P, J_PD, J_XM    = depolymerization!(zeros(n_polymers), p.depolymerization_pars, P, X)

    # forcing 
    spl_ex   = root_exudation(p.forcing_pars)
   # spl_root = root_decay(p.forcing_pars)

    # rhs
    @. du[1:n_polymers] = - J_P + J_VP + J_XP
    #@. du[1+n_polymers:n_polymers+n_monomers] = spl_ex(t)
    @. du[1+n_polymers:n_polymers+n_monomers] = spl_ex(t)  - J_D + J_PD - J_DM + J_ED + J_VD + J_XD + J_D_ads
    @. du[1+n_polymers+n_monomers:n_polymers+n_monomers+n_microbes] =  J_DE - (p.metabolism_pars.k_E - r)*E - J_E
    @. du[1+n_polymers+n_monomers+n_microbes:n_polymers+n_monomers+2*n_microbes] = r*V - J_V
    @. du[1+n_polymers+n_monomers+2*n_microbes:n_polymers+n_monomers+2*n_microbes+n_enzymes] = J_EX - J_X - J_XM + J_X_ads
    @. du[1+n_polymers+n_monomers+2*n_microbes+n_enzymes:n_polymers+n_monomers+2*n_microbes+n_enzymes+n_monomers] = J_DM - J_D_ads
    @. du[1+n_polymers+n_monomers+2*n_microbes+n_enzymes+n_monomers:n_polymers+n_monomers+2*n_microbes+2*n_enzymes+n_monomers] = J_XM - J_X_ads
    @. du[1+n_polymers+n_monomers+2*n_microbes+2*n_enzymes+n_monomers:n_polymers+n_monomers+2*n_microbes+2*n_enzymes+n_monomers+n_microbes] = rG_CO2 + rM_CO2 + rX_CO2 + J_DE_CO2

    return du
end

########################################################################################################################
