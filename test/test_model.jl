using DEBmicroTrait, Test

########################################
n_polymers = 0
n_monomers = 1
n_microbes = 1
n_enzymes  = 1
n_minerals = 0
p_set      = Setup(n_polymers, n_monomers, n_microbes, n_enzymes, n_minerals)
########################################

########################################
# split_state_batch
u     = ones(p_set.dim)
D    = u[1+n_polymers:n_polymers+n_monomers]
E    = u[1+n_polymers+n_monomers:n_polymers+n_monomers+n_microbes]
V    = u[1+n_polymers+n_monomers+n_microbes:n_polymers+n_monomers+2*n_microbes]
X    = u[1+n_polymers+n_monomers+2*n_microbes:n_polymers+n_monomers+2*n_microbes+n_enzymes]
CO2  = u[1+n_polymers+n_monomers+2*n_microbes+n_enzymes:n_polymers+n_monomers+2*n_microbes+n_enzymes+n_microbes]
########################################

########################################
# metabolism
k_E     = 0.2*ones(n_microbes)
y_EV    = 1.0*ones(n_microbes)
k_M     = 0.1*ones(n_microbes)
y_EM    = 1.0*ones(n_microbes)
α_X     = 0.0*ones(n_microbes)
y_EX    = 1.0*ones(n_microbes)
f_αX    = ones(n_enzymes)

p_met                     = MetabolismC(k_E, y_EV, k_M, y_EM, α_X, y_EX, f_αX)
r                         = DEBmicroTrait.growth!(0.0*ones(1), p_met, E, V)
x, rG_CO2, rM_CO2, rX_CO2 = DEBmicroTrait.growth_production!(r, p_met, E, V)
J_EX                      = DEBmicroTrait.enzyme_production!(x, p_met, V)
########################################

########################################
# assimilation
K_D_0        = rand(n_monomers,n_microbes)
K_D          = K_D_0
N_SB_D       = rand(n_monomers,n_microbes)
N_SB         = N_SB_D
N_C          = rand(n_monomers)
y_DE         = rand(n_monomers)
N_X          = zeros(n_monomers)

p_ass        = AssimilationC(N_SB,K_D,y_DE,N_C,N_X)
J_DE         = DEBmicroTrait.assimilation!(zeros(p_set.n_microbes), p_ass, D, V)
J_DE_CO2     = DEBmicroTrait.assimilation_production!(zeros(p_set.n_microbes), p_ass, D, V)
J_D          = DEBmicroTrait.uptake!(zeros(p_set.n_monomers), p_ass, D, V)
########################################

########################################
# turnover
γ_V_0        = ones(n_microbes)
γ_V_1        = ones(n_microbes)
γ_X          = ones(n_enzymes)
γ_D_ads      = zeros(n_monomers)
γ_X_ads      = zeros(n_enzymes)
f_ED         = ones(n_monomers)
f_VD         = ones(n_monomers)
f_VP         = ones(n_polymers)
f_V          = ones(n_microbes)
f_XD         = ones(n_monomers)
f_XP         = zeros(n_polymers)
f_X          = zeros(n_enzymes)

p_turn       = Turnover(γ_V_0,γ_V_1,γ_X,γ_D_ads,γ_X_ads,f_ED,f_VD,f_VP,f_V,f_XD,f_XP,f_X)
J_ED         = DEBmicroTrait.reserve_recycling!(zeros(p_set.n_monomers), p_turn, E)
J_X          = DEBmicroTrait.enzyme_decay!(zeros(p_set.n_enzymes), p_turn, X)
J_XD, J_XP   = DEBmicroTrait.enzyme_recycling!(zeros(p_set.n_monomers), p_turn, X)
J_V          = DEBmicroTrait.biomass_turnover!(zeros(p_set.n_microbes), p_turn, V)
J_VD, J_VP   = DEBmicroTrait.biomass_recycling!(zeros(p_set.n_monomers), p_turn, V)
J_E          = DEBmicroTrait.biomass_turnover!(zeros(p_set.n_microbes), p_turn, E)
########################################

########################################
# batch_model
du = zeros(p_set.dim)

dD = @. du[1+n_polymers:n_polymers+n_monomers] = - J_D + J_ED + J_VD + J_XD
dE = @. du[1+n_polymers+n_monomers:n_polymers+n_monomers+n_microbes] =  J_DE - (p_met.k_E - r)*E - J_E
dV = @. du[1+n_polymers+n_monomers+n_microbes:n_polymers+n_monomers+2*n_microbes] = r*V - J_V
dX = @. du[1+n_polymers+n_monomers+2*n_microbes:n_polymers+n_monomers+2*n_microbes+n_enzymes] = J_EX - J_X
dCO2 = @. du[1+n_polymers+n_monomers+2*n_microbes+n_enzymes:n_polymers+n_monomers+2*n_microbes+n_enzymes+n_microbes] = rG_CO2 + rM_CO2 + rX_CO2 + J_DE_CO2

mass_balance = @. dD + dE + dV + dX + dCO2
@test mass_balance ≈ [0.0] atol = 1e-10
########################################

########################################
# rhs function
p = Params(p_set, p_met, p_ass, nothing, p_turn)
du = zeros(p_set.dim)
du = DEBmicroTrait.batch_model!(du, u, p, 0.0)
@test sum(du) ≈ mass_balance[1] atol = 1e-10
########################################
