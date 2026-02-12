using DEBmicroTrait, Test, Distributions

########################################
n_polymers = rand(1:500)
n_monomers = rand(1:500)
n_microbes = rand(1:500)
n_enzymes  = rand(1:500)
n_minerals = 0
p_set      = Setup(n_polymers, n_monomers, n_microbes, n_enzymes, n_minerals)
########################################

########################################
# split_state_batch
u     = ones(p_set.dim)
n_polymers = p_set.n_polymers
n_monomers = p_set.n_monomers
n_microbes = p_set.n_microbes
n_enzymes  = p_set.n_enzymes

P    = u[1:n_polymers]
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
α_X     = 0.01*ones(n_microbes)
y_EX    = 1.0*ones(n_microbes)
f_αX    = rand(Dirichlet(n_enzymes,1))
mingt   = 1.0*ones(n_microbes)

p_met   = MetabolismC(k_E, y_EV, k_M, y_EM, α_X, y_EX, f_αX, mingt)
r                         = DEBmicroTrait.growth!(0.0*ones(1), p_met, E, V)
x, rG_CO2, rM_CO2, rX_CO2 = DEBmicroTrait.growth_production!(r, p_met, E, V)
J_EX                      = DEBmicroTrait.enzyme_production!(x, p_met, V)
########################################

########################################
# assimilation
N_SB         = rand(n_monomers, n_microbes)
K_D          = rand(n_monomers, n_microbes)
y_DE         = rand(n_monomers, n_microbes)
N_C          = rand(n_monomers)

p_ass        = AssimilationC(N_SB, K_D, y_DE, N_C)

J_DE         = zeros(n_microbes)
J_DE_CO2     = zeros(n_microbes)
J_D          = zeros(n_monomers)
J_DE         = DEBmicroTrait.assimilation!(J_DE, p_ass, D, V)
J_DE_CO2     = DEBmicroTrait.assimilation_production!(J_DE_CO2, p_ass, D, V)
J_D          = DEBmicroTrait.uptake!(J_D, p_ass, D, V)
########################################

# depolymerization
Χ_0          = rand(n_polymers)
V_E          = rand(n_polymers)
α_kin_P      = rand(n_polymers, n_enzymes)
K_EP_P       = rand(n_polymers, n_enzymes)

f_PD = zeros(n_polymers, n_monomers)
for i in 1:n_polymers
    f_PD[i,:] = rand(Dirichlet(n_monomers,1))
end
p_depoly     = Depolymerization(Χ_0, V_E, α_kin_P, K_EP_P, f_PD)
J_P, J_PD   = DEBmicroTrait.depolymerization!(zeros(n_polymers), p_depoly, P, X)

########################################
# turnover
γ_V_0     = ones(n_microbes)
γ_V_1     = ones(n_microbes)
γ_X       = ones(n_enzymes)
f_ED      = rand(Dirichlet(n_monomers,1))
f_VD      = rand(Dirichlet(n_monomers,1))
f_VP      = rand(Dirichlet(n_polymers,1))
f_V       = ones(n_microbes)
f_XD      = rand(Dirichlet(n_monomers,1))
f_XP      = rand(Dirichlet(n_polymers,1))
f_X       = ones(n_enzymes)

p_turn       = Turnover(γ_V_0,γ_V_1,γ_X,f_ED,f_VD,f_VP,f_V,f_XD,f_XP,f_X)
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

dP = @. du[1:n_polymers] = - J_P + J_VP + J_XP
dD = @. du[1+n_polymers:n_polymers+n_monomers] = J_PD - J_D + J_ED + J_VD + J_XD
dE = @. du[1+n_polymers+n_monomers:n_polymers+n_monomers+n_microbes] =  J_DE - (p_met.k_E - r)*E - J_E
dV = @. du[1+n_polymers+n_monomers+n_microbes:n_polymers+n_monomers+2*n_microbes] = r*V - J_V
dX = @. du[1+n_polymers+n_monomers+2*n_microbes:n_polymers+n_monomers+2*n_microbes+n_enzymes] = J_EX - J_X
dCO2 = @. du[1+n_polymers+n_monomers+2*n_microbes+n_enzymes:n_polymers+n_monomers+2*n_microbes+n_enzymes+n_microbes] = rG_CO2 + rM_CO2 + rX_CO2 + J_DE_CO2

mass_balance = sum(dP) + sum(dD) + sum(dE) + sum(dV) + sum(dX) + sum(dCO2)
@test mass_balance ≈ 0.0 atol = 1e-6
########################################

########################################
# rhs function
p = Params(p_set, p_met, p_ass, p_depoly, p_turn)
du = zeros(p_set.dim)
du = DEBmicroTrait.batch_model!(du, u, p, 0.0)
@test sum(du) ≈ mass_balance[1] atol = 1e-6
########################################