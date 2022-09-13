using DEBmicroTrait, Test

n_microbes = rand(1:1000)
n_enzymes  = rand(1:100)
n_substrates = rand(1:100)

γ_V_0     = ones(n_microbes)
γ_V_1     = ones(n_microbes)
γ_X       = ones(n_enzymes)
γ_D_ads   = zeros(n_substrates)
γ_X_ads   = zeros(n_enzymes)
f_ED      = ones(n_substrates)
f_VD      = ones(n_microbes)
f_VP      = ones(n_microbes)
f_V       = ones(n_microbes)
f_XD      = ones(n_enzymes)
f_XP      = ones(n_enzymes)
f_X       = zeros(n_enzymes)

p     = Turnover(γ_V_0,γ_V_1,γ_X,γ_D_ads,γ_X_ads,f_ED,f_VD,f_VP,f_V,f_XD,f_XP,f_X)

B     = rand(n_microbes)
J_B   = zeros(n_microbes)
J_B   = DEBmicroTrait.biomass_turnover!(J_B, p, B)
@test size(J_B,1) == n_microbes

X     = rand(n_enzymes)
J_X   = zeros(n_enzymes)
J_X   = DEBmicroTrait.enzyme_decay!(J_X, p, X)
@test size(J_X,1) == n_enzymes

E     = rand(n_microbes)
J_ED  = rand(n_microbes)
J_ED  = DEBmicroTrait.reserve_recycling!(J_ED, p, E)
@test size(J_ED,1) == n_substrates

X     = rand(n_enzymes)
J_XD  = rand(n_substrates)
J_XD  = DEBmicroTrait.enzyme_recycling!(J_XD, p, X)
@test size(J_ED,1) == n_substrates
