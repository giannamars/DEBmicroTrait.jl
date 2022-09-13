using DEBmicroTrait, Test

n_substrates = rand(1:100)
n_consumers  = rand(1:1000)

N_SB         = rand(n_substrates, n_consumers)
K_D          = rand(n_substrates, n_consumers)
y_DE         = rand(n_substrates, n_consumers)
N_C          = rand(n_substrates)
N_X          = rand(n_substrates)

p            = AssimilationC(N_SB, K_D, y_DE, N_C, N_X)

D            = rand(n_substrates)
V            = rand(n_consumers)
J_DE         = zeros(n_consumers)
J_DE_CO2     = zeros(n_consumers)
J_D          = zeros(n_substrates)
J_DE         = DEBmicroTrait.assimilation!(J_DE, p, D, V)
J_DE_CO2     = DEBmicroTrait.assimilation_production!(J_DE_CO2, p, D, V)
J_D          = DEBmicroTrait.uptake!(J_D, p, D, V)

@test size(J_DE,1) == n_consumers
@test size(J_DE_CO2,1) == n_consumers
@test size(J_D,1) == n_substrates
