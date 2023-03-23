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

# different morphologies

V_c   = [1e-20]
N_p   = 10000*ones(1,1)
rho_p = DEBmicroTrait.transporter_density(V_c, N_p)
#rho_p = 0.01*ones(1,1)
D_S   = [1e-10]
K_0   = DEBmicroTrait.specific_reference_affinity(V_c, rho_p, D_S)
Vmax  = 180.0 * N_p

affinity = Vmax./K_0