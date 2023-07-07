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

n_minerals   = rand(1:10)
M            = rand(n_minerals)
k2p_D        = 180.0*60*60*ones(n_substrates, n_consumers)
k2p_M        = ones(n_substrates, n_minerals)
k2p          = hcat(k2p_D, k2p_M)
K_D_0        = rand(n_substrates, n_consumers)
K_D_M        = rand(n_substrates, n_minerals)
K_D          = hcat(K_D_0, K_D_M)
N_SB_0       = rand(n_substrates, n_consumers)
N_SB_M       = ones(n_substrates, n_minerals)
N_SB         = hcat(N_SB_0, N_SB_M)
p            = AssimilationCM(N_SB, K_D, y_DE, N_C, N_X, M)

J_DE         = zeros(n_substrates, n_consumers)
J_DE         = DEBmicroTrait.assimilation!(J_DE, p, D, V)
J_DE_CO2     = zeros(n_substrates, n_consumers)
J_DE_CO2     = DEBmicroTrait.assimilation_production!(J_DE_CO2, p, D, V)
J_D          = zeros(n_substrates, n_consumers)
J_D, J_DM    = DEBmicroTrait.uptake!(J_D, p, D, V)

J_D
J_DM

n_substrates = size(D,1)
n_consumers  = size(V,1)
n_minerals   = size(M,1)
k2p_D        = 180.0*60*60*ones(n_substrates, n_consumers)
k2p_M        = ones(n_substrates, n_minerals)
k2p          = hcat(k2p_D, k2p_M)
K_D_0        = rand(n_substrates, n_consumers)
K_D_M        = rand(n_substrates, n_minerals)
K_D          = hcat(K_D_0, K_D_M)
N_SB_0       = rand(n_substrates, n_consumers)
N_SB_M       = ones(n_substrates, n_minerals)
N_SB         = hcat(N_SB_0, N_SB_M)
ECA          = DEBmicroTrait.ECA_kinetics!(zeros(n_substrates,n_consumers+n_minerals), D, vcat(V, M), K_D, k2p, N_SB)

J_DE = vcat(sum((1.0 .- p.y_DE).*p.N_C.*ECA[:, 1:n_consumers], dims=1)...)