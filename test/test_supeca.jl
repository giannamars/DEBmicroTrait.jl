using DEBmicroTrait, Test

n_A = rand(1:100)
n_B = rand(1:100)
n_E = rand(1:100)

F_c_A          = zeros(n_E)
F_c_B          = zeros(n_E)
F_r_A          = zeros(n_A)
F_r_B          = zeros(n_B)
SUPECA         = zeros(n_A,n_B,n_E)

A              = rand(n_A)
B              = rand(n_B)
E              = rand(n_E)

K              = rand(n_A,n_B,n_E)
N_SB           = rand(n_A,n_B,n_E)
k2p            = rand(n_A,n_B,n_E)

F_c_A = DEBmicroTrait.normalized_substrate_flux!(F_c_A, A, K[:,1,:])
@test size(F_c_A) == (n_E,)
F_c_B = DEBmicroTrait.normalized_substrate_flux!(F_c_B, B, K[1,:,:])
@test size(F_c_B) == (n_E,)
F_r_A = DEBmicroTrait.conjugate_substrate_flux!(F_r_A, E, K[:,1,:], N_SB[:,1,:])
@test size(F_r_A) == (n_A,)
F_r_B = DEBmicroTrait.conjugate_substrate_flux!(F_r_B, E, K[1,:,:], N_SB[1,:,:])
@test size(F_r_B) == (n_B,)
SUPECA = DEBmicroTrait.SUPECA_kinetics!(SUPECA, A, B, E, K, k2p, N_SB)
@test size(SUPECA) == (n_A,n_B,n_E)
