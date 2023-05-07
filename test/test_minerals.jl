using DEBmicroTrait 

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

