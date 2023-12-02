using DEBmicroTrait 

n_monomers = 1
n_microbes = 1
n_minerals = 1

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

df_metabolites          = CSV.read(joinpath(dir, "files/input/root_exudates.csv"), DataFrame, missingstring="N/A");

dG_carboxyl = -5.78
dG_methyl   = -0.88
dG_amino    = -2.12
dG_phoph    = -2.71

dG_sorption            = df_metabolites.Hydroxyl