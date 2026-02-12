using DEBmicroTrait, Test, Distributions

n_monomers = 2
n_polymers = 5
n_enzymes  = 6
P          = rand(n_polymers)
X          = rand(n_enzymes)
#
Χ_0        = [10,10,10,10,10]
ρ_P        = [1.35, 1.5, 1.4, 0.4, 1.0]
R_E        = [3.29e-9, 4e-9, 4e-9, 4.4e-9, 5e-9,3.9e-9]
V_E        = [5, 0.24, 4,0.19, 3,4]
MW_E       = [48e3, 65e3, 59.7e3, 70e3, 50e3, 52e3]
D_E        = DEBmicroTrait.aqueous_diffusivity(MW_E)
#
R_P        = 20e-3*ones(n_polymers)
#
α_kin_P    = DEBmicroTrait.binding_sites(ρ_P, R_E)
K_EP_0     = DEBmicroTrait.polymer_affinity(R_E, R_P, V_E, D_E)
#
f_PD       = zeros(n_polymers, n_monomers)
for i in 1:n_polymers
    f_PD[i,:] = rand(Dirichlet(n_monomers,1))
end
# 

p = Depolymerization(Χ_0, V_E, α_kin_P, K_EP_0, f_PD)
J_P = zeros(n_polymers)
J_P, J_PD = DEBmicroTrait.depolymerization!(J_P, p, P, X)

@test size(J_P,1) == n_polymers
@test size(J_PD,1) == n_monomers
@test sum(J_PD) ≈ sum(J_P) atol=1e-8
#######################################################################
#                     test paper values
# enzyme characteristics
R_E        = [3.29e-9]
V_E        = [5.0]
MW_E       = [48e3]
#D_E        = DEBmicroTrait.aqueous_diffusivity(MW_E)
D_E        = [1e-10]
# Polymer characteristics
ρ_P        = [1.35]
R_P        = [1].*1e-3
#
K_EP_0     = DEBmicroTrait.polymer_affinity(R_E, R_P, V_E, D_E)
@test K_EP_0 ≈ [0.244] atol=0.001
#
k_max_RP   = DEBmicroTrait.max_hydrolyis_rate(R_E, R_P, V_E)
#
R_P        = [1, 10, 20, 30, 40, 50].*1e-3
K_EP_0     = DEBmicroTrait.polymer_affinity(R_E, R_P, V_E, D_E)