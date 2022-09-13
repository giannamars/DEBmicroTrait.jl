using DEBmicroTrait, Test

n_polymers = 0
n_monomers = rand(1:100)
n_microbes = rand(1:1000)
n_enzymes  = rand(1:100)
n_minerals = 0
p_set      = Setup(n_polymers, n_monomers, n_microbes, n_enzymes, n_minerals)
@test p_set.dim == n_polymers + n_monomers + 3*n_microbes + n_enzymes
########################################

########################################
# split_state_batch
u     = ones(p_set.dim)
D    = u[1+n_polymers:n_polymers+n_monomers]
E    = u[1+n_polymers+n_monomers:n_polymers+n_monomers+n_microbes]
V    = u[1+n_polymers+n_monomers+n_microbes:n_polymers+n_monomers+2*n_microbes]
X    = u[1+n_polymers+n_monomers+2*n_microbes:n_polymers+n_monomers+2*n_microbes+n_enzymes]
CO2  = u[1+n_polymers+n_monomers+2*n_microbes+n_enzymes:n_polymers+n_monomers+2*n_microbes+n_enzymes+n_microbes]
@test size(D,1) == n_monomers
@test size(E,1) == n_microbes
@test size(V,1) == n_microbes
@test size(X,1) == n_enzymes
@test size(CO2,1) == n_microbes
