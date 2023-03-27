using DEBmicroTrait

# Interception probability
N_p = 3140.0*ones(1,1)
####################################################
# Shape = sphere (r_c)
r_c   = [1e-6]
V_c_sphere   = @. (4/3)*π*r_c^3
p_int_sphere = DEBmicroTrait.interception_probability(V_c_sphere, N_p)

####################################################
# Shape = spheroid (b~l)
l_c   = [0.5e-6]
b_c   = [0.25e-6]
V_c_spheroid   = @. (4/3)*b_c^2*l_c
p_int_spheroid = DEBmicroTrait.interception_probability(l_c, b_c, N_p)
####################################################

####################################################
# Shape = stalk (b<<<l)
l_c_stalk   = [0.5e-6]
b_c_stalk   = [5e-10]
V_c_stalk   = @. (4/3)*b_c^2*l_c
p_int_stalk = DEBmicroTrait.interception_probability(l_c_stalk, b_c_stalk, N_p,"stalk")
####################################################

# Half-saturation constant
D_S   = [1e-10]
####################################################
# Shape = sphere (r_c)
rho_p = DEBmicroTrait.transporter_density(V_c_sphere, N_p)
K_0_sphere   = DEBmicroTrait.specific_reference_affinity(V_c_sphere, rho_p, D_S)
####################################################
# Shape = spheroid (b~l)
K_0_spheroid   = DEBmicroTrait.specific_reference_affinity(l_c, b_c, N_p, D_S)
####################################################
# Shape = stalk (b<<l)
K_0_stalk   = DEBmicroTrait.specific_reference_affinity(l_c_stalk, b_c_stalk, N_p, D_S, "stalk")

# Max. specific uptake rate
min_gen_time = [4.0]
gram_stain = ["+"]
####################################################
# Shape = sphere (r_c)
N_SB_sphere  = DEBmicroTrait.transporter_density_to_monomer_uptake_sites(V_c_sphere, rho_p, min_gen_time, gram_stain)
Vmax_sphere  = 180.0*60^2 * N_SB_sphere
####################################################
# Shape = spheroid (b~l)
N_SB_spheroid  = DEBmicroTrait.transporter_density_to_monomer_uptake_sites(V_c_spheroid, N_p, min_gen_time, gram_stain, "spheroid")
Vmax_spheroid  = 180.0*60^2 * N_SB_spheroid
####################################################
# Shape = stalk (b<<l)
N_SB_stalk = DEBmicroTrait.transporter_density_to_monomer_uptake_sites(V_c_stalk, N_p, min_gen_time, gram_stain, "spheroid")
Vmax_stalk  = 180.0*60^2 * N_SB_stalk

# Affinity constant (S->0)
####################################################
# Shape = sphere (r_c)
affinity_sphere = Vmax_sphere./K_0_sphere
####################################################
# Shape = spheroid (b~l)
affinity_spheroid = Vmax_spheroid./K_0_spheroid
####################################################
# Shape = stalk (b<<l)
affinity_stalk = Vmax_stalk./K_0_stalk
