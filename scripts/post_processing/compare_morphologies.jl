using DEBmicroTrait
V_c   = [1e-18]
N_p   = 10000*ones(1,1)
rho_p = DEBmicroTrait.transporter_density(V_c, N_p)
#rho_p = 0.01*ones(1,1)
D_S   = [1e-10]
K_0   = DEBmicroTrait.specific_reference_affinity(V_c, rho_p, D_S)
min_gen_time = [4.0]
N_SB  = DEBmicroTrait.transporter_density_to_monomer_uptake_sites(V_c, rho_p, min_gen_time, ["+"])
Vmax  = 180.0*60^2 * N_SB

affinity = Vmax./K_0