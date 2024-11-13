using DEBmicroTrait

# glucose
molecular_weight = 180.156*ones(1)
D_S = DEBmicroTrait.aqueous_diffusivity(molecular_weight)

# cell
V_c = 1e-18*ones(1)
ρ_p = 0.001*ones(1,1)
N_cell = 1*ones(1)

# Angelo
pct_sand        = 28.0
pct_clay        = 27.0
theta_grav      = 0.23
rho_bulk        = 1.66 - 0.308*(1.66)^0.5
s_sat           = theta_grav*rho_bulk
ϕ_w, τ_g, τ_w, δ = DEBmicroTrait.soil_affinity_properties(pct_sand, pct_clay, s_sat)
D_eff = DEBmicroTrait.effective_diffusivity(molecular_weight, pct_sand, pct_clay, s_sat)
K_D = DEBmicroTrait.specific_affinity(V_c, ρ_p, D_S, pct_sand, pct_clay, s_sat, N_cell)

# Hopland
pct_sand        = 45.0
pct_clay        = 19.0
theta_grav      = 0.17
rho_bulk        = 1.66 - 0.308*(1.39)^0.5
s_sat           = theta_grav*rho_bulk
ϕ_w, τ_g, τ_w, δ = DEBmicroTrait.soil_affinity_properties(pct_sand, pct_clay, s_sat)
D_eff = DEBmicroTrait.effective_diffusivity(molecular_weight, pct_sand, pct_clay, s_sat)
K_D = DEBmicroTrait.specific_affinity(V_c, ρ_p, D_S, pct_sand, pct_clay, s_sat, N_cell)

# Sedgwick
pct_sand        = 38.0
pct_clay        = 34.0
theta_grav      = 0.25
rho_bulk        = 1.66 - 0.308*(1.99)^0.5
s_sat           = theta_grav*rho_bulk
ϕ_w, τ_g, τ_w, δ = DEBmicroTrait.soil_affinity_properties(pct_sand, pct_clay, s_sat)
D_eff = DEBmicroTrait.effective_diffusivity(molecular_weight, pct_sand, pct_clay, s_sat)
K_D = DEBmicroTrait.specific_affinity(V_c, ρ_p, D_S, pct_sand, pct_clay, s_sat, N_cell)