using DEBmicroTrait, Roots, Test

## MetabolismC

k_E     = 0.2*ones(1)
y_EV    = 1.0*ones(1)
k_M     = 0.1*ones(1)
y_EM    = 1.0*ones(1)
α_X     = 0.0*ones(1)
y_EX    = 1.0*ones(1)
f_αX    = 1.0*ones(1)
mingt   = 1.0*ones(1)

p       = MetabolismC(k_E, y_EV, k_M, y_EM, α_X, y_EX, f_αX,  mingt)
E       = 1.0*ones(1)
V       = 1.0*ones(1)

r       = DEBmicroTrait.growth!(0.0*ones(1), p, E, V)
r_analytic = @. (k_E*E/V - k_M*y_EM)/(E/V + y_EV)
@test r ≈ r_analytic atol=1e-3
x, rG_CO2, rM_CO2, rX_CO2 = DEBmicroTrait.growth_production!(r, p, E, V)
@test x == [0.0]
@test rG_CO2 == [0.0]
@test rX_CO2 == [0.0]


n_consumers = rand(1:1000)
k_E     = 0.2*ones(n_consumers)
y_EV    = 1.0*ones(n_consumers)
k_M     = 0.1*ones(n_consumers)
y_EM    = 1.0*ones(n_consumers)
α_X     = 0.0*ones(n_consumers)
f_αX     = 1.0*ones(n_consumers)
y_EX    = 1.0*ones(n_consumers)
mingt   = 1.0*ones(n_consumers)

p       = MetabolismC(k_E, y_EV, k_M, y_EM, α_X, y_EX, f_αX, mingt)
E       = 1.0*ones(n_consumers)
V       = 1.0*ones(n_consumers)

r       = DEBmicroTrait.growth!(0.0*ones(n_consumers), p, E, V)
@test size(r,1) == n_consumers
