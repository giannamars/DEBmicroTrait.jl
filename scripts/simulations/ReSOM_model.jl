using DEBmicroTrait
using CSV, DataFrames, Statistics 
using JLD
using OrdinaryDiffEq
using Plots

dir                     = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()
initb                   = load(joinpath(dir, "files/output/isolates_batch_init.jld"))

n_polymers = 1
n_monomers = 1
n_microbes = 1
n_enzymes  = 10
n_minerals = 1
p                       = DEBmicroTrait.init_ReSOM_model(n_polymers, n_monomers, n_microbes, n_enzymes, n_minerals)

u0                      = zeros(p.setup_pars.dim)
u0[1:n_polymers]                                                          .= 10.0
u0[1+n_polymers+n_monomers:n_polymers+n_monomers+n_microbes]              .= 0.9*8.2e-5
u0[1+n_polymers+n_monomers+n_microbes:n_polymers+n_monomers+2*n_microbes] .= 0.1*8.2e-5
u0[1+n_polymers:n_polymers+n_monomers]                                    .= 1.25

tspan             = (0.0,500.0)
prob              = ODEProblem(DEBmicroTrait.ReSOM_model!,u0,tspan,p)
sol               = solve(prob, Tsit5())

plot(sol, idxs=2)

du = zeros(p.setup_pars.dim)
du = [sum(DEBmicroTrait.ReSOM_model!(du, sol.u[i], p, 0.0)) for i in 1:size(sol.t,1)]

plot(du)