using DEBmicroTrait
using CSV, DataFrames, Statistics 
using JLD
using DifferentialEquations
using Plots

dir                     = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()

# Load media composition: Formula, Molecular weight, Medium concentration
df_metabolites = CSV.read(joinpath(dir, "files/input/1_10_R2A_medium.csv"), DataFrame, missingstring="N/A")

# Load isolate parameterization, note: assimilation parameterization depends on media composition
assimilation            = load(joinpath(dir, "files/output/isolates_assimilation_10_R2A.jld"))
enzymes                 = load(joinpath(dir, "files/output/isolates_enzymes.jld"))
maintenance             = load(joinpath(dir, "files/output/isolates_maintenance.jld"))
protein_synthesis       = load(joinpath(dir, "files/output/isolates_protein_synthesis.jld"))
turnover                = load(joinpath(dir, "files/output/isolates_turnover.jld"))
initb                   = load(joinpath(dir, "files/output/isolates_batch_init.jld"))
id_isolate = 30 # HB15
n_monomers = 43

p                 = DEBmicroTrait.init_mixed_medium(id_isolate, n_monomers, assimilation, enzymes, maintenance, protein_synthesis, turnover)
n_polymers        = p.setup_pars.n_polymers
n_monomers        = p.setup_pars.n_monomers
n_microbes        = p.setup_pars.n_microbes

u0                                                                         = zeros(p.setup_pars.dim)
u0[1+n_polymers+n_monomers:n_polymers+n_monomers+n_microbes]              .= 0.9*initb["Bio0"][id_isolate]
u0[1+n_polymers+n_monomers+n_microbes:n_polymers+n_monomers+2*n_microbes] .= 0.1*initb["Bio0"][id_isolate]
u0[1+n_polymers:n_polymers+n_monomers]                                    .= df_metabolites.Concentration


tspan             = (0.0,172.0)
prob              = ODEProblem(DEBmicroTrait.batch_model!,u0,tspan,p)
sol               = solve(prob, alg_hints=[:stiff])

E = sol[44,:]
V = sol[45,:]
Biomass = E .+ V
plot(sol.t, Biomass)

r = [DEBmicroTrait.growth!(0.0*ones(1), p.metabolism_pars, [sol[i][44]], [sol[i][45]])[1] for i in 1:size(sol.t,1)]

y_VM  = p.metabolism_pars.y_EM./p.metabolism_pars.y_EV
m_E   = E./V
j_EC  = m_E.*(p.metabolism_pars.k_E .- r)
j_EM  = p.metabolism_pars.k_M.*p.metabolism_pars.y_EM
j_EM  = 0.1
jEM   = min.(j_EC, j_EM)
jVM   = (j_EM .- jEM).*y_VM./p.metabolism_pars.y_EM

plot(sol.t, jEM)
plot!(sol.t, jVM)
plot!(sol.t, jVM)
plot!(sol.t, E)
plot!(sol.t, sol[1,:])