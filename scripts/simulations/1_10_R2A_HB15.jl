using DEBmicroTrait
using CSV, DataFrames, Statistics
using JLD
using OrdinaryDiffEq
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
sol               = solve(prob, Tsit5())

plot(sol, idxs=1:43)



