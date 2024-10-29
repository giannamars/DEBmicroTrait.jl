using DEBmicroTrait
using CSV, DataFrames, Statistics
using JLD
using Plots
using DifferentialEquations

dir             = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()
df_isolates     = CSV.read(joinpath(dir, "files/input/dement_isolates.csv"), DataFrame, missingstring="")
dropmissing!(df_isolates)

# metabolite traits
Formula  = "C6H12O6"
N_C      = [6]
Molecular_weight = [180.15588]

assimilation            = load(joinpath(dir, "files/output/dement_isolates_assimilation.jld"))
enzymes                 = load(joinpath(dir, "files/output/dement_isolates_enzymes.jld"))
maintenance             = load(joinpath(dir, "files/output/dement_isolates_maintenance.jld"))
protein_synthesis       = load(joinpath(dir, "files/output/dement_isolates_protein_synthesis.jld"))
turnover                = load(joinpath(dir, "files/output/dement_isolates_turnover.jld"))
initb                   = load(joinpath(dir, "files/output/dement_isolates_batch_init.jld"))

condition(u,t,integrator) = u[1] - 1e-5
affect!(integrator)       = terminate!(integrator)
cb                        = ContinuousCallback(condition,affect!)

BGE_tseries       = zeros(size(df_isolates,1), 1, 500)

for i in 1:size(df_isolates,1)
    for j in 1:1
        id_isolate = i
        id_monomer = j
        p                 = DEBmicroTrait.init_batch_model(id_isolate, id_monomer, assimilation, enzymes, maintenance, protein_synthesis, turnover)
        n_polymers        = p.setup_pars.n_polymers
        n_monomers        = p.setup_pars.n_monomers
        n_microbes        = p.setup_pars.n_microbes

        u0                                                                         = zeros(p.setup_pars.dim)
        u0[1+n_polymers+n_monomers:n_polymers+n_monomers+n_microbes]              .= 0.9*initb["Bio0"][id_isolate]
        u0[1+n_polymers+n_monomers+n_microbes:n_polymers+n_monomers+2*n_microbes] .= 0.1*initb["Bio0"][id_isolate]
        u0[1+n_polymers:n_polymers+n_monomers]                                    .= 1.25

        tspan             = (0.0,1000.0)
        prob              = ODEProblem(DEBmicroTrait.batch_model!,u0,tspan,p)
        try
            sol               = solve(prob, alg_hints=[:stiff], callback=cb)

            du   = zeros(p.setup_pars.dim)
            BR   = [DEBmicroTrait.batch_model!(du, sol.u[i], p, 0)[end] for i in 1:size(sol.t,1)]
            BP   = [DEBmicroTrait.batch_model!(du, sol.u[i], p, 0)[2] + DEBmicroTrait.batch_model!(du, sol.u[i], p, 0)[3] for i in 1:size(sol.t,1)]
            BGE  = @. BP/(BP + BR)
            for k in 1:length(sol.t)
                BGE_tseries[i,j,k] = BGE[k]
            end
        catch e 
            println("Error occurred: ", e)
            for k in 1:500
                BGE_tseries[i,j,k] = NaN
            end
        end
    end
end

# output BGE-growth

BGE_tseries[BGE_tseries.>1.0].=0
BGE_tseries[BGE_tseries.<=0.0].=0

BGE_median = zeros(size(df_isolates,1), 1)
for i in 1:size(df_isolates,1)
    for j in 1:1
        try
            BGE_median[i,j]  = median(filter(!iszero, BGE_tseries[i,j,:]))
        catch
            BGE_median[i,j]  = NaN
        end
    end
end

nan_count = count(isnan, BGE_median)

df_isolates[!, :CUE] .= BGE_median
CSV.write(joinpath(dir, "files/output/dement_isolates_CUE.csv"), df_isolates)

using GLM
df_CUE = dropmissing(df_isolates)

linearRegressor = lm(@formula(log(CUE) ~ log(genome_length)), dropmissing(df_CUE))
r2(linearRegressor)