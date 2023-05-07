using DEBmicroTrait
using CSV, DataFrames, Statistics
using JLD
using DifferentialEquations

dir                     = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()
df_isolates             = CSV.read(joinpath(dir, "files/input/isolates2traits.csv"), DataFrame, missingstring="N/A")
df_metabolites          = CSV.read(joinpath(dir, "files/input/root_exudates.csv"), DataFrame, missingstring="N/A")

assimilation            = load(joinpath(dir, "files/output//isolates_assimilation.jld"))
enzymes                 = load(joinpath(dir, "files/output/isolates_enzymes.jld"))
maintenance             = load(joinpath(dir, "files/output/isolates_maintenance.jld"))
protein_synthesis       = load(joinpath(dir, "files/output/isolates_protein_synthesis.jld"))
turnover                = load(joinpath(dir, "files/output/isolates_turnover.jld"))
initb                   = load(joinpath(dir, "files/output/isolates_batch_init.jld"))


condition(u,t,integrator) = u[1] - 1e-5
affect!(integrator)       = terminate!(integrator)
cb                        = ContinuousCallback(condition,affect!)

BGE_tseries       = zeros(39, 83, 500) # BGE = BP/(BP + BR)
r_tseries         = zeros(39, 83, 500) # realized growth rate

for i in 1:39
    for j in 1:83
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

        tspan             = (0.0,500.0)
        prob              = ODEProblem(DEBmicroTrait.batch_model!,u0,tspan,p)
        sol               = solve(prob, alg_hints=[:stiff], callback=cb)

        du   = zeros(p.setup_pars.dim)
        BR   = [DEBmicroTrait.batch_model!(du, sol.u[i], p, 0)[end] for i in 1:size(sol.t,1)]
        BP   = [DEBmicroTrait.batch_model!(du, sol.u[i], p, 0)[2] + DEBmicroTrait.batch_model!(du, sol.u[i], p, 0)[3] for i in 1:size(sol.t,1)]
        BGE  = @. BP/(BP + BR)
        for k in 1:length(sol.t)
            BGE_tseries[i,j,k] = BGE[k]
        end

        r    = [DEBmicroTrait.growth!(0.0*ones(1), p.metabolism_pars, [sol[i][2]], [sol[i][3]])[1] for i in 1:size(sol.t,1)]
        for k in 1:length(sol.t)
            r_tseries[i,j,k] = r[k]
        end

    end
end

BGE_tseries[BGE_tseries.>1.0].=0
BGE_tseries[BGE_tseries.<=0.0].=0
BGE_median = zeros(39,83)
for i in 1:39
    for j in 1:83
        try
            BGE_median[i,j]  = median(filter(!iszero, BGE_tseries[i,j,:]))
        catch
            BGE_median[i,j]  = NaN
        end
    end
end

