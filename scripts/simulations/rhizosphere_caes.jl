using DEBmicroTrait
using CSV, DataFrames, Statistics
using JLD, JLD2
using Plots
gr()
using DifferentialEquations


dir                     = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()
df_mags                 = CSV.read(joinpath(dir, "files/input/greenlon-H-mags2traits.csv"), DataFrame, missingstring="")
df_metabolites          = CSV.read(joinpath(dir, "files/input/greenlon-H-monomers.csv"), DataFrame, missingstring="")

assimilation            = load(joinpath(dir, "files/output/mags_assimilation.jld"))
enzymes                 = load(joinpath(dir, "files/output/mags_enzymes.jld"))
maintenance             = load(joinpath(dir, "files/output/mags_maintenance.jld"))
protein_synthesis       = load(joinpath(dir, "files/output/mags_protein_synthesis.jld"))
turnover                = load(joinpath(dir, "files/output/mags_turnover.jld"))
initb                   = load(joinpath(dir, "files/output/mags_batch_init.jld"))

condition(u,t,integrator) = u[1] - 1e-5
affect!(integrator)       = terminate!(integrator)
cb                        = ContinuousCallback(condition,affect!)


id_polymer = 5 #[cellulose, chitin, mannan, xylan, xyloglucan, glucan, protein]
id_monomer = 1


r_tseries           = zeros(size(df_mags,1), 1000)
Bio_tseries         = zeros(size(df_mags,1), 1000)
BGE_tseries         = zeros(size(df_mags,1), 1000)
J_DE_tseries        = zeros(size(df_mags,1), 1000)

for i in 1:size(df_mags,1)
    id_isolate = i
    p                 = DEBmicroTrait.init_batch_model(id_isolate, id_monomer, id_polymer, assimilation, enzymes, maintenance, protein_synthesis, turnover)
    n_polymers        = p.setup_pars.n_polymers
    n_monomers        = p.setup_pars.n_monomers
    n_microbes        = p.setup_pars.n_microbes
    n_enzymes         = p.setup_pars.n_enzymes

    u0                                                                         = zeros(p.setup_pars.dim)
    u0[1+n_polymers+n_monomers:n_polymers+n_monomers+n_microbes]              .= 0.9*initb["Bio0"][id_isolate]
    u0[1+n_polymers+n_monomers+n_microbes:n_polymers+n_monomers+2*n_microbes] .= 0.1*initb["Bio0"][id_isolate]
    u0[1+n_polymers:n_polymers+n_monomers]                                    .= 1.25
    u0[1:n_polymers]                                                          .= 10.0
    u0[1+n_polymers+n_monomers+2*n_microbes:n_polymers+n_monomers+2*n_microbes+n_enzymes] .= 0.0

    tspan             = (0.0,1000.0)
    prob              = ODEProblem(DEBmicroTrait.batch_model!,u0,tspan,p)
    sol               = solve(prob, alg_hints=[:stiff], callback=cb)

    r    = [DEBmicroTrait.growth!(0.0*ones(1), p.metabolism_pars, [sol[i][3]], [sol[i][4]])[1] for i in 1:size(sol.t,1)]
    for k in 1:length(sol.t)
        r_tseries[i,k] = r[k]
    end

    for k in 1:length(sol.t)
        Bio = sol[k][3].+sol[k][4]
        Bio_tseries[i,k] = Bio[1]
    end

    du   = zeros(p.setup_pars.dim)
    BR   = [DEBmicroTrait.batch_model!(du, sol.u[i], p, 0)[end] for i in 1:size(sol.t,1)]
    BP   = [DEBmicroTrait.batch_model!(du, sol.u[i], p, 0)[3] + DEBmicroTrait.batch_model!(du, sol.u[i], p, 0)[4] for i in 1:size(sol.t,1)]
    BGE  = @. BP/(BP + BR)
    for k in 1:length(sol.t)
        BGE_tseries[i,k] = BGE[k]
    end

    for k in 1:length(sol.t)
    J_DE  = DEBmicroTrait.assimilation!(zeros(1), p.assimilation_pars, [sol[k][2]], [sol[k][4]])
    J_DE_tseries[i,k] = J_DE[1]
    end
end

r_median        = zeros(size(r_tseries,1))
Bio_median      = zeros(size(Bio_tseries,1))
BGE_median      = zeros(size(BGE_tseries,1))
J_DE_median     = zeros(size(J_DE_tseries,1))

for i in 1:size(r_tseries,1)
    try
        r_median[i]  = median(filter(!iszero, r_tseries[i,:]))
    catch
        r_median[i,]  = NaN
    end
    try
        Bio_median[i]  = median(filter(x -> x != 0 && x < 10, Bio_tseries[i,:]))

    catch
        Bio_median[i]  = NaN
    end
    try
        BGE_median[i]  = median(filter(x -> x > 0 && x <= 1, BGE_tseries[i,:]))

    catch
        BGE_median[i]  = NaN
    end
        try
        J_DE_median[i]  = median(filter(!iszero, J_DE_tseries[i,:]))
    catch
        J_DE_median[i,]  = NaN
    end
end
 

df_mags.r_median            = r_median
df_mags.relabund_median     = Bio_median./maximum(Bio_median)
df_mags.bge_median          = BGE_median
df_mags.J_DE_median         = J_DE_median

CSV.write(joinpath(dir, "files/output/greenlon-H-xyloglucan.csv"), df_mags)

