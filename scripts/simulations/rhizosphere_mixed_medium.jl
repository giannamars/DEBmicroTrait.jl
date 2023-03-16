using DEBmicroTrait
using CSV, DataFrames, Statistics, Tables
using JLD
using OrdinaryDiffEq

dir                     = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()
df_isolates             = CSV.read(joinpath(dir, "files/input/isolates2traits.csv"), DataFrame, missingstring="N/A")
df_metabolites          = CSV.read(joinpath(dir, "files/input/root_exudates.csv"), DataFrame, missingstring="N/A")

assimilation            = load(joinpath(dir, "files/output/isolates_assimilation.jld"))
enzymes                 = load(joinpath(dir, "files/output/isolates_enzymes.jld"))
maintenance             = load(joinpath(dir, "files/output/isolates_maintenance.jld"))
protein_synthesis       = load(joinpath(dir, "files/output/isolates_protein_synthesis.jld"))
turnover                = load(joinpath(dir, "files/output/isolates_turnover.jld"))
initb                   = load(joinpath(dir, "files/output/isolates_batch_init.jld"))

initb["Bio0"][37] = 0.0001

condition(u,t,integrator) = u[83] - 1e-5
affect!(integrator)       = terminate!(integrator)
cb                        = ContinuousCallback(condition,affect!)

percent_uptake    = zeros(83, 39)
J_D_tseries       = zeros(39, 83, 500)

p                 = DEBmicroTrait.init_mixed_medium(id_isolate, 83, assimilation, enzymes, maintenance, protein_synthesis, turnover)
n_polymers        = p.setup_pars.n_polymers
n_monomers        = p.setup_pars.n_monomers
n_microbes        = p.setup_pars.n_microbes

u0                                                                         = zeros(p.setup_pars.dim)
u0[1+n_polymers+n_monomers:n_polymers+n_monomers+n_microbes]              .= 0.9*initb["Bio0"][id_isolate]
u0[1+n_polymers+n_monomers+n_microbes:n_polymers+n_monomers+2*n_microbes] .= 0.1*initb["Bio0"][id_isolate]
u0[1+n_polymers:n_polymers+n_monomers]                                    .= 1.25/n_monomers


tspan             = (0.0,30.0)
prob              = ODEProblem(DEBmicroTrait.batch_model!,u0,tspan,p)
sol               = solve(prob, Tsit5())

plot(sol, idxs=1:83)

for i in 1:39
    id_isolate = i

    p                 = DEBmicroTrait.init_mixed_medium(id_isolate, assimilation, enzymes, maintenance, protein_synthesis, turnover)
    n_polymers        = p.setup_pars.n_polymers
    n_monomers        = p.setup_pars.n_monomers
    n_microbes        = p.setup_pars.n_microbes

    u0                                                                         = zeros(p.setup_pars.dim)
    u0[1+n_polymers+n_monomers:n_polymers+n_monomers+n_microbes]              .= 0.9*initb["Bio0"][id_isolate]
    u0[1+n_polymers+n_monomers+n_microbes:n_polymers+n_monomers+2*n_microbes] .= 0.1*initb["Bio0"][id_isolate]
    u0[1+n_polymers:n_polymers+n_monomers]                                    .= 1.25/n_monomers


    tspan             = (0.0,1500.0)
    prob              = ODEProblem(DEBmicroTrait.batch_model!,u0,tspan,p)
    sol               = solve(prob, alg_hints=[:stiff], callback=cb)

    for i = 1:n_monomers
        percent_uptake[i,id_isolate] = (sol[1][i] - sol[end][i])/sol[1][i]
    end

    for k in 1:length(sol.t)
        J_D      = DEBmicroTrait.uptake!(zeros(n_monomers), p.assimilation_pars, sol[k][1:n_monomers], [sol[k][3]])
        J_D_tseries[i,:,k] = J_D[:]
    end

end

J_D_median = zeros(39,83)
for i in 1:39
    for j in 1:83
        try
            J_D_median[i,j]  = median(filter(!iszero, J_D_tseries[i,j,:]))
        catch
            J_D_median[i,j]  = NaN
        end
    end
end

CSV.write(joinpath(dir, "files/output/isolates_levins.csv"), Tables.table(J_D_median))

df_out = DataFrame()

monomer = Array{String}(undef,39,83)
for i in 1:83
     monomer[:,i] .= df_metabolites.Name[i]
end
df_out.monomer = vec(monomer)

ontology = Array{String}(undef,39,83)
for i in 1:83
     ontology[:,i] .= df_metabolites.Ontology[i]
 end
df_out.ontology = vec(ontology)

df_out.isolate = vec(repeat(df_isolates.Isolate, 83))
df_out.response= vec(repeat(df_isolates.Rhizosphere_response, 83))


reform = Array{Float64}(undef,39,83)
for i in 1:83
     reform[:,i] .= percent_uptake[i,:]
end

df_out.percent_uptake = vec(reform)

CSV.write(joinpath(dir, "files/output/isolates_mixed_medium.csv"), df_out)
