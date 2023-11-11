using BenchmarkTools
using CSV
using DEBmicroTrait
using DataFrames
using DifferentialEquations
using JLD
using JLD2
using Plots

dir                     = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()
df_isolates             = CSV.read(joinpath(dir, "files/input/isolates2traits.csv"), DataFrame, missingstring="N/A");
df_metabolites          = CSV.read(joinpath(dir, "files/input/root_exudates.csv"), DataFrame, missingstring="N/A");

assimilation            = load(joinpath(dir, "files/output/isolates_assimilation_plant.jld"));
enzymes                 = load(joinpath(dir, "files/output/isolates_enzymes.jld"));
maintenance             = load(joinpath(dir, "files/output/isolates_maintenance.jld"));
protein_synthesis       = load(joinpath(dir, "files/output/isolates_protein_synthesis.jld"));
turnover                = load(joinpath(dir, "files/output/isolates_turnover.jld"));
initb                   = load(joinpath(dir, "files/output/isolates_batch_init.jld"));
avena_ex                = load(joinpath(dir, "files/avena/avena_exudation_2006_2022.jld2"));
avena_root              = load(joinpath(dir, "files/avena/avena_root_2006_2022.jld2"));

p                       = DEBmicroTrait.init_rhizosphere_model(assimilation, enzymes, maintenance, protein_synthesis, turnover, avena_ex, avena_root);

u0                                                                                      = zeros(p.setup_pars.dim);
u0[1:p.setup_pars.n_polymers]                                                          .= 0.0;
u0[1+p.setup_pars.n_polymers+p.setup_pars.n_monomers:p.setup_pars.n_polymers+p.setup_pars.n_monomers+p.setup_pars.n_microbes]              .= 0.9.*initb["Bio0"];
u0[1+p.setup_pars.n_polymers+p.setup_pars.n_monomers+p.setup_pars.n_microbes:p.setup_pars.n_polymers+p.setup_pars.n_monomers+2*p.setup_pars.n_microbes] .= 0.1.*initb["Bio0"];
u0[1+p.setup_pars.n_polymers:p.setup_pars.n_polymers+p.setup_pars.n_monomers]                                    .= 1.25;


tspan             = (0.0, 200.0);
prob              = ODEProblem(DEBmicroTrait.rhizosphere_model!,u0,tspan,p);
sol               = solve(prob, alg_hints=[:auto]);

plot(sol, idxs=87:92)

