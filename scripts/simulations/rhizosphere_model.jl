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
u0[1:p.setup_pars.n_polymers]                                                          .= 10.0;
u0[1+p.setup_pars.n_polymers+p.setup_pars.n_monomers:p.setup_pars.n_polymers+p.setup_pars.n_monomers+p.setup_pars.n_microbes]              .= 0.9.*initb["Bio0"];
u0[1+p.setup_pars.n_polymers+p.setup_pars.n_monomers+p.setup_pars.n_microbes:p.setup_pars.n_polymers+p.setup_pars.n_monomers+2*p.setup_pars.n_microbes] .= 0.1.*initb["Bio0"];
u0[1+p.setup_pars.n_polymers:p.setup_pars.n_polymers+p.setup_pars.n_monomers]                                    .= 0.0;

 # forcing 
spl_ex   = DEBmicroTrait.root_exudation(p.forcing_pars)
ex_norm =  assimilation["LCMS"]
spl_root = DEBmicroTrait.root_decay(p.forcing_pars)

dosetimes_root = collect(1:8760)
function affect1!(integrator) 
    integrator.u[1] += spl_root(integrator.t)
end
cb1 = PresetTimeCallback(dosetimes_root, affect1!)

function affect2!(integrator) 
    integrator.u[2:7] .+= spl_ex(integrator.t).*ex_norm[:,1]
end
dosetimes_week3 = collect(0:1200)
cb2 = PresetTimeCallback(dosetimes_week3, affect2!)

function affect3!(integrator) 
    integrator.u[2:7] .+= spl_ex(integrator.t).*ex_norm[:,2]
end
dosetimes_week6 = collect(1201:2000)
cb3 = PresetTimeCallback(dosetimes_week6, affect3!)

function affect4!(integrator) 
    integrator.u[2:7] .+= spl_ex(integrator.t).*ex_norm[:,3]
end
dosetimes_week9 = collect(2001:3000)
cb4 = PresetTimeCallback(dosetimes_week9, affect4!)

function affect5!(integrator) 
    integrator.u[2:7] .+= spl_ex(integrator.t).*ex_norm[:,4]
end
dosetimes_week12 = collect(3001:4500)
cb5 = PresetTimeCallback(dosetimes_week12, affect5!)




cb = CallbackSet(cb1, cb2, cb3, cb4, cb5)

tspan             = (0.0, 8760.0);
prob              = ODEProblem(DEBmicroTrait.rhizosphere_model!,u0,tspan,p);
sol               = solve(prob, alg_hints=[:auto], callback = cb);

plot(sol, idxs=8:45)

