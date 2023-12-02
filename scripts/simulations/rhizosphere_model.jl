using BenchmarkTools
using CSV
using DEBmicroTrait
using DataFrames
using DifferentialEquations
using JLD
using JLD2
using Plots

####################################################################
# I/O 
dir                     = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()
df_isolates             = CSV.read(joinpath(dir, "files/input/isolates2traits.csv"), DataFrame, missingstring="N/A");
df_metabolites          = CSV.read(joinpath(dir, "files/input/root_exudates.csv"), DataFrame, missingstring="N/A");

####################################################################
# parameterization 
assimilation            = load(joinpath(dir, "files/output/isolates_assimilation_plant.jld"));
enzymes                 = load(joinpath(dir, "files/output/isolates_enzymes.jld"));
maintenance             = load(joinpath(dir, "files/output/isolates_maintenance.jld"));
protein_synthesis       = load(joinpath(dir, "files/output/isolates_protein_synthesis.jld"));
turnover                = load(joinpath(dir, "files/output/isolates_turnover.jld"));
initb                   = load(joinpath(dir, "files/output/isolates_batch_init.jld"));
avena_ex                = load(joinpath(dir, "files/avena/avena_exudation_2006_2022.jld2"));
avena_root              = load(joinpath(dir, "files/avena/avena_root_2006_2022.jld2"));
p                       = DEBmicroTrait.init_rhizosphere_model(assimilation, enzymes, maintenance, protein_synthesis, turnover, avena_ex, avena_root);
#save(joinpath(dir, "files/output/rhizosphere_model_p.jld"), "params", p)

####################################################################
# forcing 
month_hours = 365.25/12*24;
months = 12;
nyears  = 2009-2006;
tstop = month_hours*months*nyears;
xticks = collect(0:month_hours:tstop);
# root decay  
spl_root = DEBmicroTrait.root_decay(p.forcing_pars);
dosetimes_root = collect(1:tstop);
function affect1!(integrator) 
    integrator.u[1] += abs.(spl_root(integrator.t + integrator.t*(nyears-1)))
end;
cb1 = PresetTimeCallback(dosetimes_root, affect1!);
# exudation rate
spl_ex   = DEBmicroTrait.root_exudation(p.forcing_pars);
ex_norm =  assimilation["LCMS"];

tmp3 =Vector{Vector{Float64}}(undef, nyears);
tmp6 =Vector{Vector{Float64}}(undef, nyears);
tmp9 =Vector{Vector{Float64}}(undef, nyears);
tmp12 =Vector{Vector{Float64}}(undef, nyears);
for i in 1:nyears
    tmp3[i] = collect(xticks[1+months*(i-1)]+1:xticks[3+months*(i-1)])
    tmp6[i] = collect(xticks[3+months*(i-1)]+1:xticks[4+months*(i-1)])
    tmp9[i] = collect(xticks[4+months*(i-1)]+1:xticks[5+months*(i-1)])
    tmp12[i] = collect(xticks[5+months*(i-1)]+1:xticks[7+months*(i-1)])
end;

dosetimes_week3  = vcat(tmp3...);
dosetimes_week6  = vcat(tmp6...);
dosetimes_week9  = vcat(tmp9...);
dosetimes_week12 = vcat(tmp12...);

function affect2!(integrator) 
    integrator.u[2:7] .+= abs.(spl_ex(integrator.t + integrator.t*(nyears-1)).*ex_norm[:,1])
end;
cb2 = PresetTimeCallback(dosetimes_week3, affect2!);

function affect3!(integrator) 
    integrator.u[2:7] .+= abs.(spl_ex(integrator.t + integrator.t*(nyears-1)).*ex_norm[:,2])
end;
cb3 = PresetTimeCallback(dosetimes_week6, affect3!);

function affect4!(integrator) 
    integrator.u[2:7] .+= abs.(spl_ex(integrator.t + integrator.t*(nyears-1)).*ex_norm[:,3])
end;
cb4 = PresetTimeCallback(dosetimes_week9, affect4!);

function affect5!(integrator) 
    integrator.u[2:7] .+= abs.(spl_ex(integrator.t + integrator.t*(nyears-1)).*ex_norm[:,4])
end;
cb5 = PresetTimeCallback(dosetimes_week12, affect5!);
cb = CallbackSet(cb1, cb2, cb3, cb4, cb5);

####################################################################
# simulation
u0                                                                                      = zeros(p.setup_pars.dim);
u0[1:p.setup_pars.n_polymers]                                                          .= 10.0;
u0[1+p.setup_pars.n_polymers+p.setup_pars.n_monomers:p.setup_pars.n_polymers+p.setup_pars.n_monomers+p.setup_pars.n_microbes] .= 0.9.*median(initb["Bio0"]).*(df_isolates.week0./sum(df_isolates.week0));
u0[1+p.setup_pars.n_polymers+p.setup_pars.n_monomers+p.setup_pars.n_microbes:p.setup_pars.n_polymers+p.setup_pars.n_monomers+2*p.setup_pars.n_microbes] .= 0.1.*median(initb["Bio0"]).*(df_isolates.week0./sum(df_isolates.week0));
u0[1+p.setup_pars.n_polymers:p.setup_pars.n_polymers+p.setup_pars.n_monomers]          .= 0.0;
u0[1+p.setup_pars.n_polymers+p.setup_pars.n_monomers+2*p.setup_pars.n_microbes:p.setup_pars.n_polymers+p.setup_pars.n_monomers+2*p.setup_pars.n_microbes+p.setup_pars.n_enzymes] .= 0.0;

tspan             = (0, tstop);
prob              = ODEProblem(DEBmicroTrait.rhizosphere_model!,u0,tspan,p);
sol               = solve(prob, alg_hints=[:stiff], callback = cb);

####################################################################
# plotting
n_polymers                = p.setup_pars.n_polymers;
n_monomers                = p.setup_pars.n_monomers;
n_microbes                = p.setup_pars.n_microbes;
n_enzymes                 = p.setup_pars.n_enzymes;
n_minerals                = p.setup_pars.n_minerals;

polymersidx = 1:p.setup_pars.n_polymers;
monomersidx = 1+n_polymers:n_polymers+n_monomers;
reservesidx = 1+n_polymers+n_monomers:n_polymers+n_monomers+n_microbes;
structureidx = 1+n_polymers+n_monomers+n_microbes:n_polymers+n_monomers+2*n_microbes;
enzymesidx = 1+n_polymers+n_monomers+2*n_microbes:n_polymers+n_monomers+2*n_microbes+n_enzymes;
monomersadsidx = 1+n_polymers+n_monomers+2*n_microbes+n_enzymes:n_polymers+n_monomers+2*n_microbes+n_enzymes+n_monomers;
enzymesadsidx = 1+n_polymers+n_monomers+2*n_microbes+n_enzymes+n_monomers:n_polymers+n_monomers+2*n_microbes+2*n_enzymes+n_monomers;
co2idx = 1+n_polymers+n_monomers+2*n_microbes+2*n_enzymes+n_monomers:n_polymers+n_monomers+2*n_microbes+2*n_enzymes+n_monomers+n_microbes;

plot(sol, idxs=monomersidx)

####################################################################
# post-processing
sample_times              = floor.(Int, (vcat(dosetimes_week3[end], dosetimes_week6[end], dosetimes_week9[end], dosetimes_week12[end])));


E                         = sol[reservesidx, sample_times];
V                         = sol[structureidx, sample_times];
Bio                       = @. E+V;
rel_Bio                   = zeros(n_microbes, size(sample_times,1));
for i in 1:n_microbes
    for j in 1:size(sample_times,1)
        rel_Bio[i,j] = Bio[i,j]./sum(Bio[:,j])
    end
end;

growth_rate               = zeros(n_microbes, size(sample_times,1));
for i in 1:n_microbes
    for j in eachindex(sample_times)
        growth_rate[i,j] = DEBmicroTrait.growth!(0.0*ones(1), p.metabolism_pars, [E[i,j]], [V[i,j]])[1]
    end
end;

BR                       = zeros(n_microbes, size(sample_times, 1));
du                       = similar(u0);
for i in 1:n_microbes
    for j in eachindex(sample_times)
        BR[i,j] = DEBmicroTrait.rhizosphere_model!(du, sol.u[sample_times[j]], p, 0)[93+i]
    end
end;

BP                       = zeros(n_microbes, size(sample_times, 1));
for i in 1:n_microbes
    for j in eachindex(sample_times)
        BP[i,j] = DEBmicroTrait.rhizosphere_model!(du, sol.u[sample_times[j]], p, 0)[7+i] + DEBmicroTrait.rhizosphere_model!(du, sol.u[sample_times[j]], p, 0)[46+i]
    end
end;


k2p_D        = 180.0*60*60*ones(n_monomers, n_microbes);
k2p_M        = ones(n_monomers, n_minerals);
k2p          = hcat(k2p_D, k2p_M);
D            = sol[monomersidx, sample_times];
M            = p.assimilation_pars.M;
K_D          = p.assimilation_pars.K_D;
N_SB         = p.assimilation_pars.N_SB;
CUE          =  zeros(n_microbes, size(sample_times, 1));       
for j in eachindex(sample_times)
    ECA  = DEBmicroTrait.ECA_kinetics!(zeros(n_monomers,n_microbes+n_minerals), D[:,j], vcat(V[:,j], M), K_D, k2p, N_SB).*p.assimilation_pars.N_C
    J_D  = vcat(sum(ECA[:, 1:n_microbes], dims=1)...)
    CUE[:,j] = @. 1 - BR[:,j]/(J_D)
end;


df_out = DataFrame();
df_out.isolate = repeat(df_isolates.Isolate, size(sample_times,1));
df_out.response = repeat(df_isolates.Rhizosphere_response, size(sample_times,1));
df_out.sample   = vcat(repeat(["week3"], 39), repeat(["week6"], 39), repeat(["week9"], 39), repeat(["week12"], 39));
df_out.relabundance = vec(rel_Bio);
df_out.measurement = vcat(df_isolates.week3./sum(df_isolates.week3), df_isolates.week6./sum(df_isolates.week6), df_isolates.week9./sum(df_isolates.week9), df_isolates.week12./sum(df_isolates.week12));
df_out.rgrowth = vec(growth_rate);
df_out.BR = vec(BR);
df_out.BP = vec(BP);
df_out.CUE = vec(CUE);
df_out.biomass = vec(Bio);
CSV.write(joinpath(dir, "files/output/rhizosphere_model_2008.csv"), df_out)


