using BenchmarkTools
using CSV
using DEBmicroTrait
using DataFrames
using DifferentialEquations
using JLD
using JLD2
using Plots; gr()

####################################################################
# I/O 
dir                     = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()
df_isolates             = CSV.read(joinpath(dir, "files/input/isolates2traits.csv"), DataFrame, missingstring="N/A");
df_metabolites          = CSV.read(joinpath(dir, "files/input/root_exudates.csv"), DataFrame, missingstring="N/A");

####################################################################
# Parameterization 
assimilation            = load(joinpath(dir, "files/output/isolates_assimilation_plant.jld"));
enzymes                 = load(joinpath(dir, "files/output/isolates_enzymes.jld"));
maintenance             = load(joinpath(dir, "files/output/isolates_maintenance.jld"));
protein_synthesis       = load(joinpath(dir, "files/output/isolates_protein_synthesis.jld"));
turnover                = load(joinpath(dir, "files/output/isolates_turnover.jld"));
initb                   = load(joinpath(dir, "files/output/isolates_batch_init.jld"));
avena_ex                = load(joinpath(dir, "files/avena/avena_exudation_2006_2022.jld2"));
avena_root              = load(joinpath(dir, "files/avena/avena_root_2006_2022.jld2"));
p                       = DEBmicroTrait.init_rhizosphere_model(assimilation, enzymes, maintenance, protein_synthesis, turnover, avena_ex, avena_root);

####################################################################
# Dimensions
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


####################################################################
# Plant forcing: 2006
week_hours = 24*7;
month_hours = 365.25/12*24;
months = 12;
nyears  = 2007-2006;
tstop = month_hours*months*nyears;
week3start = 8760*(nyears-1)+1+804;
deadhr = 4214;

####################################################################
# Plant forcing: 2012
week_hours = 24*7;
month_hours = 365.25/12*24;
months = 12;
nyears  = 2013-2006;
tstop = month_hours*months*nyears;
week3start = 1+1007;
deadhr = 4578;

####################################################################
# Plant forcing: 2013
week_hours = 24*7;
month_hours = 365.25/12*24;
months = 12;
nyears  = 2014-2006;
tstop = month_hours*months*nyears;
week3start = 1+1007;
deadhr = 3740;

####################################################################
# Plant forcing: 2018
week_hours = 24*7;
month_hours = 365.25/12*24;
months = 12;
nyears  = 2019-2006;
tstop = month_hours*months*nyears;
week3start = 1+345;
deadhr = 4384;


# root decay  
spl_root = DEBmicroTrait.root_decay(p.forcing_pars);
dosetimes_root = collect(1:tstop);
function affect1!(integrator) 
    integrator.u[polymersidx] .+= abs.(spl_root(integrator.t + 8760*(nyears-1)))
end;
cb1 = PresetTimeCallback(dosetimes_root, affect1!);

# exudation rate
spl_ex   = DEBmicroTrait.root_exudation(p.forcing_pars);
ex_norm =  assimilation["LCMS"];

dosetimes_week3 = collect(week3start:week3start+week_hours*3);
dosetimes_week6 = collect(week3start+week_hours*3+1:week3start+week_hours*6);
dosetimes_week9 = collect(week3start+week_hours*6+1:week3start+week_hours*9);
dosetimes_week12 = collect(week3start+week_hours*9+1:deadhr);

function affect2!(integrator) 
    integrator.u[monomersidx] .+= abs.(spl_ex(integrator.t + 8760*(nyears-1)).*ex_norm[:,1])
end;
cb2 = PresetTimeCallback(dosetimes_week3, affect2!);

function affect3!(integrator) 
    integrator.u[monomersidx] .+= abs.(spl_ex(integrator.t + 8760*(nyears-1)).*ex_norm[:,2])
end;
cb3 = PresetTimeCallback(dosetimes_week6, affect3!);

function affect4!(integrator) 
    integrator.u[monomersidx] .+= abs.(spl_ex(integrator.t + 8760*(nyears-1)).*ex_norm[:,3])
end;
cb4 = PresetTimeCallback(dosetimes_week9, affect4!);

function affect5!(integrator) 
    integrator.u[monomersidx] .+= abs.(spl_ex(integrator.t + 8760*(nyears-1)).*ex_norm[:,4])
end;
cb5 = PresetTimeCallback(dosetimes_week12, affect5!);
cb = CallbackSet(cb1, cb2, cb3, cb4, cb5);

####################################################################
# Simulation
u0                   = zeros(p.setup_pars.dim);
u0[polymersidx]     .= 10.0;
u0[monomersidx]     .= 0.0;
u0[reservesidx]     .= 0.9.*median(initb["Bio0"]).*(df_isolates.week0./sum(df_isolates.week0));
u0[structureidx]    .= 0.1.*median(initb["Bio0"]).*(df_isolates.week0./sum(df_isolates.week0));
u0[enzymesidx]      .= 0.0;

tspan             = (0, deadhr);
prob              = ODEProblem(DEBmicroTrait.rhizosphere_model!,u0,tspan,p);
sol               = solve(prob, alg_hints=[:stiff], callback = cb);

####################################################################
# Plotting
plot(sol, idxs=monomersidx)

####################################################################
# Post-processing
sample_times              = floor.(Int, (vcat(0, dosetimes_week3[end], dosetimes_week6[end], dosetimes_week9[end], week3start+week_hours*12, deadhr)));

Bio                       = zeros(n_microbes, size(sample_times,1));
for i in 1:n_microbes
    for j in eachindex(sample_times)
        E            = sol(sample_times[j], idxs = reservesidx[i])
        V            = sol(sample_times[j], idxs = structureidx[i])
        Bio[i,j]     = @. E + V
    end
end;

rel_Bio                   = zeros(n_microbes, size(sample_times,1));
for i in 1:n_microbes
    for j in eachindex(sample_times)
        rel_Bio[i,j] = Bio[i,j]./sum(Bio[:,j])
    end
end;


growth_rate               = zeros(n_microbes, size(sample_times,1));
for i in 1:n_microbes
    for j in eachindex(sample_times)
        E            = sol(sample_times[j], idxs = reservesidx[i])
        V            = sol(sample_times[j], idxs = structureidx[i])
        growth_rate[i,j] = DEBmicroTrait.growth!(0.0*ones(1), p.metabolism_pars, [E], [V])[1]
    end
end;

BR                       = zeros(n_microbes, size(sample_times, 1));
du                       = similar(u0);
for j in eachindex(sample_times)
    BR[:,j] = DEBmicroTrait.rhizosphere_model!(du, sol(sample_times[j]), p, 0)[co2idx]
end;

BP                       = zeros(n_microbes, size(sample_times, 1));
du                       = similar(u0);
for j in eachindex(sample_times)
    BP[:,j] = DEBmicroTrait.rhizosphere_model!(du, sol(sample_times[j]), p, 0)[reservesidx] + DEBmicroTrait.rhizosphere_model!(du, sol(sample_times[j]), p, 0)[structureidx]
end


k2p_D        = 180.0*60*60*ones(n_monomers, n_microbes);
k2p_M        = ones(n_monomers, n_minerals);
k2p          = hcat(k2p_D, k2p_M);
M            = p.assimilation_pars.M;
K_D          = p.assimilation_pars.K_D;
N_SB         = p.assimilation_pars.N_SB;
CUE          =  zeros(n_microbes, size(sample_times, 1));       
for j in eachindex(sample_times)
    D           = sol(sample_times[j], idxs = monomersidx)
    V           = sol(sample_times[j], idxs = structureidx)
    ECA  = DEBmicroTrait.ECA_kinetics!(zeros(n_monomers,n_microbes+n_minerals), D, vcat(V, M), K_D, k2p, N_SB).*p.assimilation_pars.N_C
    J_D  = vcat(sum(ECA[:, 1:n_microbes], dims=1)...)
    CUE[:,j] = @. 1 - BR[:,j]/(J_D)
end; 


df_out = DataFrame();
df_out.phylum = repeat(df_isolates.Phylum, size(sample_times,1));
df_out.isolate = repeat(df_isolates.Isolate, size(sample_times,1));
df_out.response = repeat(df_isolates.Rhizosphere_response, size(sample_times,1));
df_out.sample   = vcat(repeat(["week0"], 39), repeat(["week3"], 39), repeat(["week6"], 39), repeat(["week9"], 39), repeat(["week12"], 39), repeat(["dead"], 39));
df_out.relabundance = vec(rel_Bio);
df_out.rgrowth = vec(growth_rate);
df_out.BR = vec(BR);
df_out.BP = vec(BP);
df_out.CUE = vec(CUE);
df_out.biomass = vec(Bio);
df_out.measurement = vcat(df_isolates.week0./sum(df_isolates.week0), df_isolates.week3./sum(df_isolates.week3), df_isolates.week6./sum(df_isolates.week6), df_isolates.week9./sum(df_isolates.week9), df_isolates.week12./sum(df_isolates.week12), repeat([NaN], 39));
CSV.write(joinpath(dir, "files/output/rhizosphere_model_0_2018.csv"), df_out)
