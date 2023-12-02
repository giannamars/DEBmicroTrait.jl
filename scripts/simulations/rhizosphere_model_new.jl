using DifferentialEquations

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

# root decay  
spl_root = DEBmicroTrait.root_decay(p.forcing_pars);
dosetimes_root = collect(1:tstop);
function affect1!(integrator) 
    integrator.u[polymersidx] .+= abs.(spl_root(integrator.t))
end;
cb1 = PresetTimeCallback(dosetimes_root, affect1!);

# exudation rate
spl_ex   = DEBmicroTrait.root_exudation(p.forcing_pars);
ex_norm =  assimilation["LCMS"];

dosetimes_week3 = collect(week3start:week3start+week_hours*3);
dosetimes_week6 = collect(week3start+week_hours*3+1:week3start+week_hours*6);
dosetimes_week9 = collect(week3start+week_hours*6+1:week3start+week_hours*9);
dosetimes_week12 = collect(week3start+week_hours*9+1:4214);

function affect2!(integrator) 
    integrator.u[monomersidx] .+= abs.(spl_ex(integrator.t).*ex_norm[:,1])
end;
cb2 = PresetTimeCallback(dosetimes_week3, affect2!);

function affect3!(integrator) 
    integrator.u[monomersidx] .+= abs.(spl_ex(integrator.t).*ex_norm[:,2])
end;
cb3 = PresetTimeCallback(dosetimes_week6, affect3!);

function affect4!(integrator) 
    integrator.u[monomersidx] .+= abs.(spl_ex(integrator.t).*ex_norm[:,3])
end;
cb4 = PresetTimeCallback(dosetimes_week9, affect4!);

function affect5!(integrator) 
    integrator.u[monomersidx] .+= abs.(spl_ex(integrator.t).*ex_norm[:,4])
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

plot(sol, idxs = monomersidx)
####################################################################
# Post-processing
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

df_out = DataFrame();
df_out.isolate = repeat(df_isolates.Isolate, size(sample_times,1));
df_out.response = repeat(df_isolates.Rhizosphere_response, size(sample_times,1));
df_out.sample   = vcat(repeat(["week3"], 39), repeat(["week6"], 39), repeat(["week9"], 39), repeat(["week12"], 39));
df_out.relabundance = vec(rel_Bio);
df_out.measurement = vcat(df_isolates.week3./sum(df_isolates.week3), df_isolates.week6./sum(df_isolates.week6), df_isolates.week9./sum(df_isolates.week9), df_isolates.week12./sum(df_isolates.week12));
CSV.write(joinpath(dir, "files/output/rhizosphere_model_2006.csv"), df_out)
