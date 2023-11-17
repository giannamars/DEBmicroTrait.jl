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

df_isolates.week0./maximum(df_isolates.week0)

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

# forcing 
month_hours = 365.25/12*24
months = 12
nyears  = 2007-2006
tstop = month_hours*months*nyears
xticks = collect(0:month_hours:tstop)
# root decay  
spl_root = DEBmicroTrait.root_decay(p.forcing_pars)
dosetimes_root = collect(1:tstop)
function affect1!(integrator) 
    integrator.u[1] += abs.(spl_root(integrator.t))
end
cb1 = PresetTimeCallback(dosetimes_root, affect1!)
# exudation rate
spl_ex   = DEBmicroTrait.root_exudation(p.forcing_pars)
ex_norm =  assimilation["LCMS"]

tmp3 =Vector{Vector{Float64}}(undef, nyears)
tmp6 =Vector{Vector{Float64}}(undef, nyears)
tmp9 =Vector{Vector{Float64}}(undef, nyears)
tmp12 =Vector{Vector{Float64}}(undef, nyears)
for i in 1:nyears
    tmp3[i] = collect(xticks[1+months*(i-1)]+1:xticks[3+months*(i-1)])
    tmp6[i] = collect(xticks[3+months*(i-1)]+1:xticks[4+months*(i-1)])
    tmp9[i] = collect(xticks[4+months*(i-1)]+1:xticks[5+months*(i-1)])
    tmp12[i] = collect(xticks[5+months*(i-1)]+1:xticks[7+months*(i-1)])
end

dosetimes_week3  = vcat(tmp3...)
dosetimes_week6  = vcat(tmp6...)
dosetimes_week9  = vcat(tmp9...)
dosetimes_week12 = vcat(tmp12...)


function affect2!(integrator) 
    integrator.u[2:7] .+= abs.(spl_ex(integrator.t).*ex_norm[:,1])
end
cb2 = PresetTimeCallback(dosetimes_week3, affect2!)

function affect3!(integrator) 
    integrator.u[2:7] .+= abs.(spl_ex(integrator.t).*ex_norm[:,2])
end
cb3 = PresetTimeCallback(dosetimes_week6, affect3!)

function affect4!(integrator) 
    integrator.u[2:7] .+= abs.(spl_ex(integrator.t).*ex_norm[:,3])
end
cb4 = PresetTimeCallback(dosetimes_week9, affect4!)

function affect5!(integrator) 
    integrator.u[2:7] .+= abs.(spl_ex(integrator.t).*ex_norm[:,4])
end
cb5 = PresetTimeCallback(dosetimes_week12, affect5!)
cb = CallbackSet(cb1, cb2, cb3, cb4, cb5)

# simulate
u0                                                                                      = zeros(p.setup_pars.dim);
u0[1:p.setup_pars.n_polymers]                                                          .= 10.0;
#u0[1+p.setup_pars.n_polymers+p.setup_pars.n_monomers:p.setup_pars.n_polymers+p.setup_pars.n_monomers+p.setup_pars.n_microbes]              .= 0.9.*initb["Bio0"];
#u0[1+p.setup_pars.n_polymers+p.setup_pars.n_monomers+p.setup_pars.n_microbes:p.setup_pars.n_polymers+p.setup_pars.n_monomers+2*p.setup_pars.n_microbes] .= 0.1.*initb["Bio0"];
u0[1+p.setup_pars.n_polymers:p.setup_pars.n_polymers+p.setup_pars.n_monomers]                                    .= 0.0;
#u0[1+p.setup_pars.n_polymers+p.setup_pars.n_monomers+2*p.setup_pars.n_microbes:p.setup_pars.n_polymers+p.setup_pars.n_monomers+2*p.setup_pars.n_microbes+p.setup_pars.n_enzymes] .= 10.0
u0[1+p.setup_pars.n_polymers+p.setup_pars.n_monomers+2*p.setup_pars.n_microbes:p.setup_pars.n_polymers+p.setup_pars.n_monomers+2*p.setup_pars.n_microbes+p.setup_pars.n_enzymes] .= 0.0

tspan             = (0, tstop);
prob              = ODEProblem(DEBmicroTrait.rhizosphere_model!,u0,tspan,p);
sol               = solve(prob, alg_hints=[:stiff], callback = cb);

# plotting
n_polymers                = p.setup_pars.n_polymers
n_monomers                = p.setup_pars.n_monomers
n_microbes                = p.setup_pars.n_microbes
n_enzymes                 = p.setup_pars.n_enzymes
n_minerals                = p.setup_pars.n_minerals

polymersidx = 1:p.setup_pars.n_polymers
monomersidx = 1+n_polymers:n_polymers+n_monomers
reservesidx = 1+n_polymers+n_monomers:n_polymers+n_monomers+n_microbes
structureidx = 1+n_polymers+n_monomers+n_microbes:n_polymers+n_monomers+2*n_microbes
enzymesidx = 1+n_polymers+n_monomers+2*n_microbes:n_polymers+n_monomers+2*n_microbes+n_enzymes
monomersadsidx = 1+n_polymers+n_monomers+2*n_microbes+n_enzymes:n_polymers+n_monomers+2*n_microbes+n_enzymes+n_monomers
enzymesadsidx = 1+n_polymers+n_monomers+2*n_microbes+n_enzymes+n_monomers:n_polymers+n_monomers+2*n_microbes+2*n_enzymes+n_monomers
co2idx = 1+n_polymers+n_monomers+2*n_microbes+2*n_enzymes+n_monomers:n_polymers+n_monomers+2*n_microbes+2*n_enzymes+n_monomers+n_microbes

plot(sol, idxs=structureid)

# post-processing
week3                     = 

E                         = sol[reserves, :];
V                         = sol[structure, :];
growth_rate               = zeros(n_microbes, size(sol.t, 1));
for i in 1:n_microbes
    for j in 1:size(sol.t,1)
        growth_rate[i,j] = DEBmicroTrait.growth!(0.0*ones(1), p.metabolism_pars, [E[i,j]], [V[i,j]])[1]
    end
end;

BR                       = zeros(n_microbes, size(sol.t, 1));
du                       = similar(u0)
for i in 1:n_microbes
    for j in 1:size(sol.t,1)
        BR[i,j] = DEBmicroTrait.rhizosphere_model!(du, sol.u[j], p, 0)[i]
    end
end;

BP                       = zeros(n_microbes, size(sol.t, 1));
for i in 1:n_microbes
    for j in 1:size(sol.t,1)
        BP[i,j] = DEBmicroTrait.rhizosphere_model!(du, sol.u[j], p, 0)[7+i] + DEBmicroTrait.rhizosphere_model!(du, sol.u[j], p, 0)[46+i]
    end
end;

BGE  = @. BP/(BP + BR)




#x                         = [DEBmicroTrait.growth_production!(r[i], p.metabolism_pars, E[:,i], V[:,i])[1] for i in 1:size(sol.t,1)]
#J_EX                      = [DEBmicroTrait.enzyme_production!(x[i], p.metabolism_pars, V[:,i])[1] for i in 1:size(sol.t,1)]