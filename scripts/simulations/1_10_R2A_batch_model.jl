using DEBmicroTrait
using CSV, DataFrames, Statistics
using JLD
using DifferentialEquations

dir                     = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()

# Load media composition: Formula, Molecular weight, Medium concentration
df_metabolites = CSV.read(joinpath(dir, "files/input/1_10_R2A_medium.csv"), DataFrame, missingstring="N/A")

# Load isolate parameterization, note: assimilation parameterization depends on media composition
assimilation            = load(joinpath(dir, "files/output/isolates_assimilation_10_R2A.jld"))
enzymes                 = load(joinpath(dir, "files/output/isolates_enzymes.jld"))
maintenance             = load(joinpath(dir, "files/output/isolates_maintenance.jld"))
protein_synthesis       = load(joinpath(dir, "files/output/isolates_protein_synthesis.jld"))
turnover                = load(joinpath(dir, "files/output/isolates_turnover.jld"))
initb                   = load(joinpath(dir, "files/output/isolates_batch_init.jld"))

r_tseries         = zeros(39, 500)
D_tseries         = zeros(39, 500)
E_tseries         = zeros(39, 500)
V_tseries         = zeros(39, 500)
t_tseries         = zeros(39, 500)

for i in 1:39
    id_isolate = i

    p                 = DEBmicroTrait.init_mixed_medium_r2(id_isolate, assimilation, enzymes, maintenance, protein_synthesis, turnover)
    n_polymers        = p.setup_pars.n_polymers
    n_monomers        = p.setup_pars.n_monomers
    n_microbes        = p.setup_pars.n_microbes

    u0                                                                         = zeros(p.setup_pars.dim)
    u0[1+n_polymers+n_monomers:n_polymers+n_monomers+n_microbes]              .= 0.9*initb["Bio0"][id_isolate]
    u0[1+n_polymers+n_monomers+n_microbes:n_polymers+n_monomers+2*n_microbes] .= 0.1*initb["Bio0"][id_isolate]
    u0[1+n_polymers:n_polymers+n_monomers]                                    .= df_metabolites.Concentration


    tspan             = (0.0,500.0)
    prob              = ODEProblem(DEBmicroTrait.batch_model!,u0,tspan,p)
    sol               = solve(prob, alg_hints=[:stiff])

    r    = [DEBmicroTrait.growth!(0.0*ones(1), p.metabolism_pars, [sol[i][44]], [sol[i][45]])[1] for i in 1:size(sol.t,1)]

    for k in 1:length(sol.t)
        r_tseries[i,k] = r[k]
    end

    for k in 1:length(sol.t)
        t_tseries[i,k] =  sol.t[k]
        D_tseries[i,k] =  sol[k][1]
        E_tseries[i,k] =  sol[k][44]
        V_tseries[i,k] =  sol[k][45]
    end
end


r_median        = zeros(39)
for i in 1:39
    try
        r_median[i]  = median(filter(!iszero, r_tseries[i,:]))
    catch
        r_median[i]  = NaN
    end
end


df_isolates     = CSV.read(joinpath(dir, "files/input/isolates2traits.csv"), DataFrame, missingstring="N/A")
gmax = log(2)./df_isolates.Min_gen_time

df_out = DataFrame()
df_out.rmedian = r_median
df_out.gmax = gmax


# Statistics
using GLM
linearRegressor = lm(@formula(log(rmedian) ~ log(gmax)), df_out)
r2(linearRegressor)
# test for significant difference to 1-1 line
df_out.diff = r_median-gmax
linearRegressor = lm(@formula(diff ~ rmedian), df_out)
# RMSD
RMSD = sqrt((1/(39-1)*sum((r_median - gmax).^2)))

