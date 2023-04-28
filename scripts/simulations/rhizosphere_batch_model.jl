using DEBmicroTrait
using CSV, DataFrames, Statistics
using JLD
using OrdinaryDiffEq

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
BR_tseries        = zeros(39, 83, 500) # BR = assimilation respiration + growth respiration + enzyme production respiration + maintenance respiration
BP_tseries        = zeros(39, 83, 500) # BP = E + V production
r_tseries         = zeros(39, 83, 500) # realized growth rate
x_tseries         = zeros(39, 83, 500) # constitutive enzyme production rate
rG_CO2_tseries    = zeros(39, 83, 500) # growth respiration
rM_CO2_tseries    = zeros(39, 83, 500) # maintenance respiration
rX_CO2_tseries    = zeros(39, 83, 500) # enzyme production respiration
J_EX_tseries      = zeros(39, 83, 500) # total enzyme production rate, i.e. x*V
J_DE_tseries      = zeros(39, 83, 500) # assimilation flux
J_DE_CO2_tseries  = zeros(39, 83, 500) # assimilation respiration
J_D_tseries       = zeros(39, 83, 500) # uptake rate
J_ED_tseries      = zeros(39, 83, 500) # recycling of reserve biomass E to substrates D
J_V_tseries       = zeros(39, 83, 500) # strucutral biomass turnover rate
J_E_tseries       = zeros(39, 83, 500) # reserve biomass turnover rate
t_tseries         = zeros(39, 83, 500) # time
D_tseries         = zeros(39, 83, 500) # Substrate concentration
E_tseries         = zeros(39, 83, 500) # Reverve concentration
V_tseries         = zeros(39, 83, 500) # Structural biomass concentration
X_tseries         = zeros(39, 83, 500) # Enzyme concentration
CO2_tseries       = zeros(39, 83, 500) # CO2 concentration
N_cells_tseries   = zeros(39, 83, 500) # Conversion from total biomass to number of cells
maintenance_tseries    = zeros(39, 83, 500) # maintenance respiration 

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

        tspan             = (0.0,1000.0)
        prob              = ODEProblem(DEBmicroTrait.batch_model!,u0,tspan,p)
        sol               = solve(prob, Tsit5())

        for k in 1:length(sol.t)
            t_tseries[i,j,k] = sol.t[k]
            D_tseries[i,j,k] =  sol[k][1]
            E_tseries[i,j,k] =  sol[k][2]
            V_tseries[i,j,k] =  sol[k][3]
            X_tseries[i,j,k] =  sol[k][4]
            CO2_tseries[i,j,k] =  sol[k][5]
            Bio = sol[k][2].+sol[k][3]
            N_cells_tseries[i,j,k]  = @. Bio[1]*1e-6*12.011/(initb["rhoB"]*initb["Md"])[1]
        end
        du   = zeros(p.setup_pars.dim)
        BR   = [DEBmicroTrait.batch_model!(du, sol.u[i], p, 0)[end] for i in 1:size(sol.t,1)]
        BP   = [DEBmicroTrait.batch_model!(du, sol.u[i], p, 0)[2] + DEBmicroTrait.batch_model!(du, sol.u[i], p, 0)[3] for i in 1:size(sol.t,1)]
        BGE  = @. BP/(BP + BR)
        for k in 1:length(sol.t)
            BGE_tseries[i,j,k] = BGE[k]
        end
        for k in 1:length(sol.t)
            BP_tseries[i,j,k] = BP[k]
        end
        for k in 1:length(sol.t)
            BR_tseries[i,j,k] = BR[k]
        end
        r    = [DEBmicroTrait.growth!(0.0*ones(1), p.metabolism_pars, [sol[i][2]], [sol[i][3]])[1] for i in 1:size(sol.t,1)]
        for k in 1:length(sol.t)
            r_tseries[i,j,k] = r[k]
        end

        for k in 1:length(sol.t)
            x, rG_CO2, rM_CO2, rX_CO2 = DEBmicroTrait.growth_production!(r[k], p.metabolism_pars, [sol[k][2]], [sol[k][3]])
            x_tseries[i,j,k] = x[1]
            rG_CO2_tseries[i,j,k] = rG_CO2[1]
            rM_CO2_tseries[i,j,k] = rM_CO2[1]
            rX_CO2_tseries[i,j,k] = rX_CO2[1]
            J_EX   = DEBmicroTrait.enzyme_production!(x[1], p.metabolism_pars, [sol[k][3]])
            J_EX_tseries[i,j,k] = J_EX[1]
        end

        for k in 1:length(sol.t)
            J_DE  = DEBmicroTrait.assimilation!(zeros(1), p.assimilation_pars, [sol[k][1]], [sol[k][3]])
            J_DE_tseries[i,j,k] = J_DE[1]
            J_DE_CO2 = DEBmicroTrait.assimilation_production!(zeros(1), p.assimilation_pars, [sol[k][1]], [sol[k][3]])
            J_DE_CO2_tseries[i,j,k] = J_DE_CO2[1]
            J_D      = DEBmicroTrait.uptake!(zeros(1), p.assimilation_pars, [sol[k][1]], [sol[k][3]])
            J_D_tseries[i,j,k] = J_D[1]
        end

        for k in 1:length(sol.t)
            J_ED         = DEBmicroTrait.reserve_recycling!(zeros(1), p.turnover_pars, [sol[k][2]])
            J_ED_tseries[i,j,k] = J_ED[1]
            J_V          = DEBmicroTrait.biomass_turnover!(zeros(1), p.turnover_pars, [sol[k][3]])
            J_V_tseries[i,j,k] = J_V[1]
            J_E          = DEBmicroTrait.biomass_turnover!(zeros(1), p.turnover_pars, [sol[k][2]])
            J_E_tseries[i,j,k] = J_E[1]
        end

        for k in 1:length(sol.t)
            x, rG_CO2, rM_CO2, rX_CO2 = DEBmicroTrait.growth_production!(r[k], p.metabolism_pars, [sol[k][2]], [sol[k][3]])
            maintenance_tseries[i,j,k] = rM_CO2[1]./(rG_CO2[1]+rM_CO2[1]+rX_CO2[1])
        end

    end
end



plot(t_tseries[1,2,:], J_DE_CO2_tseries[1,2,:])
plot!(t_tseries[1,2,:], rG_CO2_tseries[1,2,:])


r_median        = zeros(39,83)
x_median        = zeros(39,83)
rG_CO2_median   = zeros(39,83)
rM_CO2_median   = zeros(39,83)
rX_CO2_median   = zeros(39,83)
J_EX_median     = zeros(39,83)
J_DE_median     = zeros(39,83)
J_DE_CO2_median = zeros(39,83)
J_D_median      = zeros(39,83)
J_ED_median     = zeros(39,83)
J_V_median      = zeros(39,83)
J_E_median      = zeros(39,83)
latency         = zeros(39,83)


for i in 1:39
    for j in 1:83
        try
            r_median[i,j]  = median(filter(!iszero, r_tseries[i,j,:]))
        catch
            r_median[i,j]  = NaN
        end
        try
            x_median[i,j]  = median(filter(!iszero, x_tseries[i,j,:]))
        catch
            x_median[i,j]  = NaN
        end
        try
            rG_CO2_median[i,j]  = median(filter(!iszero, rG_CO2_tseries[i,j,:]))
        catch
            rG_CO2_median[i,j]  = NaN
        end
        try
            rM_CO2_median[i,j]  = median(filter(!iszero, rM_CO2_tseries[i,j,:]))
        catch
            rM_CO2_median[i,j]  = NaN
        end
        try
            rX_CO2_median[i,j]  = median(filter(!iszero, rX_CO2_tseries[i,j,:]))
        catch
            rX_CO2_median[i,j]  = NaN
        end
        try
            J_EX_median[i,j]  = median(filter(!iszero, J_EX_tseries[i,j,:]))
        catch
            J_EX_median[i,j]  = NaN
        end
        try
            J_DE_median[i,j]  = median(filter(!iszero, J_DE_tseries[i,j,:]))
        catch
            J_DE_median[i,j]  = NaN
        end
        try
            J_DE_CO2_median[i,j]  = median(filter(!iszero, J_DE_CO2_tseries[i,j,:]))
        catch
            J_DE_CO2_median[i,j]  = NaN
        end
        try
            J_D_median[i,j]  = median(filter(!iszero, J_D_tseries[i,j,:]))
        catch
            J_D_median[i,j]  = NaN
        end
        try
            J_ED_median[i,j]  = median(filter(!iszero, J_ED_tseries[i,j,:]))
        catch
            J_ED_median[i,j]  = NaN
        end
        try
            J_V_median[i,j]  = median(filter(!iszero, J_V_tseries[i,j,:]))
        catch
            J_V_median[i,j]  = NaN
        end
        try
            J_E_median[i,j]  = median(filter(!iszero, J_E_tseries[i,j,:]))
        catch
            J_E_median[i,j]  = NaN
        end
        try
            latency[i,j]    = t_tseries[i,j,:][findmax(rG_CO2_tseries[i,j,:])[2]]
        catch
            latency[i,j]    = NaN
        end
    end
end

latency
latency[latency .>= 1000] .= NaN


dreserve_median = zeros(39,83)
mE_tseries = E_tseries./V_tseries

for i in 1:39
    for j in 1:83
        try
            dreserve_median[i,j]  = median(filter(!isnan, mE_tseries[i,j,:]))
        catch
            dreserve_median[i,j]  = NaN
        end
    end
end

# output fluxes
df_out         = DataFrame()
df_out.rgrowth = vec(r_median)
df_out.xenzyme = vec(x_median)
df_out.rGco2   = vec(rG_CO2_median)
df_out.rMco2   = vec(rM_CO2_median)
df_out.rGco2   = vec(rG_CO2_median)
df_out.rXco2   = vec(rX_CO2_median)
df_out.jEX     = vec(J_EX_median)
df_out.jDE     = vec(J_DE_median)
df_out.jDEco2  = vec(J_DE_CO2_median)
df_out.jD      = vec(J_D_median)
df_out.jED     = vec(J_ED_median)
df_out.jV      = vec(J_V_median)
df_out.jE      = vec(J_E_median)
df_out.dreserve = vec(dreserve_median)
df_out.response = repeat(df_isolates.Rhizosphere_response, 83)
df_out.isolate = repeat(df_isolates.Isolate, 83)  # species
ontology = Array{String}(undef,39,83)
for i in 1:83
     ontology[:,i] .= df_metabolites.Ontology[i]
 end
 df_out.ontology = vec(ontology)
monomer = Array{String}(undef,39,83)
for i in 1:83
    monomer[:,i] .= df_metabolites.Name[i]
 end
df_out.monomer = vec(monomer)
df_out.latency = vec(latency)
CSV.write(joinpath(dir, "files/output/isolates_batch_model_fluxes_latency1.csv"), df_out)

# output BGE-growth

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

BP_norm_tseries = BP_tseries./N_cells_tseries
BP_norm_tseries[BP_norm_tseries.<0.0] .=NaN

BP_median = zeros(39,83)
for i in 1:39
    for j in 1:83
        try
            BP_median[i,j]  = median(filter(!isnan, BP_norm_tseries[i,j,:]))
        catch
            BP_median[i,j]  = NaN
        end
    end
end

BR_norm_tseries = BR_tseries./N_cells_tseries
BR_norm_tseries[BR_norm_tseries.<0.0] .=NaN

BR_median = zeros(39,83)
for i in 1:39
    for j in 1:83
        try
            BR_median[i,j]  = median(filter(!isnan, BR_norm_tseries[i,j,:]))
        catch
            BR_median[i,j]  = NaN
        end
    end
end

r_median = zeros(39,83)
for i in 1:39
    for j in 1:83
        try
            r_median[i,j]  = median(filter(!iszero, r_tseries[i,j,:]))
        catch
            r_median[i,j]  = NaN
        end
    end
end

df_out_bge         = DataFrame()
df_out_bge.isolate = repeat(df_isolates.Isolate, 83)  # species
df_out_bge.class = repeat(df_isolates.Class, 83)
df_out_bge.phylum = repeat(df_isolates.Phylum, 83)
df_out_bge.response = repeat(df_isolates.Rhizosphere_response, 83)
df_out_bge.mingt = repeat(df_isolates.Min_gen_time, 83)
df_out_bge.rrn = repeat(df_isolates.rRNA_genes, 83)
df_out_bge.genomesize = repeat(df_isolates.Genome_size, 83)
monomer = Array{String}(undef,39,83)
for i in 1:83
     monomer[:,i] .= df_metabolites.Name[i]
 end
df_out_bge.monomer = vec(monomer)
ontology = Array{String}(undef,39,83)
for i in 1:83
     ontology[:,i] .= df_metabolites.Ontology[i]
 end
df_out_bge.ontology = vec(ontology)
df_out_bge.BGE = vec(BGE_median)
df_out_bge.rgrowth = vec(r_median)
df_out_bge.BP = vec(BP_median)
df_out_bge.BR = vec(BR_median)
CSV.write(joinpath(dir, "files/output/isolates_batch_model_BGE.csv"), df_out_bge)

# filter by substrate preference

df_all = DataFrame()
df_all.BGE = vec(BGE_median)
df_all.BP = vec(BP_median)
df_all.BR = vec(BR_median)
df_all.rgrowth = vec(r_median)
df_all.xenzyme = vec(x_median)
df_all.rGco2   = vec(rG_CO2_median)
df_all.rMco2   = vec(rM_CO2_median)
df_all.rGco2   = vec(rG_CO2_median)
df_all.rXco2   = vec(rX_CO2_median)
df_all.jEX     = vec(J_EX_median)
df_all.jDE     = vec(J_DE_median)
df_all.jDEco2  = vec(J_DE_CO2_median)
df_all.jD      = vec(J_D_median)
df_all.jED     = vec(J_ED_median)
df_all.jV      = vec(J_V_median)
df_all.jE      = vec(J_E_median)
df_all.dreserve = vec(dreserve_median)
df_all.response = repeat(df_isolates.Rhizosphere_response, 83)
df_all.isolate = repeat(df_isolates.Isolate, 83)  # species
df_all.class = repeat(df_isolates.Class, 83)
df_all.phylum = repeat(df_isolates.Phylum, 83)
df_all.response = repeat(df_isolates.Rhizosphere_response, 83)
df_all.mingt = repeat(df_isolates.Min_gen_time, 83)
df_all.rrn = repeat(df_isolates.rRNA_genes, 83)
df_all.genomesize = repeat(df_isolates.Genome_size, 83)
df_all.kM = repeat(maintenance["kM"], 83)
df_all.alphaenz = repeat(enzymes["alpha"], 83)
df_all.gV0 = repeat(turnover["gV0"], 83)
df_all.kE = repeat(protein_synthesis["kE"], 83)
df_all.yEV = repeat(protein_synthesis["yEV"], 83)
NSB = zeros(39,83)
for i in 1:83
    NSB[:,i] = assimilation["NSB"][i,:]
end
df_all.NSB = vec(NSB)
KD = zeros(39,83)
for i in 1:83
    KD[:,i] = assimilation["KD"][i,:]
end
df_all.KD = vec(KD)
yDE = zeros(39,83)
for i in 1:83
    yDE[:,i] = assimilation["yDE"][i,:]
end
df_all.yDE = vec(yDE)

# thermodynamic efficiency
dGcox              = zeros(83)
dGcat              = zeros(83)
dGAn               = zeros(83)
λ_base             = zeros(83)
N_C                = zeros(83)
eta                = zeros(83)
chemFormBiom       = [1, 1.64, 0.2, 0.38, 0, 0, 0]

for i in 1:83
    elementstring = convert(String,df_metabolites.Formula[i])
    N_C[i] = DEBmicroTrait.extract_composition(elementstring)[1]
    out = DEBmicroTrait.get_lambda(elementstring, chemFormBiom)
    dGcox[i] = out[2][3]
    dGcat[i] = out[2][5]
    dGAn[i]  = out[2][8]
    λ_base[i]     = out[1][1]
    eta      = @. dGAn/(λ_base*dGcat)
end


monomer = Array{String}(undef,39,83)
for i in 1:83
     monomer[:,i] .= df_metabolites.Name[i]
 end
df_all.monomer = vec(monomer)
ontology = Array{String}(undef,39,83)
for i in 1:83
     ontology[:,i] .= df_metabolites.Ontology[i]
 end
df_all.ontology = vec(ontology)
eta_eff = zeros(39,83)
for i in 1:83
     eta_eff[:,i] .= eta[i]
 end
df_all.eta = vec(eta_eff)

CSV.write(joinpath(dir, "files/output/isolates_batch_model_all.csv"), df_all)
