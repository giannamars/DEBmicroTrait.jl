using DEBmicroTrait
using CSV, DataFrames, Statistics
using JLD
using DifferentialEquations

df_isolates     = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates2traits.csv", DataFrame, missingstring="N/A")
df_metabolites  = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/exudation_properties.csv", DataFrame, missingstring="N/A")

assimilation      = load("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates_assimilation.jld")
enzymes           = load("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates_enzymes.jld")
maintenance       = load("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates_maintenance.jld")
protein_synthesis = load("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates_protein_synthesis.jld")
turnover          = load("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates_turnover.jld")
initb             = load("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates_batch_init.jld")


condition(u,t,integrator) = u[1] - 1e-5
affect!(integrator)       = terminate!(integrator)
cb                        = ContinuousCallback(condition,affect!)

BGE_tseries       = zeros(39, 83, 500)
BR_tseries        = zeros(39, 83, 500)
BP_tseries        = zeros(39, 83, 500)
r_tseries         = zeros(39, 83, 500)
x_tseries         = zeros(39, 83, 500)
rG_CO2_tseries    = zeros(39, 83, 500)
rM_CO2_tseries    = zeros(39, 83, 500)
rX_CO2_tseries    = zeros(39, 83, 500)
J_EX_tseries      = zeros(39, 83, 500)
J_DE_tseries      = zeros(39, 83, 500)
J_DE_CO2_tseries  = zeros(39, 83, 500)
J_D_tseries       = zeros(39, 83, 500)
J_ED_tseries      = zeros(39, 83, 500)
J_V_tseries       = zeros(39, 83, 500)
J_E_tseries       = zeros(39, 83, 500)
t_tseries         = zeros(39, 83, 500)
D_tseries         = zeros(39, 83, 500)
E_tseries         = zeros(39, 83, 500)
V_tseries         = zeros(39, 83, 500)
X_tseries         = zeros(39, 83, 500)
CO2_tseries       = zeros(39, 83, 500)
N_cells_tseries   = zeros(39, 83, 500)
maintenance_tseries    = zeros(39, 83, 500)

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
        sol               = solve(prob, alg_hints=[:stiff], callback=cb)

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
    end
end


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
df_out.rgrowth = vec(r_median)    # growth rate
df_out.xenzyme = vec(x_median)    # enzyme rate
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
CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates_batch_model_fluxes.csv", df_out)

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
CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates_batch_model_BGE.csv", df_out_bge)
#############################################################################################################
# output training data

df_out = DataFrame()
df_out.BGE = vec(BGE_median)
df_out.rgrowth = vec(r_median)
#df_out.genomesize = repeat(df_isolates.Genome_size, 83)/1e6
#df_out.mingt = repeat(df_isolates.Min_gen_time, 83)
#df_out.rrn = repeat(df_isolates.rRNA_genes, 83)
#df_out.gramstain = repeat(df_isolates.gram_stain, 83)
#df_out.gramstain[df_out.gramstain .== "(+)"] .= "1"
#df_out.gramstain[df_out.gramstain .== "(-)"] .= "2"
df_out.kE = repeat(protein_synthesis["kE"], 83)
df_out.yEV = repeat(protein_synthesis["yEV"], 83)
df_out.aX = repeat(enzymes["alpha"], 83)
#df_out.gV = turnover["gV1"]*ones(3237)
V_cell = DEBmicroTrait.genome_size_to_cell_volume(convert(Array{Float64,1}, df_isolates.Genome_size))
SAV = DEBmicroTrait.surface_area_volume_ratio(V_cell)
#df_out.sav = repeat(SAV./1e6, 83)
yields = Array{Float64}(undef,39,83)
for i in 1:83
    yields[:,i] .= assimilation["yDE"][i,:]
end
df_out.yield = vec(yields)
rhos = Array{Float64}(undef,39,83)
for i in 1:83
    rhos[:,i] .= assimilation["rho"][i,:]
end
#df_out.rho = vec(rhos)

Vmaxs = Array{Float64}(undef,39,83)
for i in 1:83
 Vmaxs[:,i] .= 180.0*60^2*assimilation["NSB"][i,:]
end
df_out.Vmax = vec(Vmaxs)

Ks = Array{Float64}(undef,39,83)
for i in 1:83
    Ks[:,i] .= assimilation["KD"][i,:]
end
df_out.KD = vec(Ks)
df_out.BGE = vec(BGE_median)
df_out.rgrowth = vec(r_median)
df_out.response = repeat(df_isolates.Rhizosphere_response, 83)

df_out_pn = filter(x->(x.response.=="negative" || x.response.=="positive") , df_out)
df_out_u = filter(x->(x.response.=="undefined") , df_out)
CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates_batch_model_train.csv", df_out_pn)
CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates_batch_model_test.csv", df_out_u)

#############################################################################################################
df_out_levin = DataFrame(assimilation_median, :auto)
df_out_levin.response = df_isolates.Rhizosphere_response
CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/first_manuscript/files_pub/isolates_batch_model_levin.csv", df_out_levin)
