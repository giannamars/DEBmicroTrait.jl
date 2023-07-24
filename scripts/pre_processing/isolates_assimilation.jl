using DEBmicroTrait
using CSV, DataFrames, Statistics
using Roots
using JLD

########################################
# I/O
dir                     = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()
df_isolates             = CSV.read(joinpath(dir, "files/input/isolates2traits.csv"), DataFrame, missingstring="N/A")
df_metabolites          = CSV.read(joinpath(dir, "files/input/root_exudates.csv"), DataFrame, missingstring="N/A")
########################################
# metabolite traits
df_metabolites.Formula  = convert.(String, df_metabolites.Formula)
N_C                     = zeros(size(df_metabolites.Name,1))
for i in 1:size(df_metabolites.Name,1)
    elementstring       = df_metabolites.Formula[i]
    N_C[i]              = DEBmicroTrait.extract_composition(elementstring)[1]
end
########################################
# isolate traits
V_cell                  = DEBmicroTrait.genome_size_to_cell_volume(convert(Array{Float64,1}, df_isolates.Genome_size))
Min_gen_time            = df_isolates.Min_gen_time
Gram_stain              = convert(Array{String,1}, df_isolates.gram_stain)
rrn_copies              = convert(Array{Float64,1}, df_isolates.rRNA_genes)
y_EM                    = ones(size(V_cell,1))

z_sugars                = reshape(convert(Array{Float64,1}, df_isolates.z_sugars./df_isolates.Genome_size*1e6),1,39)
z_organics              = reshape(convert(Array{Float64,1}, df_isolates.z_organic_acids./df_isolates.Genome_size*1e6),1,39)
z_aminos                = reshape(convert(Array{Float64,1}, df_isolates.z_amino_acids./df_isolates.Genome_size*1e6), 1,39)
z_fattys                = reshape(convert(Array{Float64,1}, df_isolates.z_fatty_acids./df_isolates.Genome_size*1e6),1,39)
z_nucleos               = reshape(convert(Array{Float64,1}, df_isolates.z_nucleotides./df_isolates.Genome_size*1e6),1,39)
z_auxins                = reshape(convert(Array{Float64,1}, df_isolates.z_auxins./df_isolates.Genome_size*1e6),1,39)
genome_distr            = vcat(z_sugars, z_organics, z_aminos, z_fattys, z_nucleos, z_auxins)
########################################
# estimate transporter density
ρ_ps                    = zeros(size(df_metabolites.Name,1), size(V_cell,1))
y_DEs                   = zeros(size(df_metabolites.Name,1), size(V_cell,1))

for j in 1:size(df_metabolites.Name,1)
    if df_metabolites.Ontology[j] == "Sugars"
        for i in 1:size(V_cell,1)
            find_ρ(x)   = DEBmicroTrait.constrain_transporter_density_cost(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_metabolites.Formula[j])
            ρ_p         = Roots.find_zero(find_ρ, 1.0)
            closure     = genome_distr[:,i]./sum(genome_distr[:,i])
            ρ_ps[j,i]   = ρ_p[1].*closure[1]
            y_DE        = DEBmicroTrait.yield_transporter_density_cost(ρ_ps[j,i], [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_metabolites.Formula[j])
            y_DEs[j,i]  = y_DE[1]
        end
    elseif df_metabolites.Ontology[j] == "Organic acids"
        for i in 1:size(V_cell,1)
            find_ρ(x)   = DEBmicroTrait.constrain_transporter_density_cost(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_metabolites.Formula[j])
            ρ_p         = Roots.find_zero(find_ρ, 1.0)
            closure     = genome_distr[:,i]./sum(genome_distr[:,i])
            ρ_ps[j,i]   = ρ_p[1].*closure[2]
            y_DE        = DEBmicroTrait.yield_transporter_density_cost(ρ_ps[j,i], [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_metabolites.Formula[j])
            y_DEs[j,i]  = y_DE[1]
        end
    elseif df_metabolites.Ontology[j] == "Amino acids"
        for i in 1:size(V_cell,1)
            find_ρ(x)   = DEBmicroTrait.constrain_transporter_density_cost(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_metabolites.Formula[j])
            ρ_p         = Roots.find_zero(find_ρ, 1.0)
            closure     = genome_distr[:,i]./sum(genome_distr[:,i])
            ρ_ps[j,i]   = ρ_p[1].*closure[3]
            y_DE        = DEBmicroTrait.yield_transporter_density_cost(ρ_ps[j,i], [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_metabolites.Formula[j])
            y_DEs[j,i]  = y_DE[1]
        end
    elseif df_metabolites.Ontology[j] == "Fatty acids"
        for i in 1:size(V_cell,1)
            find_ρ(x)   = DEBmicroTrait.constrain_transporter_density_cost(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_metabolites.Formula[j])
            ρ_p         = Roots.find_zero(find_ρ, 1.0)
            closure     = genome_distr[:,i]./sum(genome_distr[:,i])
            ρ_ps[j,i]   = ρ_p[1].*closure[4]
            y_DE        = DEBmicroTrait.yield_transporter_density_cost(ρ_ps[j,i], [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_metabolites.Formula[j])
            y_DEs[j,i]  = y_DE[1]
        end
    elseif df_metabolites.Ontology[j] == "Nucleotides"
        for i in 1:size(V_cell,1)
            find_ρ(x)   = DEBmicroTrait.constrain_transporter_density_cost(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_metabolites.Formula[j])
            ρ_p         = Roots.find_zero(find_ρ, 1.0)
            closure     = genome_distr[:,i]./sum(genome_distr[:,i])
            ρ_ps[j,i]   = ρ_p[1].*closure[5]
            y_DE        = DEBmicroTrait.yield_transporter_density_cost(ρ_ps[j,i], [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_metabolites.Formula[j])
            y_DEs[j,i]  = y_DE[1]
        end
    else
        for i in 1:size(V_cell,1)
            find_ρ(x)   = DEBmicroTrait.constrain_transporter_density_cost(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_metabolites.Formula[j])
            ρ_p         = Roots.find_zero(find_ρ, 1.0)
            closure     = genome_distr[:,i]./sum(genome_distr[:,i])
            ρ_ps[j,i]   = ρ_p[1].*closure[6]
            y_DE        = DEBmicroTrait.yield_transporter_density_cost(ρ_ps[j,i], [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_metabolites.Formula[j])
            y_DEs[j,i]  = y_DE[1]
        end
    end
end

ρ_ps[ρ_ps.==0.0] .= 1e-12
median(ρ_ps)

########################################
# calculate uptake traits
N_SB              = DEBmicroTrait.transporter_density_to_monomer_uptake_sites(V_cell, ρ_ps, Min_gen_time, Gram_stain)
Vmax              = @. 180.0*60^2*N_SB.*N_C
#
D_S               = DEBmicroTrait.aqueous_diffusivity(df_metabolites.Molecular_weight)
K_D               = DEBmicroTrait.specific_reference_affinity(V_cell, ρ_ps, D_S)
#
a_s               = Vmax./K_D

########################################
# I/O
save(joinpath(dir, "files/output/isolates_assimilation.jld"), "rho", ρ_ps, "NSB", N_SB, "KD", K_D, "yEM", y_EM, "yDE", y_DEs, "NC", N_C)