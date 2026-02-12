using DEBmicroTrait
using CSV, DataFrames, Statistics
using Roots
using JLD

########################################
# I/O
dir                     = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()
df_mags                 = CSV.read(joinpath(dir, "files/input/greenlon-H-mags2traits.csv"), DataFrame, missingstring="")
df_metabolites          = CSV.read(joinpath(dir, "files/input/greenlon-H-monomers.csv"), DataFrame, missingstring="")
########################################
# metabolite traits
df_metabolites.Formula  = convert.(String, df_metabolites.Formula)
N_C                     = zeros(size(df_metabolites.Name,1))
for i in 1:size(df_metabolites.Name,1)
    elementstring       = df_metabolites.Formula[i]
    N_C[i]              = DEBmicroTrait.extract_composition(elementstring)[1]
end
########################################
# mag traits
genome_bp = df_mags."bin length" .* 1e6                   # bp
r_cell    = df_mags.spherical_equivalent_diameter ./ 2 .* 1e-6   # m
V_from_radius = passmissing(DEBmicroTrait.cell_radius_to_cell_volume).(r_cell)
V_from_genome = DEBmicroTrait.genome_size_to_cell_volume(genome_bp)
V_cell = coalesce.(V_from_radius, V_from_genome)
rrn_copies = DEBmicroTrait.genome_size_to_rRNA_copy_number(genome_bp)
Min_gen_time    = df_mags.mingentime
Gram_stain      = repeat(["-"], size(df_mags,1))
y_EM                    = ones(size(V_cell,1))

z_sugars                = reshape(convert(Array{Float64,1}, df_mags."monosaccharide transport"),1,192)
z_aminos                = reshape(convert(Array{Float64,1}, df_mags."free amino acids transport"), 1,192)
genome_distr            = vcat(z_sugars, z_aminos)
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
    else
        for i in 1:size(V_cell,1)
            find_ρ(x)   = DEBmicroTrait.constrain_transporter_density_cost(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_metabolites.Formula[j])
            ρ_p         = Roots.find_zero(find_ρ, 1.0)
            closure     = genome_distr[:,i]./sum(genome_distr[:,i])
            ρ_ps[j,i]   = ρ_p[1].*closure[2]
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
save("/Users/GLMarschmann/.julia/dev/DEBmicroTrait.jl/files/output/mags_assimilation.jld", "rho", ρ_ps, "NSB", N_SB, "KD", K_D, "yEM", y_EM, "yDE", y_DEs, "NC", N_C)
