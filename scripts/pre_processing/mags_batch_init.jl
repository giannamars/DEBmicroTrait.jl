using DEBmicroTrait
using CSV, DataFrames, Statistics
using JLD

########################################
# I/O
dir                     = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()
df_mags     = CSV.read(joinpath(dir, "files/input/greenlon-H-mags2traits.csv"), DataFrame, missingstring="")
########################################

########################################
genome_bp = df_mags."bin length" .* 1e6                   # bp
r_cell    = df_mags.spherical_equivalent_diameter ./ 2 .* 1e-6   # m
V_from_radius = passmissing(DEBmicroTrait.cell_radius_to_cell_volume).(r_cell)
V_from_genome = DEBmicroTrait.genome_size_to_cell_volume(genome_bp)
V_cell = coalesce.(V_from_radius, V_from_genome)
rrn_copies = DEBmicroTrait.genome_size_to_rRNA_copy_number(genome_bp)
Min_gen_time    = df_mags.mingentime
Gram_stain      = repeat(["-"], size(df_mags,1))
########################################
dry_mass        = 0.47*DEBmicroTrait.cell_volume_to_dry_mass(V_cell, gmax, Gram_stain)
ρ_bulk          = 1.0 # g/cm^3
N_cells         = 1e3
Bio_0           = N_cells*1e6*ρ_bulk*dry_mass./12.011
########################################

########################################
# I/O
save("/Users/glmarschmann/.julia/dev/DEBmicroTrait.jl/files/output/mags_batch_init.jld", "Bio0", Bio_0, "Md", dry_mass, "rhoB", ρ_bulk, "Ncells", N_cells)
########################################
