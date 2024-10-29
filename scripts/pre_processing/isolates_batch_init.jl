using DEBmicroTrait
using CSV, DataFrames, Statistics
using JLD

########################################
# I/O
dir             = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()
df_isolates     = CSV.read(joinpath(dir, "files/input/dement_isolates.csv"), DataFrame, missingstring="")
dropmissing!(df_isolates)
########################################
Genome_size     = convert(Array{Float64,1}, df_isolates.genome_length)
V_cell          = DEBmicroTrait.genome_size_to_cell_volume(Genome_size)
rrn_copies      = convert(Array{Float64,1}, df_isolates.rrncopies)
Min_gen_time    = df_isolates.mingentime
Gram_stain      = repeat(["+"], 5369)
########################################
dry_mass        = 0.47*DEBmicroTrait.cell_volume_to_dry_mass(V_cell, gmax, Gram_stain)
ρ_bulk          = 1.0 # g/cm^3
N_cells         = 1e3
Bio_0           = N_cells*1e6*ρ_bulk*dry_mass./12.011
########################################

########################################
# I/O
save("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/output/dement_isolates_batch_init.jld", "Bio0", Bio_0, "Md", dry_mass, "rhoB", ρ_bulk, "Ncells", N_cells)
########################################
