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

########################################
k_M    = DEBmicroTrait.cell_volume_to_specific_maintenance_rate(V_cell, Min_gen_time, Gram_stain)
k_M_med = median(k_M)
########################################

########################################
# I/O
save("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/output/dement_isolates_maintenance.jld", "kM", k_M)
