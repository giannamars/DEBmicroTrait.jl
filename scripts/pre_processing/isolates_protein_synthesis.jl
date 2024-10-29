using DEBmicroTrait
using CSV, DataFrames, Statistics
using JLD

########################################
# I/O
dir             = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()
df_isolates     = CSV.read(joinpath(dir, "files/input/dement_isolates.csv"), DataFrame, missingstring="")
dropmissing!(df_isolates)
########################################

########################################
Genome_size     = convert(Array{Float64,1}, df_isolates.genome_length)
V_cell          = DEBmicroTrait.genome_size_to_cell_volume(Genome_size)
rrn_copies      = convert(Array{Float64,1}, df_isolates.rrncopies)
Min_gen_time    = df_isolates.mingentime
Gram_stain      = repeat(["+"], 5369)
########################################

########################################
gmax            = log(2)./Min_gen_time
V_p             = DEBmicroTrait.cell_volume_to_protein_volume(V_cell)
V_r             = DEBmicroTrait.cell_volume_to_ribosome_volume(V_cell, gmax)
k_E             = DEBmicroTrait.translation_power(V_p, V_r, Min_gen_time)
########################################

########################################
y_EV            = DEBmicroTrait.relative_translation_efficiency_regression(rrn_copies)
########################################

########################################
# I/O
save("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/output/dement_isolates_protein_synthesis.jld", "kE", k_E, "yEV", y_EV, "mingt", Min_gen_time)
########################################
