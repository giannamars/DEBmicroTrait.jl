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
Min_gen_time    = df_isolates.mingentime
Gram_stain      = repeat(["+"], 5369)
zh              =  df_isolates.total_carbohydrate_degradation./Genome_size*1e6
α_X             =  1e-2*(df_isolates.total_carbohydrate_degradation./Genome_size*1e6)./maximum(df_isolates.total_carbohydrate_degradation./Genome_size*1e6)
########################################

########################################
# I/O
save("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/output/dement_isolates_enzymes.jld", "zh", zh, "alpha", α_X)
########################################
