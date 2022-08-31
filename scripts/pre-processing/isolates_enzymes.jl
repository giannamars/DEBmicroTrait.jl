using DEBmicroTrait
using CSV, DataFrames, Statistics
using JLD

########################################
# I/O
df_isolates     = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates2traits.csv", DataFrame, missingstring="N/A")
########################################

########################################
V_cell          = DEBmicroTrait.genome_size_to_cell_volume(convert(Array{Float64,1}, df_isolates.Genome_size))
Min_gen_time    = df_isolates.Min_gen_time
Gram_stain      = convert(Array{String,1}, df_isolates.gram_stain)
zh              =  df_isolates.z_hydrolases./df_isolates.Genome_size*1e6
α_X             =  1e-2*(df_isolates.z_hydrolases./df_isolates.Genome_size*1e6)./maximum(df_isolates.z_hydrolases./df_isolates.Genome_size*1e6)
########################################

########################################
# I/O
save("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates_enzymes.jld", "zh", zh, "alpha", α_X)
########################################
