using DEBmicroTrait
using CSV, DataFrames, Statistics
using JLD

########################################
# I/O
dir                     = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()
df_mags     = CSV.read(joinpath(dir, "files/input/greenlon-H-mags2traits.csv"), DataFrame, missingstring="")
########################################

########################################
genome_bp = df_mags."bin length" .* 1e6                     # bp
r_cell    = df_mags.spherical_equivalent_diameter ./ 2 .* 1e-6   # m
V_from_radius = passmissing(DEBmicroTrait.cell_radius_to_cell_volume).(r_cell)
V_from_genome = DEBmicroTrait.genome_size_to_cell_volume(genome_bp)
V_cell = coalesce.(V_from_radius, V_from_genome)
rrn_copies = DEBmicroTrait.genome_size_to_rRNA_copy_number(genome_bp)
Min_gen_time    = df_mags.mingentime
Gram_stain      = repeat(["-"], size(df_mags,1))

########################################
zh              =  df_mags.totalcazy./maximum(df_mags.totalcazy)
α_X             =  1e-2*zh
########################################
breakdown_cols = [
    "cellulose breakdown",
    "chitin breakdown",
    "heteromannan breakdown",
    "xylan and heteroxylan breakdown",
    "xyloglucan breakdown",
    "mixed linkage glucan breakdown",
    "protein degradation"
]

zh_class = vcat([reshape(df_mags[!, col], 1, :) for col in breakdown_cols]...)
f_αX = zh_class ./ sum(zh_class; dims = 1)   # 7 × 192 matrix

########################################
# I/O
save("/Users/glmarschmann/.julia/dev/DEBmicroTrait.jl/files/output/mags_enzymes.jld", "zh", zh, "alpha", α_X, "f_alphaX", f_αX)
########################################
