using CSV, DataFrames, Statistics

dir  = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()
# Load MAGs to traits data: granularity 1
file_path = joinpath(dir, "files/input/greenlon-H-tm1.csv") 
df_tm1 = CSV.read(file_path, DataFrame, missingstring="")
# Remove predators
predators = [
    "H1-18-all-fractions_metab_8",
    "H3-18-all-fractions_maxb_conc_78",
    "A3-18-all-fractions_metab_138",
    "A1-18-all-fractions_metab_maxbC_14"
]
filter!(row -> !(row.id in predators), df_tm1)

# Load metadata
metad = CSV.read(joinpath(dir, "files/input/Greenlon-et-al_2022_additional_supplemental_data.csv"), DataFrame, missingstring="")
rename!(metad, :me => :id)
df_tm1.id = strip.(string.(df_tm1.id))
metad.id  = strip.(string.(metad.id))
if length(unique(metad.id)) != nrow(metad)
    println("Warning: Metadata contains duplicate IDs. Deduplicating...")
    unique!(metad, :id)
end
metad_subset = select(metad, :id, ["bin.length", "afe.median", "rel.abund"])
merged_data = leftjoin(df_tm1, metad_subset, on = :id)

# Define columns
cols_binary = names(merged_data)[20:30]
cols_cont   = [names(merged_data)[2:19]; names(merged_data)[31:38]]
avg_size = mean(skipmissing(merged_data[:, "bin.length"]))
genome_sizes_mb = Float32.(coalesce.(merged_data[:, "bin.length"], avg_size) ./ 1e6)
# Convert binary columns to Float32 (0.0 or 1.0)
raw_binary_matrix = Matrix(merged_data[:, cols_binary])
X_binary_clean    = Float32.(coalesce.(raw_binary_matrix, 0) .> 0)
# Convert continuous columns to Float32, coalescing missing to 0
raw_cont_matrix = Matrix(merged_data[:, cols_cont])
X_cont_clean    = Float32.(coalesce.(raw_cont_matrix, 0))
# Normalize continuous columns by genome size
X_cont_norm = X_cont_clean ./ genome_sizes_mb
# Log-transform continuous columns
X_cont_log = Float32.(log.(X_cont_norm .+ 1e-6))
# Z-score normalization
μ = mean(X_cont_log, dims=1)
σ = std(X_cont_log, dims=1)
σ[σ .== 0] .= 1.0f0 # Prevent division by zero
X_cont_standardized = (X_cont_log .- μ) ./ σ

# Final format for Flux (features x samples)
global X_continuous = permutedims(X_cont_standardized) # Shape: (26 features, N genomes)
global X_binary     = permutedims(X_binary_clean)      # Shape: (11 features, N genomes)
global genome_sizes = genome_sizes_mb                      # Vector of length N

println("Preprocessing Complete.")
println("Continuous Shape: ", size(X_continuous))
println("Binary Shape:     ", size(X_binary))