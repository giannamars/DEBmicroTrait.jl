using DEBmicroTrait
using CSV, DataFrames, Statistics
using GLM
using Plots

df          = CSV.read("/Users/glmarschmann/Desktop/marine_sizes.csv", DataFrame, missingstring="N/A")
Genome_size = df.mean_estimated_genome_length_bp
V_cell = DEBmicroTrait.genome_size_to_cell_volume(Genome_size).*1e18
V_cell_meas = df.spherical_cell_volume_mum3
df.spherical_cell_volume_pred = V_cell

# Statistics
linearRegressor = lm(@formula(log(spherical_cell_volume_mum3) ~ log(spherical_cell_volume_pred)), df)
r2(linearRegressor)

# Statistics
df.genome_volume = df.mean_estimated_genome_length_bp*1.47e-27*1e18
linearRegressor = lm(@formula(log(spherical_cell_volume_mum3) ~ log(genome_volume)), df)
r2(linearRegressor)


# test for significant difference to 1-1 line
df.diff = df.spherical_cell_volume_pred-df.spherical_cell_volume_mum3
linearRegressor = lm(@formula(diff ~ spherical_cell_volume_pred), df)
# RMSD
RMSD = sqrt((1/(304-1)*sum((df.spherical_cell_volume_pred - df.spherical_cell_volume_mum3).^2)))

scatter(df.spherical_cell_volume_mum3, df.spherical_cell_volume_pred, xscale = :log, yscale=:log)

CSV.write("/Users/glmarschmann/Desktop/marine_sizes_pred.csv", df)

# filter
smallestbacterium = 1e-21*1e18
largestbacterium = 1e-18*1e18

df_filtered = df[(df.spherical_cell_volume_mum3 .> smallestbacterium) .& (df.spherical_cell_volume_mum3 .< largestbacterium) .& (df.spherical_cell_volume_pred .> smallestbacterium) .& (df.spherical_cell_volume_pred .< largestbacterium), :]
df_filtered.genome_volume = df_filtered.mean_estimated_genome_length_bp*1.47e-27*1e18
# Statistics
linearRegressor = lm(@formula(log(spherical_cell_volume_mum3) ~ log(spherical_cell_volume_pred)), df_filtered)
r2(linearRegressor)

linearRegressor = lm(@formula(log(spherical_cell_volume_mum3) ~ log(genome_volume)), df_filtered)
r2(linearRegressor)

# test for significant difference to 1-1 line
df_filtered.diff = df_filtered.spherical_cell_volume_pred-df_filtered.spherical_cell_volume_mum3
linearRegressor = lm(@formula(diff ~ spherical_cell_volume_pred), df_filtered)
# RMSD
RMSD = sqrt((1/(166-1)*sum((df_filtered.spherical_cell_volume_pred - df_filtered.spherical_cell_volume_mum3).^2)))

gd = groupby(df_filtered, "genus")
df_family = combine(gd, [:spherical_cell_volume_pred, :spherical_cell_volume_mum3, :genome_volume] .=> median)

# Statistics
linearRegressor = lm(@formula(log(genome_volume) ~ log(spherical_cell_volume_mum3)), df_filtered)
r2(linearRegressor)
