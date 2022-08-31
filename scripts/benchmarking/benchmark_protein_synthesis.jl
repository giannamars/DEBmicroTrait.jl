using DEBmicroTrait
using CSV, DataFrames, Statistics
using GLM

################################################
# I/O
df = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/dethlefsen_schmidt_2008.csv", DataFrame, missingstring="N/A")
gdf = groupby(df, :Species_Name)
df = combine(gdf, :Protein_Quantity_Measurement => mean, :Size_Measurement => mean, :DNA_Quantity_Measurement => mean, :RNA_Quantity_Measurement => mean, :Specific_Growth_Rate_Measurement => mean, )
df_out = DataFrame()
################################################

################################################
# Dethlefsen & Schmidt (2008): translational_power = μ*P/R
# Protein data
d_p = 1.37e6
P = df.Protein_Quantity_Measurement_mean
V_P = P/d_p
df_out.V_P = V_P
# RNA data
d_r = 1.79e6
R = df.RNA_Quantity_Measurement_mean
V_R = R/d_r
df_out.V_R = V_R
# Growth rate data
gmax = df.Specific_Growth_Rate_Measurement_mean
df_out.gmax = gmax
# Cell size data
V_cell = df.Size_Measurement_mean
df_out.V_cell = V_cell
# k_E calculation
k_E = gmax.*V_P./V_R
df_out.k_E = k_E
################################################

################################################
#Kempes et al. (2016): allometric predictions
# Cell size prediction
d_DNA = 2e6
DNA = df.DNA_Quantity_Measurement_mean
V_DNA = DNA/d_DNA
v_N = 1.47e-27
L_DNA = V_DNA/v_N
V_cell_model = DEBmicroTrait.genome_size_to_cell_volume(L_DNA)
df_out.V_cell_model = V_cell_model
# Protein prediction
V_P_model = DEBmicroTrait.cell_volume_to_protein_volume(V_cell)
df_out.V_P = V_P
df_out.V_P_model = V_P_model
# RNA prediction
V_R_model = DEBmicroTrait.cell_volume_to_ribosome_volume(V_cell, gmax)
df_out.V_R = V_R
df_out.V_R_model = V_R_model
# tRNA predictions
V_tRNA_model = DEBmicroTrait.cell_volume_to_tRNA_volume(V_cell, gmax)
df_out.tRNA_model = V_tRNA_model
# mRNA predictions
V_mRNA_model = DEBmicroTrait.cell_volume_to_mRNA_volume(V_cell, gmax)
df_out.mRNA_model = V_mRNA_model
# total RNA
V_mtR_model = V_R_model .+ V_tRNA_model
df_out.V_mtR_model = V_mtR_model
# k_E prediction
k_E_model= DEBmicroTrait.translation_power(V_P_model, V_R_model, log(2)./gmax)
df_out.k_E_model = k_E_model
################################################

################################################
# Statistics
linearRegressor = lm(@formula(log(k_E) ~ log(k_E_model)), df_out)
r2(linearRegressor)
# test for significant difference to 1-1 line
df_out.diff = k_E_model-k_E
linearRegressor = lm(@formula(diff ~ k_E_model), df_out)
# RMSD
RMSD = sqrt((1/(12-1)*sum((k_E_model - k_E).^2)))
################################################

################################################
# I/O
CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/manuscript/files/benchmark_protein_synthesis.csv", df_out)
################################################
