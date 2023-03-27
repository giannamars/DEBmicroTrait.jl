using DEBmicroTrait
using CSV, DataFrames, Statistics
using Plots

dir                     = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()
df_isolates             = CSV.read(joinpath(dir, "files/input/isolates2traits.csv"), DataFrame, missingstring="N/A")

# HB15
####################################################################################
df_HB15                 = df_isolates[df_isolates.Abbreviation .== "HB15", :]
Genome_size_HB15        = convert(Array{Float64,1}, df_HB15.Genome_size)
V_cell_HB15             = DEBmicroTrait.genome_size_to_cell_volume(Genome_size)
r_cell_HB15             = (3*V_cell_HB15/(4*pi)).^(1/3)

df_live_R2A_HB15 = CSV.read(joinpath(dir, "files/Stanford_cell_size_data/Live_R2A_15_snap_1.csv"), DataFrame, missingstring="NaN")
df_live_LB_HB15 = CSV.read(joinpath(dir, "files/Stanford_cell_size_data/Live_LB_15_snap_1.csv"), DataFrame, missingstring="NaN")
df_fixed_R2A_HB15 = CSV.read(joinpath(dir, "files/Stanford_cell_size_data/Fixed_R2A_15_snap_1.csv"), DataFrame, missingstring="NaN")
df_fixed_LB_HB15 = CSV.read(joinpath(dir, "files/Stanford_cell_size_data/Fixed_LB_15_snap_1.csv"), DataFrame, missingstring="NaN")

histogram(df_live_R2A_HB15.width)
median(skipmissing(df_fixed_R2A_HB15.width))


V_cell_live_R2A_HB15 = @. (4/3)*(df_live_R2A_HB15.width)^2*(df_live_R2A_HB15.len)
V_cell_fixed_R2A_HB15 = @. (4/3)*(df_fixed_R2A_HB15.width)^2*(df_fixed_R2A_HB15.len)


df_HA54                 = df_isolates[df_isolates.Abbreviation .== "HA54", :]
Genome_size_HA54        = convert(Array{Float64,1}, df_HA54.Genome_size)
V_cell_HA54             = DEBmicroTrait.genome_size_to_cell_volume(Genome_size_HA54)
r_cell_HA54            = (3*V_cell_HA54/(4*pi)).^(1/3)


df_live_R2A_HA54 = CSV.read(joinpath(dir, "files/Stanford_cell_size_data/Live_R2A_54_snap_1.csv"), DataFrame, missingstring="NaN")
df_live_LB_HA54 = CSV.read(joinpath(dir, "files/Stanford_cell_size_data/Live_LB_54_snap_1.csv"), DataFrame, missingstring="NaN")
df_fixed_R2A_HA54 = CSV.read(joinpath(dir, "files/Stanford_cell_size_data/Fixed_R2A_54_snap_1.csv"), DataFrame, missingstring="NaN")
df_fixed_LB_HA54 = CSV.read(joinpath(dir, "files/Stanford_cell_size_data/Fixed_LB_54_snap_1.csv"), DataFrame, missingstring="NaN")

V_cell_live_R2A_HA54 = @. (4/3)*(df_live_R2A_HA54.width)^2*(df_live_R2A_HA54.len)
V_cell_fixed_R2A_HA54 = @. (4/3)*(df_fixed_R2A_HA54.width)^2*(df_fixed_R2A_HA54.len)

# HB15
####################################################################################
df_HA57                 = df_isolates[df_isolates.Abbreviation .== "HA57", :]
Genome_size_HA57        = convert(Array{Float64,1}, df_HA57.Genome_size)
V_cell_HA57             = DEBmicroTrait.genome_size_to_cell_volume(Genome_size_HA57)
r_cell_HA57             = (3*V_cell_HA57/(4*pi)).^(1/3)

df_live_R2A_HA57 = CSV.read(joinpath(dir, "files/Stanford_cell_size_data/Live_R2A_57_snap_1.csv"), DataFrame, missingstring="NaN")
df_live_LB_HA57 = CSV.read(joinpath(dir, "files/Stanford_cell_size_data/Live_LB_57_snap_1.csv"), DataFrame, missingstring="NaN")
df_fixed_R2A_HA57 = CSV.read(joinpath(dir, "files/Stanford_cell_size_data/Fixed_R2A_57_snap_1.csv"), DataFrame, missingstring="NaN")
df_fixed_LB_HA57 = CSV.read(joinpath(dir, "files/Stanford_cell_size_data/Fixed_LB_57_snap_1.csv"), DataFrame, missingstring="NaN")

V_cell_live_R2A_HA57 = @. (4/3)*(df_live_R2A_HA57.width)^2*(df_live_R2A_HA57.len)
V_cell_fixed_R2A_HA57 = @. (4/3)*(df_fixed_R2A_HA57.width)^2*(df_fixed_R2A_HA57.len)

median(skipmissing(df_live_R2A_HA57.width))





