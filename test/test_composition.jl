using DEBmicroTrait
using CSV, DataFrames

df = CSV.read("/Users/glmarschmann/Data/BioFuels/Munson-McGee-SAG-Statistics_Class.csv", DataFrame, missingstring="N/A")

L_DNA = convert(Array{Float64,1}, df.final_assembly_length_bp)
V_cell = @. 4/3*pi*(df.estimated_diameter_mum*1e-6)^3
μ_0 = 4e7
β_B = 1.64
gmax = μ_0*V_cell.^(β_B-1)*60^2

stoich_Bacterioda = DEBmicroTrait.cell_volume_to_cell_stoichiometry(L_DNA, V_cell, gmax::Vector{Float64})

df.stoich = stoich_Bacterioda

CSV.write("/Users/glmarschmann/Data/BioFuels/Munson-McGee-SAG-Statistics_stoich.csv", df)