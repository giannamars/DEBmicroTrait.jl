using DEBmicroTrait
using CSV, DataFrames, Statistics
using GLM
using Roots

###########################################################
# Flavobacterium johnsoniae
Genome_size = [6.1e6]
V_cell = DEBmicroTrait.genome_size_to_cell_volume(Genome_size)
Min_gen_time = [log(2)./0.2]
Gram_stain = ["-"]
rrn_copies = [6.0]
y_EM = [1.0]
# Glucose
elementstr = "C6H12O6"
Molecular_weight = [180.16]

chemFormBiom       = [1, 1.64, 0.2, 0.38, 0, 0, 0]
N_C = DEBmicroTrait.extract_composition(elementstr)[1]
out = DEBmicroTrait.get_lambda(elementstr, chemFormBiom)
yield = out[4][10]/out[4][1]
y_DE = @. 1/abs(yield)

find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, V_cell, Min_gen_time, Gram_stain, rrn_copies, y_DE, N_C, y_EM)
ρ_p = Roots.find_zero(find_ρ, 1.0)
D_S               = DEBmicroTrait.aqueous_diffusivity(Molecular_weight)
K_D               = DEBmicroTrait.specific_reference_affinity(V_cell, ρ_p*ones(1,1), D_S)*1e3

###########################################################
# Escherichia coli
Genome_size = [4.64e6]
V_cell = DEBmicroTrait.genome_size_to_cell_volume(Genome_size)
Min_gen_time = [log(2)./1.23]
Gram_stain = ["-"]
rrn_copies = [7.0]
y_EM = [1.0]

# Glucose
elementstr = "C6H12O6"
Molecular_weight = [180.16]

chemFormBiom       = [1, 1.64, 0.2, 0.38, 0, 0, 0]
N_C = DEBmicroTrait.extract_composition(elementstr)[1]
out = DEBmicroTrait.get_lambda(elementstr, chemFormBiom)
yield = out[4][10]/out[4][1]
y_DE = @. 1/abs(yield)

find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, V_cell, Min_gen_time, Gram_stain, rrn_copies, y_DE, N_C, y_EM)
ρ_p = Roots.find_zero(find_ρ, 1.0)
D_S               = DEBmicroTrait.aqueous_diffusivity(Molecular_weight)
K_D               = DEBmicroTrait.specific_reference_affinity(V_cell, ρ_p*ones(1,1), D_S)*1e3

###########################################################
# Lacticaseibacillus casei
# Lactose
Genome_size = [2.95e6]
V_cell = DEBmicroTrait.genome_size_to_cell_volume(Genome_size)
rrn_copies = [5.0]
Gram_stain = ["+"]
gmax = DEBmicroTrait.gmax_regression(rrn_copies)
V_p = DEBmicroTrait.cell_volume_to_protein_volume(V_cell)
d_p = 1.37e6
mp = d_p*V_p*12.01
Vmax = 117*1e-9./(1e3*mp)/60
Min_gen_time = [log(2)./Vmax[1]]
y_EM = [1.0]

elementstr = "C12H22O11"
Molecular_weight = [342.3]

chemFormBiom       = [1, 1.64, 0.2, 0.38, 0, 0, 0]
N_C = DEBmicroTrait.extract_composition(elementstr)[1]
out = DEBmicroTrait.get_lambda(elementstr, chemFormBiom)
yield = out[4][10]/out[4][1]
y_DE = @. 1/abs(yield)

find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, V_cell, Min_gen_time, Gram_stain, rrn_copies, y_DE, N_C, y_EM)
ρ_p = Roots.find_zero(find_ρ, 1.0)
D_S               = DEBmicroTrait.aqueous_diffusivity(Molecular_weight)
K_D               = DEBmicroTrait.specific_reference_affinity(V_cell, ρ_p*ones(1,1), D_S)*1e3

###########################################################
# Corynebacterium
Genome_size = [2.47e6]
V_cell = DEBmicroTrait.genome_size_to_cell_volume(Genome_size)
Min_gen_time = [log(2)./0.15]
Gram_stain = ["+"]
rrn_copies = [4.0]
y_EM = [1.0]

# Glucose
elementstr = "C6H12O6"
Molecular_weight = [180.16]

chemFormBiom       = [1, 1.64, 0.2, 0.38, 0, 0, 0]
N_C = DEBmicroTrait.extract_composition(elementstr)[1]
out = DEBmicroTrait.get_lambda(elementstr, chemFormBiom)
yield = out[4][10]/out[4][1]
y_DE = @. 1/abs(yield)

find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, V_cell, Min_gen_time, Gram_stain, rrn_copies, y_DE, N_C, y_EM)
ρ_p = Roots.find_zero(find_ρ, 1.0)
N_SB              = DEBmicroTrait.transporter_density_to_monomer_uptake_sites(V_cell, ρ_p*ones(1,1), Min_gen_time, Gram_stain)
Vmax              = @. 180.0*60^2*N_SB.*N_C
D_S               = DEBmicroTrait.aqueous_diffusivity(Molecular_weight)
K_D               = DEBmicroTrait.specific_reference_affinity(V_cell, ρ_p*ones(1,1), D_S)*1e3

###########################################################
# Thiobacillus
Genome_size = [2.85e6]
V_cell = DEBmicroTrait.genome_size_to_cell_volume(Genome_size)
Gram_stain = ["-"]
rrn_copies = [2.0]
gmax = DEBmicroTrait.gmax_regression(rrn_copies)
dw = DEBmicroTrait.cell_volume_to_dry_mass(V_cell, gmax, Gram_stain)*12.01
Vmax = 74.8*1e-9./(1e3*dw)/60
Min_gen_time = [log(2)./Vmax[1]]
y_EM = [1.0]

# Fructose
elementstr = "C6H12O6"
Molecular_weight = [180.16]

chemFormBiom       = [1, 1.64, 0.2, 0.38, 0, 0, 0]
N_C = DEBmicroTrait.extract_composition(elementstr)[1]
out = DEBmicroTrait.get_lambda(elementstr, chemFormBiom)
yield = out[4][10]/out[4][1]
y_DE = @. 1/abs(yield)

find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, V_cell, Min_gen_time, Gram_stain, rrn_copies, y_DE, N_C, y_EM)
ρ_p = Roots.find_zero(find_ρ, 1.0)
D_S               = DEBmicroTrait.aqueous_diffusivity(Molecular_weight)
K_D               = DEBmicroTrait.specific_reference_affinity(V_cell, ρ_p*ones(1,1), D_S)*1e3

# Ribose
elementstr = "C5H10O5"
Molecular_weight = [150.13]

Vmax = 30.8*1e-9./(1e3*dw)/60
Min_gen_time = [log(2)./Vmax[1]]
chemFormBiom       = [1, 1.64, 0.2, 0.38, 0, 0, 0]
N_C = DEBmicroTrait.extract_composition(elementstr)[1]
out = DEBmicroTrait.get_lambda(elementstr, chemFormBiom)
yield = out[4][10]/out[4][1]
y_DE = @. 1/abs(yield)

find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, V_cell, Min_gen_time, Gram_stain, rrn_copies, y_DE, N_C, y_EM)
ρ_p = Roots.find_zero(find_ρ, 1.0)
D_S               = DEBmicroTrait.aqueous_diffusivity(Molecular_weight)
K_D               = DEBmicroTrait.specific_reference_affinity(V_cell, ρ_p*ones(1,1), D_S)*1e3


###########################################################
# Brevibacterium linens
# Tyrosine
Genome_size = [3.81e6]
V_cell = DEBmicroTrait.genome_size_to_cell_volume(Genome_size)
rrn_copies = [4.0]
Gram_stain = ["+"]
gmax = DEBmicroTrait.gmax_regression(rrn_copies)
dw = DEBmicroTrait.cell_volume_to_dry_mass(V_cell, gmax, Gram_stain)*12.01
Vmax = 0.3*1e-9./dw/60
Min_gen_time = [log(2)./Vmax[1]]
y_EM = [1.0]

elementstr = "C9H11NO3"
Molecular_weight = [181.191]

chemFormBiom       = [1, 1.64, 0.2, 0.38, 0, 0, 0]
N_C = DEBmicroTrait.extract_composition(elementstr)[1]
out = DEBmicroTrait.get_lambda(elementstr, chemFormBiom)
yield = out[4][10]/out[4][1]
y_DE = @. 1/abs(yield)

find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, V_cell, Min_gen_time, Gram_stain, rrn_copies, y_DE, N_C, y_EM)
ρ_p = Roots.find_zero(find_ρ, 1.0)
D_S               = DEBmicroTrait.aqueous_diffusivity(Molecular_weight)
K_D               = DEBmicroTrait.specific_reference_affinity(V_cell, ρ_p*ones(1,1), D_S)*1e3

# Phenylalanine
Vmax = 0.5*1e-9./dw/60
Min_gen_time = [log(2)./Vmax[1]]

elementstr = "C9H11NO2"
Molecular_weight = [165.192]

chemFormBiom       = [1, 1.64, 0.2, 0.38, 0, 0, 0]
N_C = DEBmicroTrait.extract_composition(elementstr)[1]
out = DEBmicroTrait.get_lambda(elementstr, chemFormBiom)
yield = out[4][10]/out[4][1]
y_DE = @. 1/abs(yield)

find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, V_cell, Min_gen_time, Gram_stain, rrn_copies, y_DE, N_C, y_EM)
ρ_p = Roots.find_zero(find_ρ, 1.0)
D_S               = DEBmicroTrait.aqueous_diffusivity(Molecular_weight)
K_D               = DEBmicroTrait.specific_reference_affinity(V_cell, ρ_p*ones(1,1), D_S)*1e3

# Tryptophan
Vmax = 0.1*1e-9./dw/60
Min_gen_time = [log(2)./Vmax[1]]

elementstr = "C11H12N2O2"
Molecular_weight = [204.229]

chemFormBiom       = [1, 1.64, 0.2, 0.38, 0, 0, 0]
N_C = DEBmicroTrait.extract_composition(elementstr)[1]
out = DEBmicroTrait.get_lambda(elementstr, chemFormBiom)
yield = out[4][10]/out[4][1]
y_DE = @. 1/abs(yield)

find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, V_cell, Min_gen_time, Gram_stain, rrn_copies, y_DE, N_C, y_EM)
ρ_p = Roots.find_zero(find_ρ, 1.0)
D_S               = DEBmicroTrait.aqueous_diffusivity(Molecular_weight)
K_D               = DEBmicroTrait.specific_reference_affinity(V_cell, ρ_p*ones(1,1), D_S)*1e3

###########################################################
# Streptococcus thermophylus
# isoleucine

Genome_size = [1.82e6]
V_cell = DEBmicroTrait.genome_size_to_cell_volume(Genome_size)
rrn_copies = [6.0]
Gram_stain = ["+"]
gmax = DEBmicroTrait.gmax_regression(rrn_copies)
dw = DEBmicroTrait.cell_volume_to_dry_mass(V_cell, gmax, Gram_stain)*12.011
Vmax = 0.027*1e-9./dw/60
Min_gen_time = [log(2)./Vmax[1]]
y_EM = [1.0]

elementstr = "C6H13NO2"
Molecular_weight = [131.17]
chemFormBiom       = [1, 1.64, 0.2, 0.38, 0, 0, 0]
N_C = DEBmicroTrait.extract_composition(elementstr)[1]
out = DEBmicroTrait.get_lambda(elementstr, chemFormBiom)
yield = out[4][10]/out[4][1]
y_DE = @. 1/abs(yield)

find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, V_cell, Min_gen_time, Gram_stain, rrn_copies, y_DE, N_C, y_EM)
ρ_p = Roots.find_zero(find_ρ, 1.0)
D_S               = DEBmicroTrait.aqueous_diffusivity(Molecular_weight)
K_D               = DEBmicroTrait.specific_reference_affinity(V_cell, ρ_p*ones(1,1), D_S)*1e3

# L-valine
elementstr = "C5H11NO2"
Molecular_weight = [117.151]
Vmax = 0.034*1e-9./dw/60
Min_gen_time = [log(2)./Vmax[1]]
chemFormBiom       = [1, 1.64, 0.2, 0.38, 0, 0, 0]
N_C = DEBmicroTrait.extract_composition(elementstr)[1]
out = DEBmicroTrait.get_lambda(elementstr, chemFormBiom)
yield = out[4][10]/out[4][1]
y_DE = @. 1/abs(yield)

find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, V_cell, Min_gen_time, Gram_stain, rrn_copies, y_DE, N_C, y_EM)
ρ_p = Roots.find_zero(find_ρ, 1.0)
D_S               = DEBmicroTrait.aqueous_diffusivity(Molecular_weight)
K_D               = DEBmicroTrait.specific_reference_affinity(V_cell, ρ_p*ones(1,1), D_S)*1e3

###########################################################
# Escherichia coli
Genome_size = [4.64e6]
V_cell = DEBmicroTrait.genome_size_to_cell_volume(Genome_size)

Gram_stain = ["-"]
rrn_copies = [7.0]
gmax = DEBmicroTrait.gmax_regression(rrn_copies)
nd = DEBmicroTrait.cell_volume_to_cell_number_density(V_cell, gmax, Gram_stain)
Vmax = 125*1e-12./(4.2e8./nd*1e-23)*60^2
Min_gen_time = [log(2)./Vmax[1]]
y_EM = [1.0]

# Glycerol-3-phosphate
elementstr = "C3H9O6P"
Molecular_weight = [172.073]

chemFormBiom       = [1, 1.64, 0.2, 0.38, 0, 0, 0]
N_C = DEBmicroTrait.extract_composition(elementstr)[1]
out = DEBmicroTrait.get_lambda(elementstr, chemFormBiom)
yield = out[4][10]/out[4][1]
y_DE = @. 1/abs(yield)

find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, V_cell, Min_gen_time, Gram_stain, rrn_copies, y_DE, N_C, y_EM)
ρ_p = Roots.find_zero(find_ρ, 1.0)
N_SB              = DEBmicroTrait.transporter_density_to_monomer_uptake_sites(V_cell, ρ_p*ones(1,1), Min_gen_time, Gram_stain)
Vmax              = @. 180.0*60^2*N_SB.*N_C
D_S               = DEBmicroTrait.aqueous_diffusivity(Molecular_weight)
K_D               = DEBmicroTrait.specific_reference_affinity(V_cell, ρ_p*ones(1,1), D_S)*1e3

###########################################################
# Rhizobium leguminosarum
Genome_size = [7.6e6]
V_cell = DEBmicroTrait.genome_size_to_cell_volume(Genome_size)
Gram_stain = ["-"]
rrn_copies = [3.0]
gmax = DEBmicroTrait.gmax_regression(rrn_copies)
dw = DEBmicroTrait.cell_volume_to_dry_mass(V_cell, gmax, Gram_stain)
Vmax = 0.16*1e-9./dw/12.01/60
Min_gen_time = [log(2)./Vmax[1]]
y_EM = [1.0]

# Succinate
elementstr = "C4H6O4"
Molecular_weight = [118.09]

chemFormBiom       = [1, 1.64, 0.2, 0.38, 0, 0, 0]
N_C = DEBmicroTrait.extract_composition(elementstr)[1]
out = DEBmicroTrait.get_lambda(elementstr, chemFormBiom)
yield = out[4][10]/out[4][1]
y_DE = @. 1/abs(yield)

find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, V_cell, Min_gen_time, Gram_stain, rrn_copies, y_DE, N_C, y_EM)
ρ_p = Roots.find_zero(find_ρ, 1.0)
N_SB              = DEBmicroTrait.transporter_density_to_monomer_uptake_sites(V_cell, ρ_p*ones(1,1), Min_gen_time, Gram_stain)
Vmax              = @. 180.0*60^2*N_SB.*N_C
D_S               = DEBmicroTrait.aqueous_diffusivity(Molecular_weight)
K_D               = DEBmicroTrait.specific_reference_affinity(V_cell, ρ_p*ones(1,1), D_S)*1e3

###########################################################

###########################################################
# Bradyrhizobium japonicum
Genome_size = [9.32e6]
V_cell = DEBmicroTrait.genome_size_to_cell_volume(Genome_size)
Gram_stain = ["-"]
rrn_copies = [2.0]
V_p = DEBmicroTrait.cell_volume_to_protein_volume(V_cell)
d_p = 1.37e6
mp = d_p*V_p*12.01
Vmax = 3.3*1e-9./(1e3*mp)/60
Min_gen_time = [log(2)./Vmax[1]]
y_EM = [1.0]

# Succinate
elementstr = "C4H6O4"
Molecular_weight = [118.09]

chemFormBiom       = [1, 1.64, 0.2, 0.38, 0, 0, 0]
N_C = DEBmicroTrait.extract_composition(elementstr)[1]
out = DEBmicroTrait.get_lambda(elementstr, chemFormBiom)
yield = out[4][10]/out[4][1]
y_DE = @. 1/abs(yield)

find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, V_cell, Min_gen_time, Gram_stain, rrn_copies, y_DE, N_C, y_EM)
ρ_p = Roots.find_zero(find_ρ, 1.0)
D_S               = DEBmicroTrait.aqueous_diffusivity(Molecular_weight)
K_D               = DEBmicroTrait.specific_reference_affinity(V_cell, ρ_p*ones(1,1), D_S)*1e3

###########################################################
# Pseudomonas sp.
Genome_size = [6.22e6]
V_cell = DEBmicroTrait.genome_size_to_cell_volume(Genome_size)
Gram_stain = ["-"]
rrn_copies = [6.0]
gmax = DEBmicroTrait.gmax_regression(rrn_copies)
dw = DEBmicroTrait.cell_volume_to_dry_mass(V_cell, gmax, Gram_stain)*12.011
Vmax = 159*1e-9./(1e3*dw)/60
Min_gen_time = [log(2)./Vmax[1]]
y_EM = [1.0]

# Toluene
elementstr = "C7H8"
Molecular_weight = [92.141]

chemFormBiom       = [1, 1.64, 0.2, 0.38, 0, 0, 0]
N_C = DEBmicroTrait.extract_composition(elementstr)[1]
out = DEBmicroTrait.get_lambda(elementstr, chemFormBiom)
yield = out[4][10]/out[4][1]
y_DE = @. 1/abs(yield)

find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, V_cell, Min_gen_time, Gram_stain, rrn_copies, y_DE, N_C, y_EM)
ρ_p = Roots.find_zero(find_ρ, 1.0)
D_S               = DEBmicroTrait.aqueous_diffusivity(Molecular_weight)
K_D               = DEBmicroTrait.specific_reference_affinity(V_cell, ρ_p*ones(1,1), D_S)*1e3

# Methanol
Vmax = 67*1e-9./(1e3*dw)/60
Min_gen_time = [log(2)./Vmax[1]]
y_EM = [1.0]
elementstr = "CH4O"
Molecular_weight = [32.04]

chemFormBiom       = [1, 1.64, 0.2, 0.38, 0, 0, 0]
N_C = DEBmicroTrait.extract_composition(elementstr)[1]
out = DEBmicroTrait.get_lambda(elementstr, chemFormBiom)
yield = out[4][10]/out[4][1]
y_DE = @. abs(yield)

find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, V_cell, Min_gen_time, Gram_stain, rrn_copies, y_DE, N_C, y_EM)
ρ_p = Roots.find_zero(find_ρ, 1.0)
N_SB              = DEBmicroTrait.transporter_density_to_monomer_uptake_sites(V_cell, ρ_p*ones(1,1), Min_gen_time, Gram_stain)
Vmax              = @. 180.0*60^2*N_SB.*N_C
D_S               = DEBmicroTrait.aqueous_diffusivity(Molecular_weight)
K_D               = DEBmicroTrait.specific_reference_affinity(V_cell, ρ_p*ones(1,1), D_S)*1e3


################################################
# I/O
df = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/benchmark_uptake_kinetics.csv", DataFrame, missingstring="N/A")

################################################
# Statistics
linearRegressor = lm(@formula(log(Kt) ~ log(Ktmodel)), df)
r2(linearRegressor)

# test for significant difference to 1-1 line
df.diff = df.Ktmodel-df.Kt
linearRegressor = lm(@formula(diff ~ Ktmodel), df)
# RMSD
RMSD = sqrt((1/(13-1)*sum((df.Ktmodel - df.Kt).^2)))
################################################
