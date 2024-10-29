using DEBmicroTrait
using CSV, DataFrames, Statistics
using Roots
using JLD

########################################
# I/O
dir             = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()
df_isolates     = CSV.read(joinpath(dir, "files/input/dement_isolates.csv"), DataFrame, missingstring="")
dropmissing!(df_isolates)
########################################
# metabolite traits
Formula  = "C6H12O6"
N_C      = [6]
Molecular_weight = [180.15588]
########################################
# isolate traits
Genome_size     = convert(Array{Float64,1}, df_isolates.genome_length)
V_cell          = DEBmicroTrait.genome_size_to_cell_volume(Genome_size)
rrn_copies      = convert(Array{Float64,1}, df_isolates.rrncopies)
Min_gen_time    = df_isolates.mingentime
Gram_stain      = repeat(["+"], 5369)
y_EM                    = ones(size(V_cell,1))

z_transporters = reshape(convert(Array{Float64,1}, df_isolates.total_transporters./Genome_size*1e6),1, size(V_cell,1))
genome_distr = z_transporters
########################################
# estimate transporter density
ρ_ps                    = zeros(1, size(V_cell,1))
y_DEs                   = zeros(1, size(V_cell,1))

for i in 1:size(V_cell,1)
    find_ρ(x)   = DEBmicroTrait.constrain_transporter_density_cost(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], Formula)
    ρ_p         = Roots.find_zero(find_ρ, 0.01)
    ρ_p         = isnan(ρ_p) ? 0.0 : ρ_p 
    closure     = genome_distr[:,i]./sum(genome_distr[:,i])
    ρ_ps[1,i]     = ρ_p[1].*closure[1]
    y_DE        = DEBmicroTrait.yield_transporter_density_cost(ρ_ps[i], [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], Formula)
    y_DEs[1,i]  = y_DE[1]
 end

ρ_ps = replace(ρ_ps, NaN => 0.0)
ρ_ps[ρ_ps.==0.0] .= 1e-12
median(ρ_ps)

########################################
# calculate uptake traits
N_SB              = DEBmicroTrait.transporter_density_to_monomer_uptake_sites(V_cell, ρ_ps, Min_gen_time, Gram_stain)
Vmax              = @. 180.0*60^2*N_SB.*N_C
#
D_S               = DEBmicroTrait.aqueous_diffusivity(Molecular_weight)
K_D               = DEBmicroTrait.specific_reference_affinity(V_cell, ρ_ps, D_S)
K_D = replace(K_D, NaN => 0.01)
#
a_s               = Vmax./K_D

########################################
# I/O
save("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/output/dement_isolates_assimilation.jld", "rho", ρ_ps, "NSB", N_SB, "KD", K_D, "yEM", y_EM, "yDE", y_DEs, "NC", N_C)
