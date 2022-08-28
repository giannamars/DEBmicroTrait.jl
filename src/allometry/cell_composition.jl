abstract type AbstractIsolate end

"""
    CellComposition structure

Structure that contains the cellular composition of isolates.
"""
@units @description mutable struct IsolateComposition <: AbstractIsolate
    Cell_volume::Vector{Float64}        |  m^3   | "Cell volume"
    Protein_volume::Vector{Float64}     |  m^3   | "Protein volume"
    Ribosome_volume::Vector{Float64}    |  m^3   | "Ribosome volume"
    Envelope_volume::Vector{Float64}    |  m^3   | "Cell envelope volume"
    CNP::Matrix{Float64}                | _      | "Cell stoichiometry"
    YEV::Matrix{Float64}                | mol/mol| "Yield of structure on CNP-reserve"
end

function IsolateComposition(Genome_size, Min_gen_time, Gram_stain)
    V_cell = genome_size_to_cell_volume(collect(skipmissing(Genome_size)))
    gmax   = log(2)./collect(skipmissing(Min_gen_time))
    V_protein = cell_volume_to_protein_volume(V_cell)
    V_r = cell_volume_to_ribosome_volume(V_cell, gmax)
    V_env = cell_volume_to_envelope_volume(V_cell, collect(skipmissing(Gram_stain)))
    CNP = zeros(3, size(V_cell,1))
    for i in 1:size(V_cell,1)
        CNP[:,i] = cell_volume_to_stoichiometry([V_cell[i]], [gmax[i]])
    end
    YEV = zeros(3, size(V_cell,1))
    for i in 1:size(V_cell,1)
        YEV[:,i] = cell_volume_to_yields([V_cell[i]], [gmax[i]])
    end
    return IsolateComposition(V_cell, V_protein, V_r, V_env, CNP, YEV)
end


function genome_size_to_cell_volume(L_DNA::Vector{Float64})
    # Kempes et al. (2016), Eq. 8
    v_N = 1.47e-27
    V_DNA = v_N*L_DNA
    D_0 = 3e-17
    β_D = 0.21
    V_c = (V_DNA/D_0).^(1/β_D) # m^3
end

function surface_area_volume_ratio(V_c::Vector{Float64})
    A_c = @. π^(1/3)*(6*V_c)^(2/3)
    SAV = @. A_c/V_c
end

function cell_volume_to_genome_size(V_c::Vector{Float64})
    # Kempes et al. (2016), Eq. 8
    D_0 = 3e-17
    β_D = 0.21
    V_DNA = D_0*V_c.^(β_D) # m^3
    v_N = 1.47e-27
    L_DNA = V_DNA./v_N
end

function RC_to_cell_volume(R_C::Vector{Float64})
    # Kempes et al. (2016), Eq. 8
    D_0 = 9.58
    β_D = 0.66
    V_c = 1e-18(R_C/D_0).^(1/β_D) # m^3
end

function cell_volume_to_protein_volume(V_c::Vector{Float64})
    # Kempes et al. (2016), Eq. 12
    P_0 = 3.42e-7
    β_p = 0.70
    V_p = P_0*V_c.^β_p # m^3
end

function cell_volume_to_ribosome_volume(V_c::Vector{Float64}, gmax::Vector{Float64})
    # Kempes et al. (2016), Eq. 13
    l_p = 975
    V_p = cell_volume_to_protein_volume(V_c)
    m_p = 5.81e-20
    d_p = 1.37e6
    v_p = m_p/d_p
    N_p = V_p/v_p
    ϕ = 6.20e-5
    η = 6.20e-5
    l_r = 4566
    r_r = 63
    μ = gmax./60^2
    # μ_0 = 4e7
    # β_B = 1.64
    # μ = μ_0*V_c.^(β_B-1)
    N_r = @. l_p*N_p*(ϕ/μ + 1)/(r_r/μ - l_r*(η/μ + 1))
    v_r = 3.04e-24
    V_r = N_r*v_r # m^3
end

function cell_volume_to_mRNA_volume(V_c::Vector{Float64}, gmax::Vector{Float64})
    # Kempes et al. (2016), Eq. 18
    V_r = cell_volume_to_ribosome_volume(V_c, gmax./60^2)
    v_r = 3.04e-24
    N_r = V_r/v_r
    v_mRNA = 1.43e-24
    n_mRNA = 1.08
    V_mRNA = v_mRNA*n_mRNA*N_r # m^3
end

function cell_volume_to_tRNA_volume(V_c::Vector{Float64}, gmax::Vector{Float64})
    # Kempes et al. (2016), Eq. 17
    V_r = cell_volume_to_ribosome_volume(V_c, gmax./60^2)
    v_r = 3.04e-24
    N_r = V_r/v_r
    v_tRNA = 3.10e-26
    n_tRNA = 9.3
    V_tRNA = v_tRNA*n_tRNA*N_r # m^3
end

function cell_volume_to_envelope_volume(V_c::Vector{Float64}, Gram_stain::Vector{String})
    # Kempes et al. (2016), Eq. 19
    r_c = (3*V_c/(4*π)).^(1/3)
    p_p = 0.149
    V_env = zeros(size(r_c,1))
    for i in 1:size(r_c, 1)
        if Gram_stain[i] == "(-)"
            r_env = 40.5e-9      # 27.5 if gram +
            V_env[i] = (V_c[i] - 4/3*π*(r_c[i] - r_env).^3)*(1-p_p) # m^3
        else
            r_env = 27.5e-9
            V_env[i] = (V_c[i] - 4/3*π*(r_c[i] - r_env).^3)*(1-p_p) # m^3
        end
    end
    return V_env
end

function cell_volume_to_dry_mass(V_c::Vector{Float64}, gmax::Vector{Float64}, Gram_stain::Vector{String})
    v_0 = 3.0e-17
    β_D = 0.21
    V_DNA = v_0*V_c.^β_D
    V_p = cell_volume_to_protein_volume(V_c)
    V_r = cell_volume_to_ribosome_volume(V_c, gmax)
    V_mRNA = cell_volume_to_mRNA_volume(V_c, gmax)
    V_tRNA = cell_volume_to_tRNA_volume(V_c, gmax)
    V_env = cell_volume_to_envelope_volume(V_c, Gram_stain)
    d_DNA = 2e6
    d_p = 1.37e6
    d_r = 1.79e6
    d_mem = 1.05e6
    d_RNA = 2e6
    m_d = d_DNA*V_DNA + d_p*V_p + d_r*V_r + d_mem*V_env + d_RNA*V_mRNA + d_RNA*V_tRNA  # g
end

function cell_volume_to_cell_number_density(V_c::Vector{Float64}, gmax::Vector{Float64}, Gram_stain::Vector{String})
    # assuming 47% of dry mass is C (BNID 100649)
    m_d = cell_volume_to_dry_mass(V_c, gmax, Gram_stain)
    M_C = 12.0111
    N_A = 6.022e23
    λ_B = M_C./(0.47*m_d*N_A)    # mol cells/mol C
end

function cell_volume_to_cellular_density(V_c::Vector{Float64}, gmax::Vector{Float64}, Gram_stain::Vector{String})
    # Kempes et al. (2016), Eq. S53
    v_0 = 3.0e-17
    β_D = 0.21
    V_DNA = v_0*V_c.^β_D
    V_p = cell_volume_to_protein_volume(V_c)
    V_r = cell_volume_to_ribosome_volume(V_c, gmax)
    V_mRNA = cell_volume_to_mRNA_volume(V_c, gmax)
    V_tRNA = cell_volume_to_tRNA_volume(V_c, gmax)
    V_env = cell_volume_to_envelope_volume(V_c, Gram_stain)
    V_w = V_c - (V_DNA + V_p + V_r + V_env + V_tRNA + V_mRNA)
    d_DNA = 2e6
    d_p = 1.37e6
    d_r = 1.79e6
    d_mem = 1.05e6
    d_RNA = 2e6
    d_w = 1e6
    d_c = @. (d_DNA*V_DNA + d_p*V_p + d_r*V_r + d_mem*V_env + d_RNA*V_mRNA + d_RNA*V_tRNA + d_w*V_w)/V_c    #g/m^3
end


function genome_size_to_rRNA_copy_number(L_DNA::Vector{Float64})
    # Roller et al. (2016) - Supp. Fig. 1
    rRNA  = @. 2^(L_DNA./1e6 - 2^2.8)/0.66 + 1.0
    return round.(rRNA)
end
