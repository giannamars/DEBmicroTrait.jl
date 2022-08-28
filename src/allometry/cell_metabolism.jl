function cell_volume_to_specific_maintenance_rate(V_c::Vector{Float64}, Min_gen_time::Vector{Float64}, Gram_stain::Vector{String})
    # Lynch & Marinov (2015), Eq. 1a
    gmax = log(2)./Min_gen_time
    n_0 = cell_volume_to_cell_number_density(V_c, gmax, Gram_stain)
    V_cell = V_c./1e-18         #  V_cell in μm^3
    E_M = 0.39*V_cell.^0.88     # [1e9 ATP/(cell h)]
    E_G = 1e9*E_M/26            # 26 ATP per glucose molecule: [glucose/(cell h)]
    E_C = E_G*6*24              # [C atoms/(cell/d)]
    mol_C = E_C/(6.022e23)      # [mol C/(cell d)]
    k_M = mol_C./n_0./24        # 1/h
end

function translation_power(Protein_volume::Vector{Float64}, Ribosome_volume::Vector{Float64}, Min_gen_time::Vector{Float64})
    V_p = Protein_volume
    V_r = Ribosome_volume
    gmax = log(2)./Min_gen_time
    transl_power = @. gmax*V_p/V_r
    return max.(transl_power, 1.01*gmax)
end

function relative_translation_efficiency_regression(rrn_copies::Vector{Float64})
    log2_rna = log2.(rrn_copies)
    y = @. 9.5 - 1.22*log2_rna
    return 9.6./y
end

function gmax_regression(rrn_copies::Vector{Float64})
    log2_rna = log2.(rrn_copies)
    y = @. log2(0.06) + 1.03*log2_rna
    return 2.0.^y
end

function enzyme_allocation(α_base, hydrolase_distr)
  genome_distr_norm = hydrolase_distr./sum(hydrolase_distr)
  closure = α_base.*genome_distr_norm
end


function constrain_transporter_density(ρ_p, V_c, Min_gen_time, Gram_stain, rrn_copies, y_DE, N_C, y_EM)
    N_SB = transporter_density_to_monomer_uptake_sites(V_c, ρ_p*ones(1,1), Min_gen_time, Gram_stain)
    k2p = 180.0*60^2
    j_EA_m = (1.0 .- y_DE).*N_C*k2p.*N_SB
    k_M = cell_volume_to_specific_maintenance_rate(V_c, Min_gen_time, Gram_stain)

    gmax   = log(2)./Min_gen_time
    V_p = cell_volume_to_protein_volume(V_c)
    V_r = cell_volume_to_ribosome_volume(V_c, gmax)

    k_E = translation_power(V_p, V_r, Min_gen_time)
    y_EV = relative_translation_efficiency_regression(rrn_copies)

    r_m = @. (j_EA_m-k_M*y_EM)/(j_EA_m/k_E+y_EV)
    return r_m[1] - gmax[1]
end


function constrain_transporter_density_cost(ρ_p, V_c, Min_gen_time, Gram_stain, rrn_copies, y_EM, elementstr)
    chemFormBiom       = [1, 1.64, 0.2, 0.38, 0, 0, 0]
    N_C = DEBmicroTrait.extract_composition(elementstr)[1]
    out = DEBmicroTrait.get_lambda(elementstr, chemFormBiom, V_c, Min_gen_time, Gram_stain, ρ_p)
    yield = out[4][10]/out[4][1]
    y_DE = @. 1/abs(yield)
    N_SB = transporter_density_to_monomer_uptake_sites(V_c, ρ_p*ones(1,1), Min_gen_time, Gram_stain)
    k2p = 180.0*60^2
    j_EA_m = (1.0 .- y_DE)*N_C*k2p*N_SB
    k_M = cell_volume_to_specific_maintenance_rate(V_c, Min_gen_time, Gram_stain)

    gmax   = log(2)./Min_gen_time
    V_p = cell_volume_to_protein_volume(V_c)
    V_r = cell_volume_to_ribosome_volume(V_c, gmax)

    k_E = translation_power(V_p, V_r, Min_gen_time)
    y_EV = relative_translation_efficiency_regression(rrn_copies)

    r_m = @. (j_EA_m-k_M*y_EM)/(j_EA_m/k_E+y_EV)
    return r_m[1] - gmax[1]
end

function yield_transporter_density_cost(ρ_p, V_c, Min_gen_time, Gram_stain, rrn_copies, y_EM, elementstr)
    chemFormBiom       = [1, 1.64, 0.2, 0.38, 0, 0, 0]
    N_C = DEBmicroTrait.extract_composition(elementstr)[1]
    out = DEBmicroTrait.get_lambda(elementstr, chemFormBiom, V_c, Min_gen_time, Gram_stain, ρ_p)
    yield = out[4][10]/out[4][1]
    y_DE = @. 1/abs(yield)
end


function check_growth_rate(ρ_p, V_c, Min_gen_time, Gram_stain, rrn_copies, y_DE, N_C, y_EM)
    N_SB = transporter_density_to_monomer_uptake_sites(V_c, ρ_p*ones(1,1), Min_gen_time, Gram_stain)
    k2p = 180.0*60^2
    j_EA_m = (1.0 .- y_DE)*N_C*k2p*N_SB
    #j_EA_m = (1.0/y_DE .- 1.0)*N_C*k2p*N_SB
    k_M = cell_volume_to_specific_maintenance_rate(V_c, Min_gen_time, Gram_stain)

    gmax   = log(2)./Min_gen_time
    V_p = cell_volume_to_protein_volume(V_c)
    V_r = cell_volume_to_ribosome_volume(V_c, gmax)

    k_E = translation_power(V_p, V_r, Min_gen_time)
    y_EV = relative_translation_efficiency_regression(rrn_copies)

    r_m = @. (j_EA_m-k_M*y_EM)/(j_EA_m/k_E+y_EV)
    return r_m[1]
end


function steady_state_reserve_density(ρ_p, V_c, Min_gen_time, Gram_stain, y_DE, N_C)
    # ρ_p[1,:], y_DE[1], N_C[1], i.e. loop over monomers
    ρ_p             = reshape(ρ_p, 1, size(V_cell,1))
    N_SB            = transporter_density_to_monomer_uptake_sites(V_cell, ρ_p, Min_gen_time, Gram_stain)
    k2p             = 180.0*60^2
    j_EA_m          = @. (1.0 .- y_DE)*N_C*k2p*N_SB
    gmax            = log(2)./Min_gen_time
    V_p             = cell_volume_to_protein_volume(V_cell)
    V_r             = cell_volume_to_ribosome_volume(V_cell, gmax)

    k_E             = translation_power(V_p, V_r, Min_gen_time)
    m_E_max         = j_EA_m[1,:]./k_E
end

function rate_yield_trade_off(ρ_p, α, V_c, Min_gen_time, Gram_stain, rrn_copies, y_DE, N_C, y_EM)
    N_SB            = transporter_density_to_monomer_uptake_sites(V_c, ρ_p, Min_gen_time, Gram_stain)
    k2p             = 180.0*60^2
    j_EA_m          = (1.0 .- y_DE)*N_C*k2p*N_SB
    k_M             = cell_volume_to_specific_maintenance_rate(V_c, Min_gen_time, Gram_stain)
    gmax            = log(2)./Min_gen_time
    V_p             = cell_volume_to_protein_volume(V_c)
    V_r             = cell_volume_to_ribosome_volume(V_c, gmax)
    k_E             = translation_power(V_p, V_r, Min_gen_time)
    y_EV            = relative_translation_efficiency_regression(rrn_copies)
    r_m             = @. (j_EA_m-k_M*y_EM)/(j_EA_m/k_E+(1+α)*y_EV)

    rtmp = convert(Array{Float64,1}, LinRange(1e-5, 1e1, 100000))
    filter!(x->x<=gmax[1],rtmp)
    r = rtmp
    m_E = @. (k_M*y_EM + (1+α)*r*y_EV)/(k_E-r)
    j_EA = @. m_E*k_E
    yield_coeff = @. y_EV*r./j_EA
    yield = yield_coeff[findall(x->x <= 1.0, m_E)]
    rate = r[findall(x->x <= 1.0, m_E)]
    return yield, convert(Array{Float64,1}, rate)
end
