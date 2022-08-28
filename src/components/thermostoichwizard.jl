function extract_composition(elementstring::String)
    # (C,H,N,O,S,P)
    regexes = [r"[C]\d*", r"[H]\d*", r"[N]\d*", r"[O]\d*", r"[S]\d*", r"[P]\d*"]
    chemical_indices = Vector{Int64}(undef, 6)
    for i in eachindex(regexes)
        m = collect(m.match for m in eachmatch(regexes[i], elementstring))
        if size(m,1) == 0
            chemical_indices[i] = 0
        else
            try
                chemical_indices[i] = parse(Int64, split(m[1], r"[A-Z]")[2])
            catch e
                chemical_indices[i] = 1
            end
        end
    end
    return chemical_indices
end

function get_stoich_electron_donor(elementstring::String)
    chemical_indices = extract_composition(elementstring)
    a = chemical_indices[1] #C
    b = chemical_indices[2] #H
    c = chemical_indices[3] #N
    d = chemical_indices[4] #O
    e = chemical_indices[5] #S
    f = chemical_indices[6] #P
    z = 0
    return [-1,  -(3*a+4*e-d), a, c, e, f, 5*a+b-4*c-2*d+7*e-f, -z+4*a+b-3*c-2*d+5*e-2*f, 0, 0]
end


function get_stoich_electron_acceptor()
    stoich_electron_acceptor = zeros(10)
    stoich_electron_acceptor[9] = -1  # oxygen
    stoich_electron_acceptor[7] = -4  # h+
    stoich_electron_acceptor[8] = -4  # e-
    stoich_electron_acceptor[2] =  2  # h2o
    return stoich_electron_acceptor
end


function get_stoich_catabolic_reaction(elementstring::String)
    stoich_electron_donor = get_stoich_electron_donor(elementstring)
    stoich_electron_acceptor = get_stoich_electron_acceptor()
    yEd = stoich_electron_donor[8]
    yEa = stoich_electron_acceptor[8]
    stoich_cat_rxns = stoich_electron_donor-(yEd/yEa)*stoich_electron_acceptor
    return stoich_cat_rxns
end


function get_stoich_anabolic_reaction(elementstring::String, chemFormBiom)
    chemical_indices = extract_composition(elementstring)
    stoich_electron_donor = get_stoich_electron_donor(elementstring)
    stoich_electron_acceptor = get_stoich_electron_acceptor()
    a = chemical_indices[1] #C

    #chemFormBiom = [1, 1.8, 0.2, 0.5, 0, 0, 0]  # C H_1.8 N_0.2 O_0.5 made variable
    aB = chemFormBiom[1]
    bB = chemFormBiom[2]
    cB = chemFormBiom[3]
    dB = chemFormBiom[4]
    eB = chemFormBiom[5]
    fB = chemFormBiom[6]
    zB = chemFormBiom[7]

    ySource = -1
    yH2o = -(3*aB+4*eB-dB)
    yHco3 = aB
    yNh4 = cB
    yHpo4 = eB
    yHs = fB
    yH = 5*aB+bB-4*cB-2*dB+7*eB-fB
    yE = -zB+4*aB+bB-3*cB-2*dB+5*eB-2*fB

    stoichAnStarB = [ySource,yH2o,yHco3,yNh4,yHpo4,yHs,yH,yE,0,0]
    stoichAnStarB = -stoichAnStarB
    stoichAnStarB[end] = stoichAnStarB[1]
    stoichAnStarB[1] = 0

    # O2
    stoichAnStar_O2 = stoichAnStarB+(1/a)*stoich_electron_donor
    yEana = stoichAnStar_O2[8]
    if yEana > 0.0
        yEa = stoich_electron_acceptor[8]
        stoichAn_O2 = stoichAnStar_O2-yEana/yEa*stoich_electron_acceptor
    elseif yEana < 0.0
        yEd = stoich_electron_donor[8]
        stoichAn_O2 = stoichAnStar_O2-yEana/yEd*stoich_electron_donor
    else
        stoichAn_O2 = stoichAnStar_O2
    end

    # HCO3-
    yEd = stoich_electron_donor[8]
    yEa = stoichAnStarB[8]
    stoichAn_HCO3 = stoich_electron_donor-(yEd/yEa)*stoichAnStarB
    stoichAn_HCO3 = stoichAn_HCO3/stoichAn_HCO3[10]

    return stoichAn_O2, stoichAn_HCO3
end

function get_lambda(elementstring::String, chemFormBiom)
    chemical_indices = extract_composition(elementstring)
    stoich_electron_donor = get_stoich_electron_donor(elementstring)
    stoich_cat_rxns = get_stoich_catabolic_reaction(elementstring)
    stoich_anabolic_O2, stoich_anabolic_HCO3 = get_stoich_anabolic_reaction(elementstring, chemFormBiom)

    a = chemical_indices[1] #C
    b = chemical_indices[2] #H
    c = chemical_indices[3] #N
    d = chemical_indices[4] #O
    e = chemical_indices[5] #S
    f = chemical_indices[6] #P
    z = 0

    ne = -z+4*a+b-3*c-2*d+5*e-2*f  # number of electrons transferred in D
    nosc = -(ne/a) + 4
    delGcox0PerC = 60.3-28.5*nosc  # kJ/C-mol
    delGcox0 = delGcox0PerC*a*abs(stoich_electron_donor[1])
    # - estimate delGf0 for electron donor
    delGf0_D_zero = 0
    delGf0_zero = [delGf0_D_zero, -237.2, -586.8, -79.3, -1096.1, 12.1, 0, 0, 16.4, -67]
    delGcox0_zero = dot(delGf0_zero, stoich_electron_donor)
    delGf0_D_est = (delGcox0-delGcox0_zero)/stoich_electron_donor[1]
    # - finally, delGf0
    delGf0 = delGf0_zero
    delGf0[1] = delGf0_D_est

    # - standard delG at pH=0
    delGcat0 = dot(delGf0, stoich_cat_rxns)
    delGan0_O2 = dot(delGf0, stoich_anabolic_O2)
    delGan0_HCO3 = dot(delGf0, stoich_anabolic_HCO3)

    # - stadard delG at pH=7
    R = 0.008314  # kJ/(K.mol)
    T = 298  # K
    iProton = 7  # [eD,h2o,hco3-,nh4+,hpo4**2-,hs-,h+,e-,eA,biom]
    delGcox = delGcox0+R*T*stoich_electron_donor[iProton]*log(1e-7)
    delGcat = delGcat0+R*T*stoich_cat_rxns[iProton]*log(1e-7)
    delGan_O2 = delGan0_O2+R*T*stoich_anabolic_O2[iProton]*log(1e-7)
    delGan_HCO3 = delGan0_HCO3+R*T*stoich_anabolic_HCO3[iProton]*log(1e-7)


    # The Thermodynamic Electron Equivalents Model (TEEM)
    # --------
    eta = 0.43
    delGsyn = 200 # kJ/(mol.X)

    if delGan_O2 < 0
        m_O2 = 1
    else
        m_O2 = -1
    end

    if delGan_HCO3 < 0
        m_HCO3 = 1
    else
        m_HCO3 = -1
    end

    lambda_O2 = (delGan_O2*eta^m_O2+delGsyn)/(-delGcat*eta)
    lambda_HCO3 = (delGan_HCO3*eta^m_HCO3+delGsyn)/(-delGcat*eta)

    if lambda_O2 > 0
        stoichMet_O2 = lambda_O2*stoich_cat_rxns+stoich_anabolic_O2
    else
        stoichMet_O2 = stoich_anabolic_O2
    end

    if lambda_HCO3 > 0
        stoichMet_HCO3 = lambda_HCO3*stoich_cat_rxns+stoich_anabolic_HCO3
    else
        stoichMet_HCO3 = stoich_anabolic_HCO3
    end

    delGdis_O2 = dot(delGf0,stoichMet_O2) + R*T*stoichMet_O2[iProton]*log(1e-7)
    delGdis_HCO3 = dot(delGf0,stoichMet_HCO3) + R*T*stoichMet_HCO3[iProton]*log(1e-7)

    return [lambda_O2, lambda_HCO3], [delGcox0PerC,delGcox0,delGcox,delGcat0,delGcat,delGan0_O2,delGan0_HCO3, delGan_O2,delGan_HCO3,delGdis_O2,delGdis_HCO3], [m_O2, m_HCO3], stoichMet_O2, stoichMet_HCO3
end


function get_lambda(elementstring::String, chemFormBiom, V_c, Min_gen_time, Gram_stain, ρ_p)
    chemical_indices = extract_composition(elementstring)
    stoich_electron_donor = get_stoich_electron_donor(elementstring)
    stoich_cat_rxns = get_stoich_catabolic_reaction(elementstring)
    stoich_anabolic_O2, stoich_anabolic_HCO3 = get_stoich_anabolic_reaction(elementstring, chemFormBiom)

    a = chemical_indices[1] #C
    b = chemical_indices[2] #H
    c = chemical_indices[3] #N
    d = chemical_indices[4] #O
    e = chemical_indices[5] #S
    f = chemical_indices[6] #P
    z = 0

    ne = -z+4*a+b-3*c-2*d+5*e-2*f  # number of electrons transferred in D
    nosc = -(ne/a) + 4
    delGcox0PerC = 60.3-28.5*nosc  # kJ/C-mol
    delGcox0 = delGcox0PerC*a*abs(stoich_electron_donor[1])
    # - estimate delGf0 for electron donor
    delGf0_D_zero = 0
    delGf0_zero = [delGf0_D_zero, -237.2, -586.8, -79.3, -1096.1, 12.1, 0, 0, 16.4, -67]
    delGcox0_zero = dot(delGf0_zero, stoich_electron_donor)
    delGf0_D_est = (delGcox0-delGcox0_zero)/stoich_electron_donor[1]
    # - finally, delGf0
    delGf0 = delGf0_zero
    delGf0[1] = delGf0_D_est

    # - standard delG at pH=0
    delGcat0 = dot(delGf0, stoich_cat_rxns)
    delGan0_O2 = dot(delGf0, stoich_anabolic_O2)
    delGan0_HCO3 = dot(delGf0, stoich_anabolic_HCO3)

    # - stadard delG at pH=7
    R = 0.008314  # kJ/(K.mol)
    T = 298  # K
    iProton = 7  # [eD,h2o,hco3-,nh4+,hpo4**2-,hs-,h+,e-,eA,biom]
    delGcox = delGcox0+R*T*stoich_electron_donor[iProton]*log(1e-7)
    delGcat = delGcat0+R*T*stoich_cat_rxns[iProton]*log(1e-7)
    delGan_O2 = delGan0_O2+R*T*stoich_anabolic_O2[iProton]*log(1e-7)
    delGan_HCO3 = delGan0_HCO3+R*T*stoich_anabolic_HCO3[iProton]*log(1e-7)


    # The Thermodynamic Electron Equivalents Model (TEEM)
    # --------
    eta = 0.43
    #
    r_p = 1e-9
    N_A = 6.022e23
    r_c = (3*V_c/(4*pi)).^(1/3)   # cell radius
    gmax = log(2)./Min_gen_time
    N_p = @. 4*ρ_p*r_c^2/r_p^2
    N_SB = transporter_density_to_monomer_uptake_sites(V_c, ρ_p*ones(1,1), Min_gen_time, Gram_stain)[1]
    delG_porter = @. 1.09e-19*60^2*gmax*N_SB*N_A/1e3*6 # kJ/mol.X

    delGsyn = 4.18 + delG_porter[1]  # kJ/(mol.X)
    #delGsyn = 200.0 + delG_porter[1]  # kJ/(mol.X)

    if delGan_O2 < 0
        m_O2 = 1
    else
        m_O2 = -1
    end

    if delGan_HCO3 < 0
        m_HCO3 = 1
    else
        m_HCO3 = -1
    end

    lambda_O2 = (delGan_O2*eta^m_O2+delGsyn)/(-delGcat*eta)
    lambda_HCO3 = (delGan_HCO3*eta^m_HCO3+delGsyn)/(-delGcat*eta)

    if lambda_O2 > 0
        stoichMet_O2 = lambda_O2*stoich_cat_rxns+stoich_anabolic_O2
    else
        stoichMet_O2 = stoich_anabolic_O2
    end

    if lambda_HCO3 > 0
        stoichMet_HCO3 = lambda_HCO3*stoich_cat_rxns+stoich_anabolic_HCO3
    else
        stoichMet_HCO3 = stoich_anabolic_HCO3
    end

    delGdis_O2 = dot(delGf0,stoichMet_O2) + R*T*stoichMet_O2[iProton]*log(1e-7)
    delGdis_HCO3 = dot(delGf0,stoichMet_HCO3) + R*T*stoichMet_HCO3[iProton]*log(1e-7)

    return [lambda_O2, lambda_HCO3], [delGcox0PerC,delGcox0,delGcox,delGcat0,delGcat,delGan0_O2,delGan0_HCO3, delGan_O2,delGan_HCO3,delGdis_O2,delGdis_HCO3], [m_O2, m_HCO3], stoichMet_O2, stoichMet_HCO3
end
