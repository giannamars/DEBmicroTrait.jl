function transporter_density_to_monomer_uptake_sites(V_c::Vector{Float64}, ρ_closure::Matrix{Float64}, Min_gen_time::Vector{Float64}, Gram_stain::Vector{String})
    gmax = log(2)./Min_gen_time
    n_0 = cell_volume_to_cell_number_density(V_c, gmax, Gram_stain)
    r_p = 1e-9                    # 1 μm
    r_c = (3*V_c/(4*pi)).^(1/3)
    N_SB = zeros(size(ρ_closure,1), size(ρ_closure,2))
    for i in 1:size(V_c,1)
        N_SB[:,i] = n_0[i]*4*r_c[i]^2/(r_p^2)*ρ_closure[:,i]/12.011
    end
    N_SB[N_SB.<1e-12].=0.0
    return N_SB
end

function transporter_density_to_monomer_uptake_sites(V_c::Vector{Float64}, N_p::Matrix{Float64}, Min_gen_time::Vector{Float64}, Gram_stain::Vector{String}, t::String)
    gmax = log(2)./Min_gen_time
    n_0 = cell_volume_to_cell_number_density(V_c, gmax, Gram_stain)
    N_SB = zeros(size(N_p,1), size(N_p,2))
    for i in 1:size(V_c,1)
        N_SB[:,i] = n_0[i]*N_p[:,i]/12.011
    end
    N_SB[N_SB.<1e-12].=0.0
    return N_SB
end

function transporters_closure(ρ_p::Vector{Float64}, genome_distr::Matrix{Float64})
    closure   = zeros(size(genome_distr))
    for i in 1:size(closure,2)
        closure[:,i] = ρ_p[i]*genome_distr[:,i]./sum(genome_distr[:,i])
    end
    closure[closure.==0.0].=1e-8
    return closure
end

function transporter_density(V_c::Array{Float64,1}, N_p::Array{Float64,2})
    r_p = 1e-9                    # 1 μm
    r_c = (3*V_c/(4*pi)).^(1/3)
    ρ_p = zeros(size(N_p,1), size(N_p,2))
    for i in 1:size(V_c,1)
        ρ_p[:,i] = N_p[:,i]*r_p^2/(4*r_c[i]^2)
    end
    return ρ_p
end


function interception_probability(V_c::Array{Float64,1}, N_p::Array{Float64,2})
    r_p = 1e-9                    
    r_c = (3*V_c/(4*pi)).^(1/3)
    p_int = zeros(size(N_p,1), size(N_p,2))
    for i in 1:size(V_c,1)
        p_int[:,i] = @. (N_p[:,i]*r_p[i])/(N_p[:,i]*r_p + π*r_c[i])
    end
    return p_int
end

function interception_probability(l_c::Array{Float64,1}, b_c::Array{Float64,1}, N_p::Array{Float64,2})
    r_p = 1e-9
    A_c = @. 2π*b_c*(b_c+l_c)
    
    β_c = zeros(size(N_p,1), size(N_p,2))
    for i in 1:size(l_c,1)
        β_c[:,i] = @. 1 - N_p[:,i]*π*r_p^2/A_c[i]
    end

    P_s = @. 1 - 2b_c*r_p/(l_c*(b_c + l_c))

    p_int = zeros(size(N_p,1), size(N_p,2))
    for i in 1:size(l_c,1)
        p_int[:,i] =  @. (1-β_c[:,i])*P_s[i]/(1-β_c[:,i]*P_s[i])
    end
    return p_int
end

function interception_probability(l_c::Array{Float64,1}, b_c::Array{Float64,1}, N_p::Array{Float64,2}, t::String)
    r_p = 1e-9
    A_c = @. π^2*l_c*b_c
    
    β_c = zeros(size(N_p,1), size(N_p,2))
    for i in 1:size(l_c,1)
        β_c[:,i] = @. 1 - N_p[:,i]*π*r_p^2/A_c[i]
    end

    P_s = @. 1 - 4r_p/(π*b_c)*1/log(2l_c/b_c)

    p_int = zeros(size(N_p,1), size(N_p,2))
    for i in 1:size(l_c,1)
        p_int[:,i] =  @. (1-β_c[:,i])*P_s[i]/(1-β_c[:,i]*P_s[i])
    end
    return p_int
end


function specific_reference_affinity(V_c::Vector{Float64}, ρ_p::Matrix{Float64}, D_S::Vector{Float64})
    # N_p (number of transporters) -> ρ_p = N_p*4*pi*r_p^2/(4*pi*r_c^2) (transporter density on the cell surface)
    V_S = 180.0                   # k2p, [s]
    N_A = 6.022e23                # Avogrado constant [1/mol]
    r_p = 1e-9                    # porter radius, 1 μm
    r_c = (3*V_c/(4*pi)).^(1/3)   # cell radius
    #
    p_inv = zeros(size(ρ_p,1), size(ρ_p,2))
    for i in 1:size(ρ_p,2)
        p_inv[:,i] = @. (4*ρ_p[:,i]*r_c[i]^2/r_p + pi*r_c[i])/(4*ρ_p[:,i]*r_c[i]^2/r_p)
    end
    #
    K_SC_0 = zeros(size(ρ_p,1), size(ρ_p,2))
    for i in 1:size(ρ_p,2)
        K_SC_0[:,i] = @. ((V_S*ρ_p[:,i]*r_c[i])/(pi*D_S*r_p^2*N_A))*p_inv[:,i]
    end
    #
    any(x->x==true, isnan.(K_SC_0)) ? throw(DomainError("NaN in DEBmicroTrait.ECA_kinetics!"))  : return K_SC_0.*10
end

function specific_reference_affinity(l_c::Array{Float64,1}, b_c::Array{Float64,1}, N_p::Matrix{Float64}, D_S::Vector{Float64})
    # N_p (number of transporters) -> ρ_p = N_p*4*pi*r_p^2/(4*pi*r_c^2) (transporter density on the cell surface)
    V_S = 180.0                   # k2p, [s]
    N_A = 6.022e23                # Avogrado constant [1/mol]
    #
    p_inv = (interception_probability(l_c, b_c, N_p))^(-1)
    #
    e_c   = @. (1-b_c^2/l_c^2)^(1/2)
    l_eff = @. (e_c/tanh(e_c)^(-1))*l_c
    I_max = @. 4π*D_S*l_eff
    #
    K_SC_0 = zeros(size(N_p,1), size(N_p,2))
    for i in 1:size(N_p,2)
        K_SC_0[:,i] = @. ((V_S*N_p[:,i])/(I_max[i]*N_A))*p_inv[:,i]
    end
    #
    any(x->x==true, isnan.(K_SC_0)) ? throw(DomainError("NaN in DEBmicroTrait.ECA_kinetics!"))  : return K_SC_0.*10
end

function specific_reference_affinity(l_c::Array{Float64,1}, b_c::Array{Float64,1}, N_p::Matrix{Float64}, D_S::Vector{Float64}, t::String)
    # N_p (number of transporters) -> ρ_p = N_p*4*pi*r_p^2/(4*pi*r_c^2) (transporter density on the cell surface)
    V_S = 180.0                   # k2p, [s]
    N_A = 6.022e23                # Avogrado constant [1/mol]
    #
    p_inv = (interception_probability(l_c, b_c, N_p, t))^(-1)
    #
    e_c   = @. (1-b_c^2/l_c^2)^(1/2)
    l_eff = @. (e_c/tanh(e_c)^(-1))*l_c
    I_max = @. 4π*D_S*l_eff
    #
    K_SC_0 = zeros(size(N_p,1), size(N_p,2))
    for i in 1:size(N_p,2)
        K_SC_0[:,i] = @. ((V_S*N_p[:,i])/(I_max[i]*N_A))*p_inv[:,i]
    end
    #
    any(x->x==true, isnan.(K_SC_0)) ? throw(DomainError("NaN in DEBmicroTrait.ECA_kinetics!"))  : return K_SC_0.*10
end

function reynolds_number(u_c::Array{Float64,1}, l_c::Array{Float64,1}, ν_kin::Float64)
    # ν_kin = kinematic viscosity, u_c = swimming velocity, l_c = characteristic length
    Re = @. u_c*l_c/ν_kin
    return Re
end