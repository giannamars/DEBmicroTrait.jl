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
