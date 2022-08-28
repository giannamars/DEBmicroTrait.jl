function max_specific_death_rate(Min_gen_time::Vector{Float64})
    # Biselli et al., 2020
    gmax     = log(2)./Min_gen_time
    γ_V_0    = 0.23.*exp.(0.88.*gmax)
end
