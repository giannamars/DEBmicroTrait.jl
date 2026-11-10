const HYSOMI_BIOPAR_NAMES = (:mu, :K_L, :Y_mu, :S_tr, :r_react, :r_deact, :K_C_H, :vmax, :fm, :Y_m, :m_max, :beta, :tau)
const HYSOMI_SOILPAR_NAMES = (:k_a, :k_d, :C_S_max, :I_C_L, :I_C_H)
const HYSOMI_STATE_NAMES = (:B_a, :B_i, :C_l_L, :C_H, :C_s_L, :CO2)

const HYSOMI_BIOPAR_RANGES = [
    0.01 10.0;
    1.0 500.0;
    0.20 0.995;
    1e-7 0.3;
    0.01 100.0;
    0.01 100.0;
    10e-6 10.0;
    1e-5 1.0;
    0.2 0.9;
    0.25 0.995;
    1e-4 2.0;
    1e-4 1.0;
    1e-4 1.0;
]

const HYSOMI_SOILPAR_RANGES = [
    0.01 100.0;
    0.01 100.0;
    0.01 10.0;
    0.0 10.0;
    0.0 10.0;
]

struct HysomiParameters{B,S}
    biopar::B
    soilpar::S
end

HysomiParameters(values::NamedTuple) = HysomiParameters(values.biopar, values.soilpar)

function build_hysomi_parameters(biopar::NamedTuple, soilpar::NamedTuple)
    return HysomiParameters(biopar, soilpar)
end

function hysomi_biopar_vector(parameters::HysomiParameters)
    biopar = parameters.biopar
    return [getproperty(biopar, name) for name in HYSOMI_BIOPAR_NAMES]
end

function hysomi_soilpar_vector(parameters::HysomiParameters)
    soilpar = parameters.soilpar
    return [getproperty(soilpar, name) for name in HYSOMI_SOILPAR_NAMES]
end

function hysomi_parameter_vector(parameters::HysomiParameters)
    return vcat(hysomi_biopar_vector(parameters), hysomi_soilpar_vector(parameters))
end

struct HysomiModel{L,P}
    state_layout::L
    parameters::P
end

build_hysomi_model(state_layout::StateLayout, parameters::HysomiParameters) = HysomiModel(state_layout, parameters)

parameter_vector(model::HysomiModel) = hysomi_parameter_vector(model.parameters)

function hysomi_state_layout()
    return build_state_layout(collect(HYSOMI_STATE_NAMES), fill(1, length(HYSOMI_STATE_NAMES)))
end

function hysomi_rhs!(du, u, p, t)
    mu, K_L, Y_mu, S_tr, r_react, r_deact, K_C_H, vmax, fm, Y_m, m_max, beta, tau, k_a, k_d, C_S_max, I_C_L, I_C_H = p

    B_a = u[1]
    B_i = u[2]
    C_l_L = u[3]
    C_H = u[4]
    C_s_L = u[5]
    CO2 = u[6]

    activation = inv(one(eltype(u)) + exp((S_tr - C_l_L) / tau * S_tr))
    deactivation = one(eltype(u)) - activation

    r_a = r_react * B_i * activation
    r_d = r_deact * B_a * deactivation

    r_m_a_B = m_max * B_a
    r_m_a_M = ((m_max * C_l_L * K_L) / (m_max + C_l_L * K_L)) * B_a
    r_m_i_B = m_max * B_i * beta
    r_m_i_M = ((m_max * C_l_L * K_L) / (m_max + C_l_L * K_L)) * B_a * beta

    r_m_total = (r_m_a_B - r_m_a_M) + (r_m_i_B - r_m_i_M)
    r_m_sum = r_m_a_M + r_m_i_M
    growth = mu * B_a * (K_L * C_l_L / (mu + K_L * C_l_L))
    high_molecular_uptake = vmax * B_a * (C_H / (K_C_H + C_H))

    du[1] = growth + r_a - r_d - (1 / Y_m) * (r_m_a_B - r_m_a_M)
    du[2] = r_d - r_a - (1 / Y_m) * (r_m_i_B - r_m_i_M)
    du[3] = I_C_L + (1 - fm) * r_m_total - (1 / Y_mu) * growth - r_m_sum + high_molecular_uptake - k_a * C_l_L * (C_S_max - C_s_L) + k_d * C_s_L
    du[4] = I_C_H + fm * r_m_total - high_molecular_uptake
    du[5] = k_a * C_l_L * (C_S_max - C_s_L) - k_d * C_s_L
    du[6] = ((1 - Y_mu) / Y_mu) * growth + ((1 - Y_m) / Y_m) * r_m_total + r_m_sum

    return du
end

function hysomi_observables(u, model::HysomiModel)
    return copy(u)
end

hysomi_mass_balance(u, model::HysomiModel) = sum(u)
