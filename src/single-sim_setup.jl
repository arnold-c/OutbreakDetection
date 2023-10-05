using DrWatson
@quickactivate "OutbreakDetection"

using JLD2
using OutbreakDetection

singlesim_states_p = StateParameters(;
    N = 500_000,
    s_prop = 0.1,
    e_prop = 0.00,
    i_prop = 0.00
)

singlesim_time_p = SimTimeParameters(;
    tmin = 0.0, tmax = 365.0 * 100, tstep = 1.0
)

const LATENT_PER_DAYS = 8
const DUR_INF_DAYS = 5
const R_0 = 10.0
const LIFE_EXPECTANCY_YEARS = 62.5
const VACCINATION_COVERAGE = 0.0

const SIGMA = 1 / LATENT_PER_DAYS
const GAMMA = 1 / DUR_INF_DAYS
const MU = 1 / (LIFE_EXPECTANCY_YEARS * 365)
const ANNUAL_BIRTHS_PER_K = 1000 / LIFE_EXPECTANCY_YEARS

# beta_mean is the average transmission rate
const BETA_MEAN = calculate_beta(R_0, GAMMA, MU, 1, singlesim_states_p.init_states.N)

# Adjust the scale of the seasonal variation in infectivity i.e. beta_force scales the amplitude of cosine function
const BETA_FORCE = 0.2
const EPSILON = calculate_import_rate(MU, R_0, singlesim_states_p.init_states.N)

singlesim_dynamics_p = DynamicsParameters(
    BETA_MEAN,
    BETA_FORCE,
    SIGMA,
    GAMMA,
    MU,
    ANNUAL_BIRTHS_PER_K,
    EPSILON,
    R_0,
    VACCINATION_COVERAGE
)

jldsave("data/singlesim/single-sim_setup.jld2"; singlesim_states_p, singlesim_time_p, singlesim_dynamics_p)
