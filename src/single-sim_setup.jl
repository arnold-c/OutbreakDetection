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

singlesim_dynamics_p = DynamicsParameters(
    BETA_MEAN,
    BETA_FORCE,
    SIGMA,
    GAMMA,
    MU,
    ANNUAL_BIRTHS_PER_K,
    EPSILON,
    R0,
    VACCINATION_COVERAGE,
)

jldsave(
    "data/singlesim/single-sim_setup.jld2";
    singlesim_states_p,
    singlesim_time_p,
    singlesim_dynamics_p,
)
