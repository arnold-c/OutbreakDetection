#%%
using StatsBase: StatsBase
using OutbreakDetection: create_schematic_simulation,
    plot_schematic
using OutbreakDetectionUtils: StateParameters, DynamicsParameters, IndividualTestSpecification, OutbreakSpecification, OutbreakDetectionSpecification

#%%
states_p = StateParameters(
    ; N = 500_000,
    s_prop = 0.05,
    e_prop = 0.0,
    i_prop = 0.0
)

dynamics_p = DynamicsParameters(
    BETA_MEAN,
    BETA_FORCE,
    cos,
    SIGMA,
    GAMMA,
    MU,
    ANNUAL_BIRTHS_PER_K,
    EPSILON,
    R0,
    VACCINATION_COVERAGE,
)

noise_states_p = StateParameters(
    ; N = 500_000,
    s_prop = 0.15,
    e_prop = 0.0,
    i_prop = 0.0
)

noise_dynamics_p = DynamicsParameters(
    500_000,
    27,
    0.2,
    1 / 7,
    1 / 14,
    5.0,
    0.65;
    seasonality = sin,
)

test_specification = IndividualTestSpecification(
    0.85, 0.85, 0
)

time_p = SimTimeParameters(;
    tmin = 0.0, tmax = 365.0 * 20, tstep = 1.0
)

outbreak_specification = OutbreakSpecification(5, 30, 500)

movingavg_window = 20

outbreak_detection_specification = OutbreakDetectionSpecification(
    8,
    movingavg_window,
    1.0,
    0.75,
    "movingavg",
)

#%%
#%%
inc_vec, outbreak_status, outbreak_bounds, noise_vec, movingavg_testpositives, alertstatus_vec, alert_bounds = create_schematic_simulation(
    states_p,
    dynamics_p,
    noise_states_p,
    noise_dynamics_p,
    test_specification,
    time_p;
    seed = 12345,
    outbreak_specification = outbreak_specification,
    outbreak_detection_specification = outbreak_detection_specification,
    noise_scaling = 15,
    shift_noise = -100,
);

#%%
schematic_with_shade_fig = plot_schematic(
    inc_vec,
    outbreak_status,
    outbreak_bounds[:, 1:2],
    outbreak_specification,
    noise_vec,
    movingavg_testpositives,
    alertstatus_vec,
    alert_bounds,
    outbreak_detection_specification.alert_threshold;
    time_p = time_p,
    shade_alert_outbreak_overlap = true,
    xlims = (5, 13),
    measlesalpha = 0.4,
    testalpha = 0.4,
    outbreakcolormap = [
        "#5E5C6C",
        "#CE6F58",
    ],
    alertcolormap = [
        "#5E5C6C",
        "#2B3465",
    ],
)

save(
    plotsdir("schematic-plot.svg"),
    schematic_with_shade_fig,
)
