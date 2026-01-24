#%%
using DrWatson
@quickactivate "OutbreakDetection"

using StatsBase: StatsBase

using OutbreakDetection: create_schematic_simulation,
    plot_schematic

using OutbreakDetectionCore
using Dates: Day
using CairoMakie: save

#%%
time_parameters = SimTimeParameters(;
    tmin = 0.0,
    tmax = 365.0 * 100.0,
    tstep = 1.0,
)

births_per_k_pop = 27.0

common_disease_dynamics_parameters = CommonDiseaseDynamicsParameters(;
    births_per_k_pop = births_per_k_pop,
    nsims = 1,
)

#%%
measles_R0 = 16.0
measles_initial_s_prop = 0.05

# Set initial population proportions so that initial Reff = 0.8
measles_population_state_parameters = StateParameters(
    500_000,
    Dict(
        :s_prop => measles_initial_s_prop,
        :e_prop => 0.0,
        :i_prop => 0.0,
        :r_prop => 1.0 - measles_initial_s_prop,
    ),
)

# @guerraBasicReproductionNumber2017 @gastanaduyMeasles2019
measles_dynamics_parameters = TargetDiseaseDynamicsParameters(;
    R_0 = measles_R0,
    latent_period = Day(10.0),
    infectious_duration = Day(8.0),
    beta_force = 0.2,
    seasonality = SeasonalityFunction(CosineSeasonality()),
    min_vaccination_coverage = 0.8,
    max_vaccination_coverage = 0.8,
)


# Choose Rubella-like parameters for the dynamic noise in the measles simulations
# @papadopoulosEstimatesBasicReproduction2022 @RubellaCDCYellow.
rubella_dynamical_noise_parameters = DynamicalNoiseParameters(;
    R_0 = 5.0,
    latent_period = Day(7.0),
    infectious_duration = Day(14.0),
    correlation = "in-phase",
    poisson_component = 0.15,
)

dynamics_spec = DynamicsParameterSpecification(
    measles_population_state_parameters,
    measles_dynamics_parameters,
    common_disease_dynamics_parameters
)

dynamics_parameters = DynamicsParameters(dynamics_spec; seed = 12345)

# Create noise specification with fixed vaccination coverage
noise_spec = DynamicalNoiseSpecification(
    rubella_dynamical_noise_parameters,
    0.6  # vaccination_coverage
)
test_specification = IndividualTestSpecification(
    0.85, 0.85, 0
)

outbreak_specification = OutbreakSpecification(5, 30, 500)

# Manually construct OutbreakDetectionSpecification with all fields
alert_threshold = 8
percent_visit_clinic = 1.0
percent_clinic_tested = 0.75
alert_method = OutbreakDetectionCore.AlertMethod(OutbreakDetectionCore.MovingAverage(28))
percent_tested = percent_visit_clinic * percent_clinic_tested
dirpath = joinpath(
    "alertmethod_MovingAverage",
    "alertthreshold_$(alert_threshold)",
    "moveavglag_$(alert_method.window)",
    "perc_visit_clinic_$(percent_visit_clinic)",
    "perc_clinic_tested_$(percent_clinic_tested)"
)

outbreak_detection_specification = OutbreakDetectionSpecification(
    alert_threshold,
    alert_method.window,
    percent_visit_clinic,
    percent_clinic_tested,
    percent_tested,
    alert_method,
    dirpath,
)

#%%
inc_vec, outbreak_status, outbreak_bounds, noise_vec, movingavg_testpositives, alertstatus_vec, alert_bounds = create_schematic_simulation(
    measles_population_state_parameters,
    dynamics_parameters,
    dynamics_spec,
    noise_spec,
    test_specification,
    time_parameters;
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
    time_p = time_parameters,
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
    DrWatson.plotsdir("schematic-plot.svg"),
    schematic_with_shade_fig,
)
