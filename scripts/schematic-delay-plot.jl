#%%
using DrWatson
@quickactivate "OutbreakDetection"

using StatsBase: StatsBase

using OutbreakDetection: create_schematic_simulation, supplement_plots, trim_schematic_window, plot_test_delay_panel!

using OutbreakDetectionCore
using Dates: Day
using CairoMakie: Axis, Figure, hidexdecorations!, hlines!, lines!,
    linkxaxes!, save, text!, vspan!, ylims!
using Try

mkpath(supplement_plots())

include(
    DrWatson.scriptsdir("plotting-setup.jl")
)


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
    0.1  # vaccination_coverage
)
outbreak_specification = OutbreakSpecification(5, 30, 500)

# Manually construct OutbreakDetectionSpecification with all fields
alert_threshold = 6
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
delay_xlims = (7, 9)
perfect_tests = [
    ("Perfect test, 0-day delay", IndividualTestSpecification(1.0, 1.0, 0)),
    ("Perfect test, 14-day delay", IndividualTestSpecification(1.0, 1.0, 14)),
]

delay_panels = map(perfect_tests) do (label, test_specification)
    _, _, outbreak_bounds, _, movingavg_testpositives, alertstatus_vec,
        alert_bounds = create_schematic_simulation(
        measles_population_state_parameters,
        dynamics_parameters,
        dynamics_spec,
        noise_spec,
        test_specification,
        time_parameters;
        seed = 12345,
        outbreak_specification = outbreak_specification,
        outbreak_detection_specification = outbreak_detection_specification,
        noise_scaling = 1,
        shift_noise = 0,
    )

    return (
        label = label,
        panel_data = trim_schematic_window(
            movingavg_testpositives,
            outbreak_bounds,
            alert_bounds,
            time_parameters;
            xlims = delay_xlims,
        ),
        alert_status = alertstatus_vec,
    )
end

max_test_positives = maximum(
    maximum(panel.panel_data.values) for panel in delay_panels
)

delay_fig = Figure(size = (1300, 900))
delay_axes = [
    Axis(delay_fig[row, 1]) for row in eachindex(delay_panels)
]

for (index, (ax, panel)) in enumerate(zip(delay_axes, delay_panels))
    plot_test_delay_panel!(
        ax,
        panel.label,
        panel.panel_data,
        panel.alert_status,
        time_parameters,
        outbreak_detection_specification.alert_threshold;
        alertcolormap = [
            "#5E5C6C",
            "#2B3465",
        ],
        outbreakcolormap = [
            "#5E5C6C",
            "#CE6F58",
        ],
        measlesalpha = 0.4,
        testalpha = 0.4,
        linewidth = 5,
        ylabelsize = 28,
        thresholdfontsize = 26,
        xticklabelsize = 26,
        yticklabelsize = 26,
        show_xlabel = index == length(delay_axes),
    )
    ylims!(ax, 0, max_test_positives * 1.05)
end

linkxaxes!(delay_axes...)

save(
    supplement_plots("schematic-plot-test-delay.svg"),
    delay_fig,
)

save(
    supplement_plots("schematic-plot-test-delay.eps"),
    delay_fig,
    size = (2250, 1600),
    pt_per_unit = 72 / 300,
)
