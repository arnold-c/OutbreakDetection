#%%
using OutbreakDetectionCore
using Dates
using StructArrays

#%%
time_parameters = SimTimeParameters(;
    tmin = 0.0,
    tmax = 365.0 * 20.0,
    tstep = 1.0,
)

births_per_k_pop = 27.0

burnin_target_Reff = 0.9

common_disease_dynamics_parameters = CommonDiseaseDynamicsParameters(;
    births_per_k_pop = births_per_k_pop,
    nsims = 10,
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
    min_vaccination_coverage = 0.7,
    max_vaccination_coverage = 0.9,
)

# Choose Rubella-like parameters for the dynamic noise in the measles simulations
# @papadopoulosEstimatesBasicReproduction2022 @RubellaCDCYellow.
rubella_dynamical_noise_spec = DynamicalNoiseSpecification(;
    R_0 = 5.0,
    latent_period = 7.0,
    duration_infection = 14.0,
    correlation = "in-phase",
    poisson_component = 0.15,
)

#%%
measles_ensemble_specification = create_ensemble_specifications(
    "measles",
    measles_population_state_parameters,
    time_parameters,
    measles_dynamics_parameters,
    common_disease_dynamics_parameters,
    rubella_dynamical_noise_spec,
)

#%%
# noise_level_vec = [1.0, 2.0, 4.0, 8.0]
noise_level_vec = [1.0]
# noise_type_description_vec = [:static, :dynamic]
noise_type_description_vec = [:static]

test_specification_vec = [
    IndividualTestSpecification(val, val, lag) for
        (val, lag) in zip((0.85, 0.9, 1.0, 1.0), (fill(0, 3)..., 14))
]
percent_tested_vec = collect(0.1:0.1:1.0)

#%%
alert_method_vec = AlertMethod[AlertMethod(MovingAverage())]
accuracy_metric_vec = AccuracyMetric[AccuracyMetric(BalancedAccuracy())]
threshold_bounds = (; lower = 0.0, upper = 20.0)

#%%
specification_vecs = ScenarioSpecificationVecs(;
    ensemble_specification_vec = [measles_ensemble_specification],
    noise_level_vec = noise_level_vec,
    noise_type_description_vec = noise_type_description_vec,
    test_specification_vec = test_specification_vec,
    percent_tested_vec = percent_tested_vec,
    alert_method_vec = alert_method_vec,
    accuracy_metric_vec = accuracy_metric_vec,
    threshold_bounds_vec = [threshold_bounds],
)

#%%
create_scenarios_structvector(specification_vecs)

#%%
optimized_threshold_results = run_scenario_threshold_optimization(
    specification_vecs;
    force = true,
    save_results = true,
    save_checkpoints = true,
    disable_time_check = true,
)
