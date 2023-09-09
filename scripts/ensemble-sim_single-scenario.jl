#%%
using DrWatson
@quickactivate "OutbreakDetection"

using UnPack

# include("ensemble-sim.jl")
include("../src/OutbreakDetection.jl")
using .OutbreakDetection

#%%
ensemble_single_scenario_spec = ScenarioSpecification(
    EnsembleSpecification(
        ("seasonal-infectivity-import", "tau-leaping"),
        StateParameters(
            500_000,
            Dict(
                :s_prop => 0.1,
                :e_prop => 0.01,
                :i_prop => 0.01,
                :r_prop => 0.88,
            ),
        ),
        DynamicsParameters(500_000, 10, 0.2),
        SimTimeParameters(; tmin = 0.0, tmax = 365.0 * 100, tstep = 1.0),
        1_000,
    ),
    OutbreakSpecification(5, 30, 500),
    create_static_NoiseSpecification(
        [10.0],
        SimTimeParameters(; tmin = 0.0, tmax = 365.0 * 100, tstep = 1.0),
        0.0,
        0.1,
        1_000,
    ),
    OutbreakDetectionSpecification(10, 7, 0.3, 0.3, 3),
    IndividualTestSpecification(0.8, 0.8),
)

#%%
ensemble_single_scenario_sol = get_ensemble_file(
    "solution", ensemble_single_scenario_spec.ensemble_specification
)

ensemble_single_scenario_quantiles = get_ensemble_file(
    "95", ensemble_single_scenario_spec.ensemble_specification
)

ensemble_single_scenario_detection = get_scenario_file(
    "scenario", ensemble_single_scenario_spec
)

#%%
@unpack ensemble_seir_arr,
dynamics_parameters,
time_parameters = ensemble_single_scenario_sol

@unpack ensemble_seir_summary, caption = ensemble_single_scenario_quantiles

@unpack OT_Chars, incarr, testarr = ensemble_single_scenario_detection
