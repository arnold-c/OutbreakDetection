#%%
using DrWatson
@quickactivate "OutbreakDetection"

using ProgressMeter
using FLoops
using NaNMath: NaNMath

using OutbreakDetection

includet(srcdir("makie-plotting-setup.jl"))

#%%
sensitivity_vec = collect(0.8:0.2:1.0)
specificity_vec = collect(0.8:0.2:1.0)
detectthreshold_vec = collect(4:2:14)

#%%
ensemble_scenario_spec_vec = Vector{ScenarioSpecification}(
    undef,
    length(sensitivity_vec) * length(detectthreshold_vec) +
    length(detectthreshold_vec),
)
ensemble_chars_vec = Vector(
    undef, length(ensemble_scenario_spec_vec)
)

#%%
ensemble_specification = EnsembleSpecification(
    ("seasonal-infectivity-import", "tau-leaping"),
    StateParameters(
        500_000,
        Dict(
            :s_prop => 0.1,
            :e_prop => 0.0,
            :i_prop => 0.0,
            :r_prop => 0.9,
        ),
    ),
    DynamicsParameters(500_000, 10, 0.2),
    SimTimeParameters(;
        tmin = 0.0, tmax = 365.0 * 100, tstep = 1.0
    ),
    100,
)
noise_specification = NoiseSpecification("poisson", 1.0)
outbreak_specification = OutbreakSpecification(5, 30, 500)

moving_avg_detection_lag = 7
test_result_lag = 0

#%%
percent_visit_clinic = collect(0.2:0.2:1.0)
