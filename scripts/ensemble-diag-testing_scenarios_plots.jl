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
alertthreshold_vec = [collect(4:2:14)..., collect(18:4:30)...]

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
    DynamicsParameters(500_000, 10, 0.2; vaccination_coverage = 0.0),
    SimTimeParameters(;
        tmin = 0.0, tmax = 365.0 * 100, tstep = 1.0
    ),
    100,
)
noise_specification = NoiseSpecification("poisson", 1.0)
outbreak_specification = OutbreakSpecification(5, 30, 500)

moving_avg_detection_lag = 7
test_result_lag = 0
percent_visit_clinic = 0.6
percent_clinic_tested_vec = collect(0.2:0.2:1.0)

threshold_comparison_params = (
    sensitivity_vec = sensitivity_vec,
    specificity_vec = specificity_vec,
    alertthreshold_vec = alertthreshold_vec,
    ensemble_specification = ensemble_specification,
    noise_specification = noise_specification,
    outbreak_specification = outbreak_specification,
    moving_avg_detection_lag = moving_avg_detection_lag,
    test_result_lag = test_result_lag,
    percent_visit_clinic = percent_visit_clinic,
)

#%%
@showprogress for percent_clinic_tested in percent_clinic_tested_vec
    plot_all_threshold_comparisons(
        percent_clinic_tested, threshold_comparison_params
    )
end
