#%%
using DrWatson
@quickactivate "OutbreakDetection"

using OutbreakDetectionUtils:
    IndividualTestSpecification, DynamicalNoiseSpecification
using StatsBase: mean


include(projectdir("manuscript", "optimal-thresholds-loading.jl"));

if false
    include("optimal-thresholds-loading.jl")
end

#%%
mean_elisa_0d_delays = map(
    test_chars -> mean(vcat(getproperty(test_chars, :detectiondelays)...)),
    filter(
        chars ->
            chars.individual_test_specification ==
            IndividualTestSpecification(1.0, 1.0, 0),
        dynamical_noise_optimal_solutions,
    ).outbreak_threshold_chars,
)

mean_rdt_90_8x_dynamical_delays = map(
    test_chars -> mean(vcat(getproperty(test_chars, :detectiondelays)...)),
    filter(
        chars ->
            chars.individual_test_specification ==
            IndividualTestSpecification(0.9, 0.9, 0),
        dynamical_noise_optimal_solutions,
    ).outbreak_threshold_chars,
)

mean_rdt_85_8x_dynamical_delays = map(
    test_chars -> mean(vcat(getproperty(test_chars, :detectiondelays)...)),
    filter(
        chars ->
            chars.individual_test_specification ==
            IndividualTestSpecification(0.85, 0.85, 0),
        dynamical_noise_optimal_solutions,
    ).outbreak_threshold_chars,
)

mapreduce(
    vcat,
    (
        ("ELISA", mean_elisa_0d_delays),
        ("90%", mean_rdt_90_8x_dynamical_delays),
        ("85%", mean_rdt_85_8x_dynamical_delays),
    ),
) do (label, mean_delays_vec)
    Dict(label => round.(extrema(mean_delays_vec); digits = 1))
end

#%%
mean_elisa_0d_alert_duration = map(
    test_chars -> mean(vcat(getproperty(test_chars, :alert_duration_vec)...)),
    filter(
        chars ->
            chars.individual_test_specification ==
            IndividualTestSpecification(1.0, 1.0, 0),
        dynamical_noise_optimal_solutions,
    ).outbreak_threshold_chars,
)

mean_rdt_90_8x_dynamical_alert_duration = map(
    test_chars -> mean(vcat(getproperty(test_chars, :alert_duration_vec)...)),
    filter(
        chars ->
            chars.individual_test_specification ==
            IndividualTestSpecification(0.9, 0.9, 0),
        dynamical_noise_optimal_solutions,
    ).outbreak_threshold_chars,
)

mean_rdt_85_8x_dynamical_alert_duration = map(
    test_chars -> mean(vcat(getproperty(test_chars, :alert_duration_vec)...)),
    filter(
        chars ->
            chars.individual_test_specification ==
            IndividualTestSpecification(0.85, 0.85, 0),
        dynamical_noise_optimal_solutions,
    ).outbreak_threshold_chars,
)

mapreduce(
    vcat,
    (
        ("ELISA", mean_elisa_0d_alert_duration),
        ("90%", mean_rdt_90_8x_dynamical_alert_duration),
        ("85%", mean_rdt_85_8x_dynamical_alert_duration),
    ),
) do (label, mean_alert_duration_vec)
    Dict(label => round.(extrema(mean_alert_duration_vec); digits = 1))
end
