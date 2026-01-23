export OptimizationResult

"""
    OptimizationResult

Complete results from EWS hyperparameter optimization including scenario details and performance.

This struct provides a comprehensive record of an optimization run, storing both the complete
scenario configuration used during optimization and the optimal hyperparameter values and
performance metrics discovered. It serves as the primary output format for optimization
functions and enables reproducible analysis of optimization results.

The result includes all parameters needed to recreate the optimization scenario, making it
suitable for result storage, comparison across different scenarios, and further analysis.
This is also used when new scenarios are prepared to be run, checking if the results already
exist, and if they do, load the existing results to avoid redundant computations.

# Fields
"""
Base.@kwdef struct OptimizationResult
    ensemble_specification::EnsembleSpecification
    noise_level::Float64
    noise_type_description::Symbol
    test_specification::IndividualTestSpecification
    percent_tested::Float64
    alert_method::AlertMethod
    accuracy_metric::AccuracyMetric
    threshold_bounds::@NamedTuple{lower::Float64, upper::Float64}
    outbreak_specification::OutbreakSpecification
    optimal_threshold::Float64
    accuracies::Vector{Float64}
    proportion_alerts_correct::Vector{Float64}
    proportion_outbreaks_detected::Vector{Float64}
    detection_delays::Vector{Int64}
    unavoidable_cases::Vector{Int64}
    alert_durations::Vector{Int64}
    outbreak_durations::Vector{Int64}
    proportion_timeseries_in_alert::Vector{Float64}
    proportion_timeseries_in_outbreak::Vector{Float64}
end
