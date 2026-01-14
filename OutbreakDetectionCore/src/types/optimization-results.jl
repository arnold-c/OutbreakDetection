export OptimizationResult,
    AlertClassificationResults,
    OptimizedValues

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
    optimal_threshold::Float64
    accuracy::Float64
    proportion_outbreaks_detected::Float64
    proportion_alerts_correct::Float64
    mean_detection_delay::Float64
end

"""
    AlertClassificationResults

Stores binary classification results from outbreak detection analysis.

Contains the confusion matrix components and total counts for evaluating detection performance
across all outbreaks in each simulation.

# Fields
- `true_positives::Float64`: Number of outbreaks correctly identified by alert criteria
- `false_positives::Float64`: Number of alerts that do not correspond to an outbreak
- `false_negatives::Float64`: Number of outbreaks not identified by alert criteria
- `nsims::Int64`: Total number of simulations

# Notes
The classification counts are stored as Float64 to support weighted or fractional classifications,
while the total simulation counts remain as integers. This struct serves as an intermediate
representation for calculating performance metrics like sensitivity, specificity, and accuracy.
"""
Base.@kwdef struct AlertClassificationResults
    true_positives::Float64
    false_positives::Float64
    false_negatives::Float64
    nsims::Int64
end

"""
    OptimizedValues

Optimal hyperparameter values and performance metrics from alert threshold optimization.

This struct stores the results of threshold optimization for outbreak detection, containing
both the optimal threshold value found during optimization and the corresponding performance
metrics achieved with that threshold.

# Fields
- `alert_threshold::Float64`: Optimal alert threshold determination
- `accuracy::Float64`: Overall classification accuracy achieved with optimal parameters (0.0 to 1.0)
- `proportion_outbreaks_detected::Float64`: True positive rate (sensitivity) given the optimal threshold
- `proportion_alerts_correct::Float64`: Proportion of alerts that are associated with an outbreak (Positive predictive value) given the optimal threshold
- `mean_detection_delay::Float64`: The mean delay in days between the start of the outbreak and when the alert is triggered

# Example
"""
Base.@kwdef struct OptimizedValues
    alert_threshold::Float64
    accuracy::Float64
    proportion_outbreaks_detected::Float64
    proportion_alerts_correct::Float64
    mean_detection_delay::Float64
end
