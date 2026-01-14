export ScenarioSpecificationVecs,
    OptimizationScenario


"""
    ScenarioSpecificationVecs
"""
Base.@kwdef struct ScenarioSpecificationVecs
    ensemble_specification_vec::Vector{EnsembleSpecification}
    noise_level_vec::Vector{Float64}
    noise_type_description_vec::Vector{Symbol}
    test_specification_vec::Vector{IndividualTestSpecification}
    percent_tested_vec::Vector{Float64}
    alert_method_vec::Vector{AlertMethod}
    accuracy_metric_vec::Vector{AccuracyMetric}
    threshold_bounds_vec::Vector{@NamedTuple{lower::Float64, upper::Float64}}
end

"""
	OptimizationScenario
"""
Base.@kwdef struct OptimizationScenario
    ensemble_specification::EnsembleSpecification
    noise_level::Float64
    noise_type_description::Symbol
    test_specification::IndividualTestSpecification
    percent_tested::Float64
    alert_method::AlertMethod
    accuracy_metric::AccuracyMetric
    threshold_bounds::@NamedTuple{lower::Float64, upper::Float64}

    function OptimizationScenario(
            ensemble_specification,
            noise_level,
            noise_type_description,
            test_specification,
            percent_tested,
            alert_method,
            accuracy_metric,
            threshold_bounds,
        )
        @assert noise_type_description in [:static, :dynamic] "The noise type must be either :static or :dynamic. Received $noise_type_description"
        @assert threshold_bounds.lower < threshold_bounds.upper
        @assert threshold_bounds.lower >= 0.0
        return new(
            ensemble_specification,
            noise_level,
            noise_type_description,
            test_specification,
            percent_tested,
            alert_method,
            accuracy_metric,
            threshold_bounds
        )
    end
end
