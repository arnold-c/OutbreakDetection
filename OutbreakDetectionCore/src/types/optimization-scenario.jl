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
    function OptimizationScenario(
            ensemble_specification,
            noise_level,
            noise_type_description,
            test_specification,
            percent_tested,
            alert_method,
        )
        @assert noise_type_description in [:static, :dynamic] "The noise type must be either :static or :dynamic. Received $noise_type_description"
        return new(
            ensemble_specification,
            noise_level,
            noise_type_description,
            test_specification,
            percent_tested,
            alert_method,
        )
    end
end
