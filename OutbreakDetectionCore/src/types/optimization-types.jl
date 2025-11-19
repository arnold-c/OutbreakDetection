export OptimalThresholdCharacteristics, OptimizationMethods, MSO,
    AbstractOptimizationMethods

"""
    OptimalThresholdCharacteristics

Characteristics of the optimal detection threshold for a scenario.

# Fields
- `outbreak_threshold_chars::StructVector{OutbreakThresholdChars}`: Threshold characteristics
- `individual_test_specification::IndividualTestSpecification`: Test specification
- `noise_specification::NoiseSpecification`: Noise specification
- `percent_clinic_tested::AbstractFloat`: Proportion tested at clinic
- `alert_threshold::Real`: Optimal alert threshold
- `accuracy::AbstractFloat`: Detection accuracy at optimal threshold
"""
struct OptimalThresholdCharacteristics{
        T1 <: StructVector{<:OutbreakThresholdChars},
        T2 <: IndividualTestSpecification,
        T3 <: NoiseSpecification,
        T4 <: AbstractFloat,
        TReal <: Real,
    }
    outbreak_threshold_chars::T1
    individual_test_specification::T2
    noise_specification::T3
    percent_clinic_tested::T4
    alert_threshold::TReal
    accuracy::T4
end

"""
    AbstractOptimizationMethods

Abstract type for optimization methods.
"""
abstract type AbstractOptimizationMethods end

"""
    MSO

Multistart optimization method.
"""
struct MSO end

"""
    OptimizationMethods

Sum type for available optimization methods.
"""
LightSumTypes.@sumtype OptimizationMethods(MSO) <: AbstractOptimizationMethods
