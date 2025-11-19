export ScenarioSpecification

"""
    ScenarioSpecification

Complete specification for a detection scenario.

# Fields
- `ensemble_specification::EnsembleSpecification`: Ensemble parameters
- `outbreak_specification::OutbreakSpecification`: Outbreak definition
- `noise_specification::NoiseSpecification`: Noise characteristics
- `outbreak_detection_specification::OutbreakDetectionSpecification`: Detection parameters
- `individual_test_specification::IndividualTestSpecification`: Test characteristics
- `dirpath::AbstractString`: Directory path for output

# Constructor
    ScenarioSpecification(ensemble_specification, outbreak_specification,
                         noise_specification, outbreak_detection_specification,
                         individual_test_specification)

Creates a `ScenarioSpecification` with automatically generated directory path.
"""
struct ScenarioSpecification{
        T1 <: EnsembleSpecification,
        T2 <: OutbreakSpecification,
        T3 <: NoiseSpecification,
        T4 <: OutbreakDetectionSpecification,
        T5 <: IndividualTestSpecification,
        T6 <: AbstractString,
    }
    ensemble_specification::T1
    outbreak_specification::T2
    noise_specification::T3
    outbreak_detection_specification::T4
    individual_test_specification::T5
    dirpath::T6
end

function ScenarioSpecification(
        ensemble_specification::EnsembleSpecification,
        outbreak_specification::OutbreakSpecification,
        noise_specification::NoiseSpecification,
        outbreak_detection_specification::OutbreakDetectionSpecification,
        individual_test_specification::IndividualTestSpecification,
    )
    dirpath = joinpath(
        ensemble_specification.dirpath,
        outbreak_specification.dirpath,
        getdirpath(noise_specification),
        outbreak_detection_specification.dirpath,
        "testsens_$(individual_test_specification.sensitivity)",
        "testspec_$(individual_test_specification.specificity)",
        "testlag_$(individual_test_specification.test_result_lag)",
    )

    return ScenarioSpecification(
        ensemble_specification,
        outbreak_specification,
        noise_specification,
        outbreak_detection_specification,
        individual_test_specification,
        dirpath,
    )
end
