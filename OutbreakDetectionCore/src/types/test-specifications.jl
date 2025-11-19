export IndividualTestSpecification, TestPositivity,
    get_test_description, table_test_type, plot_test_description

"""
    IndividualTestSpecification

Specification for individual diagnostic test characteristics.

# Fields
- `sensitivity::AbstractFloat`: Test sensitivity (true positive rate)
- `specificity::AbstractFloat`: Test specificity (true negative rate)
- `test_result_lag::Integer`: Lag in days for test results
"""
struct IndividualTestSpecification{T1 <: AbstractFloat, T2 <: Integer}
    sensitivity::T1
    specificity::T1
    test_result_lag::T2
end

"""
    get_test_description(test_specification::IndividualTestSpecification)

Get a human-readable description of the test specification.
"""
function get_test_description(test_specification::IndividualTestSpecification)
    description = Match.@match test_specification begin
        IndividualTestSpecification(1.0, 0.0, 0) => "Clinical case definition"
        IndividualTestSpecification(0.98, 0.98, x::Int) || IndividualTestSpecification(0.95, 0.98, x::Int) => "ELISA-like ($(test_specification.sensitivity * 100)% sens, $(test_specification.specificity * 100)% spec, $(x) day lag)"
        IndividualTestSpecification(1.0, 1.0, x::Int) => "Perfect test ($(x) day lag)"
        IndividualTestSpecification(x::AbstractFloat, x::AbstractFloat, 0) where {x < 1.0} => "RDT-like ($(test_specification.sensitivity * 100)% sens/spec)"
    end
    return description
end

"""
    table_test_type(sensitivity, specificity, test_lag)

Get a table-formatted test type description.
"""
function table_test_type(sensitivity, specificity, test_lag)
    return Match.@match (sensitivity, specificity, test_lag) begin
        (1.0, 0.0, 0) => "Clinical Case Definition"
        (x::AbstractFloat, x::AbstractFloat, 0) where {x < 1.0} => "Imperfect Test ($(Int64(round(sensitivity * 100; digits = 0)))%)"
        (1.0, 1.0, x::Int) => "Perfect Test"
    end
end

"""
    plot_test_description(test_specification::IndividualTestSpecification)

Get a plot-formatted test description.
"""
function plot_test_description(test_specification::IndividualTestSpecification)
    return Match.@match test_specification begin
        IndividualTestSpecification(1.0, 0.0, 0) => "Clinical Case Definition"
        IndividualTestSpecification(x::AbstractFloat, x::AbstractFloat, 0) where {x < 1.0} => "Imperfect Test ($(Int64(round(x * 100; digits = 0)))%)"
        IndividualTestSpecification(1.0, 1.0, x::Int) => "Perfect Test ($x-day lag)"
    end
end

"""
    TestPositivity

Test positivity rates over different time windows.

# Fields
- `one_day::AbstractArray{AbstractFloat}`: 1-day test positivity
- `seven_day::AbstractArray{AbstractFloat}`: 7-day test positivity
- `fourteen_day::AbstractArray{AbstractFloat}`: 14-day test positivity
- `thirty_day::AbstractArray{AbstractFloat}`: 30-day test positivity
"""
struct TestPositivity{T1 <: AbstractArray{<:AbstractFloat}}
    one_day::T1
    seven_day::T1
    fourteen_day::T1
    thirty_day::T1
end

function TestPositivity(true_positive_vec, total_test_vec, alert_vec)
    return TestPositivity(
        calculate_test_positivity(
            true_positive_vec, total_test_vec, alert_vec, 1
        ),
        calculate_test_positivity(
            true_positive_vec, total_test_vec, alert_vec, 7
        ),
        calculate_test_positivity(
            true_positive_vec, total_test_vec, alert_vec, 14
        ),
        calculate_test_positivity(
            true_positive_vec, total_test_vec, alert_vec, 30
        ),
    )
end
