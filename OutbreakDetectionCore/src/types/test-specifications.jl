export IndividualTestSpecification,
    TestPositivity,
    get_test_description,
    table_test_type,
    plot_test_description

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
    sens = test_specification.sensitivity
    spec = test_specification.specificity
    lag = test_specification.test_result_lag

    description = if sens == 1.0 && spec == 0.0 && lag == 0
        "Clinical case definition"
    elseif (sens == 0.98 && spec == 0.98) || (sens == 0.95 && spec == 0.98)
        "ELISA-like ($(sens * 100)% sens, $(spec * 100)% spec, $(lag) day lag)"
    elseif sens == 1.0 && spec == 1.0
        "Perfect test ($(lag) day lag)"
    elseif sens == spec && sens < 1.0 && lag == 0
        "RDT-like ($(sens * 100)% sens/spec)"
    else
        error("Unknown test specification pattern: sensitivity=$(sens), specificity=$(spec), lag=$(lag)")
    end

    return description
end

"""
    table_test_type(sensitivity, specificity, test_lag)

Get a table-formatted test type description.
"""
function table_test_type(
        sensitivity::Float64,
        specificity::Float64,
        test_lag::Int64
    )
    return if sensitivity == 1.0 && specificity == 0.0 && test_lag == 0
        "Clinical Case Definition"
    elseif sensitivity == specificity && sensitivity < 1.0 && test_lag == 0
        "Imperfect Test ($(Int64(round(sensitivity * 100; digits = 0)))%)"
    elseif sensitivity == 1.0 && specificity == 1.0
        "Perfect Test"
    else
        error("Unknown test type pattern: sensitivity=$(sensitivity), specificity=$(specificity), test_lag=$(test_lag)")
    end
end

"""
    plot_test_description(test_specification::IndividualTestSpecification)

Get a plot-formatted test description.
"""
function plot_test_description(test_specification::IndividualTestSpecification)
    sens = test_specification.sensitivity
    spec = test_specification.specificity
    lag = test_specification.test_result_lag

    return if sens == 1.0 && spec == 0.0 && lag == 0
        "Clinical Case Definition"
    elseif sens == spec && sens < 1.0 && lag == 0
        "Imperfect Test ($(Int64(round(sens * 100; digits = 0)))%)"
    elseif sens == 1.0 && spec == 1.0
        "Perfect Test ($lag-day lag)"
    else
        error("Unknown test specification pattern: sensitivity=$(sens), specificity=$(spec), lag=$(lag)")
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
