export IndividualTestSpecification

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
