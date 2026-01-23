export sort_test_specifications

"""
    sort_test_specifications(test_specifications)

Sort test specifications first by test_result_lag (descending), then by
(sensitivity, specificity) tuple (ascending).

# Arguments
- `test_specifications`: Collection of test specifications to sort

# Returns
- Sorted collection of test specifications

# Example
```julia
specs = [
    IndividualTestSpecification(sensitivity=0.9, specificity=0.95, test_result_lag=1),
    IndividualTestSpecification(sensitivity=0.8, specificity=0.90, test_result_lag=2),
]
sorted_specs = sort_test_specifications(specs)
```
"""
function sort_test_specifications(test_specifications)
    return sort(
        sort(
            test_specifications;
            by = t -> t.test_result_lag,
            rev = true,
        );
        by = t -> (t.sensitivity, t.specificity),
        rev = false,
    )
end
