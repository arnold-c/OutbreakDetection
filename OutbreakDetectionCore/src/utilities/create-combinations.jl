export create_combinations_vec

"""
    create_combinations_vec(
    	custom_function,
    	combinations;
		init = custom_function[]
    )

Generate all combinations of input parameters and apply a custom function to each.

This utility function creates the Cartesian product of all input parameter combinations,
applies a user-defined function to each combination, and concatenates the results into
a single vector.

# Arguments
- `custom_function`: Function to apply to each combination. Should accept the elements
  of a combination as separate arguments and return a vector or array-like object.
- `combinations`: Iterable collection of parameter values to combine. Each element
  should be an iterable (e.g., vector, range) of values for one parameter.
- `init`: Initial value for the reduction operation (default: `custom_function[]`).
  Should be an empty instance of the expected return type.

# Returns
- Vector containing concatenated results from applying `custom_function` to all
  combinations of input parameters.

# Notes
- Uses `Iterators.product` to generate all combinations efficiently.
- Results are combined using `vcat`, so `custom_function` should return vector-like objects.
- TODO: make type stable

# Examples
```julia
# Generate parameter combinations for a simulation
function create_params(R0, beta_force)
    return [DynamicsParameters(R_0=R0, beta_force=beta_force)]
end

R0_values = [12.0, 16.0, 20.0]
beta_values = [0.1, 0.2, 0.3]

params = create_combinations_vec(
    create_params,
    (R0_values, beta_values)
)
# Returns vector with 9 parameter combinations
```

# See Also
- `Iterators.product`: Creates Cartesian product of iterables
- `mapreduce`: Applies function and reduces results
"""
function create_combinations_vec(
        custom_function,
        combinations;
        init = custom_function[]
    )
    combs = Iterators.product(combinations...)

    # TODO: make type stable
    return mapreduce(combination -> custom_function(combination...), vcat, combs; init = init)
end
