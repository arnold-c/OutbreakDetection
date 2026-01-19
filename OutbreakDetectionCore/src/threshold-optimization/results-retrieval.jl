export get_most_recent_hyperparam_filepath

"""
    get_most_recent_hyperparam_filepath(
        filename_base,
        filedir,
    ) -> Try.Result{String}

Retrieve the filepath of the most recent hyperparameter optimization file.

This function searches a directory for JLD2 files matching a specific filename pattern
and returns the path to the most recent file based on the datetime prefix in the filename.
It uses the `Try.jl` error handling pattern to safely handle cases where no matching
files are found or datetime parsing fails.

The function expects files to follow the naming convention:
`<datetime>_<filename_base>.jld2` where `<datetime>` is a parseable `DateTime` string.

# Arguments

  - `filename_base`: Base name of the optimization file (without datetime prefix or extension)
  - `filedir`: Directory path containing the optimization files

# Returns

  - `Try.Ok{String}`: Path to the most recent optimization file if found
  - `Try.Err{String}`: Error message if no valid optimization files are found

# Algorithm

  1. Validates that the directory exists (assertion)
  2. Reads all files in the directory
  3. Filters files matching the pattern `*<filename_base>.jld2`
  4. Extracts datetime prefixes from matching filenames
  5. Parses datetime strings and filters out invalid entries
  6. Sorts valid datetimes and selects the most recent
  7. Constructs and returns the full filepath

# Example

```julia
# Search for most recent hyperparameter file
result = get_most_recent_hyperparam_filepath(
    "threshold_optimization",
    "/path/to/results"
)

# Handle the result
if Try.isok(result)
    filepath = Try.unwrap(result)
    data = JLD2.load(filepath)
else
    error_msg = Try.unwrap_err(result)
    @warn "Failed to find optimization file: \$error_msg"
end
```

# Error Cases

Returns `Try.Err` in the following situations:
  - No files exist in the directory
  - No files match the `<filename_base>.jld2` pattern
  - All matching files have unparseable datetime prefixes

# See Also

  - [`Try.Ok`](@ref): Success wrapper for valid results
  - [`Try.Err`](@ref): Error wrapper for failure cases
  - [`Try.unwrap`](@ref): Extract value from `Try.Ok`
  - [`Try.unwrap_err`](@ref): Extract error message from `Try.Err`
"""
function get_most_recent_hyperparam_filepath(
        filename_base,
        filedir,
    )
    @assert isdir(filedir)
    optimization_files = readdir(filedir)

    if length(optimization_files) == 0
        return Try.Err("No optimization files found.")
    end

    filter_regex = Regex("(.*)$(filename_base)\\.jld2\$")

    filtered_optimization_files = filter(
        f -> contains(f, filter_regex),
        optimization_files,
    )

    if length(filtered_optimization_files) == 0
        return Try.Err("No optimization files found.")
    end

    filtered_optimization_datetimes = Vector{Union{Try.Ok, Try.Err}}(
        undef, length(filtered_optimization_files)
    )

    for (i, f) in pairs(filtered_optimization_files)
        matches = match(filter_regex, f)
        if isnothing(matches)
            filtered_optimization_datetimes[i] = Try.Err(
                "No matches for filename $(f)"
            )
            continue
        end

        filtered_optimization_datetimes[i] = Try.Ok(
            tryparse(
                Dates.DateTime,
                strip(
                    matches[1],
                    '_',
                ),
            ),
        )
    end

    filtered_optimization_datetimes = filter(
        Try.isok, filtered_optimization_datetimes
    )

    if length(filtered_optimization_datetimes) == 0
        return Try.Err("No optimization files found.")
    end

    most_recent_optimization_datetime = sort(
        Try.unwrap.(filtered_optimization_datetimes)
    )[end]

    most_recent_filepath = joinpath(
        filedir,
        string(most_recent_optimization_datetime) *
            "_$(filename_base).jld2",
    )
    return Try.Ok(most_recent_filepath)
end
