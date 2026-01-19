export reshape_optim_df_to_matrix, sort_noise_specifications

"""
    reshape_optim_df_to_matrix(optim_df)

Reshape optimization DataFrame into a matrix organized by noise specifications.
"""
function reshape_optim_df_to_matrix(
        optim_df::T1
    ) where {T1 <: DataFrames.DataFrame}
    noise_spec_vec = unique(optim_df.noise_spec)
    unique_noise_descriptions = unique(get_noise_description.(noise_spec_vec))
    num_noise_descriptions = length(unique_noise_descriptions)

    shape_noise_specifications =
        map(unique_noise_descriptions) do noise_description
        filter(
            noise_spec ->
            noise_description == get_noise_description(noise_spec),
            noise_spec_vec,
        ) |>
            v -> sort_noise_specifications(convert(Vector{typeof(v[1])}, v))
    end

    num_shape_noise_specifications = length(shape_noise_specifications[1])

    @assert length(vcat(shape_noise_specifications...)) ==
        num_noise_descriptions * num_shape_noise_specifications

    optim_arr = Array{Any}(
        undef,
        num_noise_descriptions,
        num_shape_noise_specifications,
    )

    unique_test_specifications = sort(
        sort(
            unique(optim_df.test_spec);
            by = t -> t.test_result_lag,
            rev = true,
        );
        by = t -> (t.sensitivity, t.specificity),
        rev = false,
    )

    for i in eachindex(unique_noise_descriptions)
        for (j, noise_spec) in pairs(shape_noise_specifications[i])
            subset_df = DataFrames.subset(
                optim_df,
                :noise_spec => DataFrames.ByRow(==(noise_spec)),
            )

            test_order = mapreduce(
                t -> findall(isequal(t).(subset_df.test_spec)),
                vcat,
                unique_test_specifications,
            )

            subset_df = subset_df[test_order, :]

            optim_arr[i, j] = StructArrays.StructVector(
                OptimalThresholdCharacteristics.(
                    subset_df.OT_chars,
                    subset_df.test_spec,
                    subset_df.noise_spec,
                    getfield.(
                        subset_df.outbreak_detection_spec,
                        :percent_clinic_tested,
                    ),
                    subset_df.optimal_threshold,
                    subset_df.optimal_accuracy,
                ),
            )
        end
    end

    return optim_arr
end

"""
    sort_noise_specifications(v::AbstractVector{PoissonNoiseSpecification})

Sort Poisson noise specifications by noise magnitude.
"""
function sort_noise_specifications(
        v::AbstractVector{T}
    ) where {T <: PoissonNoise}
    noise_magnitudes = getproperty.(v, :noise_mean_scaling)

    noise_orders = sortperm(noise_magnitudes)

    return v[noise_orders]
end

"""
    sort_noise_specifications(v::AbstractVector{DynamicalNoiseSpecification})

Sort dynamical noise specifications by vaccination coverage (descending).
"""
function sort_noise_specifications(
        v::AbstractVector{T}
    ) where {T <: DynamicalNoise}
    noise_magnitudes = getproperty.(v, :vaccination_coverage)

    noise_orders = sortperm(noise_magnitudes; rev = true)

    return v[noise_orders]
end
