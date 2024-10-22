using StatsBase: StatsBase

function create_inc_infec_arr(
    ensemble_inc_vecs, outbreak_specification::OutbreakSpecification
)
    ensemble_inc_arr = zeros(
        Int64, size(ensemble_inc_vecs, 1), 3, size(ensemble_inc_vecs, 2)
    )

    ensemble_thresholds_vec = Vector{Array{Int64,2}}(
        undef, size(ensemble_inc_vecs, 2)
    )

    create_inc_infec_arr!(
        ensemble_inc_arr,
        ensemble_thresholds_vec,
        ensemble_inc_vecs,
        outbreak_specification.outbreak_threshold,
        outbreak_specification.minimum_outbreak_duration,
        outbreak_specification.minimum_outbreak_size,
    )

    return ensemble_inc_arr, ensemble_thresholds_vec
end

function create_inc_infec_arr!(
    ensemble_inc_arr, ensemble_thresholds_vec, ensemble_inc_vecs,
    outbreakthreshold, minoutbreakdur,
    minoutbreaksize,
)
    @inbounds for sim in axes(ensemble_inc_vecs, 2)
        convert_svec_to_matrix!(
            @view(ensemble_inc_arr[:, 1, sim]),
            @view(ensemble_inc_vecs[:, sim])
        )

        ensemble_inc_arr[:, 2, sim] .=
            @view(ensemble_inc_arr[:, 1, sim]) .>= outbreakthreshold

        abovethresholdrle = StatsBase.rle(@view(ensemble_inc_arr[:, 2, sim]))

        outbreak_thresholds = calculate_outbreak_thresholds(
            abovethresholdrle; ncols = 5
        )

        classify_all_outbreaks!(
            @view(ensemble_inc_arr[:, 3, sim]),
            outbreak_thresholds,
            @view(ensemble_inc_arr[:, 1, sim]),
            minoutbreakdur,
            minoutbreaksize,
        )

        ensemble_thresholds_vec[sim] = filter_only_outbreaks(
            outbreak_thresholds
        )
    end
    return nothing
end

function calculate_outbreak_thresholds(outbreakrle; ncols = 5)
    @assert ncols >= 3

    # Calculate upper and lower indices of consecutive days of infection
    outbreakaccum = accumulate(+, outbreakrle[2])
    upperbound_indices = findall(isequal(1), outbreakrle[1])

    outbreak_thresholds = zeros(Int64, length(upperbound_indices), ncols)

    @inbounds outbreak_thresholds[:, 2] .= @view(
        outbreakaccum[upperbound_indices]
    )
    map!(
        x -> x - 1 == 0 ? 1 : outbreakaccum[x - 1] + 1,
        @view(outbreak_thresholds[:, 1]),
        upperbound_indices,
    )

    return outbreak_thresholds
end

function calculate_outbreak_duration!(outbreak_thresholds)
    outbreak_thresholds[:, 3] .=
        calculate_outbreak_duration.(
            outbreak_thresholds[:, 1], outbreak_thresholds[:, 2]
        )
    return nothing
end

function calculate_outbreak_duration(outbreak_thresholds_row)
    return calculate_outbreak_duration(
        outbreak_thresholds_row[1], outbreak_thresholds_row[2]
    )
end

function calculate_outbreak_duration(lower_time, upper_time)
    return upper_time - lower_time + 1
end

function classify_all_outbreaks!(
    outbreakstatus_vec,
    all_thresholds_arr,
    incidence_vec,
    minoutbreakdur,
    minoutbreaksize,
)
    for (row, (lower, upper, outbreakdur)) in
        pairs(eachrow(all_thresholds_arr[:, 1:3]))
        all_thresholds_arr[row, 3] = calculate_outbreak_duration(lower, upper)

        all_thresholds_arr[row, 4] = calculate_outbreak_size(
            incidence_vec, lower, upper
        )

        all_thresholds_arr[row, 5] = classify_outbreak(
            outbreakdur,
            minoutbreakdur,
            all_thresholds_arr[row, 4],
            minoutbreaksize,
        )

        @view(outbreakstatus_vec[lower:upper]) .= @view(
            all_thresholds_arr[row, 5]
        )
    end

    return nothing
end

function calculate_outbreak_size(incidence_vec, lower_time, upper_time)
    return sum(@view(incidence_vec[lower_time:upper_time]))
end

function filter_only_outbreaks(all_thresholds_arr)
    return @view(
        all_thresholds_arr[(all_thresholds_arr[:, 5] .== 1), 1:4]
    )
end

function calculate_period_sum(incvec)
    return sum(incvec)
end

function classify_outbreak(
    lower_time,
    upper_time,
    minoutbreakdur,
    periodsumvec,
    minoutbreaksize,
)
    return classify_outbreak(
        calculate_outbreak_duration(lower_time, upper_time),
        minoutbreakdur,
        periodsumvec,
        minoutbreaksize,
    )
end

function classify_outbreak(
    outbreakdur, minoutbreakdur, periodsumvec, minoutbreaksize
)
    if outbreakdur >= minoutbreakdur && periodsumvec >= minoutbreaksize
        return 1
    end
    return 0
end
