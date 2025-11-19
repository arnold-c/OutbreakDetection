using Match: Match

export create_testing_arrs, create_testing_arrs!

"""
    create_testing_arrs(incarr, noisearr, outbreak_detect_spec, 
                        individual_test_spec)

Create testing arrays for outbreak detection simulation.

Returns testarr and test_movingavg_arr.
"""
function create_testing_arrs(
        incarr,
        noisearr,
        outbreak_detect_spec::OutbreakDetectionSpecification,
        individual_test_spec::IndividualTestSpecification,
    )
    testarr = zeros(Int64, size(incarr, 1), 7, size(incarr, 3))
    test_movingavg_arr = zeros(Int64, size(incarr, 1), size(incarr, 3))

    create_testing_arrs!(
        testarr,
        test_movingavg_arr,
        incarr,
        noisearr,
        outbreak_detect_spec.alert_method.method_name,
        outbreak_detect_spec.alert_threshold,
        outbreak_detect_spec.moving_average_lag,
        outbreak_detect_spec.percent_tested,
        individual_test_spec.test_result_lag,
        individual_test_spec.sensitivity,
        individual_test_spec.specificity,
    )

    return testarr, test_movingavg_arr
end

"""
    create_testing_arrs!(testarr, test_movingavg_arr, incarr, noisearr, 
                         alert_method, alertthreshold, moveavglag, 
                         perc_tested, testlag, testsens, testspec)

In-place creation of testing arrays.
"""
function create_testing_arrs!(
        testarr,
        test_movingavg_arr,
        incarr,
        noisearr,
        alert_method,
        alertthreshold,
        moveavglag,
        perc_tested,
        testlag,
        testsens,
        testspec,
    )
    tlength = size(testarr, 1)

    for sim in axes(incarr, 3)
        # Number of infectious individuals tested
        calculate_tested!(
            @view(testarr[:, 1, sim]), @view(incarr[:, 1, sim]), perc_tested
        )

        # Number of noise individuals tested
        calculate_tested!(
            @view(testarr[:, 2, sim]), @view(noisearr[:, sim]), perc_tested
        )

        # Number of test positive INFECTED individuals
        calculate_true_positives!(
            @view(testarr[:, 3, sim]),
            @view(testarr[:, 1, sim]),
            tlength,
            testlag,
            testsens,
        )

        # Number of test positive NOISE individuals
        calculate_noise_positives!(
            @view(testarr[:, 4, sim]),
            @view(testarr[:, 2, sim]),
            tlength,
            testlag,
            testspec,
        )

        # Number of test positive TOTAL individuals
        @. testarr[:, 5, sim] =
            @view(testarr[:, 3, sim]) + @view(testarr[:, 4, sim])

        # Calculate moving average of TOTAL test positives
        calculate_movingavg!(
            @view(test_movingavg_arr[:, sim]),
            @view(testarr[:, 5, sim]),
            moveavglag,
        )

        detectoutbreak_args = Match.@match alert_method begin
            "movingavg" => (
                @view(testarr[:, 6, sim]),
                @view(test_movingavg_arr[:, sim]),
                alertthreshold,
            )
            "dailythreshold_movingavg" => (
                @view(testarr[:, 6, sim]),
                @view(testarr[:, 5, sim]),
                @view(test_movingavg_arr[:, sim]),
                alertthreshold,
            )
        end

        # TOTAL Test positive individuals trigger outbreak response
        detectoutbreak!(detectoutbreak_args...)

        # Triggered outbreak equal to actual outbreak status
        @. testarr[:, 7, sim] =
            @view(testarr[:, 6, sim]) == @view(incarr[:, 3, sim])
    end

    return nothing
end
