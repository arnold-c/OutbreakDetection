using GLMakie

function compare_optimal_thresholds_chars_plot(
    optimal_thresholds_vec,
    plottingchars;
    plotdirpath = plotsdir("ensemble/optimal-thresholds"),
    plotformat = "png",
    force = false,
    kwargs...,
)
    unique_percent_clinic_tested = unique(
        optimal_thresholds_vec.percent_clinic_tested
    )

    mkpath(plotdirpath)

    for percent_clinic_tested in unique_percent_clinic_tested
        plotname = "compare-outbreak_clinic-tested-$(percent_clinic_tested)_best-thresholds"
        plotpath = joinpath(plotdirpath, "$plotname.$plotformat")

        if !isfile(plotpath) || force
            optimal_thresholds_chars = filter(
                optimal_thresholds ->
                    optimal_thresholds.percent_clinic_tested ==
                    percent_clinic_tested ||
                        (
                            optimal_thresholds.percent_clinic_tested ==
                            1.0 &&
                            optimal_thresholds.individual_test_specification ==
                            CLINICAL_CASE_TEST_SPEC
                        ),
                optimal_thresholds_vec,
            )

            plot = create_optimal_thresholds_chars_plot(
                optimal_thresholds_chars,
                plottingchars;
                kwargs...,
            )

            save(
                joinpath(
                    plotdirpath,
                    "compare-outbreak_clinic-tested-$(percent_clinic_tested)_best-thresholds.png",
                ),
                plot;
                size = (2200, 1600),
            )

            Makie.empty!(plot)

            @info "Created optimal thresholds plot for % clinic tested $(percent_clinic_tested)"

            continue
        end

        GC.gc(true)

        @info "Optimal thresholds plot for % clinic tested $(percent_clinic_tested) already exists. Skipping."
    end

    return nothing
end

function create_optimal_thresholds_chars_plot(
    optimal_thresholds_chars,
    plottingchars;
    kwargs...,
)
    number_tests = length(optimal_thresholds_chars)
    number_plotting_chars = length(plottingchars)
    midpoint_plotting_chars = Int64(
        round(number_plotting_chars / 2; digits = 0)
    )

    sort!(
        optimal_thresholds_chars;
        by = threshold ->
            threshold.individual_test_specification.specificity,
    )

    fig = Figure()

    Label(
        fig[
            1,
            (midpoint_plotting_chars + 1):(midpoint_plotting_chars + 2),
            Top(),
        ],
        "Perc Clinic Tested: $(optimal_thresholds_chars[end].percent_clinic_tested), Noise: $(getdirpath(optimal_thresholds_chars[end].noise_specification))",
    )

    thresholdschars_structarr =
        optimal_thresholds_chars.outbreak_threshold_chars

    for (column, chartuple) in pairs(plottingchars)
        x = column + 1

        thresholdschars_vec =
            getproperty.(thresholdschars_structarr, chartuple.char)

        charvecs = reduce(vcat, thresholdschars_vec)

        if !haskey(chartuple, :bins)
            if !haskey(chartuple, :binwidth)
                @error "The metric $(chartuple.char) wasn't provided with bins or a binwidth"
                break
            end
            bins = calculate_bins(charvecs, chartuple.binwidth)
        else
            bins .= chartuple.bins
        end

        if !haskey(chartuple, :label)
            label = :none
        else
            label = chartuple.label
        end

        for (row, optimal_thresholds) in pairs(optimal_thresholds_chars)
            y = row + 1

            gl = fig[y, x] = GridLayout()
            ax = Axis(gl[1, 1]; xlabel = label)

            hist!(
                ax,
                reduce(vcat, thresholdschars_vec[row]);
                bins = bins,
                color = chartuple.color,
            )

            Label(
                fig[y, 1, Left()],
                "Sens: $(optimal_thresholds.individual_test_specification.sensitivity), Spec: $(optimal_thresholds.individual_test_specification.specificity),\nLag: $(optimal_thresholds.individual_test_specification.test_result_lag), Threshold: $(optimal_thresholds.alert_threshold)";
                rotation = pi / 2,
            )

            if row < number_tests
                hidexdecorations!(ax; ticklabels = false, ticks = false)
            end
        end
    end

    rowsize!(fig.layout, 1, 5)
    colsize!(fig.layout, 1, 7)

    return fig
end

function compare_optimal_thresholds_test_chars_plot(
    optimal_thresholds_vec,
    plottingchars;
    plotdirpath = plotsdir("ensemble/optimal-thresholds"),
    plotformat = "png",
    force = false,
    kwargs...,
)
    unique_tests = unique(
        optimal_thresholds_vec.individual_test_specification
    )

    filter!(x -> !in(x, CLINICAL_TEST_SPECS), unique_tests)

    mkpath(plotdirpath)

    for test_specification in unique_tests
        plotname = "compare-outbreak_clinic-test-specification_sens-$(test_specification.sensitivity)_spec-$(test_specification.specificity)_lag-$(test_specification.test_result_lag)_best-thresholds"
        plotpath = joinpath(plotdirpath, "$plotname.$plotformat")

        if !isfile(plotpath) || force
            optimal_thresholds_chars = filter(
                optimal_thresholds ->
                    optimal_thresholds.individual_test_specification ==
                    test_specification,
                optimal_thresholds_vec,
            )

            plot = create_optimal_thresholds_test_chars_plot(
                optimal_thresholds_chars,
                plottingchars;
                kwargs...,
            )

            save(
                joinpath(
                    plotdirpath,
                    "compare-outbreak_clinic-test-specification_sens-$(test_specification.sensitivity)_spec-$(test_specification.specificity)_lag-$(test_specification.test_result_lag)_best-thresholds.png",
                ),
                plot;
                size = (2200, 1600),
            )

            Makie.empty!(plot)

            @info "Created optimal thresholds plot for test specification $(test_specification.sensitivity)_$(test_specification.specificity)_$(test_specification.test_result_lag)"
            continue
        end

        GC.gc(true)

        @info "Optimal thresholds plot for test specification $(test_specification.sensitivity)_$(test_specification.specificity)_$(test_specification.test_result_lag) already exists. Skipping."
    end

    return nothing
end

function create_optimal_thresholds_test_chars_plot(
    optimal_thresholds_chars,
    plottingchars;
    kwargs...,
)
    number_clinic_testing_rates = length(optimal_thresholds_chars)
    number_plotting_chars = length(plottingchars)
    midpoint_plotting_chars = Int64(
        round(number_plotting_chars / 2; digits = 0)
    )

    sort!(
        optimal_thresholds_chars;
        by = threshold ->
            threshold.percent_clinic_tested,
    )

    fig = Figure()

    Label(
        fig[
            1,
            (midpoint_plotting_chars + 1):(midpoint_plotting_chars + 2),
            Top(),
        ],
        "Test characteristics: Sensitivity $(optimal_thresholds_chars[end].individual_test_specification.sensitivity), Specificity $(optimal_thresholds_chars[end].individual_test_specification.specificity), Lag $(optimal_thresholds_chars[end].individual_test_specification.test_result_lag), Noise: $(getdirpath(optimal_thresholds_chars[end].noise_specification))",
    )

    thresholdschars_structarr =
        optimal_thresholds_chars.outbreak_threshold_chars

    for (column, chartuple) in pairs(plottingchars)
        x = column + 1

        thresholdschars_vec =
            getproperty.(thresholdschars_structarr, chartuple.char)

        charvecs = reduce(vcat, thresholdschars_vec)

        if !haskey(chartuple, :bins)
            if !haskey(chartuple, :binwidth)
                @error "The metric $(chartuple.char) wasn't provided with bins or a binwidth"
                break
            end
            bins = calculate_bins(charvecs, chartuple.binwidth)
        else
            bins .= chartuple.bins
        end

        if !haskey(chartuple, :label)
            label = :none
        else
            label = chartuple.label
        end

        for (row, optimal_thresholds) in pairs(optimal_thresholds_chars)
            y = row + 1

            gl = fig[y, x] = GridLayout()
            ax = Axis(gl[1, 1]; xlabel = label)

            hist!(
                ax,
                reduce(vcat, thresholdschars_vec[row]);
                bins = bins,
                color = chartuple.color,
            )

            Label(
                fig[y, 1, Left()],
                "% Clinic Tested: $(optimal_thresholds.percent_clinic_tested), Threshold: $(optimal_thresholds.alert_threshold)";
                rotation = pi / 2,
            )

            if row < number_clinic_testing_rates
                hidexdecorations!(ax; ticklabels = false, ticks = false)
            end
        end
    end

    rowsize!(fig.layout, 1, 5)
    colsize!(fig.layout, 1, 7)

    return fig
end
