export line_plot

lineplot_colors = [
    "#56B4E9"
    "#E69F00"
    repeat(["#483248"], 2)...
]

# Main plotting functions

function line_plot(
        results::StructVector{OutbreakDetectionCore.OptimizationResult};
        outcome = :accuracies,
        alert_method = OutbreakDetectionCore.AlertMethod(OutbreakDetectionCore.MovingAverage(7)),
        accuracy_metric = OutbreakDetectionCore.AccuracyMetric(OutbreakDetectionCore.BalancedAccuracy()),
        threshold_bounds = (; lower = 0.0, upper = 20.0),
        plotdirpath = DrWatson.plotsdir(),
        plotname = "line_accuracy_plot",
        plotformat = "png",
        size = (1300, 800),
        colors = lineplot_colors,
        alpha = 0.5,
        markersize = 15,
        xlabel = "Proportion Of Infected Individuals Tested",
        ylabel = "Outbreak Detection\nAccuracy",
        facet_fontsize = 24,
        legendsize = 28,
        xlabelsize = 28,
        ylabelsize = 28,
        show_x_facet_label = true,
        show_y_facet_label = true,
        ylims = (nothing, nothing),
        hidedecorations = (true, true),
        hlines = nothing,
        nbanks = 1,
        legend_rowsize = Makie.Relative(0.05),
        xlabel_rowsize = Makie.Relative(0.03),
        force = true,
        save_plot = true,
        dots = false,
        percentiles = [0.1, 0.9],
        kwargs...,
    )
    mkpath(plotdirpath)
    plotpath = joinpath(plotdirpath, "$plotname.$plotformat")

    local_colors = colors

    if !isfile(plotpath) || force
        filtered_results = filter(
            r -> r.alert_method == alert_method && r.accuracy_metric == accuracy_metric && r.threshold_bounds == threshold_bounds,
            results
        )
        # Reshape results into matrix structure
        result_matrix, unique_noise_types = reshape_optimization_results_to_matrix(filtered_results)

        fig = Figure()

        num_noise_types = length(unique_noise_types)

        # # Get unique test specifications from first cell in sorted order
        unique_test_specifications = get_unique_test_specifications_in_sorted_order(
            result_matrix[1, 1].test_specification
        )

        for i in axes(result_matrix, 1)
            noise_type = unique_noise_types[i]
            label_noise_description = get_noise_label(noise_type)

            for j in axes(result_matrix, 2)
                cell_results = result_matrix[i, j]
                noise_level = cell_results[1].noise_level

                _line_plot(
                    fig,
                    noise_level,
                    noise_type,
                    unique_test_specifications,
                    cell_results,
                    i,
                    j;
                    outcome = outcome,
                    num_noise_descriptions = num_noise_types,
                    colors = colors,
                    alpha = alpha,
                    markersize = markersize,
                    ylabel = "$label_noise_description\n" * ylabel,
                    ylabelsize = ylabelsize,
                    ylims = ylims,
                    hidedecorations = hidedecorations,
                    facet_fontsize = facet_fontsize,
                    show_x_facet_label = show_x_facet_label,
                    hlines = hlines,
                    dots = dots,
                    percentiles = percentiles,
                    kwargs...,
                )
            end
        end

        if show_y_facet_label
            map(enumerate(unique_noise_types)) do (i, noise_type)
                Box(fig[i, 0]; color = :lightgray, strokevisible = false)
                Label(
                    fig[i, 0],
                    get_noise_label(noise_type);
                    fontsize = 16,
                    rotation = pi / 2,
                    padding = (0, 0, 0, 0),
                    valign = :center,
                    tellheight = false,
                )
            end
            colsize!(fig.layout, 0, Relative(0.03))
        end

        Legend(
            fig[0, :],
            map(
                enumerate(reverse(unique_test_specifications))
            ) do (i, test_spec)
                # Calculate the original index (before reversal) to get correct color
                original_idx = length(unique_test_specifications) - i + 1

                if outcome == :alert_threshold
                    markerstyle =
                        test_spec.test_result_lag == 0 ? :circle : :utriangle
                    return [
                        MarkerElement(;
                            color = colors[original_idx],
                            marker = markerstyle,
                            markersize = markersize * 2,
                        ),
                    ]
                else
                    linestyle = test_spec.test_result_lag == 0 ? :solid : :dot
                    return [
                        PolyElement(; color = (colors[original_idx], 0.3)),
                        LineElement(;
                            color = colors[original_idx], linestyle = linestyle
                        ),
                    ]
                end
            end,
            OutbreakDetectionCore.plot_test_description.(reverse(unique_test_specifications));
            labelsize = legendsize,
            orientation = :horizontal,
            nbanks = nbanks,
        )
        xlabel_position = num_noise_types + 1
        Label(
            fig[xlabel_position, :],
            xlabel;
            fontsize = xlabelsize,
            font = :bold,
        )
        rowsize!(fig.layout, 0, legend_rowsize)
        rowsize!(fig.layout, xlabel_position, xlabel_rowsize)

        if save_plot
            Makie.save(plotpath, fig; size = size)
        end
        return fig
    end

    return nothing
end

function _line_plot(
        fig,
        noise_level,
        noise_type,
        unique_test_specifications,
        cell_results::StructVector{OutbreakDetectionCore.OptimizationResult},
        i,
        j;
        outcome = :accuracies,
        colors = lineplot_colors,
        alpha = 0.5,
        markersize = 15,
        ylabel = "Accuracy",
        ylabelsize = 28,
        show_x_facet_label = true,
        facet_fontsize = 24,
        num_noise_descriptions = 1,
        ylims = (nothing, nothing),
        hidedecorations = (true, true),
        hlines = nothing,
        dots = false,
        percentiles = [0.1, 0.9],
        kwargs...,
    )
    kwargs_dict = Dict{Symbol, Any}(kwargs)

    # Organize data by test specification
    data_by_test = map(unique_test_specifications) do test_spec
        # Filter results for this test specification
        matching_indices = findall(cell_results.test_specification) do ts
            ts.sensitivity == test_spec.sensitivity &&
                ts.specificity == test_spec.specificity &&
                ts.test_result_lag == test_spec.test_result_lag
        end

        matching_results = cell_results[matching_indices]

        # Extract x values (percent tested)
        x_values = matching_results.percent_tested

        # Extract and compute statistics for outcome
        if outcome == :alert_threshold
            # For alert threshold, just use the optimal threshold value
            y_values = matching_results.optimal_threshold
            return (
                test_spec = test_spec,
                x = x_values,
                y = y_values,
                is_scalar = true,
            )
        else
            # For other outcomes, compute statistics from vectors
            outcome_data = extract_outcome_values.(matching_results, outcome)

            # Determine if this is nested data (Vector{Vector{T}})
            is_nested = !isempty(outcome_data) &&
                outcome_data[1] isa AbstractVector{<:AbstractVector}

            # Compute statistics for each result
            stats = if is_nested
                compute_nested_summary_statistics.(outcome_data; percentiles = percentiles)
            else
                compute_summary_statistics.(outcome_data; percentiles = percentiles)
            end

            y_mean = [s.mean for s in stats]
            y_lower = [s.percentiles[1] for s in stats]
            y_upper = [s.percentiles[2] for s in stats]

            return (
                test_spec = test_spec,
                x = x_values,
                y_mean = y_mean,
                y_lower = y_lower,
                y_upper = y_upper,
                is_scalar = false,
            )
        end
    end

    if show_x_facet_label && i == 1
        # Create facet label with noise level
        x_facet_label = L"\Lambda(%$(Int64(round(noise_level; digits = 0))))"
        kwargs_dict[:x_facet_label] = x_facet_label
    end

    gl = fig[i, j] = GridLayout()

    if i != num_noise_descriptions
        xlabel = ""
    end

    if j != 1
        ylabel = ""
    end

    _line_plot_facet(
        gl,
        data_by_test;
        outcome = outcome,
        colors = colors,
        alpha = alpha,
        markersize = markersize,
        ylabel = ylabel,
        ylabelsize = ylabelsize,
        facet_fontsize = facet_fontsize,
        ylims = ylims,
        hlines = hlines,
        kwargs_dict...,
    )

    if ylims != (nothing, nothing) && hidedecorations[2] && j != 1
        hideydecorations!(contents(gl)[1])
    end

    if hidedecorations[1] && i < num_noise_descriptions
        hidexdecorations!(contents(gl)[1])
    end

    return nothing
end

function _line_plot_facet(
        gl,
        data_by_test;
        outcome = :accuracy,
        colors = lineplot_colors,
        alpha = 0.5,
        markersize = 15,
        xlabel = "Proportion Tested",
        ylabel = "Accuracy",
        xlabelsize = 28,
        ylabelsize = 28,
        facet_fontsize = 24,
        ylims = (nothing, nothing),
        hlines = nothing,
        kwargs...,
    )
    kwargs_dict = Dict(kwargs)

    ypos = haskey(kwargs_dict, :x_facet_label) ? 2 : 1
    xpos = 1

    ax = Axis(
        gl[ypos, xpos];
        ylabel = ylabel,
        ylabelsize = ylabelsize,
    )

    if haskey(kwargs_dict, :x_facet_label)
        Box(gl[1, xpos]; color = :lightgray, strokevisible = false)
        Label(
            gl[1, xpos],
            kwargs_dict[:x_facet_label];
            fontsize = facet_fontsize,
            padding = (0, 0, 0, 0),
            valign = :bottom,
            tellwidth = false,
        )
        rowsize!(gl, 2, Relative(0.9))
    end

    for (i, test_data) in pairs(data_by_test)
        test_spec = test_data.test_spec

        if outcome == :alert_threshold
            # Plot scatter with lines for alert threshold
            if test_spec.test_result_lag == 0
                markerstyle = :circle
                msize = markersize
            else
                markerstyle = :utriangle
                msize = markersize * 1.5
            end

            scatterlines!(
                ax,
                test_data.x,
                test_data.y;
                color = (colors[i], alpha),
                strokecolor = colors[i],
                marker = markerstyle,
                markersize = msize,
            )
        else
            # Plot lines with confidence bands for other outcomes
            linestyle = test_spec.test_result_lag == 0 ? :solid : :dash

            lines!(
                ax,
                test_data.x,
                test_data.y_mean;
                color = colors[i],
                linestyle = linestyle,
            )

            band!(
                ax,
                test_data.x,
                test_data.y_lower,
                test_data.y_upper;
                color = (colors[i], alpha),
            )
        end
    end

    if !isnothing(hlines)
        hlines!(ax, hlines; color = :black, linewidth = 1)
    end

    ylims!(ax, ylims)

    return nothing
end
