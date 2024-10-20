using DataFrames
using DrWatson: DrWatson
using StatsBase: StatsBase
using Match: Match
using StructArrays

lineplot_colors = [
    "#56B4E9"
    "#E69F00"
    repeat(["#000000"], 2)...
]

function line_plot(
    noise_spec_vec,
    ensemble_percent_clinic_tested_vec,
    optimal_threshold_test_spec_vec,
    optimal_threshold_core_params;
    outcome = :accuracy,
    plotdirpath = DrWatson.plotsdir(),
    plotname = "line_accuracy_plot",
    plotformat = "png",
    size = (2200, 1200),
    colors = lineplot_colors,
    xlabel = "Proportion Tested",
    ylabel = "Outbreak Detection\nAccuracy",
    facet_fontsize = 24,
    labelsize = 30,
    show_x_facet_label = true,
    show_y_facet_label = true,
    ylims = (nothing, nothing),
    hidedecorations = (true, true),
    clinical_hline = true,
    hlines = nothing,
    force = false,
    save_plot = true,
    kwargs...,
)
    mkpath(plotdirpath)
    plotpath = joinpath(plotdirpath, "$plotname.$plotformat")

    if !isfile(plotpath) || force
        optimal_threshold_characteristics = collect_OptimalThresholdCharacteristics(
            noise_spec_vec,
            ensemble_percent_clinic_tested_vec,
            optimal_threshold_test_spec_vec,
            optimal_threshold_core_params;
            clinical_hline = clinical_hline
        )

        return line_plot(
            optimal_threshold_characteristics;
            outcome = outcome,
            plotdirpath = plotdirpath,
            plotname = plotname,
            plotformat = plotformat,
            size = size,
            colors = colors,
            xlabel = xlabel,
            ylabel = ylabel,
            facet_fontsize = facet_fontsize,
            labelsize = labelsize,
            show_x_facet_label = show_x_facet_label,
            show_y_facet_label = show_y_facet_label,
            ylims = ylims,
            hidedecorations = hidedecorations,
            clinical_hline = clinical_hline,
            hlines = hlines,
            force = force,
            save_plot = save_plot,
            kwargs...
        )
    end

    return nothing
end

function line_plot(
    optimal_thresholds_chars_array;
    outcome = :accuracy,
    plotdirpath = DrWatson.plotsdir(),
    plotname = "line_accuracy_plot",
    plotformat = "png",
    size = (2200, 1200),
    colors = lineplot_colors,
    xlabel = "Proportion Tested",
    ylabel = "Outbreak Detection\nAccuracy",
    facet_fontsize = 24,
    labelsize = 30,
    show_x_facet_label = true,
    show_y_facet_label = true,
    ylims = (nothing, nothing),
    hidedecorations = (true, true),
    clinical_hline = true,
    hlines = nothing,
    force = false,
    save_plot = true,
    kwargs...,
)
    mkpath(plotdirpath)
    plotpath = joinpath(plotdirpath, "$plotname.$plotformat")

    if !isfile(plotpath) || force
        fig = Figure()

        unique_noise_descriptions = unique(
            map(
                char -> char[1].noise_specification,
                optimal_thresholds_chars_array,
            ),
        )
        num_noise_descriptions = length(unique_noise_descriptions)
        unique_test_specifications = unique(
            optimal_thresholds_chars_array[1, 1].individual_test_specification
        )

        for i in axes(optimal_thresholds_chars_array, 1)
            noise_description = get_noise_description(
                optimal_thresholds_chars_array[i, 1].noise_specification[1]
            )
            label_noise_description = Match.@match noise_description begin
                "poisson" => "Poisson Noise"
                "dynamical, in-phase" => "Dynamical Noise: In-Phase"
                _ => "Other Noise"
            end

            for j in axes(optimal_thresholds_chars_array, 2)
                noise_spec = optimal_thresholds_chars_array[i, j].noise_specification[1]

                _line_plot(
                    fig,
                    noise_spec,
                    unique_test_specifications,
                    optimal_thresholds_chars_array[i, j],
                    i,
                    j;
                    outcome = outcome,
                    num_noise_descriptions = num_noise_descriptions,
                    colors = colors,
                    xlabel = xlabel,
                    ylabel = "$label_noise_description\n" * ylabel,
                    ylims = ylims,
                    hidedecorations = hidedecorations,
                    facet_fontsize = facet_fontsize,
                    show_x_facet_label = show_x_facet_label,
                    hlines = hlines,
                    kwargs...,
                )
            end
        end

        if show_y_facet_label
            map(enumerate(unique_noise_descriptions)) do (i, noise_description)
                Box(fig[i, 0]; color = :lightgray, strokevisible = false)
                Label(
                    fig[i, 0],
                    titlecase(noise_description);
                    fontsize = 16,
                    rotation = pi / 2,
                    padding = (0, 0, 0, 0),
                    valign = :center,
                    tellheight = false,
                )
            end
            colsize!(fig.layout, 0, Relative(0.03))
        end
        if clinical_hline
            push!(colors, "green")
        end
        rg = r"\((.*)(\% .*\))"
        Legend(
            fig[0, :],
            map(enumerate(unique_test_specifications)) do (i, test_spec)
                linestyle = test_spec.test_result_lag == 0 ? :solid : :dot
                return [
                    PolyElement(; color = (colors[i], 0.3)),
                    LineElement(; color = colors[i], linestyle = linestyle),
                ]
            end,
            map(
                test_description ->
                    replace.(
                        test_description,
                        rg =>
                            s -> parse(Float64, match(rg, s).captures[1]) / 100,
                    ),
                get_test_description.(unique_test_specifications),
            );
            labelsize = labelsize,
            orientation = :horizontal,
        )
        rowsize!(fig.layout, 0, Relative(0.03))

        if save_plot
            Makie.save(plotpath, fig; size = size)
        end
        return fig
    end

    return nothing
end

function collect_OptimalThresholdCharacteristics(
    noise_spec_vec,
    ensemble_percent_clinic_tested_vec,
    optimal_threshold_test_spec_vec,
    optimal_threshold_core_params;
    clinical_hline = false,
)
    noise_descriptions = get_noise_description.(noise_spec_vec)
    unique_noise_descriptions = unique(noise_descriptions)

    shape_noise_specifications =
        map(unique_noise_descriptions) do noise_description
            filter(
                noise_spec ->
                    noise_description == get_noise_description(noise_spec),
                noise_spec_vec,
            )
        end

    @assert length(vcat(shape_noise_specifications...)) ==
        length(unique_noise_descriptions) *
            length(shape_noise_specifications[1])

    if clinical_hline
        optimal_threshold_test_spec_vec = vcat(
            optimal_threshold_test_spec_vec, CLINICAL_CASE_TEST_SPEC
        )
    end

    optimal_thresholds_vecs = Array{
        StructArray{OptimalThresholdCharacteristics}
    }(
        undef,
        length(unique_noise_descriptions),
        length(shape_noise_specifications[1]),
    )

    for (i, noise_description) in pairs(unique_noise_descriptions)
        label_noise_description = Match.@match noise_description begin
            "poisson" => "Poisson Noise"
            "dynamical, in-phase" => "Dynamical Noise: In-Phase"
            _ => "Other Noise"
        end

        for (j, noise_spec) in pairs(shape_noise_specifications[i])
            optimal_threshold_comparison_params = (
                noise_specification = noise_spec,
                optimal_threshold_core_params...,
            )

            optimal_thresholds_vecs[i, j] = calculate_OptimalThresholdCharacteristics(
                ensemble_percent_clinic_tested_vec,
                optimal_threshold_test_spec_vec,
                optimal_threshold_comparison_params,
            )
        end
    end
    return optimal_thresholds_vecs
end

function _line_plot(
    fig,
    noise_spec,
    unique_test_specifications,
    optimal_thresholds_vec,
    i,
    j;
    outcome = :accuracy,
    colors = lineplot_colors,
    xlabel = "Proportion Tested",
    ylabel = "Accuracy",
    show_x_facet_label = true,
    facet_fontsize = 24,
    num_noise_descriptions = 1,
    ylims = (nothing, nothing),
    hidedecorations = (true, true),
    hlines = nothing,
    kwargs...,
)
    kwargs_dict = Dict{Symbol,Any}(kwargs)

    long_df = create_optimal_threshold_summary_df(
        optimal_thresholds_vec,
        outcome;
        percentiles = [0.1, 0.9],
    )

    select!(
        long_df,
        Cols(
            :percent_clinic_tested,
            :sensitivity,
            :specificity,
            :test_lag,
            string(outcome) * "_mean",
            x -> endswith(x, "th"),
        ),
    )

    if show_x_facet_label && i == 1
        x_facet_label = "$(Int64(round(
            StatsBase.mean(
                optimal_thresholds_vec[1].outbreak_threshold_chars.mean_noise_incidence_ratio
            );
            digits = 0,
        ))):1 Noise:Signal Ratio"

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
        noise_spec,
        unique_test_specifications,
        long_df;
        outcome = outcome,
        colors = colors,
        xlabel = xlabel,
        ylabel = ylabel,
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
    noise_spec,
    unique_test_specifications,
    long_df;
    outcome = :accuracy,
    colors = lineplot_colors,
    xlabel = "Proportion Tested",
    ylabel = "Accuracy",
    facet_fontsize = 24,
    ylims = (nothing, nothing),
    hlines = nothing,
    kwargs...,
)
    kwargs_dict = Dict(kwargs)

    ypos = haskey(kwargs_dict, :x_facet_label) ? 2 : 1
    xpos = 1

    outcome_mean = string(outcome) * "_mean"
    outcome_10th = string(outcome) * "_10th"
    outcome_90th = string(outcome) * "_90th"

    ax = Axis(
        gl[ypos, xpos];
        xlabel = xlabel,
        ylabel = ylabel,
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

    for (i, test) in pairs(unique_test_specifications)
        subsetted_df = DataFrames.subset(
            long_df,
            :sensitivity =>
                x -> x .== test.sensitivity,
            :specificity =>
                x -> x .== test.specificity,
            :test_lag => x -> x .== test.test_result_lag,
        )

        if test == CLINICAL_CASE_TEST_SPEC
            hlines!(ax, subsetted_df[!, outcome_mean][1]; color = colors[end])
            continue
        end

        linestyle = test.test_result_lag == 0 ? :solid : :dash

        string(
            subsetted_df.sensitivity[1], ", ", subsetted_df.test_lag[1]
        )
        lines!(
            ax,
            subsetted_df.percent_clinic_tested,
            subsetted_df[!, outcome_mean];
            color = colors[i],
            linestyle = linestyle,
        )

        band!(
            ax,
            subsetted_df.percent_clinic_tested,
            subsetted_df[!, outcome_10th],
            subsetted_df[!, outcome_90th];
            color = colors[i],
            alpha = 0.3,
        )

        ylims!(ax, ylims)
    end

    if !isnothing(hlines)
        hlines!(ax, hlines; color = :black, linewidth = 1)
    end

    return nothing
end
