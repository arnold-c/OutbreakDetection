#%%
using DrWatson
@quickactivate "OutbreakDetection"

using StatsBase

using OutbreakDetection: N_MISSED_OUTBREAKS_COLOR,
    PERC_OUTBREAKS_DETECTED_COLOR, N_MISSED_OUTBREAKS_COLOR, N_ALERTS_COLOR

using OutbreakDetectionUtils

include(projectdir("manuscript", "scripts", "plotting-setup.jl"))

#%%
states_p = StateParameters(;
    N = 500_000, s_prop = 0.05, e_prop = 0.00, i_prop = 0.00
)

dynamics_p = DynamicsParameters(
    BETA_MEAN,
    BETA_FORCE,
    cos,
    SIGMA,
    GAMMA,
    MU,
    ANNUAL_BIRTHS_PER_K,
    EPSILON,
    R0,
    VACCINATION_COVERAGE,
)

noise_states_p = StateParameters(;
    N = 500_000, s_prop = 0.15, e_prop = 0.00, i_prop = 0.00
)

noise_dynamics_p = DynamicsParameters(
    500_000,
    27,
    0.2,
    1 / 7,
    1 / 14,
    5.0,
    0.65;
    seasonality = sin,
)

test_specification = IndividualTestSpecification(
    0.85, 0.85, 0
)

time_p = SimTimeParameters(;
    tmin = 0.0, tmax = 365.0 * 20, tstep = 1.0
)

outbreak_specification = OutbreakSpecification(5, 30, 500)

movingavg_window = 20

outbreak_detection_specification = OutbreakDetectionSpecification(
    5,
    movingavg_window,
    1.0,
    0.75,
    "movingavg",
)

#%%
function get_outbreak_status(
    inc_vec, outbreak_specification
)
    abovethreshold_vec = vec(
        inc_vec .>= outbreak_specification.outbreak_threshold
    )

    abovethresholdrle = rle(abovethreshold_vec)

    all_outbreak_bounds = calculate_outbreak_thresholds(
        abovethresholdrle
    )

    classify_all_outbreaks!(
        abovethreshold_vec,
        all_outbreak_bounds,
        inc_vec,
        outbreak_specification.minimum_outbreak_duration,
        outbreak_specification.minimum_outbreak_size,
    )

    outbreak_bounds = filter_only_outbreaks(
        all_outbreak_bounds
    )

    outbreak_status = zeros(Int64, length(inc_vec))

    for (lower, upper) in eachrow(outbreak_bounds)
        # @show lower, upper
        outbreak_status[lower:upper] .= 1
    end
    return outbreak_status, outbreak_bounds
end

#%%
function shift_vec(invec, shift::T) where {T<:Integer}
    if shift == 0
        return invec
    end
    outvec = zeros(eltype(invec), length(invec))
    if shift > 0
        for i in (shift + 1):lastindex(invec)
            outvec[i] = invec[i - shift]
        end
    end
    if shift < 0
        for i in 1:(lastindex(invec) + shift)
            outvec[i] = invec[i - shift]
        end
    end
    return outvec
end

#%%
function create_schematic_simulation(
    states_p, dynamics_p, noise_states_p, noise_dynamics_p, test_specification,
    time_p; seed = 1234,
    outbreak_specification = OutbreakSpecification(5, 30, 500),
    outbreak_detection_specification = OutbreakDetectionSpecification(
        5, 7, 1.0, 0.5, "movingavg"
    ),
    noise_scaling = 10,
    shift_noise = 0,
)
    inc_sv = seir_mod(
        states_p.init_states,
        dynamics_p,
        time_p;
        seed = seed,
    )[2]

    inc_vec =
        round.(
            calculate_movingavg(
                vec(convert_svec_to_matrix(inc_sv)),
                movingavg_window,
            )
        )

    outbreak_status, outbreak_bounds = get_outbreak_status(
        inc_vec, outbreak_specification
    )

    noise_sv = seir_mod(
        noise_states_p.init_states,
        noise_dynamics_p,
        time_p;
        seed = seed,
    )[2]

    noise_vec_tmp = calculate_movingavg(
        vec(convert_svec_to_matrix(noise_sv)) .*
        noise_scaling, movingavg_window
    )

    noise_vec = shift_vec(noise_vec_tmp, shift_noise)

    perc_tested =
        outbreak_detection_specification.percent_visit_clinic *
        outbreak_detection_specification.percent_clinic_tested

    true_positives = calculate_positives(
        calculate_true_positives!,
        inc_vec .* perc_tested,
        time_p.tlength,
        test_specification.test_result_lag,
        test_specification.sensitivity,
    )

    false_positives = calculate_positives(
        calculate_noise_positives!,
        noise_vec .* perc_tested,
        time_p.tlength,
        test_specification.test_result_lag,
        test_specification.specificity,
    )

    test_positives = true_positives .+ false_positives

    # Calculate moving average of TOTAL test positives
    movingavg_testpositives = calculate_movingavg(
        test_positives,
        outbreak_detection_specification.moving_average_lag,
    )

    alertstatus_vec = detectoutbreak(
        movingavg_testpositives,
        outbreak_detection_specification.alert_threshold,
    )
    alert_bounds = calculate_outbreak_thresholds(
        rle(alertstatus_vec .> 0); ncols = 3
    )
    OutbreakDetectionUtils.calculate_outbreak_duration!(alert_bounds)

    return inc_vec,
    outbreak_status,
    outbreak_bounds,
    noise_vec,
    movingavg_testpositives,
    alertstatus_vec,
    alert_bounds
end

#%%
function plot_schematic(
    inc_vec,
    outbreakstatus_vec,
    outbreak_bounds;
    time_p = SimTimeParameters(; tmin = 0.0, tmax = 365.0 * 10, tstep = 1.0),
    outbreakcolormap = [
        N_MISSED_OUTBREAKS_COLOR, PERC_OUTBREAKS_DETECTED_COLOR
    ],
    shade_alert_outbreak_overlap = false,
    measlesalpha = 0.5,
    kwargs...,
)
    kwargs_dict = Dict(kwargs)

    times = collect(time_p.trange)

    if haskey(kwargs_dict, :xlims)
        lower = maximum([1, kwargs_dict[:xlims][1] * 365])
        upper = minimum([Int64(time_p.tlength), kwargs_dict[:xlims][2] * 365])
        times = times[lower:upper]
        inc_vec = inc_vec[lower:upper]
        outbreakstatus_vec = outbreakstatus_vec[lower:upper]

        outbreak_bounds = outbreak_bounds[
            (outbreak_bounds[:, 1] .>= lower) .& (outbreak_bounds[:, 2] .<= upper),
            :,
        ]
    end
    outbreak_bounds_vec = vec(outbreak_bounds[:, 1:2])

    fig = Figure()
    incga = fig[1, 1] = GridLayout()
    incax = Axis(incga[1, 1]; ylabel = "Measles Incidence")

    if shade_alert_outbreak_overlap
        if !isempty(outbreak_bounds_vec)
            vspan!(
                incax,
                outbreak_bounds[:, 1],
                outbreak_bounds[:, 2];
                color = (outbreakcolormap[2], measlesalpha),
            )
        end
    end

    lines!(
        incax,
        times,
        inc_vec;
        color = outbreakstatus_vec,
        colormap = outbreakcolormap,
        linewidth = 3,
    )

    hidexdecorations!(incax)

    return fig
end

#%%
function plot_schematic(
    inc_vec,
    outbreakstatus_vec,
    outbreak_bounds,
    alertstatus_vec,
    alert_bounds,
    alertthreshold;
    time_p = SimTimeParameters(; tmin = 0.0, tmax = 365.0 * 10, tstep = 1.0),
    outbreakcolormap = [
        N_MISSED_OUTBREAKS_COLOR, PERC_OUTBREAKS_DETECTED_COLOR
    ],
    alertcolormap = [
        N_MISSED_OUTBREAKS_COLOR, N_ALERTS_COLOR
    ],
    shade_alert_outbreak_overlap = false,
    measlesalpha = 0.5,
    testalpha = 0.5,
    kwargs...)
    kwargs_dict = Dict(kwargs)

    times = collect(time_p.trange)

    if haskey(kwargs_dict, :xlims)
        lower = maximum([1, kwargs_dict[:xlims][1] * 365])
        upper = minimum([Int64(time_p.tlength), kwargs_dict[:xlims][2] * 365])
        times = times[lower:upper]
        inc_vec = inc_vec[lower:upper]
        outbreakstatus_vec = outbreakstatus_vec[lower:upper]
        alertstatus_vec = alertstatus_vec[lower:upper]

        outbreak_bounds = outbreak_bounds[
            (outbreak_bounds[:, 1] .>= lower) .& (outbreak_bounds[:, 2] .<= upper),
            :,
        ]

        alert_bounds = alert_bounds[
            (alert_bounds[:, 1] .>= lower) .& (alert_bounds[:, 2] .<= upper), :,
        ]
    end
    outbreak_bounds_vec = vec(outbreak_bounds[:, 1:2])
    alert_bounds_vec = vec(alert_bounds[:, 1:2])

    fig = Figure()
    incga = fig[1, 1] = GridLayout()
    incax = Axis(incga[1, 1]; ylabel = "Measles Incidence")

    if shade_alert_outbreak_overlap
        if !isempty(alert_bounds_vec)
            vspan!(
                incax,
                alert_bounds[:, 1],
                alert_bounds[:, 2];
                color = (alertcolormap[2], testalpha),
            )
        end
        if !isempty(outbreak_bounds_vec)
            vspan!(
                incax,
                outbreak_bounds[:, 1],
                outbreak_bounds[:, 2];
                color = (outbreakcolormap[2], measlesalpha),
            )
        end
    end

    lines!(
        incax,
        times,
        inc_vec;
        color = alertstatus_vec,
        colormap = alertcolormap,
        linewidth = 3,
    )

    hlines!(
        incax,
        alertthreshold;
        color = :black,
        linewidth = 2,
        linestyle = :dash,
    )

    text!(
        incax,
        times[1] + 1,
        alertthreshold + 0.5;
        text = "T = $alertthreshold",
        justification = :left,
    )

    hidexdecorations!(incax)

    return fig
end

#%%
inc_vec, outbreak_status, outbreak_bounds, noise_vec, movingavg_testpositives, alertstatus_vec, alert_bounds = create_schematic_simulation(
    states_p,
    dynamics_p,
    noise_states_p,
    noise_dynamics_p,
    test_specification,
    time_p;
    seed = 12345,
    outbreak_specification = outbreak_specification,
    outbreak_detection_specification = outbreak_detection_specification,
    noise_scaling = 15,
    shift_noise = -100,
);

#%%
incidence_schematic = plot_schematic(
    inc_vec,
    outbreak_status,
    outbreak_bounds[:, 1:2];
    time_p = time_p,
    shade_alert_outbreak_overlap = true,
    xlims = (5, 13),
    measlesalpha = 0.4,
    outbreakcolormap = [
        N_MISSED_OUTBREAKS_COLOR,
        "#CE6F58",
    ],
)

save(
    plotsdir("schematic-plot_incidence.svg"),
    incidence_schematic;
    size = (1300, 800 / 3),
)

#%%
measles_incidence_threshold = 3

inc_alertstatus_vec = detectoutbreak(
    inc_vec,
    measles_incidence_threshold,
)
inc_alert_bounds = calculate_outbreak_thresholds(
    rle(inc_alertstatus_vec .> 0); ncols = 3
)
OutbreakDetectionUtils.calculate_outbreak_duration!(inc_alert_bounds)

incidence_threshold_schematic = plot_schematic(
    inc_vec,
    outbreak_status,
    outbreak_bounds[:, 1:2],
    inc_alertstatus_vec,
    inc_alert_bounds,
    measles_incidence_threshold;
    time_p = time_p,
    shade_alert_outbreak_overlap = true,
    xlims = (5, 13),
    measlesalpha = 0.4,
    testalpha = 0.4,
    outbreakcolormap = [
        N_MISSED_OUTBREAKS_COLOR,
        "#CE6F58",
    ],
    alertcolormap = [
        N_MISSED_OUTBREAKS_COLOR,
        "#2B3465",
    ],
)

save(
    plotsdir("schematic-plot_incidence-threshold.svg"),
    incidence_threshold_schematic;
    size = (1300, 800 / 3),
)

#%%

inc_alertstatus_vec = detectoutbreak(
    inc_vec,
    6,
)
inc_alert_bounds = calculate_outbreak_thresholds(
    rle(inc_alertstatus_vec .> 0); ncols = 3
)
OutbreakDetectionUtils.calculate_outbreak_duration!(inc_alert_bounds)

incidence_threshold_schematic_6 = plot_schematic(
    inc_vec,
    outbreak_status,
    outbreak_bounds[:, 1:2],
    inc_alertstatus_vec,
    inc_alert_bounds,
    6;
    time_p = time_p,
    shade_alert_outbreak_overlap = true,
    xlims = (5, 13),
    measlesalpha = 0.4,
    testalpha = 0.4,
    outbreakcolormap = [
        N_MISSED_OUTBREAKS_COLOR,
        "#CE6F58",
    ],
    alertcolormap = [
        N_MISSED_OUTBREAKS_COLOR,
        "#2B3465",
    ],
)

save(
    plotsdir("schematic-plot_incidence-threshold-6.svg"),
    incidence_threshold_schematic_6;
    size = (1300, 800 / 3),
)

#%%
function plot_schematic(
    inc_vec,
    outbreakstatus_vec,
    outbreak_bounds,
    outbreak_specification,
    noise_vec,
    testpositive_vec,
    alertstatus_vec,
    alert_bounds,
    alertthreshold;
    time_p = SimTimeParameters(; tmin = 0.0, tmax = 365.0 * 10, tstep = 1.0),
    outbreakcolormap = [
        N_MISSED_OUTBREAKS_COLOR, PERC_OUTBREAKS_DETECTED_COLOR
    ],
    alertcolormap = [
        N_MISSED_OUTBREAKS_COLOR, N_ALERTS_COLOR
    ],
    shade_alert_outbreak_overlap = false,
    measlesalpha = 0.5,
    testalpha = 0.5,
    kwargs...,
)
    kwargs_dict = Dict(kwargs)

    times = collect(time_p.trange)

    if haskey(kwargs_dict, :xlims)
        lower = maximum([1, kwargs_dict[:xlims][1] * 365])
        upper = minimum([Int64(time_p.tlength), kwargs_dict[:xlims][2] * 365])
        times = times[lower:upper]
        inc_vec = inc_vec[lower:upper]
        outbreakstatus_vec = outbreakstatus_vec[lower:upper]
        noise_vec = noise_vec[lower:upper]
        testpositive_vec = testpositive_vec[lower:upper]
        alertstatus_vec = alertstatus_vec[lower:upper]

        outbreak_bounds = outbreak_bounds[
            (outbreak_bounds[:, 1] .>= lower) .& (outbreak_bounds[:, 2] .<= upper),
            :,
        ]

        alert_bounds = alert_bounds[
            (alert_bounds[:, 1] .>= lower) .& (alert_bounds[:, 2] .<= upper), :,
        ]
    end
    outbreak_bounds_vec = vec(outbreak_bounds[:, 1:2])
    alert_bounds_vec = vec(alert_bounds[:, 1:2])

    fig = Figure()
    incga = fig[1, 1] = GridLayout()
    noisega = fig[2, 1] = GridLayout()
    testga = fig[3, 1] = GridLayout()
    incax = Axis(incga[1, 1]; ylabel = "Measles Incidence")
    noiseax = Axis(noisega[1, 1]; ylabel = "Noise Incidence")
    testax = Axis(testga[1, 1]; xlabel = "Time", ylabel = "Test Positives")

    if shade_alert_outbreak_overlap
        if !isempty(outbreak_bounds_vec)
            map(
                ax -> vspan!(
                    ax,
                    outbreak_bounds[:, 1],
                    outbreak_bounds[:, 2];
                    color = (outbreakcolormap[2], measlesalpha),
                ),
                [incax, noiseax, testax],
            )
        end

        if !isempty(alert_bounds_vec)
            vspan!(
                testax,
                alert_bounds[:, 1],
                alert_bounds[:, 2];
                color = (alertcolormap[2], testalpha),
            )
        end
    end

    lines!(
        incax,
        times,
        inc_vec;
        color = outbreakstatus_vec,
        colormap = outbreakcolormap,
        linewidth = 3,
    )

    lines!(
        noiseax,
        times,
        noise_vec;
        color = N_MISSED_OUTBREAKS_COLOR,
        linewidth = 3,
    )

    lines!(
        testax,
        times,
        testpositive_vec;
        color = alertstatus_vec,
        colormap = alertcolormap,
        linewidth = 3,
    )

    hlines!(
        testax,
        alertthreshold;
        color = :black,
        linewidth = 2,
        linestyle = :dash,
    )

    text!(
        testax,
        times[1] + 1,
        alertthreshold + 0.5;
        text = "T = $alertthreshold",
        justification = :left,
    )

    map(
        ax -> hidexdecorations!(ax),
        [
            noiseax,
            incax,
            testax,
        ],
    )

    linkxaxes!(noiseax, incax, testax)

    return fig
end

#%%
function plot_schematic(
    inc_vec,
    outbreakstatus_vec,
    outbreak_bounds,
    noise_vec,
    testpositive_vec;
    time_p = SimTimeParameters(; tmin = 0.0, tmax = 365.0 * 10, tstep = 1.0),
    outbreakcolormap = [
        N_MISSED_OUTBREAKS_COLOR, PERC_OUTBREAKS_DETECTED_COLOR
    ],
    shade_alert_outbreak_overlap = false,
    measlesalpha = 0.5,
    kwargs...,
)
    kwargs_dict = Dict(kwargs)

    times = collect(time_p.trange)

    if haskey(kwargs_dict, :xlims)
        lower = maximum([1, kwargs_dict[:xlims][1] * 365])
        upper = minimum([Int64(time_p.tlength), kwargs_dict[:xlims][2] * 365])
        times = times[lower:upper]
        inc_vec = inc_vec[lower:upper]
        outbreakstatus_vec = outbreakstatus_vec[lower:upper]
        noise_vec = noise_vec[lower:upper]
        testpositive_vec = testpositive_vec[lower:upper]

        outbreak_bounds = outbreak_bounds[
            (outbreak_bounds[:, 1] .>= lower) .& (outbreak_bounds[:, 2] .<= upper),
            :,
        ]
    end
    outbreak_bounds_vec = vec(outbreak_bounds[:, 1:2])

    fig = Figure()
    incga = fig[1, 1] = GridLayout()
    noisega = fig[2, 1] = GridLayout()
    testga = fig[3, 1] = GridLayout()
    incax = Axis(incga[1, 1]; ylabel = "Measles Incidence")
    noiseax = Axis(noisega[1, 1]; ylabel = "Noise Incidence")
    testax = Axis(testga[1, 1]; xlabel = "Time", ylabel = "Test Positives")

    if shade_alert_outbreak_overlap
        if !isempty(outbreak_bounds_vec)
            map(
                ax -> vspan!(
                    ax,
                    outbreak_bounds[:, 1],
                    outbreak_bounds[:, 2];
                    color = (outbreakcolormap[2], measlesalpha),
                ),
                [incax, noiseax, testax],
            )
        end
    end

    lines!(
        incax,
        times,
        inc_vec;
        color = outbreakstatus_vec,
        colormap = outbreakcolormap,
        linewidth = 3,
    )

    lines!(
        noiseax,
        times,
        noise_vec;
        color = N_MISSED_OUTBREAKS_COLOR,
        linewidth = 3,
    )

    lines!(
        testax,
        times,
        testpositive_vec;
        color = N_MISSED_OUTBREAKS_COLOR,
        linewidth = 3,
    )

    map(
        ax -> hidexdecorations!(ax),
        [
            noiseax,
            incax,
            testax,
        ],
    )

    linkxaxes!(noiseax, incax, testax)

    return fig
end

#%%
schematic_with_shade_fig = plot_schematic(
    inc_vec,
    outbreak_status,
    outbreak_bounds[:, 1:2],
    outbreak_specification,
    noise_vec,
    movingavg_testpositives,
    alertstatus_vec,
    alert_bounds,
    outbreak_detection_specification.alert_threshold;
    time_p = time_p,
    shade_alert_outbreak_overlap = true,
    xlims = (5, 13),
    measlesalpha = 0.4,
    testalpha = 0.4,
    outbreakcolormap = [
        N_MISSED_OUTBREAKS_COLOR,
        "#CE6F58",
    ],
    alertcolormap = [
        N_MISSED_OUTBREAKS_COLOR,
        "#2B3465",
    ],
)

save(
    plotsdir("schematic-plot_threshold-5.svg"),
    schematic_with_shade_fig,
)

#%%
test_6_1_alertstatus_vec = detectoutbreak(
    movingavg_testpositives,
    6.2,
)
test_6_1_alert_bounds = calculate_outbreak_thresholds(
    rle(test_6_1_alertstatus_vec .> 0); ncols = 3
)

schematic_with_shade_fig = plot_schematic(
    inc_vec,
    outbreak_status,
    outbreak_bounds[:, 1:2],
    outbreak_specification,
    noise_vec,
    movingavg_testpositives,
    test_6_1_alertstatus_vec,
    test_6_1_alert_bounds,
    6;
    time_p = time_p,
    shade_alert_outbreak_overlap = true,
    xlims = (5, 13),
    measlesalpha = 0.4,
    testalpha = 0.4,
    outbreakcolormap = [
        N_MISSED_OUTBREAKS_COLOR,
        "#CE6F58",
    ],
    alertcolormap = [
        N_MISSED_OUTBREAKS_COLOR,
        "#2B3465",
    ],
)

save(
    plotsdir("schematic-plot_threshold-6.svg"),
    schematic_with_shade_fig,
)

#%%
schematic_no_shade_fig = plot_schematic(
    inc_vec,
    outbreak_status,
    outbreak_bounds[:, 1:2],
    noise_vec,
    movingavg_testpositives;
    time_p = time_p,
    shade_alert_outbreak_overlap = true,
    xlims = (5, 13),
    measlesalpha = 0.4,
    testalpha = 0.4,
    outbreakcolormap = [
        N_MISSED_OUTBREAKS_COLOR,
        "#CE6F58",
    ],
)

save(
    plotsdir("schematic-plot_no-threshold.svg"),
    schematic_no_shade_fig,
)

#%%
test_5_1_alertstatus_vec = detectoutbreak(
    movingavg_testpositives,
    5.105,
)
test_5_1_alert_bounds = calculate_outbreak_thresholds(
    rle(test_5_1_alertstatus_vec .> 0); ncols = 3
)

schematic_with_shade_fig = plot_schematic(
    inc_vec,
    outbreak_status,
    outbreak_bounds[:, 1:2],
    outbreak_specification,
    noise_vec,
    movingavg_testpositives,
    test_5_1_alertstatus_vec,
    test_5_1_alert_bounds,
    5;
    time_p = time_p,
    shade_alert_outbreak_overlap = true,
    xlims = (5, 13),
    measlesalpha = 0.4,
    testalpha = 0.4,
    outbreakcolormap = [
        N_MISSED_OUTBREAKS_COLOR,
        "#CE6F58",
    ],
    alertcolormap = [
        N_MISSED_OUTBREAKS_COLOR,
        "#2B3465",
    ],
)

save(
    plotsdir("schematic-plot_threshold-5.svg"),
    schematic_with_shade_fig,
)

#%%
test_3_alertstatus_vec = detectoutbreak(
    movingavg_testpositives,
    measles_incidence_threshold,
)
test_3_alert_bounds = calculate_outbreak_thresholds(
    rle(test_3_alertstatus_vec .> 0); ncols = 3
)

schematic_with_shade_fig = plot_schematic(
    inc_vec,
    outbreak_status,
    outbreak_bounds[:, 1:2],
    outbreak_specification,
    noise_vec,
    movingavg_testpositives,
    test_3_alertstatus_vec,
    test_3_alert_bounds,
    3;
    time_p = time_p,
    shade_alert_outbreak_overlap = true,
    xlims = (5, 13),
    measlesalpha = 0.4,
    testalpha = 0.4,
    outbreakcolormap = [
        N_MISSED_OUTBREAKS_COLOR,
        "#CE6F58",
    ],
    alertcolormap = [
        N_MISSED_OUTBREAKS_COLOR,
        "#2B3465",
    ],
)

save(
    plotsdir("schematic-plot_threshold-3.svg"),
    schematic_with_shade_fig,
)

#%%
function inc_noise_schematic(
    inc_vec,
    outbreakstatus_vec,
    outbreak_bounds,
    noise_vec;
    time_p = SimTimeParameters(; tmin = 0.0, tmax = 365.0 * 10, tstep = 1.0),
    noisealpha = 0.2,
    measlesalpha = 0.5,
    outbreakcolormap = [
        (N_MISSED_OUTBREAKS_COLOR, measlesalpha),
        (PERC_OUTBREAKS_DETECTED_COLOR, measlesalpha),
    ],
    noise_color = :darkred,
    shade_alert_outbreak_overlap = false,
    size = (1300, 800 / 3),
    kwargs...,
)
    kwargs_dict = Dict(kwargs)

    times = collect(time_p.trange)

    if haskey(kwargs_dict, :xlims)
        lower = maximum([1, kwargs_dict[:xlims][1] * 365])
        upper = minimum([Int64(time_p.tlength), kwargs_dict[:xlims][2] * 365])
        times = times[lower:upper]
        inc_vec = inc_vec[lower:upper]
        outbreakstatus_vec = outbreakstatus_vec[lower:upper]
        noise_vec = noise_vec[lower:upper]

        outbreak_bounds = outbreak_bounds[
            (outbreak_bounds[:, 1] .>= lower) .& (outbreak_bounds[:, 2] .<= upper),
            :,
        ]
    end
    outbreak_bounds_vec = vec(outbreak_bounds[:, 1:2])

    fig = Figure(; size = size)
    incga = fig[1, 1] = GridLayout()
    incax = Axis(incga[1, 1]; ylabel = "Incidence")

    if shade_alert_outbreak_overlap
        if !isempty(outbreak_bounds_vec)
            vspan!(
                incax,
                outbreak_bounds[:, 1],
                outbreak_bounds[:, 2];
                color = (outbreakcolormap[2], measlesalpha),
            )
        end
    end

    band!(
        incax,
        times,
        repeat([0.0], length(times)),
        noise_vec;
        color = (noise_color, noisealpha),
    )

    lines!(
        incax,
        times,
        inc_vec;
        color = outbreakstatus_vec,
        colormap = outbreakcolormap,
        linewidth = 3,
    )

    if haskey(kwargs_dict, :ylims)
        ylims!(incax, kwargs_dict[:ylims])
    end

    # band!(
    #     incax,
    #     times,
    #     repeat([0.0], length(times)),
    #     inc_vec;
    #     color = outbreakstatus_vec,
    #     colormap = outbreakcolormap,
    #     alpha = measlesalpha,
    # )
    #
    # lines!(
    #     incax,
    #     times,
    #     noise_vec;
    #     color = noise_color,
    #     linewidth = 3,
    # )

    hidexdecorations!(incax)

    return fig
end

#%%
dynamical_noise = create_noise_arr(
    DynamicalNoiseSpecification(
        "dynamical",
        5.0,
        7,
        14,
        "in-phase",
        0.15,
        0.85,
    ),
    inc_vec;
    ensemble_specification = ensemble_specification,
    seed = 1234,
)[1][
    :, 1
];
mean(dynamical_noise)

dynamical_noise = round.(
    calculate_movingavg(
        dynamical_noise,
        movingavg_window,
    )
)

dynamical_noise_schematic = inc_noise_schematic(
    inc_vec,
    outbreak_status,
    outbreak_bounds[:, 1:2],
    dynamical_noise;
    time_p = time_p,
    xlims = (5, 13),
    ylims = (0, 20),
    size = (1300, 800 / 3),
    shade_alert_outbreak_overlap = false,
    measlesalpha = 1,
    noisealpha = 0.4,
)

save(
    plotsdir("dynamical-noise-schematic.svg"),
    dynamical_noise_schematic;
)

#%%
poisson_noise = create_noise_arr(
    PoissonNoiseSpecification("poisson", 2.0),
    inc_vec;
    seed = 1234,
)[1]

poisson_noise = round.(
    calculate_movingavg(
        poisson_noise,
        movingavg_window,
    )
)

poisson_noise_schematic = inc_noise_schematic(
    inc_vec,
    outbreak_status,
    outbreak_bounds[:, 1:2],
    poisson_noise;
    time_p = time_p,
    xlims = (5, 13),
    ylims = (0, 20),
    size = (1300, 800 / 3),
    shade_alert_outbreak_overlap = false,
    measlesalpha = 1.0,
    noisealpha = 0.4,
)

save(
    plotsdir("poisson-noise-schematic.svg"),
    poisson_noise_schematic;
)
