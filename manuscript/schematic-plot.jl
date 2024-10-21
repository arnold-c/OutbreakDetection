using DrWatson
@quickactivate "OutbreakDetection"

using GLMakie
using StatsBase

using OutbreakDetection
using OutbreakDetectionUtils

include(srcdir("makie-plotting-setup.jl"))

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
    8,
    movingavg_window,
    1.0,
    0.75,
    "movingavg",
)

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
        rle(alertstatus_vec .> 0); ncols = 2
    )

    return inc_vec,
    outbreak_status,
    outbreak_bounds,
    noise_vec,
    movingavg_testpositives,
    alertstatus_vec,
    alert_bounds
end

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
    kwargs...,
)
    kwargs_dict = Dict(kwargs)

    times = collect(time_p.trange)

    if haskey(kwargs_dict, :xlims)
        lower = kwargs_dict[:xlims][1] * 365
        upper = kwargs_dict[:xlims][2] * 365
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
        outbreak_bounds_vec = vec(outbreak_bounds)

        alert_bounds = alert_bounds[
            (alert_bounds[:, 1] .>= lower) .& (alert_bounds[:, 2] .<= upper), :,
        ]
        alert_bounds_vec = vec(alert_bounds)
    end

    fig = Figure()
    noisega = fig[1, 1] = GridLayout()
    incga = fig[2, 1] = GridLayout()
    testga = fig[3, 1] = GridLayout()
    noiseax = Axis(noisega[1, 1]; ylabel = "Noise")
    incax = Axis(incga[1, 1]; ylabel = "Incidence")
    testax = Axis(testga[1, 1]; xlabel = "Time", ylabel = "Test Positives")

    lines!(
        noiseax,
        times,
        noise_vec;
        color = :black,
        linewidth = 3,
    )

    lines!(
        incax,
        times,
        inc_vec;
        color = outbreakstatus_vec,
        colormap = outbreakcolormap,
        linewidth = 3,
    )

    hlines!(
        incax,
        outbreak_specification.outbreak_threshold;
        color = :black,
        linewidth = 2,
        linestyle = :dash,
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

    if shade_alert_outbreak_overlap
        if !isempty(outbreak_bounds_vec)
            vspan!(
                incax,
                outbreak_bounds[:, 1],
                outbreak_bounds[:, 2];
                color = (outbreakcolormap[2], 0.2),
            )

            vlines!(testax,
                outbreak_bounds_vec;
                color = outbreakcolormap[2],
                linewidth = 3,
            )
        end

        if !isempty(alert_bounds)
            vspan!(
                testax,
                alert_bounds[:, 1],
                alert_bounds[:, 2];
                color = (alertcolormap[2], 0.2),
            )

            vlines!(
                incax,
                alert_bounds_vec;
                color = alertcolormap[2],
                linewidth = 3,
            )
        end
    end

    for (label, layout) in
        zip(["a", "b", "c"], [noisega[1, 1], incga[1, 1], testga[1, 1]])
        Label(layout[1, 1, TopLeft()], label;
            fontsize = 30,
            font = :bold,
            padding = (0, 0, 20, 0),
            halign = :right)
    end

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

schematic_with_shade_fig = plot_schematic(
    inc_vec,
    outbreak_status,
    outbreak_bounds[:, 1:2],
    outbreak_specification,
    noise_vec,
    movingavg_testpositives,
    alertstatus_vec,
    alert_bounds[
        (@view(alert_bounds[:, 2]) .- @view(alert_bounds[:, 1]) .> 30), :,
    ],
    outbreak_detection_specification.alert_threshold; time_p = time_p,
    shade_alert_outbreak_overlap = true,
    xlims = (5, 13),
)
