#%%
using DrWatson
@quickactivate "OutbreakDetection"
using Revise: includet

using GLMakie
using StatsBase

using OutbreakDetection

includet(srcdir("makie-plotting-setup.jl"))

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
    seasonality = sin
)

test_specification = IndividualTestSpecification(0.85, 0.85, 0)

time_p = SimTimeParameters(; tmin = 0.0, tmax = 365.0 * 20, tstep = 1.0)

outbreak_specification = OutbreakSpecification(5, 30, 500)

outbreak_detection_specification = OutbreakDetectionSpecification(
    5,
    7,
    1.0,
    0.75,
    "movingavg"
)

function get_outbreak_status(
    inc_vec, outbreak_specification
)
    abovethreshold_vec = vec(
        inc_vec .>= outbreak_specification.outbreak_threshold
    )

    abovethresholdrle = rle(abovethreshold_vec)

    all_outbreak_thresholds = calculate_outbreak_thresholds(
        abovethresholdrle; ncols = 4
    )

    OutbreakDetection.classify_all_outbreaks!(
        inc_vec,
        abovethreshold_vec,
        all_outbreak_thresholds,
        outbreak_specification.minimum_outbreak_duration,
        outbreak_specification.minimum_outbreak_size,
    )

    outbreak_thresholds = OutbreakDetection.filter_only_outbreaks(
        all_outbreak_thresholds
    )

    outbreak_status = zeros(Int64, length(inc_vec))

    for (lower, upper) in eachrow(outbreak_thresholds)
        # @show lower, upper
        outbreak_status[lower:upper] .= 1
    end
    return outbreak_status
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
)
    inc_sv = seir_mod(
        states_p.init_states,
        dynamics_p,
        time_p;
        seed = seed,
    )[2]

    inc_vec = vec(convert_svec_to_matrix(inc_sv))

    outbreak_status = get_outbreak_status(
        inc_vec, outbreak_specification
    )

    noise_sv = seir_mod(
        noise_states_p.init_states,
        noise_dynamics_p,
        time_p;
        seed = seed
    )[2]

    noise_vec = calculate_movingavg(
        vec(convert_svec_to_matrix(noise_sv)) .* noise_scaling, 7
    )

    perc_tested =
        outbreak_detection_specification.percent_visit_clinic *
        outbreak_detection_specification.percent_clinic_tested

    true_positives = OutbreakDetection.calculate_positives(
        calculate_true_positives!,
        inc_vec .* perc_tested,
        time_p.tlength,
        test_specification.test_result_lag,
        test_specification.sensitivity,
    )

    false_positives = OutbreakDetection.calculate_positives(
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
        outbreak_detection_specification.moving_average_lag
    )

    alertstatus_vec = OutbreakDetection.detectoutbreak(
        movingavg_testpositives,
        outbreak_detection_specification.alert_threshold,
    )

    return inc_vec,
    outbreak_status, noise_vec, movingavg_testpositives,
    alertstatus_vec
end

inc_vec, outbreak_status, noise_vec, movingavg_testpositives, alertstatus_vec = create_schematic_simulation(
    states_p,
    dynamics_p,
    noise_states_p,
    noise_dynamics_p,
    test_specification,
    time_p;
    seed = 12345,
    outbreak_detection_specification = outbreak_detection_specification,
    noise_scaling = 10,
);

function plot_schematic(
    inc_vec,
    outbreakstatus_vec,
    outbreak_specification,
    noise_vec,
    testpositive_vec,
    alertstatus_vec,
    alertthreshold;
    time_p = SimTimeParameters(; tmin = 0.0, tmax = 365.0 * 10, tstep = 1.0),
    outbreakcolormap = [
        N_MISSED_OUTBREAKS_COLOR, PERC_OUTBREAKS_DETECTED_COLOR
    ],
    alertcolormap = [
        N_MISSED_OUTBREAKS_COLOR, N_ALERTS_COLOR
    ],
    years = false,
    kwargs...,
)
    if years == true
        times = time_p.trange / 365
    else
        times = time_p.trange
    end

    kwargs_dict = Dict(kwargs)

    if haskey(kwargs_dict, :xlims)
        lower = kwargs_dict[:xlims][1] * 365
        upper = kwargs_dict[:xlims][2] * 365
        times = times[lower:upper]
        inc_vec = inc_vec[lower:upper]
        outbreakstatus_vec = outbreakstatus_vec[lower:upper]
        noise_vec = noise_vec[lower:upper]
        testpositive_vec = testpositive_vec[lower:upper]
        alertstatus_vec = alertstatus_vec[lower:upper]
    end

    fig = Figure()
    noiseax = Axis(fig[1, 1]; ylabel = "Noise")
    incax = Axis(fig[2, 1]; ylabel = "Incidence")
    testax = Axis(fig[3, 1]; xlabel = "Time", ylabel = "Test Positives")

    lines!(
        noiseax,
        times,
        noise_vec;
        color = :black,
        linewidth = 2.5,
    )

    lines!(
        incax,
        times,
        inc_vec;
        color = outbreakstatus_vec,
        colormap = outbreakcolormap,
        linewidth = 2.5,
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
        linewidth = 2.5,
    )

    hlines!(
        testax,
        alertthreshold;
        color = :black,
        linewidth = 2,
        linestyle = :dash
    )

    map(
        ax -> hidexdecorations!(ax),
        [noiseax, incax, testax],
    )

    linkxaxes!(noiseax, incax, testax)

    return fig
end

plot_schematic(
    inc_vec, outbreak_status, outbreak_specification, noise_vec,
    movingavg_testpositives,
    alertstatus_vec,
    outbreak_detection_specification.alert_threshold;
    time_p = time_p,
    years = true,
    xlims = (5, 13),
)
