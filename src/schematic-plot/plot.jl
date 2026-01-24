export plot_schematic

"""
    plot_schematic(
        optimization_result::OutbreakDetectionCore.OptimizationResult;
        seed = 1234,
        noise_scaling = 10.0,
        shift_noise = 0,
        movingavg_window = 20,
        kwargs...
    )

Create a schematic plot from an OptimizationResult by regenerating the simulation.

This method extracts parameters from an OptimizationResult, regenerates the simulation
data, and creates a schematic visualization showing the outbreak detection process.

# Arguments
- `optimization_result::OptimizationResult`: Result from optimization containing all
  necessary parameters for simulation

# Keyword Arguments
- `seed`: Random seed for reproducibility (default: 1234)
- `noise_scaling`: Multiplicative factor for noise intensity (default: 10.0)
- `shift_noise`: Number of time steps to shift noise signal (default: 0)
- `movingavg_window`: Window size for smoothing incidence data (default: 20)
- `kwargs...`: Additional keyword arguments passed to the base plot_schematic function
  (e.g., xlims, shade_alert_outbreak_overlap, measlesalpha, testalpha, colormaps)

# Returns
- `Figure`: Makie figure object containing the schematic plot

# Examples
```julia
# From a StructVector of results, plot a specific result
results::StructVector{OptimizationResult} = load_optimization_results(...)
fig = plot_schematic(results[1]; xlims = (5, 13), shade_alert_outbreak_overlap = true)

# Plot with custom noise parameters
fig = plot_schematic(
    results[1];
    seed = 42,
    noise_scaling = 15.0,
    shift_noise = -100,
    shade_alert_outbreak_overlap = true
)

# Filter results and plot
filtered_results = results[results.noise_level .== 10.0]
fig = plot_schematic(filtered_results[1])
```

# Notes
- This method regenerates simulation data from the parameters stored in OptimizationResult
- The simulation uses the optimal threshold found during optimization
- Noise parameters (noise_scaling, shift_noise) are not stored in OptimizationResult,
  so they must be provided as keyword arguments if different from defaults
- The percent_visit_clinic parameter is assumed to be 1.0 as it's not stored in
  OptimizationResult
"""
function plot_schematic(
        optimization_result::OutbreakDetectionCore.OptimizationResult;
        seed = 1234,
        noise_scaling = 10.0,
        shift_noise = 0,
        movingavg_window = 20,
        kwargs...
    )
    # Extract parameters from OptimizationResult
    ensemble_spec = optimization_result.ensemble_specification
    test_spec = optimization_result.test_specification
    outbreak_spec = optimization_result.outbreak_specification
    alert_method = optimization_result.alert_method
    optimal_threshold = optimization_result.optimal_threshold

    # Get time parameters from ensemble specification
    time_p = ensemble_spec.time_p

    # Extract state and dynamics parameters for main epidemic
    states_p = ensemble_spec.target_disease_states_p
    dynamics_p = ensemble_spec.target_disease_dynamics_p

    # Extract noise parameters
    noise_spec = ensemble_spec.noise_specification
    noise_states_p = noise_spec.noise_states_p
    noise_dynamics_p = noise_spec.noise_dynamics_p

    # Calculate percent tested
    percent_tested = optimization_result.percent_tested

    # Create outbreak detection specification from alert method and optimal threshold
    # Extract moving average lag from alert method
    moving_avg_lag = if alert_method.method isa OutbreakDetectionCore.MovingAverage
        alert_method.method.lag
    else
        movingavg_window  # fallback to default
    end

    outbreak_detection_spec = OutbreakDetectionCore.OutbreakDetectionSpecification(
        optimal_threshold,
        moving_avg_lag,
        1.0,  # percent_visit_clinic (assumed, as not stored in OptimizationResult)
        percent_tested,
        "movingavg"
    )

    # Generate simulation data
    inc_vec, outbreak_status, outbreak_bounds, noise_vec, movingavg_testpositives,
        alertstatus_vec, alert_bounds = create_schematic_simulation(
        states_p,
        dynamics_p,
        noise_states_p,
        noise_dynamics_p,
        test_spec,
        time_p;
        seed = seed,
        outbreak_specification = outbreak_spec,
        outbreak_detection_specification = outbreak_detection_spec,
        noise_scaling = noise_scaling,
        shift_noise = shift_noise,
        movingavg_window = movingavg_window
    )

    # Call the base plot_schematic function
    return plot_schematic(
        inc_vec,
        outbreak_status,
        outbreak_bounds[:, 1:2],
        outbreak_spec,
        noise_vec,
        movingavg_testpositives,
        alertstatus_vec,
        alert_bounds,
        optimal_threshold;
        time_p = time_p,
        kwargs...
    )
end

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
            "#5E5C6C",
            "#CE6F58",
        ],
        alertcolormap = [
            "#F4A157",
            "#2B3465",
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

        if !isempty(alert_bounds)
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

    hlines!(
        incax,
        outbreak_specification.outbreak_threshold;
        color = :black,
        linewidth = 2,
        linestyle = :dash,
    )

    lines!(
        noiseax,
        times,
        noise_vec;
        color = outbreakcolormap[1],
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

    for (label, layout) in
        zip(["a", "b", "c"], [incga[1, 1], noisega[1, 1], testga[1, 1]])
        Label(
            layout[1, 1, TopLeft()], label;
            fontsize = 30,
            font = :bold,
            padding = (0, 0, 20, 0),
            halign = :right
        )
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
