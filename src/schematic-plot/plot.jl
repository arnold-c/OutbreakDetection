export plot_schematic, visualize_timeseries

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
    time_p = ensemble_spec.time_parameters

    # Extract state and dynamics parameters for main epidemic
    states_p = ensemble_spec.state_parameters
    dynamics_spec = ensemble_spec.dynamics_parameter_specification
    dynamics_p = OutbreakDetectionCore.DynamicsParameters(
        dynamics_spec;
        seed = seed
    )

    # Extract noise parameters
    noise_spec = OutbreakDetectionCore.DynamicalNoiseSpecification(
        ensemble_spec.dynamical_noise_params,
        0.6  # vaccination_coverage
    )


    # Calculate percent tested
    percent_tested = optimization_result.percent_tested

    # Create outbreak detection specification from alert method and optimal threshold
    # Extract moving average lag from alert method
    moving_avg_lag = if alert_method isa OutbreakDetectionCore.MovingAverage
        alert_method.window
    else
        movingavg_window  # fallback to default
    end

    outbreak_detection_spec = OutbreakDetectionCore.OutbreakDetectionSpecification(
        optimal_threshold,
        moving_avg_lag,
        1.0,  # percent_visit_clinic (assumed, as not stored in OptimizationResult)
        percent_tested,
        percent_tested,
        alert_method,
        ""
    )

    # Generate simulation data
    inc_vec, outbreak_status, outbreak_bounds, noise_vec, movingavg_testpositives,
        alertstatus_vec, alert_bounds = create_schematic_simulation(
        states_p,
        dynamics_p,
        dynamics_spec,
        noise_spec,
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

"""
    visualize_timeseries(
        optimization_result::OutbreakDetectionCore.OptimizationResult;
        sim_number = 1,
        seed = 1234,
        noise_scaling = 10.0,
        shift_noise = 0,
        movingavg_window = 20,
        kwargs...
    )

Visualize the actual simulation timeseries for an optimized result.

This function recreates the ensemble simulation and noise vectors from an OptimizationResult,
allowing visualization of specific simulation runs. Unlike `plot_schematic`, which creates
a single representative simulation, this function recreates the full ensemble and displays
a specific simulation number.

# Arguments
- `optimization_result::OptimizationResult`: Result from optimization containing all
  necessary parameters for simulation

# Keyword Arguments
- `sim_number`: Which simulation from the ensemble to visualize (default: 1)
- `seed`: Random seed for reproducibility (default: 1234)
- `noise_scaling`: Multiplicative factor for noise intensity (default: 10.0)
- `shift_noise`: Number of time steps to shift noise signal (default: 0)
- `movingavg_window`: Window size for smoothing incidence data (default: 20)
- `kwargs...`: Additional keyword arguments passed to the base plot_schematic function
  (e.g., xlims, shade_alert_outbreak_overlap, measlesalpha, testalpha, colormaps)

# Returns
- `Figure`: Makie figure object containing the timeseries visualization

# Examples
```julia
# Load optimization results
optimized_filedir = OutbreakDetectionCore.outdir("ensemble", "threshold-optimization")
optimized_threshold_results = Try.unwrap(
    OutbreakDetectionCore.load_previous_optimization_results_structvector(
        optimized_filedir,
        "threshold-optimization.jld2",
        joinpath(optimized_filedir, "checkpoints"),
    )
)

# Visualize first simulation of first result
fig = visualize_timeseries(optimized_threshold_results[1])

# Visualize 5th simulation with custom parameters
fig = visualize_timeseries(
    optimized_threshold_results[1];
    sim_number = 5,
    xlims = (5, 13),
    shade_alert_outbreak_overlap = true
)
```

# Notes
- This function regenerates the full ensemble simulation from parameters stored in OptimizationResult
- The simulation uses the optimal threshold found during optimization
- Noise parameters (noise_scaling, shift_noise) are not stored in OptimizationResult,
  so they must be provided as keyword arguments if different from defaults
"""
function visualize_timeseries(
        optimization_result::OutbreakDetectionCore.OptimizationResult;
        sim_number = 1,
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
    time_p = ensemble_spec.time_parameters

    # Extract state and dynamics parameters for main epidemic
    states_p = ensemble_spec.state_parameters
    dynamics_spec = ensemble_spec.dynamics_parameter_specification

    # Extract noise parameters
    noise_spec = OutbreakDetectionCore.DynamicalNoiseSpecification(
        ensemble_spec.dynamical_noise_params,
        0.6  # vaccination_coverage
    )

    # Calculate percent tested
    percent_tested = optimization_result.percent_tested

    # Create outbreak detection specification from alert method and optimal threshold
    # Extract moving average lag from alert method
    moving_avg_lag = if alert_method isa OutbreakDetectionCore.MovingAverage
        alert_method.window
    else
        movingavg_window  # fallback to default
    end

    outbreak_detection_spec = OutbreakDetectionCore.OutbreakDetectionSpecification(
        optimal_threshold,
        moving_avg_lag,
        1.0,  # percent_visit_clinic (assumed, as not stored in OptimizationResult)
        percent_tested,
        percent_tested,
        alert_method,
        ""
    )

    # Generate ensemble simulation
    ensemble_results = OutbreakDetectionCore.generate_single_ensemble(
        ensemble_spec;
        seed = seed
    )

    # Check that sim_number is valid
    if sim_number < 1 || sim_number > length(ensemble_results)
        error("sim_number must be between 1 and $(length(ensemble_results))")
    end

    # Extract incidence from the specified simulation
    inc_vec_raw = ensemble_results[sim_number].incidence

    # Apply moving average smoothing
    inc_vec = round.(
        OutbreakDetectionCore.calculate_movingavg(
            inc_vec_raw,
            movingavg_window,
        )
    )

    # Get outbreak status
    outbreak_status, outbreak_bounds = get_outbreak_status(
        inc_vec, outbreak_spec
    )

    # Generate noise using proper framework
    # Convert DynamicalNoiseSpecification back to DynamicalNoiseParameters for ensemble spec
    vac_cov = noise_spec.vaccination_coverage
    noise_params = OutbreakDetectionCore.DynamicalNoiseParameters(;
        R_0 = noise_spec.R_0,
        latent_period = Dates.Day(round(Int, noise_spec.latent_period)),
        infectious_duration = Dates.Day(round(Int, noise_spec.infectious_duration)),
        correlation = noise_spec.correlation,
        poisson_component = noise_spec.poisson_component,
        vaccination_bounds = [vac_cov, vac_cov]  # Equal bounds for fixed vaccination coverage
    )

    # Create a minimal ensemble specification for noise generation
    noise_ensemble_spec = OutbreakDetectionCore.EnsembleSpecification(
        "timeseries_viz",
        states_p,
        time_p,
        dynamics_spec,
        noise_params,
        ensemble_spec.nsims,  # Generate same number of noise simulations
        ""  # dirpath
    )

    # Recreate noise dynamics
    noise_dynamics_spec = OutbreakDetectionCore.recreate_noise_dynamics_parameter_specification(
        noise_spec,
        noise_ensemble_spec
    )

    noise_dynamics_p = OutbreakDetectionCore.DynamicsParameters(;
        beta_mean = noise_dynamics_spec.beta_mean,
        beta_force = noise_dynamics_spec.beta_force,
        seasonality = noise_dynamics_spec.seasonality,
        sigma = noise_dynamics_spec.sigma,
        gamma = noise_dynamics_spec.gamma,
        mu = noise_dynamics_spec.mu,
        annual_births_per_k = noise_dynamics_spec.annual_births_per_k,
        epsilon = noise_dynamics_spec.epsilon,
        R_0 = noise_dynamics_spec.R_0,
        vaccination_coverage = noise_spec.vaccination_coverage
    )

    # Generate noise
    noise_result = OutbreakDetectionCore.create_noise_vecs(
        noise_spec,
        noise_ensemble_spec,
        noise_dynamics_p;
        seed = seed
    )

    # Extract noise incidence for the specified simulation
    noise_vec_raw = LightSumTypes.variant(noise_result).incidence[sim_number]

    noise_vec_tmp = OutbreakDetectionCore.calculate_movingavg(
        noise_vec_raw .* noise_scaling,
        movingavg_window
    )

    noise_vec = shift_vec(noise_vec_tmp, shift_noise)

    perc_tested =
        outbreak_detection_spec.percent_visit_clinic *
        outbreak_detection_spec.percent_clinic_tested

    # Calculate number tested (convert to Int)
    tested_inc_vec = round.(Int, inc_vec .* perc_tested)
    tested_noise_vec = round.(Int, noise_vec .* perc_tested)

    # Allocate output vectors
    true_positives = Vector{Int64}(undef, time_p.tlength)
    false_positives = Vector{Int64}(undef, time_p.tlength)

    # Calculate true positives from incidence
    OutbreakDetectionCore.calculate_true_positives_vec!(
        true_positives,
        tested_inc_vec,
        time_p.tlength,
        test_spec,
    )

    # Calculate false positives from noise
    OutbreakDetectionCore.calculate_false_positives_vec!(
        false_positives,
        tested_noise_vec,
        time_p.tlength,
        test_spec,
    )

    test_positives = true_positives .+ false_positives

    # Calculate moving average of TOTAL test positives
    movingavg_testpositives = OutbreakDetectionCore.calculate_movingavg(
        test_positives,
        outbreak_detection_spec.moving_average_lag,
    )

    # Detect alerts where moving average exceeds threshold
    alertstatus_vec = movingavg_testpositives .>= outbreak_detection_spec.alert_threshold
    alert_rle = StatsBase.rle(alertstatus_vec)
    alert_thresholds = OutbreakDetectionCore.calculate_above_threshold_bounds(alert_rle)

    # Create a matrix with [lower, upper, duration] for compatibility with plot function
    alert_bounds = hcat(
        alert_thresholds.lower_bounds,
        alert_thresholds.upper_bounds,
        alert_thresholds.duration
    )

    # Call the base plot_schematic function
    fig = plot_schematic(
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

    # Add description with noise type, noise level, and test characteristics
    noise_type = optimization_result.noise_type_description
    noise_lvl = optimization_result.noise_level
    sens = test_spec.sensitivity
    spec = test_spec.specificity
    lag = test_spec.test_result_lag

    description_text = "Noise: $(noise_type), Level: $(noise_lvl) | Test: Sens=$(sens), Spec=$(spec), Lag=$(lag)d"

    Label(
        fig[0, 1],
        description_text;
        fontsize = 12,
        halign = :center,
        tellwidth = false
    )

    return fig
end
