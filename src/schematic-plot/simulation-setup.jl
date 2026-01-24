export get_outbreak_status,
    shift_vec,
    create_schematic_simulation

"""
    get_outbreak_status(inc_vec, outbreak_specification)

Identifies outbreak periods in an incidence time series based on specified
outbreak criteria.

Determines which time periods constitute outbreaks by applying threshold-based
detection, run-length encoding, and classification based on minimum duration
and size requirements. Returns both a binary outbreak status vector and the
temporal bounds of detected outbreaks.

# Arguments
- `inc_vec`: Vector of incidence values over time
- `outbreak_specification`: OutbreakSpecification containing outbreak_threshold,
  minimum_outbreak_duration, and minimum_outbreak_size

# Returns
- `outbreak_status`: Binary vector (true/false) indicating outbreak status at each
  time point
- `outbreak_bounds`: Matrix with rows representing [start, end] indices of each
  detected outbreak
"""
function get_outbreak_status(
        inc_vec,
        outbreak_specification
    )
    abovethreshold_vec = vec(
        inc_vec .>= outbreak_specification.outbreak_threshold
    )

    abovethresholdrle = StatsBase.rle(abovethreshold_vec)

    threshold_bounds = OutbreakDetectionCore.calculate_above_threshold_bounds(
        abovethresholdrle
    )

    outbreak_thresholds = OutbreakDetectionCore.classify_all_outbreaks(
        inc_vec,
        threshold_bounds,
        outbreak_specification.minimum_outbreak_duration,
        outbreak_specification.minimum_outbreak_size,
    )

    outbreak_status = zeros(Bool, length(inc_vec))

    for i in eachindex(outbreak_thresholds.lower_bounds)
        lower = outbreak_thresholds.lower_bounds[i]
        upper = outbreak_thresholds.upper_bounds[i]
        outbreak_status[lower:upper] .= true
    end

    # Create a matrix with [lower, upper] for compatibility with plot function
    outbreak_bounds = hcat(
        outbreak_thresholds.lower_bounds,
        outbreak_thresholds.upper_bounds
    )

    return outbreak_status, outbreak_bounds
end

"""
    shift_vec(invec, shift::T) where {T <: Integer}

Shifts a vector by a specified number of positions, padding with zeros.

Positive shift values move elements forward in time (right), while negative
values move elements backward in time (left). The resulting vector maintains
the same length as the input, with zeros filling the vacated positions.

# Arguments
- `invec`: Input vector to be shifted
- `shift`: Integer number of positions to shift. Positive values shift right
  (forward in time), negative values shift left (backward in time), zero
  returns the original vector

# Returns
- `outvec`: Shifted vector of the same length and element type as `invec`,
  with zeros in positions that have no corresponding input value
"""
function shift_vec(invec, shift::T) where {T <: Integer}
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

"""
    create_schematic_simulation(states_p, dynamics_p, dynamics_spec,
        noise_spec, test_specification, time_p; seed = 1234,
        outbreak_specification = OutbreakSpecification(5, 30, 500),
        outbreak_detection_specification = OutbreakDetectionSpecification(...),
        noise_scaling = 10, shift_noise = 0, movingavg_window = 20)

Creates a complete outbreak detection simulation for schematic visualization.

Generates synthetic epidemic data by running SEIR models for both true
infections and noise, applies diagnostic testing with specified sensitivity and
specificity, calculates moving averages, and performs outbreak detection. This
function is designed to produce all necessary data for creating schematic plots
that illustrate the outbreak detection process.

# Arguments
- `states_p`: StateParameters containing initial SEIR compartment states for
  the main epidemic
- `dynamics_p`: DynamicsParameters for the main epidemic transmission dynamics
- `dynamics_spec`: DynamicsParameterSpecification for the main epidemic
- `noise_spec`: DynamicalNoiseSpecification for the noise process
- `test_specification`: IndividualTestSpecification containing test
  characteristics (sensitivity, specificity, result lag)
- `time_p`: Time parameters specifying simulation duration and time step

# Keyword Arguments
- `seed`: Random seed for reproducibility (default: 1234)
- `outbreak_specification`: OutbreakSpecification defining outbreak criteria
  (default: threshold=5, duration=30, size=500)
- `outbreak_detection_specification`: OutbreakDetectionSpecification defining
  detection parameters
- `noise_scaling`: Multiplicative factor for noise intensity (default: 10)
- `shift_noise`: Number of time steps to shift noise signal (default: 0)
- `movingavg_window`: Window size for smoothing incidence data (default: 20)

# Returns
A tuple containing:
- `inc_vec`: Smoothed incidence vector from the main epidemic
- `outbreak_status`: Binary vector indicating true outbreak periods
- `outbreak_bounds`: Matrix of [start, end] indices for each outbreak
- `noise_vec`: Shifted and scaled noise incidence vector
- `movingavg_testpositives`: Moving average of total test positives
- `alertstatus_vec`: Binary vector indicating when alerts are triggered
- `alert_bounds`: Matrix of [start, end, duration] for each alert period
"""
function create_schematic_simulation(
        states_p,
        dynamics_p,
        dynamics_spec,
        noise_spec,
        test_specification,
        time_p;
        seed = 1234,
        outbreak_specification = OutbreakDetectionCore.OutbreakSpecification(5, 30, 500),
        outbreak_detection_specification = nothing,
        noise_scaling = 10,
        shift_noise = 0,
        movingavg_window = 20
    )
    # Run main epidemic simulation
    seir_run = OutbreakDetectionCore.seir_mod(
        StaticArrays.SVector(states_p.init_states),
        dynamics_p,
        time_p;
        seed = seed,
    )
    inc_vec = seir_run.incidence

    inc_vec =
        round.(
        OutbreakDetectionCore.calculate_movingavg(
            inc_vec,
            movingavg_window,
        )
    )

    outbreak_status, outbreak_bounds = get_outbreak_status(
        inc_vec, outbreak_specification
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
    ensemble_spec = OutbreakDetectionCore.EnsembleSpecification(
        "schematic",
        states_p,
        time_p,
        dynamics_spec,
        noise_params,
        1,  # nsims
        ""  # dirpath
    )

    # Recreate noise dynamics
    noise_dynamics_spec = OutbreakDetectionCore.recreate_noise_dynamics_parameter_specification(
        noise_spec,
        ensemble_spec
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
        ensemble_spec,
        noise_dynamics_p;
        seed = seed
    )

    # Extract noise incidence (first simulation)
    noise_vec_raw = LightSumTypes.variant(noise_result).incidence[1]

    noise_vec_tmp = OutbreakDetectionCore.calculate_movingavg(
        noise_vec_raw .* noise_scaling,
        movingavg_window
    )

    noise_vec = shift_vec(noise_vec_tmp, shift_noise)

    perc_tested =
        outbreak_detection_specification.percent_visit_clinic *
        outbreak_detection_specification.percent_clinic_tested

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
        test_specification,
    )

    # Calculate false positives from noise
    OutbreakDetectionCore.calculate_false_positives_vec!(
        false_positives,
        tested_noise_vec,
        time_p.tlength,
        test_specification,
    )

    test_positives = true_positives .+ false_positives

    # Calculate moving average of TOTAL test positives
    movingavg_testpositives = OutbreakDetectionCore.calculate_movingavg(
        test_positives,
        outbreak_detection_specification.moving_average_lag,
    )

    # Detect alerts where moving average exceeds threshold
    alertstatus_vec = movingavg_testpositives .>= outbreak_detection_specification.alert_threshold
    alert_rle = StatsBase.rle(alertstatus_vec)
    alert_thresholds = OutbreakDetectionCore.calculate_above_threshold_bounds(alert_rle)

    # Create a matrix with [lower, upper, duration] for compatibility with plot function
    alert_bounds = hcat(
        alert_thresholds.lower_bounds,
        alert_thresholds.upper_bounds,
        alert_thresholds.duration
    )

    return inc_vec,
        outbreak_status,
        outbreak_bounds,
        noise_vec,
        movingavg_testpositives,
        alertstatus_vec,
        alert_bounds
end
