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

    all_outbreak_bounds = OutbreakDetectionCore.calculate_outbreak_thresholds(
        abovethresholdrle
    )

    OutbreakDetectionCore.classify_all_outbreaks!(
        abovethreshold_vec,
        all_outbreak_bounds,
        inc_vec,
        outbreak_specification.minimum_outbreak_duration,
        outbreak_specification.minimum_outbreak_size,
    )

    outbreak_bounds = OutbreakDetectionCore.filter_only_outbreaks(
        all_outbreak_bounds
    )

    outbreak_status = zeros(Bool, length(inc_vec))

    for (lower, upper) in eachrow(outbreak_bounds)
        # @show lower, upper
        outbreak_status[lower:upper] .= true
    end
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
    create_schematic_simulation(states_p, dynamics_p, noise_states_p,
        noise_dynamics_p, test_specification, time_p; seed = 1234,
        outbreak_specification = OutbreakSpecification(5, 30, 500),
        outbreak_detection_specification = OutbreakDetectionSpecification(
            5, 7, 1.0, 0.5, "movingavg"),
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
- `noise_states_p`: StateParameters for the noise/background process
- `noise_dynamics_p`: DynamicsParameters for the noise process
- `test_specification`: IndividualTestSpecification containing test
  characteristics (sensitivity, specificity, result lag)
- `time_p`: Time parameters specifying simulation duration and time step

# Keyword Arguments
- `seed`: Random seed for reproducibility (default: 1234)
- `outbreak_specification`: OutbreakSpecification defining outbreak criteria
  (default: threshold=5, duration=30, size=500)
- `outbreak_detection_specification`: OutbreakDetectionSpecification defining
  detection parameters (default: threshold=5, lag=7, clinic_visit=1.0,
  clinic_tested=0.5, method="movingavg")
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
        noise_states_p,
        noise_dynamics_p,
        test_specification,
        time_p;
        seed = 1234,
        outbreak_specification = OutbreakDetectionCore.OutbreakSpecification(5, 30, 500),
        outbreak_detection_specification = OutbreakDetectionCore.OutbreakDetectionSpecification(5, 7, 1.0, 0.5, "movingavg"),
        noise_scaling = 10,
        shift_noise = 0,
        movingavg_window = 20
    )
    inc_sv = OutbreakDetectionCore.seir_mod(
        states_p.init_states,
        dynamics_p,
        time_p;
        seed = seed,
    )[2]

    inc_vec =
        round.(
        OutbreakDetectionCore.calculate_movingavg(
            vec(convert_svec_to_matrix(inc_sv)),
            movingavg_window,
        )
    )

    outbreak_status, outbreak_bounds = get_outbreak_status(
        inc_vec, outbreak_specification
    )

    noise_sv = OutbreakDetectionCore.seir_mod(
        noise_states_p.init_states,
        noise_dynamics_p,
        time_p;
        seed = seed,
    )[2]

    noise_vec_tmp = OutbreakDetectionCore.calculate_movingavg(
        vec(OutbreakDetectionCore.convert_svec_to_matrix(noise_sv)) .*
            noise_scaling, movingavg_window
    )

    noise_vec = shift_vec(noise_vec_tmp, shift_noise)

    perc_tested =
        outbreak_detection_specification.percent_visit_clinic *
        outbreak_detection_specification.percent_clinic_tested

    true_positives = OutbreakDetectionCore.calculate_positives(
        calculate_true_positives!,
        inc_vec .* perc_tested,
        time_p.tlength,
        test_specification.test_result_lag,
        test_specification.sensitivity,
    )

    false_positives = OutbreakDetectionCore.calculate_positives(
        calculate_noise_positives!,
        noise_vec .* perc_tested,
        time_p.tlength,
        test_specification.test_result_lag,
        test_specification.specificity,
    )

    test_positives = true_positives .+ false_positives

    # Calculate moving average of TOTAL test positives
    movingavg_testpositives = OutbreakDetectionCore.calculate_movingavg(
        test_positives,
        outbreak_detection_specification.moving_average_lag,
    )

    alertstatus_vec = OutbreakDetectionCore.detectoutbreak(
        movingavg_testpositives,
        outbreak_detection_specification.alert_threshold,
    )
    alert_bounds = OutbreakDetectionCore.calculate_outbreak_thresholds(
        rle(alertstatus_vec); ncols = 3
    )
    OutbreakDetectionCore.calculate_outbreak_duration!(alert_bounds)

    return inc_vec,
        outbreak_status,
        outbreak_bounds,
        noise_vec,
        movingavg_testpositives,
        alertstatus_vec,
        alert_bounds
end
