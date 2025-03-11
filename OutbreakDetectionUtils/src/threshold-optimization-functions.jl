using StatsBase: StatsBase
using NaNMath: NaNMath
using UnPack: @unpack
using Optim: Optim
using DataFrames: DataFrames

function setup_optimization(ensemble_param_dict)
    UnPack.@unpack ensemble_spec,
    seed,
    outbreak_spec = ensemble_param_dict

    UnPack.@unpack state_parameters,
    dynamics_parameters, time_parameters,
    nsims =
        ensemble_spec

    UnPack.@unpack tstep, tlength, trange = time_parameters

    ensemble_seir_vecs = Array{typeof(state_parameters.init_states),2}(
        undef,
        tlength,
        nsims,
    )

    ensemble_inc_vecs = Array{typeof(StaticArrays.SVector(0)),2}(
        undef,
        tlength,
        nsims,
    )

    ensemble_beta_arr = zeros(Float64, tlength)

    for sim in axes(ensemble_inc_vecs, 2)
        run_seed = seed + (sim - 1)

        seir_mod!(
            @view(ensemble_seir_vecs[:, sim]),
            @view(ensemble_inc_vecs[:, sim]),
            ensemble_beta_arr,
            state_parameters.init_states,
            dynamics_parameters,
            time_parameters;
            seed = run_seed,
        )
    end

    ensemble_inc_arr, ensemble_thresholds_vec = create_inc_infec_arr(
        ensemble_inc_vecs, outbreak_spec
    )

    return ensemble_inc_arr, ensemble_thresholds_vec
end

function run_optimization(
    OT_chars_param_dict;
    guess_threshold = 2.0,
    threshold_lower_bound = 0.5,
    threshold_upper_bound = 50.0,
    optim_method = Optim.Brent(),
    optim_options = Optim.Options(),
    maxiters = 100000,
    maxtime = 1000.0,
)
    UnPack.@unpack scenario_spec, ensemble_inc_arr, thresholds_vec, seed =
        OT_chars_param_dict
    UnPack.@unpack noise_specification,
    outbreak_detection_specification,
    individual_test_specification = scenario_spec

    noise_array = create_noise_arr(
        noise_specification,
        ensemble_inc_arr;
        ensemble_specification = scenario_spec.ensemble_specification,
        seed = seed,
    )[1]

    obj_inputs = (;
        ensemble_inc_arr,
        noise_array,
        outbreak_detection_specification,
        individual_test_specification,
        thresholds_vec,
    )

    optim_sol = Optim.optimize(
        t -> objective_function(t, obj_inputs),
        guess_threshold,
        threshold_lower_bound,
        threshold_upper_bound,
        optim_method,
    )

    return sol
end

function objective_function(
    alert_threshold_vec,
    inputs,
)
    @assert length(alert_threshold_vec) == 1

    @unpack ensemble_inc_arr,
    noise_array,
    outbreak_detection_specification,
    individual_test_specification,
    thresholds_vec,
    output_df = inputs

    outbreak_detection_specification = OutbreakDetectionSpecification(
        alert_threshold_vec[1],
        outbreak_detection_specification.moving_average_lag,
        outbreak_detection_specification.percent_visit_clinic,
        outbreak_detection_specification.percent_clinic_tested,
        outbreak_detection_specification.alert_method.method_name,
    )

    testarr = create_testing_arrs(
        ensemble_inc_arr,
        noise_array,
        outbreak_detection_specification,
        individual_test_specification,
    )[1]

    objective = calculate_ensemble_objective_metric(
        testarr, ensemble_inc_arr, thresholds_vec
    )

    # DataFrames.push!(
    # 	output_df,
    # 	("threshold" = alert_threshold_vec[1], "loss" = objective),
    # 	cols = :union
    # )
    #
    return objective
end

function calculate_ensemble_objective_metric(
    testarr, infecarr, thresholds_vec
)
    mean_accuracy = map(axes(infecarr, 3)) do sim
        dailychars = calculate_daily_detection_characteristics(
            @view(testarr[:, 6, sim]), @view(infecarr[:, 3, sim])
        )
        alertrle = StatsBase.rle(@view(testarr[:, 6, sim]))
        outbreakbounds = thresholds_vec[sim]
        alertbounds = calculate_outbreak_thresholds(alertrle; ncols = 3)
        # calculate_outbreak_duration!(alertbounds)

        accuracy = calculate_outbreak_detection_accuracy(
            outbreakbounds, alertbounds
        )
    end

    return 1 - NaNMath.mean(mean_accuracy)
end

function calculate_outbreak_detection_accuracy(
    outbreakbounds, alertbounds
)
    filtered_matched_bounds = match_outbreak_detection_bounds(
        outbreakbounds, alertbounds
    )[1]

    noutbreaks = size(outbreakbounds, 1)
    nalerts = size(alertbounds, 1)

    n_true_outbreaks_detected = length(
        Set(@view(filtered_matched_bounds[:, 1]))
    )
    n_correct_alerts = size(filtered_matched_bounds, 1)

    perc_true_outbreaks_detected = n_true_outbreaks_detected / noutbreaks
    perc_alerts_correct = n_correct_alerts / nalerts # c.f. PPV

    accuracy = NaNMath.mean([perc_true_outbreaks_detected, perc_alerts_correct])

    return accuracy
end
