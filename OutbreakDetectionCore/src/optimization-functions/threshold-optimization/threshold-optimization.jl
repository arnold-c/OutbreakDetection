export threshold_optimization

function threshold_optimization(
        scenario::OptimizationScenario,
        test_positives_container::TestPositiveContainer,
        opt_params::ThresholdOptimizationParameters,
    )
    # Create tracker instance for this scenario
    tracker = OptimizationTracker()

    # Create objective function closure that updates tracker
    objective = params -> multistart_objective_function(
        params,
        scenario,
        tracker
    )

    # Setup multistart problem
    problem = MultistartOptimization.MinimizationProblem(
        objective,
        scenario.threshold_bounds.lower,
        scenario.threshold_bounds.upper
    )

    # Configure local optimization method
    local_method = MultistartOptimization.NLopt_local_method(
        local_algorithm;
        xtol_rel = opt_params.xtol_rel,
        xtol_abs = opt_params.xtol_abs,
        maxeval = opt_params.maxeval,
    )

    # Configure multistart method (TikTak uses Sobol sequences)
    multistart_method = MultistartOptimization.TikTak(opt_params.n_sobol_points)

    # Run optimization
    MultistartOptimization.multistart_minimization(
        multistart_method,
        local_method,
        problem
    )
    return OptimizationResult(
        # From scenario
        ensemble_specification = scenario.ensemble_specification,
        noise_level = scenario.noise_level,
        noise_type_description = scenario.noise_type_description,
        test_specification = scenario.test_specification,
        percent_tested = scenario.percent_tested,
        alert_method = scenario.alert_method,
        accuracy_metric = scenario.accuracy_metric,
        # From optimized values
        optimal_threshold = tracker.optimal_threshold,
        accuracy = tracker.best_accuracy,
        proportion_outbreaks_detected = tracker.proportion_outbreaks_detected,
        proportion_alerts_correct = tracker.proportion_alerts_correct,
        mean_detection_delay = tracker.mean_detection_delay
    )
end
