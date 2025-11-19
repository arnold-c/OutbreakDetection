module OutbreakDetectionUtils

# Utilities
include("./utilities/DrWatson-helpers.jl")

# Types
include("./types/time-parameters.jl")
include("./types/state-parameters.jl")
include("./types/dynamics-parameters.jl")
include("./types/ensemble-specifications.jl")
include("./types/outbreak-specifications.jl")
include("./types/detection-specifications.jl")
include("./types/test-specifications.jl")
include("./types/noise-specifications.jl")
include("./types/scenario-specifications.jl")
include("./types/optimization-types.jl")

# Simulation (needed before constants)
include("./simulation/transmission-functions.jl")

# Constants (depend on transmission functions)
include("./constants/dynamics-constants.jl")
include("./constants/test-constants.jl")

# Simulation (continued)
include("./simulation/seir-model.jl")
include("./simulation/ensemble-simulation.jl")
include("./simulation/ensemble-outbreak-detection.jl")

# Noise
include("./noise/noise-generation.jl")

# Utilities (data processing)
include("./utilities/cleaning-functions.jl")
include("./utilities/collect-thresholds-vec_functions.jl")
include("./utilities/create-combinations.jl")

# Detection
include("./detection/outbreak-thresholds.jl")
include("./detection/outbreak-classification.jl")
include("./detection/moving-average.jl")
include("./detection/alert-detection.jl")
include("./detection/detection-characteristics.jl")

# Diagnostic Testing
include("./diagnostic-testing/calculate-num-tested.jl")
include("./diagnostic-testing/calculate-num-positive.jl")
include("./diagnostic-testing/calculate-test-positivity.jl")
include("./diagnostic-testing/create-test-arrays.jl")

include("optimal-threshold-functions.jl")
export calculate_optimal_threshold, calculate_OptimalThresholdCharacteristics,
    calculate_optimal_threshold_summaries,
    create_optimal_thresholds_df, create_wide_optimal_thresholds_df,
    create_and_save_xlsx_optimal_threshold_summaries,
    create_optimal_threshold_summary_df,
    create_wide_optimal_threshold_summary_df,
    create_all_wide_optimal_threshold_summary_dfs,
    save_xlsx_optimal_threshold_summaries,
    create_and_save_xlsx_optimal_threshold_summaries

include("threshold-optimization-functions.jl")
export run_optimization,
    setup_optimization,
    objective_function,
    calculate_ensemble_objective_metric,
    calculate_outbreak_detection_accuracy,
    optimization_wrapper

include("scenario-optimizations.jl")
export run_scenario_optimizations,
    check_missing_scenario_optimizations,
    run_missing_scenario_optimizations!,
    get_most_recent_optimization_filepath,
    reshape_optim_df_to_matrix,
    sort_noise_specifications

end
