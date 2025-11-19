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

# Optimal Thresholds
include("./optimal-thresholds/threshold-calculation.jl")
include("./optimal-thresholds/threshold-summaries.jl")
include("./optimal-thresholds/results-dataframes.jl")
include("./optimal-thresholds/results-export.jl")

# Threshold Optimization
include("./optimization-functions/common/objective-function.jl")
include("./optimization-functions/threshold-optimization/optimization-setup.jl")
include("./optimization-functions/threshold-optimization/optimization-wrapper.jl")

# Scenario Optimization
include("./optimization-functions/scenario-optimization/scenario-creation.jl")
include("./optimization-functions/scenario-optimization/scenario-confirmation.jl")
include("./optimization-functions/scenario-optimization/results-loading.jl")
include("./optimization-functions/scenario-optimization/results-reshaping.jl")
include("./optimization-functions/scenario-optimization/optimization-wrapper.jl")

end
