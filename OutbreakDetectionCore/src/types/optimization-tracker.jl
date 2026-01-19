export OptimizationTracker

"""
    OptimizationTracker

Mutable struct to track the best solution and its metrics during optimization.
"""
Base.@kwdef mutable struct OptimizationTracker
    best_loss::Float64 = Inf
    optimal_threshold::Float64 = 0.0
    best_accuracy::Float64 = 0.0
    proportion_outbreaks_detected::Float64 = 0.0
    proportion_alerts_correct::Float64 = 0.0
    mean_detection_delay::Float64 = Inf
end
