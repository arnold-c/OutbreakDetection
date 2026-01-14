export AccuracyMetric,
    F1,
    BalancedAccuracy

abstract type AbstractAccuracyMetric end

struct F1 <: AbstractAccuracyMetric end
struct BalancedAccuracy <: AbstractAccuracyMetric end

LightSumTypes.@sumtype AccuracyMetric(F1, BalancedAccuracy) <: AbstractAccuracyMetric
