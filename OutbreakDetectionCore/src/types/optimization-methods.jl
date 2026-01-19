export OptimizationMethods,
    MSO,
    AbstractOptimizationMethods

abstract type AbstractOptimizationMethods end

"""
    MSO

Multistart optimization method.
"""
struct MSO end

"""
    OptimizationMethods

Sum type for available optimization methods.
"""
LightSumTypes.@sumtype OptimizationMethods(MSO) <: AbstractOptimizationMethods
