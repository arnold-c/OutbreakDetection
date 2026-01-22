export NoiseVaccinationOptimizationParameters,
    ThresholdOptimizationParameters

"""
    NoiseVaccinationOptimizationParameters

Parameters for optimizing vaccination coverage in dynamical noise simulations.

This type configures the multistart optimization algorithm used to find the
vaccination coverage that produces a target noise level. It uses Sobol
sequences for global search initialization and NLopt for local optimization.

# Fields
- `n_sobol_points::Int64`: Number of Sobol sequence points for multistart initialization (default: 100)
- `local_algorithm`: NLopt algorithm for local optimization (default: NLopt.LN_BOBYQA)
- `maxeval::Int64`: Maximum function evaluations per local optimization (default: 1000)
- `xtol_rel::Float64`: Relative tolerance on parameters (default: 1.0e-3)
- `xtol_abs::Float64`: Absolute tolerance on parameters (default: 1.0e-3)
- `atol::Float64`: Absolute difference tolerance threshold for convergence (default: 1.0e-2)

# Algorithm Details
The optimization uses a two-stage approach:
1. **Global search**: Sobol sequence generates `n_sobol_points` initial points
   uniformly distributed across the vaccination bounds
2. **Local refinement**: Each Sobol point is used as a starting point for
   local optimization with the specified NLopt algorithm

# Convergence Criteria
Optimization is considered successful if:
- NLopt returns a success code (SUCCESS, XTOL_REACHED, FTOL_REACHED, or STOPVAL_REACHED)
- Absolute difference between achieved and target noise is less than `atol`

# Examples
```julia
# Default parameters (good for most cases)
opt_params = NoiseVaccinationOptimizationParameters()

# Fast optimization for testing
opt_params = NoiseVaccinationOptimizationParameters(
    n_sobol_points = 10,
    maxeval = 100,
    atol = 0.5
)

# High-precision optimization
opt_params = NoiseVaccinationOptimizationParameters(
    n_sobol_points = 200,
    maxeval = 2000,
    xtol_rel = 1.0e-4,
    xtol_abs = 1.0e-4,
    atol = 1.0e-3
)

# Use with optimization
result = optimize_dynamic_noise_params(
    target_scaling,
    mean_target_incidence,
    noise_spec,
    ensemble_spec,
    base_dynamics,
    opt_params
)
```

# Performance Notes
- More Sobol points increase global search coverage but slow optimization
- Tighter tolerances increase precision but require more evaluations
- BOBYQA (default) is derivative-free and works well for noisy objectives
- Consider relaxing tolerances for initial exploration

# See Also
- [`optimize_dynamic_noise_params`](@ref): Main optimization function
- [`DynamicalNoiseParameters`](@ref): Noise specification with bounds
"""
Base.@kwdef struct NoiseVaccinationOptimizationParameters
    n_sobol_points::Int64 = 100
    local_algorithm = NLopt.LN_BOBYQA
    maxeval::Int64 = 1000
    xtol_rel::Float64 = 1.0e-3
    xtol_abs::Float64 = 1.0e-3
    atol::Float64 = 1.0e-2
end

"""
    ThresholdOptimizationParameters

Parameters for optimizing detection thresholds in outbreak detection algorithms.

This type configures the multistart optimization algorithm used to find the
optimal threshold that minimizes a loss function (e.g., balancing false positives
and false negatives). It uses Sobol sequences for global search initialization
and NLopt for local optimization.

# Fields
- `n_sobol_points::Int64`: Number of Sobol sequence points for multistart initialization (default: 100)
- `local_algorithm`: NLopt algorithm for local optimization (default: NLopt.LN_BOBYQA)
- `maxeval::Int64`: Maximum function evaluations per local optimization (default: 1000)
- `xtol_rel::Float64`: Relative tolerance on parameters (default: 1.0e-3)
- `xtol_abs::Float64`: Absolute tolerance on parameters (default: 1.0e-3)
- `atol::Float64`: Absolute difference tolerance threshold for convergence (default: 1.0e-2)

# Algorithm Details
The optimization uses a two-stage approach:
1. **Global search**: Sobol sequence generates `n_sobol_points` initial points
   uniformly distributed across the threshold bounds
2. **Local refinement**: Each Sobol point is used as a starting point for
   local optimization with the specified NLopt algorithm

# Convergence Criteria
Optimization is considered successful if:
- NLopt returns a success code (SUCCESS, XTOL_REACHED, FTOL_REACHED, or STOPVAL_REACHED)
- Absolute difference between achieved and target loss is less than `atol`

# Examples
```julia
# Default parameters (good for most cases)
opt_params = ThresholdOptimizationParameters()

# Fast optimization for testing
opt_params = ThresholdOptimizationParameters(
    n_sobol_points = 10,
    maxeval = 100,
    atol = 0.5
)

# High-precision optimization
opt_params = ThresholdOptimizationParameters(
    n_sobol_points = 200,
    maxeval = 2000,
    xtol_rel = 1.0e-4,
    xtol_abs = 1.0e-4,
    atol = 1.0e-3
)
```

# Performance Notes
- More Sobol points increase global search coverage but slow optimization
- Tighter tolerances increase precision but require more evaluations
- BOBYQA (default) is derivative-free and works well for noisy objectives
- Consider relaxing tolerances for initial exploration

# See Also
- [`OptimizationTracker`](@ref): Tracks best solution during optimization
- [`NoiseVaccinationOptimizationParameters`](@ref): Similar parameters for noise optimization
"""
Base.@kwdef struct ThresholdOptimizationParameters
    n_sobol_points::Int64 = 100
    local_algorithm = NLopt.LN_BOBYQA
    maxeval::Int64 = 1000
    xtol_rel::Float64 = 1.0e-3
    xtol_abs::Float64 = 1.0e-3
    atol::Float64 = 1.0e-2
end
