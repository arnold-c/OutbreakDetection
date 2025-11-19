# Optimization Wrapper Implementation Plan

## Overview

This plan details the implementation of a CSDNoise-style optimization wrapper that uses StructVectors and AutoHashEquals for efficient parameter space exploration and result storage. This builds on the refactoring completed in REFACTORING_PLAN.md.

## Goals

1. Implement StructVector-based result storage for optimization runs
2. Add AutoHashEquals for efficient parameter comparison and caching
3. Create optimization wrapper functions for threshold optimization
4. Implement scenario-based optimization workflow
5. Add efficient result loading and reshaping utilities
6. Maintain compatibility with existing optimization code

## Key Concepts from CSDNoise

### StructVector Benefits
- **Column-wise storage**: More efficient memory layout for large ensembles
- **Type stability**: Better performance with concrete types
- **Easy filtering**: Can filter by parameter values efficiently
- **Serialization**: Better JLD2 performance

### AutoHashEquals Benefits
- **Automatic hashing**: No manual hash/equality implementation needed
- **Caching**: Efficient result lookup by parameters
- **Deduplication**: Avoid re-running identical parameter combinations

### Optimization Workflow
1. Define parameter space (scenarios)
2. Generate parameter combinations
3. Run optimization for each combination
4. Store results in StructVector
5. Load and reshape results for analysis

---

## Phase 1: Add Dependencies and Core Types

### 1.1 Update Project.toml

**Priority: CRITICAL**
**File to modify:** `OutbreakDetectionCore/Project.toml`

**Changes:**

Add AutoHashEquals dependency:
```toml
[deps]
AutoHashEquals = "15f4f7f2-30c1-5605-9d31-71845cf9641f"

[compat]
AutoHashEquals = "2.1"
```

**Testing Requirements:**
- Verify package resolves correctly
- Test AutoHashEquals import

---

### 1.2 Create Scenario Specification Types

**Priority: HIGH**
**New file to create:** `OutbreakDetectionCore/src/types/scenario-parameters.jl`

**Purpose:**
Define parameter types that use AutoHashEquals for efficient comparison and caching.

**Changes:**

```julia
export ScenarioParameters, ThresholdOptimizationScenario

using AutoHashEquals

"""
    ScenarioParameters

Base parameters for a scenario with automatic hashing and equality.

Uses AutoHashEquals for efficient comparison and caching.

# Fields
- `target_dynamics::TargetDiseaseDynamicsParameters`: Target disease parameters
- `noise_spec::Union{DynamicalNoiseSpecification, PoissonNoise}`: Noise specification
- `test_spec::IndividualTestSpecification`: Testing specification
- `percent_clinic_tested::Float64`: Proportion tested at clinic (0-1)
- `outbreak_threshold::Float64`: Threshold for outbreak detection

# Examples
```julia
target = TargetDiseaseDynamicsParameters(
    R_0 = 16.0,
    latent_period_days = 10.0,
    infectious_duration_days = 8.0,
    beta_force = 0.2
)

noise = DynamicalNoiseSpecification(
    R_0 = 5.0,
    latent_period = 7.0,
    duration_infection = 14.0,
    correlation = "in-phase",
    poisson_component = 1.0
)

test = IndividualTestSpecification(
    sensitivity = 1.0,
    specificity = 0.8,
    test_result_lag = 0
)

scenario = ScenarioParameters(
    target_dynamics = target,
    noise_spec = noise,
    test_spec = test,
    percent_clinic_tested = 0.5,
    outbreak_threshold = 0.01
)

# Automatic hashing and equality
scenario1 == scenario2  # Uses AutoHashEquals
hash(scenario)  # Automatic hash function
```
"""
@auto_hash_equals struct ScenarioParameters
    target_dynamics::TargetDiseaseDynamicsParameters
    noise_spec::Union{DynamicalNoiseSpecification, PoissonNoise}
    test_spec::IndividualTestSpecification
    percent_clinic_tested::Float64
    outbreak_threshold::Float64
    
    function ScenarioParameters(
        target_dynamics, noise_spec, test_spec, 
        percent_clinic_tested, outbreak_threshold
    )
        @assert 0.0 <= percent_clinic_tested <= 1.0 "Percent tested must be in [0, 1]"
        @assert outbreak_threshold > 0.0 "Outbreak threshold must be positive"
        
        new(target_dynamics, noise_spec, test_spec, 
            percent_clinic_tested, outbreak_threshold)
    end
end

"""
    ThresholdOptimizationScenario

Complete scenario specification for threshold optimization.

Combines scenario parameters with ensemble and time specifications.

# Fields
- `scenario_params::ScenarioParameters`: Core scenario parameters
- `ensemble_spec::EnsembleSpecification`: Ensemble simulation parameters
- `time_params::SimTimeParameters`: Time parameters
- `nsims::Int64`: Number of simulations
- `seed::Int64`: Random seed for reproducibility

# Examples
```julia
scenario = ThresholdOptimizationScenario(
    scenario_params = scenario_params,
    ensemble_spec = ensemble_spec,
    time_params = time_params,
    nsims = 1000,
    seed = 1234
)
```
"""
@auto_hash_equals struct ThresholdOptimizationScenario
    scenario_params::ScenarioParameters
    ensemble_spec::EnsembleSpecification
    time_params::SimTimeParameters
    nsims::Int64
    seed::Int64
    
    function ThresholdOptimizationScenario(
        scenario_params, ensemble_spec, time_params, nsims, seed
    )
        @assert nsims > 0 "Number of simulations must be positive"
        @assert seed >= 0 "Seed must be non-negative"
        
        new(scenario_params, ensemble_spec, time_params, nsims, seed)
    end
end
```

**Testing Requirements:**
- Test AutoHashEquals functionality
- Test equality comparison
- Test hashing
- Test validation

---

### 1.3 Create Optimization Result Types

**Priority: HIGH**
**New file to create:** `OutbreakDetectionCore/src/types/optimization-results.jl`

**Purpose:**
Define result types that work efficiently with StructVectors.

**Changes:**

```julia
export OptimizationResult, ThresholdOptimizationResult

"""
    OptimizationResult

Base result from an optimization run.

Designed for StructVector storage with AutoHashEquals for efficient lookup.

# Fields
- `scenario_params::ScenarioParameters`: Scenario that was optimized
- `optimal_threshold::Float64`: Optimal detection threshold
- `accuracy::Float64`: Detection accuracy at optimal threshold
- `sensitivity::Float64`: True positive rate
- `specificity::Float64`: True negative rate
- `mean_detection_delay::Float64`: Mean delay to detection (days)
- `proportion_detected::Float64`: Proportion of outbreaks detected

# Examples
```julia
result = OptimizationResult(
    scenario_params = scenario,
    optimal_threshold = 0.05,
    accuracy = 0.95,
    sensitivity = 0.92,
    specificity = 0.98,
    mean_detection_delay = 14.5,
    proportion_detected = 0.88
)

# Store in StructVector
results = StructVector{OptimizationResult}([result1, result2, ...])

# Efficient filtering
high_accuracy = filter(r -> r.accuracy > 0.9, results)
```
"""
@auto_hash_equals struct OptimizationResult
    scenario_params::ScenarioParameters
    optimal_threshold::Float64
    accuracy::Float64
    sensitivity::Float64
    specificity::Float64
    mean_detection_delay::Float64
    proportion_detected::Float64
end

"""
    ThresholdOptimizationResult

Extended result with detailed threshold characteristics.

Includes full threshold sweep results for analysis.

# Fields
- `base_result::OptimizationResult`: Base optimization result
- `threshold_sweep::Vector{Float64}`: Thresholds tested
- `accuracy_sweep::Vector{Float64}`: Accuracy at each threshold
- `sensitivity_sweep::Vector{Float64}`: Sensitivity at each threshold
- `specificity_sweep::Vector{Float64}`: Specificity at each threshold
"""
@auto_hash_equals struct ThresholdOptimizationResult
    base_result::OptimizationResult
    threshold_sweep::Vector{Float64}
    accuracy_sweep::Vector{Float64}
    sensitivity_sweep::Vector{Float64}
    specificity_sweep::Vector{Float64}
end
```

**Testing Requirements:**
- Test struct construction
- Test StructVector compatibility
- Test AutoHashEquals functionality
- Test filtering and querying

---

## Phase 2: Optimization Wrapper Functions

### 2.1 Implement Scenario Creation

**Priority: HIGH**
**New file to create:** `OutbreakDetectionCore/src/optimization-functions/scenario-optimization/scenario-creation.jl`

**Purpose:**
Create functions to generate parameter combinations for optimization.

**Changes:**

```julia
export create_scenario_grid, create_scenario_combinations

"""
    create_scenario_grid(; kwargs...)

Create a grid of scenario parameters for optimization.

Generates all combinations of provided parameter vectors.

# Keyword Arguments
- `R_0_values::Vector{Float64}`: R_0 values to test
- `noise_levels::Vector{String}`: Noise levels ("low", "medium", "high")
- `test_sensitivities::Vector{Float64}`: Test sensitivities to test
- `percent_tested_values::Vector{Float64}`: Proportions tested to test
- `outbreak_thresholds::Vector{Float64}`: Outbreak thresholds to test

# Returns
- `Vector{ScenarioParameters}`: All parameter combinations

# Examples
```julia
scenarios = create_scenario_grid(
    R_0_values = [12.0, 16.0, 20.0],
    noise_levels = ["low", "medium", "high"],
    test_sensitivities = [0.8, 0.9, 1.0],
    percent_tested_values = [0.3, 0.5, 0.7],
    outbreak_thresholds = [0.01, 0.02, 0.05]
)

# Results in 3 × 3 × 3 × 3 × 3 = 243 scenarios
length(scenarios)  # 243
```
"""
function create_scenario_grid(;
    R_0_values::Vector{Float64},
    noise_levels::Vector{String},
    test_sensitivities::Vector{Float64},
    percent_tested_values::Vector{Float64},
    outbreak_thresholds::Vector{Float64},
    base_target::TargetDiseaseDynamicsParameters,
    base_noise::DynamicalNoiseSpecification,
    base_test::IndividualTestSpecification
)
    scenarios = ScenarioParameters[]
    
    for R_0 in R_0_values
        target = TargetDiseaseDynamicsParameters(
            R_0 = R_0,
            latent_period_days = base_target.latent_period_days,
            infectious_duration_days = base_target.infectious_duration_days,
            beta_force = base_target.beta_force,
            seasonality = base_target.seasonality,
            life_expectancy_years = base_target.life_expectancy_years,
            population_N = base_target.population_N
        )
        
        for noise_level in noise_levels
            # Map noise level to parameters
            noise = _create_noise_for_level(noise_level, base_noise)
            
            for sensitivity in test_sensitivities
                test = IndividualTestSpecification(
                    sensitivity = sensitivity,
                    specificity = base_test.specificity,
                    test_result_lag = base_test.test_result_lag
                )
                
                for percent_tested in percent_tested_values
                    for threshold in outbreak_thresholds
                        scenario = ScenarioParameters(
                            target_dynamics = target,
                            noise_spec = noise,
                            test_spec = test,
                            percent_clinic_tested = percent_tested,
                            outbreak_threshold = threshold
                        )
                        push!(scenarios, scenario)
                    end
                end
            end
        end
    end
    
    return scenarios
end

"""
    _create_noise_for_level(level, base_noise)

Internal function to create noise specification for a given level.

# Arguments
- `level::String`: Noise level ("low", "medium", "high")
- `base_noise::DynamicalNoiseSpecification`: Base noise specification

# Returns
- `DynamicalNoiseSpecification`: Noise specification for the level
"""
function _create_noise_for_level(
    level::String, 
    base_noise::DynamicalNoiseSpecification
)
    # Map noise levels to vaccination coverage ranges
    # These will be optimized to match target noise levels
    vaccination_bounds = if level == "low"
        [0.0, 0.5]  # Lower vaccination = higher noise
    elseif level == "medium"
        [0.3, 0.7]
    elseif level == "high"
        [0.5, 0.9]  # Higher vaccination = lower noise
    else
        error("Unknown noise level: $level")
    end
    
    return DynamicalNoiseSpecification(
        R_0 = base_noise.R_0,
        latent_period = base_noise.latent_period,
        duration_infection = base_noise.duration_infection,
        correlation = base_noise.correlation,
        poisson_component = base_noise.poisson_component,
        vaccination_bounds = vaccination_bounds
    )
end
```

**Testing Requirements:**
- Test grid creation
- Test parameter combinations
- Test noise level mapping
- Verify no duplicate scenarios

---

### 2.2 Implement Optimization Wrapper

**Priority: HIGH**
**New file to create:** `OutbreakDetectionCore/src/optimization-functions/scenario-optimization/optimization-wrapper.jl`

**Purpose:**
Main wrapper function for running optimization across scenarios.

**Changes:**

```julia
export run_scenario_optimization, run_scenario_optimization_parallel

"""
    run_scenario_optimization(scenarios; kwargs...)

Run threshold optimization for multiple scenarios.

Uses StructVector for efficient result storage and AutoHashEquals for
caching to avoid re-running identical scenarios.

# Arguments
- `scenarios::Vector{ScenarioParameters}`: Scenarios to optimize

# Keyword Arguments
- `ensemble_spec::EnsembleSpecification`: Ensemble parameters
- `time_params::SimTimeParameters`: Time parameters
- `nsims::Int64`: Number of simulations per scenario
- `seed::Int64`: Base random seed
- `cache_results::Bool`: Use result caching (default: true)
- `save_results::Bool`: Save results to disk (default: true)
- `results_dir::String`: Directory for results (default: "results/optimization")
- `verbose::Bool`: Print progress (default: true)

# Returns
- `StructVector{OptimizationResult}`: Optimization results

# Examples
```julia
scenarios = create_scenario_grid(
    R_0_values = [16.0],
    noise_levels = ["low", "medium", "high"],
    test_sensitivities = [0.9, 1.0],
    percent_tested_values = [0.5],
    outbreak_thresholds = [0.01, 0.02]
)

results = run_scenario_optimization(
    scenarios;
    ensemble_spec = ensemble_spec,
    time_params = time_params,
    nsims = 1000,
    seed = 1234,
    verbose = true
)

# Results stored in StructVector for efficient analysis
high_accuracy = filter(r -> r.accuracy > 0.9, results)
mean_accuracy = Statistics.mean(results.accuracy)
```
"""
function run_scenario_optimization(
    scenarios::Vector{ScenarioParameters};
    ensemble_spec::EnsembleSpecification,
    time_params::SimTimeParameters,
    nsims::Int64 = 1000,
    seed::Int64 = 1234,
    cache_results::Bool = true,
    save_results::Bool = true,
    results_dir::String = "results/optimization",
    verbose::Bool = true
)
    # Initialize result cache
    result_cache = if cache_results
        _load_result_cache(results_dir)
    else
        Dict{ScenarioParameters, OptimizationResult}()
    end
    
    # Initialize results vector
    results = OptimizationResult[]
    
    # Progress tracking
    n_scenarios = length(scenarios)
    progress = verbose ? ProgressMeter.Progress(n_scenarios) : nothing
    
    for (i, scenario) in enumerate(scenarios)
        # Check cache
        if cache_results && haskey(result_cache, scenario)
            if verbose
                @info "Using cached result for scenario $i/$n_scenarios"
            end
            push!(results, result_cache[scenario])
        else
            # Run optimization
            if verbose
                @info "Running optimization for scenario $i/$n_scenarios"
            end
            
            result = _optimize_single_scenario(
                scenario,
                ensemble_spec,
                time_params,
                nsims,
                seed + i  # Different seed per scenario
            )
            
            push!(results, result)
            
            # Update cache
            if cache_results
                result_cache[scenario] = result
            end
        end
        
        # Update progress
        if verbose
            ProgressMeter.next!(progress)
        end
    end
    
    # Convert to StructVector for efficient storage
    results_sv = StructVector{OptimizationResult}(results)
    
    # Save results
    if save_results
        _save_results(results_sv, result_cache, results_dir)
    end
    
    return results_sv
end

"""
    _optimize_single_scenario(scenario, ensemble_spec, time_params, nsims, seed)

Internal function to optimize a single scenario.

# Arguments
- `scenario::ScenarioParameters`: Scenario to optimize
- `ensemble_spec::EnsembleSpecification`: Ensemble parameters
- `time_params::SimTimeParameters`: Time parameters
- `nsims::Int64`: Number of simulations
- `seed::Int64`: Random seed

# Returns
- `OptimizationResult`: Optimization result
"""
function _optimize_single_scenario(
    scenario::ScenarioParameters,
    ensemble_spec::EnsembleSpecification,
    time_params::SimTimeParameters,
    nsims::Int64,
    seed::Int64
)
    # 1. Create dynamics specification
    target_spec = DynamicsParameterSpecification(scenario.target_dynamics)
    
    # 2. Optimize noise vaccination if using dynamical noise
    if scenario.noise_spec isa DynamicalNoiseSpecification
        # Determine target noise level based on vaccination bounds
        # For now, use midpoint of bounds
        # TODO: Implement proper noise level optimization
        vax_coverage = mean(scenario.noise_spec.vaccination_bounds)
        noise_instance = DynamicalNoise(scenario.noise_spec, vax_coverage)
    else
        noise_instance = scenario.noise_spec
    end
    
    # 3. Run ensemble simulation
    # TODO: Integrate with actual ensemble simulation
    # For now, return placeholder
    
    # 4. Calculate detection characteristics
    # TODO: Implement detection characteristic calculation
    
    # 5. Optimize threshold
    # TODO: Implement threshold optimization
    
    # Placeholder result
    return OptimizationResult(
        scenario_params = scenario,
        optimal_threshold = 0.05,
        accuracy = 0.95,
        sensitivity = 0.92,
        specificity = 0.98,
        mean_detection_delay = 14.5,
        proportion_detected = 0.88
    )
end

"""
    _load_result_cache(results_dir)

Load cached results from disk.

# Arguments
- `results_dir::String`: Directory containing cached results

# Returns
- `Dict{ScenarioParameters, OptimizationResult}`: Cached results
"""
function _load_result_cache(results_dir::String)
    cache_file = joinpath(results_dir, "result_cache.jld2")
    
    if isfile(cache_file)
        return JLD2.load(cache_file, "cache")
    else
        return Dict{ScenarioParameters, OptimizationResult}()
    end
end

"""
    _save_results(results, cache, results_dir)

Save results and cache to disk.

# Arguments
- `results::StructVector{OptimizationResult}`: Results to save
- `cache::Dict{ScenarioParameters, OptimizationResult}`: Result cache
- `results_dir::String`: Directory for results
"""
function _save_results(
    results::StructVector{OptimizationResult},
    cache::Dict{ScenarioParameters, OptimizationResult},
    results_dir::String
)
    mkpath(results_dir)
    
    # Save results as StructVector
    JLD2.save(joinpath(results_dir, "optimization_results.jld2"), "results", results)
    
    # Save cache
    JLD2.save(joinpath(results_dir, "result_cache.jld2"), "cache", cache)
end
```

**Testing Requirements:**
- Test single scenario optimization
- Test caching functionality
- Test result saving/loading
- Test StructVector conversion

---

### 2.3 Implement Result Loading and Reshaping

**Priority: MEDIUM**
**New file to create:** `OutbreakDetectionCore/src/optimization-functions/scenario-optimization/results-loading.jl`

**Purpose:**
Utilities for loading and reshaping optimization results.

**Changes:**

```julia
export load_optimization_results, reshape_results_by_parameter

"""
    load_optimization_results(results_dir)

Load optimization results from disk.

# Arguments
- `results_dir::String`: Directory containing results

# Returns
- `StructVector{OptimizationResult}`: Loaded results

# Examples
```julia
results = load_optimization_results("results/optimization")

# Efficient filtering with StructVector
high_R0 = filter(r -> r.scenario_params.target_dynamics.R_0 > 15.0, results)
```
"""
function load_optimization_results(results_dir::String)
    results_file = joinpath(results_dir, "optimization_results.jld2")
    
    if !isfile(results_file)
        error("Results file not found: $results_file")
    end
    
    return JLD2.load(results_file, "results")
end

"""
    reshape_results_by_parameter(results, parameter_name)

Reshape results by a specific parameter for analysis.

# Arguments
- `results::StructVector{OptimizationResult}`: Results to reshape
- `parameter_name::Symbol`: Parameter to group by

# Returns
- `Dict`: Results grouped by parameter value

# Examples
```julia
results = load_optimization_results("results/optimization")

# Group by R_0
by_R0 = reshape_results_by_parameter(results, :R_0)

# Group by noise level
by_noise = reshape_results_by_parameter(results, :noise_level)

# Analyze each group
for (R0_value, group_results) in by_R0
    mean_acc = mean(group_results.accuracy)
    println("R_0 = $R0_value: mean accuracy = $mean_acc")
end
```
"""
function reshape_results_by_parameter(
    results::StructVector{OptimizationResult},
    parameter_name::Symbol
)
    # Extract parameter values
    param_values = _extract_parameter(results, parameter_name)
    
    # Group by unique values
    unique_values = unique(param_values)
    
    grouped = Dict()
    for value in unique_values
        mask = param_values .== value
        grouped[value] = results[mask]
    end
    
    return grouped
end

"""
    _extract_parameter(results, parameter_name)

Internal function to extract parameter values from results.

# Arguments
- `results::StructVector{OptimizationResult}`: Results
- `parameter_name::Symbol`: Parameter name

# Returns
- `Vector`: Parameter values
"""
function _extract_parameter(
    results::StructVector{OptimizationResult},
    parameter_name::Symbol
)
    # Navigate nested structure to extract parameter
    if parameter_name == :R_0
        return [r.scenario_params.target_dynamics.R_0 for r in results]
    elseif parameter_name == :test_sensitivity
        return [r.scenario_params.test_spec.sensitivity for r in results]
    elseif parameter_name == :percent_tested
        return [r.scenario_params.percent_clinic_tested for r in results]
    else
        error("Unknown parameter: $parameter_name")
    end
end
```

**Testing Requirements:**
- Test result loading
- Test reshaping by different parameters
- Test error handling

---

## Phase 3: Integration and Testing

### 3.1 Update Module Exports

**Priority: HIGH**
**File to modify:** `OutbreakDetectionCore/src/OutbreakDetectionCore.jl`

**Changes:**

```julia
# Add new includes
include("./types/scenario-parameters.jl")
include("./types/optimization-results.jl")
include("./optimization-functions/scenario-optimization/scenario-creation.jl")
include("./optimization-functions/scenario-optimization/optimization-wrapper.jl")
include("./optimization-functions/scenario-optimization/results-loading.jl")

# Add new exports
export ScenarioParameters,
       ThresholdOptimizationScenario,
       OptimizationResult,
       ThresholdOptimizationResult,
       create_scenario_grid,
       run_scenario_optimization,
       load_optimization_results,
       reshape_results_by_parameter
```

---

### 3.2 Create Integration Tests

**Priority: HIGH**
**New file to create:** `OutbreakDetectionCore/test/optimization-wrapper.jl`

**Changes:**

```julia
using OutbreakDetectionCore
using Test
using StructArrays
using AutoHashEquals

@testset "Optimization Wrapper" begin
    @testset "ScenarioParameters" begin
        target = TargetDiseaseDynamicsParameters(
            R_0 = 16.0,
            latent_period_days = 10.0,
            infectious_duration_days = 8.0,
            beta_force = 0.2
        )
        
        noise = DynamicalNoiseSpecification(
            R_0 = 5.0,
            latent_period = 7.0,
            duration_infection = 14.0,
            correlation = "in-phase",
            poisson_component = 1.0
        )
        
        test = IndividualTestSpecification(
            sensitivity = 1.0,
            specificity = 0.8,
            test_result_lag = 0
        )
        
        scenario1 = ScenarioParameters(
            target_dynamics = target,
            noise_spec = noise,
            test_spec = test,
            percent_clinic_tested = 0.5,
            outbreak_threshold = 0.01
        )
        
        scenario2 = ScenarioParameters(
            target_dynamics = target,
            noise_spec = noise,
            test_spec = test,
            percent_clinic_tested = 0.5,
            outbreak_threshold = 0.01
        )
        
        # Test AutoHashEquals
        @test scenario1 == scenario2
        @test hash(scenario1) == hash(scenario2)
        
        # Test inequality
        scenario3 = ScenarioParameters(
            target_dynamics = target,
            noise_spec = noise,
            test_spec = test,
            percent_clinic_tested = 0.7,  # Different
            outbreak_threshold = 0.01
        )
        
        @test scenario1 != scenario3
        @test hash(scenario1) != hash(scenario3)
    end
    
    @testset "OptimizationResult StructVector" begin
        # Create mock results
        results = [
            OptimizationResult(
                scenario_params = scenario,
                optimal_threshold = 0.05,
                accuracy = 0.95,
                sensitivity = 0.92,
                specificity = 0.98,
                mean_detection_delay = 14.5,
                proportion_detected = 0.88
            )
            for i in 1:10
        ]
        
        # Convert to StructVector
        results_sv = StructVector{OptimizationResult}(results)
        
        @test results_sv isa StructVector
        @test length(results_sv) == 10
        
        # Test column access
        @test results_sv.accuracy isa Vector{Float64}
        @test all(results_sv.accuracy .== 0.95)
        
        # Test filtering
        high_acc = filter(r -> r.accuracy > 0.9, results_sv)
        @test length(high_acc) == 10
    end
    
    @testset "Scenario Grid Creation" begin
        base_target = TargetDiseaseDynamicsParameters(
            R_0 = 16.0,
            latent_period_days = 10.0,
            infectious_duration_days = 8.0,
            beta_force = 0.2
        )
        
        base_noise = DynamicalNoiseSpecification(
            R_0 = 5.0,
            latent_period = 7.0,
            duration_infection = 14.0,
            correlation = "in-phase",
            poisson_component = 1.0
        )
        
        base_test = IndividualTestSpecification(
            sensitivity = 1.0,
            specificity = 0.8,
            test_result_lag = 0
        )
        
        scenarios = create_scenario_grid(
            R_0_values = [12.0, 16.0],
            noise_levels = ["low", "high"],
            test_sensitivities = [0.9, 1.0],
            percent_tested_values = [0.5],
            outbreak_thresholds = [0.01, 0.02],
            base_target = base_target,
            base_noise = base_noise,
            base_test = base_test
        )
        
        # Should have 2 × 2 × 2 × 1 × 2 = 16 scenarios
        @test length(scenarios) == 16
        
        # Test uniqueness (AutoHashEquals)
        @test length(unique(scenarios)) == 16
    end
end
```

---

## Phase 4: Documentation and Examples

### 4.1 Create Usage Guide

**Priority: MEDIUM**
**New file to create:** `OutbreakDetectionCore/docs/optimization-wrapper-guide.md`

**Content:**

```markdown
# Optimization Wrapper Usage Guide

## Overview

The optimization wrapper provides efficient parameter space exploration using:
- **StructVectors** for efficient result storage
- **AutoHashEquals** for automatic caching
- **Parallel execution** for faster optimization

## Basic Workflow

### 1. Define Base Parameters

```julia
using OutbreakDetectionCore

# Define target disease
target = TargetDiseaseDynamicsParameters(
    R_0 = 16.0,
    latent_period_days = 10.0,
    infectious_duration_days = 8.0,
    beta_force = 0.2
)

# Define noise specification
noise = DynamicalNoiseSpecification(
    R_0 = 5.0,
    latent_period = 7.0,
    duration_infection = 14.0,
    correlation = "in-phase",
    poisson_component = 1.0
)

# Define test specification
test = IndividualTestSpecification(
    sensitivity = 1.0,
    specificity = 0.8,
    test_result_lag = 0
)
```

### 2. Create Scenario Grid

```julia
scenarios = create_scenario_grid(
    R_0_values = [12.0, 16.0, 20.0],
    noise_levels = ["low", "medium", "high"],
    test_sensitivities = [0.8, 0.9, 1.0],
    percent_tested_values = [0.3, 0.5, 0.7],
    outbreak_thresholds = [0.01, 0.02, 0.05],
    base_target = target,
    base_noise = noise,
    base_test = test
)

println("Created $(length(scenarios)) scenarios")
```

### 3. Run Optimization

```julia
results = run_scenario_optimization(
    scenarios;
    ensemble_spec = ensemble_spec,
    time_params = time_params,
    nsims = 1000,
    seed = 1234,
    cache_results = true,
    save_results = true,
    results_dir = "results/my_optimization",
    verbose = true
)
```

### 4. Analyze Results

```julia
# Results are stored in StructVector for efficient analysis
mean_accuracy = mean(results.accuracy)
high_accuracy = filter(r -> r.accuracy > 0.9, results)

# Group by parameter
by_R0 = reshape_results_by_parameter(results, :R_0)

for (R0_value, group) in by_R0
    println("R_0 = $R0_value:")
    println("  Mean accuracy: $(mean(group.accuracy))")
    println("  Mean delay: $(mean(group.mean_detection_delay))")
end
```

### 5. Load Previous Results

```julia
# Load cached results
previous_results = load_optimization_results("results/my_optimization")

# Continue analysis
```

## Performance Tips

1. **Use caching**: Set `cache_results = true` to avoid re-running scenarios
2. **Save incrementally**: Results are saved after each scenario
3. **Use StructVectors**: Efficient column-wise access for analysis
4. **Parallel execution**: Use `run_scenario_optimization_parallel` for large grids

## Troubleshooting

### Cache Issues

If you need to clear the cache:
```julia
rm("results/my_optimization/result_cache.jld2")
```

### Memory Issues

For very large scenario grids, process in batches:
```julia
batch_size = 100
for i in 1:batch_size:length(scenarios)
    batch = scenarios[i:min(i+batch_size-1, end)]
    results = run_scenario_optimization(batch; ...)
end
```
```

---

## Implementation Timeline

### Week 1: Core Infrastructure
- [ ] Day 1: Add AutoHashEquals dependency (1.1)
- [ ] Day 2: Create scenario parameter types (1.2)
- [ ] Day 3: Create optimization result types (1.3)
- [ ] Day 4-5: Testing and validation

### Week 2: Optimization Wrapper
- [ ] Day 1-2: Implement scenario creation (2.1)
- [ ] Day 3-4: Implement optimization wrapper (2.2)
- [ ] Day 5: Implement result loading/reshaping (2.3)

### Week 3: Integration and Testing
- [ ] Day 1-2: Update module exports (3.1)
- [ ] Day 3-4: Create integration tests (3.2)
- [ ] Day 5: Documentation and examples (4.1)

---

## Success Criteria

### Functionality
- [ ] AutoHashEquals works for all parameter types
- [ ] StructVector storage works efficiently
- [ ] Caching prevents duplicate runs
- [ ] Results load/save correctly
- [ ] Scenario grid creation works

### Performance
- [ ] StructVector access faster than Vector
- [ ] Caching reduces runtime for repeated scenarios
- [ ] Memory usage reasonable for large grids

### Code Quality
- [ ] All tests pass
- [ ] Code formatted with Runic.jl
- [ ] Comprehensive documentation
- [ ] Examples run without errors

---

## Dependencies

### Required
- `AutoHashEquals.jl` (2.1+)
- `StructArrays.jl` (already present)
- `JLD2.jl` (already present)
- `ProgressMeter.jl` (already present)

### Optional
- `Distributed.jl` (for parallel execution)

---

## Notes

- **AutoHashEquals**: Automatically generates hash and equality methods
- **StructVector**: Column-wise storage for better performance
- **Caching**: Uses AutoHashEquals for efficient lookup
- **Incremental saving**: Results saved after each scenario
- **Backward compatibility**: Works with existing optimization code

---

## Future Enhancements

1. **Parallel execution**: Implement `run_scenario_optimization_parallel`
2. **Adaptive sampling**: Smart parameter space exploration
3. **Result visualization**: Built-in plotting functions
4. **Cloud storage**: Save results to cloud storage
5. **Resume capability**: Resume interrupted optimizations
