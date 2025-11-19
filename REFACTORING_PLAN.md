# OutbreakDetection Refactoring Plan

## Overview
This plan details the structural refactoring of OutbreakDetection to align with the improved patterns from CSDNoise, including the addition of dynamical noise vaccination optimization capabilities.

## Goals
1. Remove global constants in favor of explicit keyword arguments
2. Adopt hierarchical parameter specification pattern
3. Add comprehensive struct validation and documentation
4. Implement dynamical noise vaccination optimization
5. Improve type design with sum types and clear separation of concerns
6. **Migrate to StructVectors for efficient ensemble storage**
7. **Use LightSumTypes for type-safe seasonality and noise specifications**
8. **Add endemic equilibrium calculation for noise simulations**
9. Maintain backward compatibility where possible during transition

## Key Additions from CSDNoise

### StructVector Integration
- All result types (`SEIRRun`, `DynamicalNoiseRun`, `PoissonNoiseRun`) designed for StructVector compatibility
- Efficient column-wise storage for large ensembles
- Better memory layout and performance

### Type Safety with Sum Types
- `SeasonalityFunction`: Type-safe cosine/sine seasonality
- `NoiseSpecification`: Type-safe Poisson/Dynamical noise variants
- Pattern matching for cleaner code

### Endemic Equilibrium
- Calculate equilibrium proportions for SEIR model with vaccination
- Used to set initial states for noise simulations
- Handles edge cases (R_eff ≤ 1) with `Try.jl`

### Simplified Noise Workflow
Unlike CSDNoise, OutbreakDetection does **not** use:
- Burnin periods
- Target Reff calculations
- Post-burnin vaccination ranges

Instead, vaccination coverage is determined **solely by optimization** to match target noise levels.

---

## Phase 1: Core Type System Refactoring

### 1.1 Dynamics Parameters Restructuring

**Priority: CRITICAL**
**Files to modify:**
- `OutbreakDetectionCore/src/types/dynamics-parameters.jl`
- `OutbreakDetectionCore/src/constants/dynamics-constants.jl` (DELETE)

**Changes:**

#### Step 1.1.1: Create New Type Hierarchy
Create three-tier hierarchy similar to CSDNoise:

```julia
# First, define seasonality sum types (like CSDNoise)
abstract type AbstractSeasonalityFunction end
struct CosineSeasonality end
struct SineSeasonality end
LightSumTypes.@sumtype SeasonalityFunction(
    CosineSeasonality,
    SineSeasonality
) <: AbstractSeasonalityFunction

# Tier 1: User-facing high-level specification
Base.@kwdef struct TargetDiseaseDynamicsParameters
    R_0::Float64
    latent_period_days::Float64
    infectious_duration_days::Float64
    beta_force::Float64
    seasonality::SeasonalityFunction = SeasonalityFunction(CosineSeasonality())  # Default to cosine
    life_expectancy_years::Float64 = 62.5
    population_N::Int64 = 500_000
end

# Tier 2: Derived specification with calculated values
Base.@kwdef struct DynamicsParameterSpecification
    beta_mean::Float64
    beta_force::Float64
    seasonality::SeasonalityFunction
    sigma::Float64
    gamma::Float64
    mu::Float64
    annual_births_per_k::Float64
    epsilon::Float64
    R_0::Float64
    population_N::Int64
end

# Tier 3: Concrete instance for simulation
Base.@kwdef struct DynamicsParameters
    beta_mean::Float64
    beta_force::Float64
    seasonality::SeasonalityFunction
    sigma::Float64
    gamma::Float64
    mu::Float64
    annual_births_per_k::Float64
    epsilon::Float64
    R_0::Float64
    vaccination_coverage::Float64
end
```

#### Step 1.1.2: Add Constructors with Validation

```julia
# Constructor: TargetDiseaseDynamicsParameters → DynamicsParameterSpecification
function DynamicsParameterSpecification(target::TargetDiseaseDynamicsParameters)
    @assert target.R_0 > 0 "R_0 must be positive"
    @assert target.latent_period_days > 0 "Latent period must be positive"
    @assert target.infectious_duration_days > 0 "Infectious duration must be positive"
    @assert 0 <= target.beta_force <= 1 "Beta force must be in [0, 1]"
    @assert target.life_expectancy_years > 0 "Life expectancy must be positive"
    @assert target.population_N > 0 "Population must be positive"

    sigma = 1.0 / target.latent_period_days
    gamma = 1.0 / target.infectious_duration_days
    mu = 1.0 / (target.life_expectancy_years * 365.0)
    annual_births_per_k = 1000.0 / target.life_expectancy_years
    beta_mean = calculate_beta(target.R_0, gamma, mu, 1, target.population_N)
    epsilon = calculate_import_rate(mu, target.R_0, target.population_N)

    return DynamicsParameterSpecification(
        beta_mean = beta_mean,
        beta_force = target.beta_force,
        seasonality = target.seasonality,
        sigma = sigma,
        gamma = gamma,
        mu = mu,
        annual_births_per_k = annual_births_per_k,
        epsilon = epsilon,
        R_0 = target.R_0,
        population_N = target.population_N
    )
end

# Constructor: DynamicsParameterSpecification → DynamicsParameters
function DynamicsParameters(
    spec::DynamicsParameterSpecification;
    vaccination_coverage::Float64 = 0.8
)
    @assert 0.0 <= vaccination_coverage <= 1.0 "Vaccination coverage must be in [0, 1]"

    return DynamicsParameters(
        beta_mean = spec.beta_mean,
        beta_force = spec.beta_force,
        seasonality = spec.seasonality,
        sigma = spec.sigma,
        gamma = spec.gamma,
        mu = spec.mu,
        annual_births_per_k = spec.annual_births_per_k,
        epsilon = spec.epsilon,
        R_0 = spec.R_0,
        vaccination_coverage = vaccination_coverage
    )
end
```

#### Step 1.1.3: Add Comprehensive Documentation
Add detailed docstrings with examples for each struct and constructor.

#### Step 1.1.4: Delete Constants File
Remove `OutbreakDetectionCore/src/constants/dynamics-constants.jl` entirely.

**Testing Requirements:**
- Test validation assertions (invalid inputs should error)
- Test constructor chain: Target → Specification → Parameters
- Test backward compatibility with existing code
- Run JET.jl type stability tests

---

### 1.2 Time Parameters Simplification

**Priority: HIGH**
**Files to modify:**
- `OutbreakDetectionCore/src/types/time-parameters.jl`

**Changes:**

#### Step 1.2.1: Simplify Struct Definition
Remove parametric types, add validation:

```julia
Base.@kwdef struct SimTimeParameters
    tmin::Float64 = 0.0
    tmax::Float64 = 365.0 * 100.0
    tstep::Float64 = 1.0
    trange::StepRangeLen{Float64, Float64, Float64, Int64}
    tspan::Tuple{Float64, Float64}
    tlength::Int64

    function SimTimeParameters(tmin, tmax, tstep, trange, tspan, tlength)
        @assert tmax > tmin + tstep "tmax must be greater than tmin + tstep"
        @assert tstep > 0.0 "tstep must be positive"
        return new(tmin, tmax, tstep, trange, tspan, tlength)
    end
end

# Keyword constructor
function SimTimeParameters(;
    tmin::Float64 = 0.0,
    tmax::Float64 = 365.0 * 100.0,
    tstep::Float64 = 1.0
)
    trange = tmin:tstep:tmax
    tspan = (tmin, tmax)
    tlength = length(trange)

    return SimTimeParameters(tmin, tmax, tstep, trange, tspan, tlength)
end
```

**Testing Requirements:**
- Test validation (invalid time ranges should error)
- Test keyword constructor with defaults
- Verify existing code compatibility

---

### 1.3 Noise Specifications Refactoring

**Priority: HIGH**
**Files to modify:**
- `OutbreakDetectionCore/src/types/noise-specifications.jl`

**Changes:**

#### Step 1.3.1: Add LightSumTypes Dependency
Add to `Project.toml`:
```toml
LightSumTypes = "0.1"
```

#### Step 1.3.2: Restructure Noise Types
Separate specification from instance:

```julia
# Abstract types for sum type variants
abstract type AbstractNoiseType end
struct PoissonNoiseType end
struct DynamicalNoiseType end

LightSumTypes.@sumtype NoiseType(PoissonNoiseType, DynamicalNoiseType) <: AbstractNoiseType

# Specification: defines parameters and bounds
Base.@kwdef struct DynamicalNoiseSpecification
    R_0::Float64
    latent_period::Float64
    duration_infection::Float64
    correlation::String  # "in-phase", "out-of-phase", "none"
    poisson_component::Float64
    vaccination_bounds::Vector{Float64} = [0.0, 1.0]

    function DynamicalNoiseSpecification(R_0, latent_period, duration_infection,
                                         correlation, poisson_component, vaccination_bounds)
        @assert R_0 > 0 "R_0 must be positive"
        @assert latent_period > 0 "Latent period must be positive"
        @assert duration_infection > 0 "Infectious duration must be positive"
        @assert correlation in ["in-phase", "out-of-phase", "none"] "Invalid correlation type"
        @assert poisson_component >= 0 "Poisson component must be non-negative"
        @assert length(vaccination_bounds) == 2 "Vaccination bounds must have 2 elements"
        @assert 0 <= vaccination_bounds[1] < vaccination_bounds[2] <= 1 "Invalid vaccination bounds"

        return new(R_0, latent_period, duration_infection, correlation,
                   poisson_component, vaccination_bounds)
    end
end

# Instance: concrete noise with specific vaccination coverage
Base.@kwdef struct DynamicalNoise
    R_0::Float64
    latent_period::Float64
    duration_infection::Float64
    correlation::String
    poisson_component::Float64
    vaccination_coverage::Float64
end

# Constructor: Specification → Instance
function DynamicalNoise(spec::DynamicalNoiseSpecification, vaccination_coverage::Float64)
    @assert 0.0 <= vaccination_coverage <= 1.0 "Vaccination coverage must be in [0, 1]"

    return DynamicalNoise(
        R_0 = spec.R_0,
        latent_period = spec.latent_period,
        duration_infection = spec.duration_infection,
        correlation = spec.correlation,
        poisson_component = spec.poisson_component,
        vaccination_coverage = vaccination_coverage
    )
end

# Poisson noise (simpler, no optimization needed)
Base.@kwdef struct PoissonNoise
    noise_mean_scaling::Float64
end

# Sum type for polymorphic noise handling
abstract type AbstractNoiseSpecification end
LightSumTypes.@sumtype NoiseSpecification(PoissonNoise, DynamicalNoise) <: AbstractNoiseSpecification
```

**Testing Requirements:**
- Test validation for all noise types
- Test sum type pattern matching
- Test specification → instance conversion

---

### 1.4 Test Specifications Enhancement

**Priority: MEDIUM**
**Files to modify:**
- `OutbreakDetectionCore/src/types/test-specifications.jl`
- `OutbreakDetectionCore/src/constants/test-constants.jl` (KEEP but simplify)

**Changes:**

#### Step 1.4.1: Add Keyword Constructor
```julia
Base.@kwdef struct IndividualTestSpecification
    sensitivity::Float64
    specificity::Float64
    test_result_lag::Int64

    function IndividualTestSpecification(sensitivity, specificity, test_result_lag)
        @assert 0.0 <= sensitivity <= 1.0 "Sensitivity must be in [0, 1]"
        @assert 0.0 <= specificity <= 1.0 "Specificity must be in [0, 1]"
        @assert test_result_lag >= 0 "Test result lag must be non-negative"
        return new(sensitivity, specificity, test_result_lag)
    end
end
```

#### Step 1.4.2: Keep Constants File
The test constants are legitimate reusable specifications, so keep them:
```julia
# OutbreakDetectionCore/src/constants/test-constants.jl
const CLINICAL_CASE_TEST_SPEC = IndividualTestSpecification(
    sensitivity = 1.0,
    specificity = 0.0,
    test_result_lag = 0
)

const EPI_LINKED_CASE_TEST_SPEC = IndividualTestSpecification(
    sensitivity = 1.0,
    specificity = 0.8,
    test_result_lag = 0
)

const CLINICAL_TEST_SPECS = (CLINICAL_CASE_TEST_SPEC, EPI_LINKED_CASE_TEST_SPEC)
```

**Testing Requirements:**
- Test validation
- Test keyword constructor
- Verify constants still work

---

### 1.5 Ensemble Specifications Simplification

**Priority: MEDIUM**
**Files to modify:**
- `OutbreakDetectionCore/src/types/ensemble-specifications.jl`

**Changes:**

#### Step 1.5.1: Simplify Struct
Remove parametric types:

```julia
Base.@kwdef struct EnsembleSpecification
    modeltypes::Tuple
    state_parameters::StateParameters
    dynamics_parameters::DynamicsParameters
    time_parameters::SimTimeParameters
    nsims::Int64
    dirpath::String
end

# Constructor with automatic dirpath generation
function EnsembleSpecification(
    modeltypes::Tuple,
    state_parameters::StateParameters,
    dynamics_parameters::DynamicsParameters,
    time_parameters::SimTimeParameters,
    nsims::Int64
)
    @assert nsims > 0 "Number of simulations must be positive"

    dirpath = outdir(
        "ensemble",
        modeltypes...,
        "N_$(state_parameters.init_states.N)",
        "r_$(state_parameters.init_state_props.r_prop)",
        "nsims_$(nsims)",
        "R0_$(dynamics_parameters.R_0)",
        "latent_period_$(round(1 / dynamics_parameters.sigma; digits = 2))",
        "infectious_period_$(round(1 / dynamics_parameters.gamma; digits = 2))",
        "vaccination_coverage_$(dynamics_parameters.vaccination_coverage)",
        "births_per_k_$(dynamics_parameters.annual_births_per_k)",
        "beta_force_$(dynamics_parameters.beta_force)",
        "tmax_$(time_parameters.tmax)",
        "tstep_$(time_parameters.tstep)"
    )

    return EnsembleSpecification(
        modeltypes,
        state_parameters,
        dynamics_parameters,
        time_parameters,
        nsims,
        dirpath
    )
end
```

**Testing Requirements:**
- Test validation
- Test dirpath generation
- Verify backward compatibility

---

### 1.6 Add Simulation Result Types (StructVectors)

**Priority: HIGH**  
**New files to create:**
- `OutbreakDetectionCore/src/types/simulation-results.jl`

**Changes:**

#### Step 1.6.1: Define SEIR Result Types
Port from CSDNoise with StructVector support:

```julia
# OutbreakDetectionCore/src/types/simulation-results.jl
export SEIRRun, NoiseRun, DynamicalNoiseRun, PoissonNoiseRun

"""
    SEIRRun

Results from a single SEIR epidemic simulation run.

Designed to work efficiently with StructVectors for ensemble storage.

# Fields
- `states::Vector{StaticArrays.SVector{5, Int64}}`: Time series of [S, E, I, R, N]
- `incidence::Vector{Int64}`: Time series of new infections
- `Reff::Vector{Float64}`: Time series of effective reproduction number
"""
Base.@kwdef struct SEIRRun
    states::Vector{StaticArrays.SVector{5, Int64}}
    incidence::Vector{Int64}
    Reff::Vector{Float64}
end

"""
    DynamicalNoiseRun

Results from noise simulations with both Poisson and dynamical noise.

# Fields
- `incidence::Vector{Vector{Int64}}`: Collection of noisy incidence time series
- `Reff::Vector{Vector{Float64}}`: Collection of noisy Reff time series
- `mean_noise::Float64`: Overall mean noise level
- `mean_poisson_noise::Float64`: Mean Poisson noise component
- `mean_dynamic_noise::Float64`: Mean dynamical noise component
"""
Base.@kwdef struct DynamicalNoiseRun
    incidence::Vector{Vector{Int64}}
    Reff::Vector{Vector{Float64}}
    mean_noise::Float64
    mean_poisson_noise::Float64
    mean_dynamic_noise::Float64
end

"""
    PoissonNoiseRun

Results from noise simulations with only Poisson noise.

# Fields
- `incidence::Vector{Vector{Int64}}`: Collection of noisy incidence time series
- `mean_noise::Float64`: Mean Poisson noise level
"""
Base.@kwdef struct PoissonNoiseRun
    incidence::Vector{Vector{Int64}}
    mean_noise::Float64
end

# Sum type for noise runs
abstract type AbstractNoiseRun end
LightSumTypes.@sumtype NoiseRun(DynamicalNoiseRun, PoissonNoiseRun) <: AbstractNoiseRun
```

**Dependencies to add:**
- `StructArrays.jl` for StructVector support
- `StaticArrays.jl` for efficient state vectors

**Testing Requirements:**
- Test struct construction
- Test StructVector compatibility
- Test sum type pattern matching

---

### 1.7 Add Endemic Equilibrium Calculation

**Priority: HIGH**  
**New files to create:**
- `OutbreakDetectionCore/src/simulation/endemic-equilibrium.jl`

**Changes:**

#### Step 1.7.1: Implement Endemic Equilibrium Function
Port from CSDNoise:

```julia
# OutbreakDetectionCore/src/simulation/endemic-equilibrium.jl
export calculate_endemic_equilibrium_proportions

"""
Calculate endemic equilibrium proportions for SEIR model with vaccination.

# Arguments
- `dynamics_params::Union{DynamicsParameters, DynamicsParameterSpecification}`: Disease dynamics
- `vaccination_coverage::Float64`: Vaccination coverage (0-1)

# Returns
- `Try.Ok((s_prop, e_prop, i_prop, r_prop))` on success
- `Try.Err(message)` if equilibrium doesn't exist (R_eff ≤ 1)
"""
function calculate_endemic_equilibrium_proportions(
    dynamics_params::Union{DynamicsParameters, DynamicsParameterSpecification},
    vaccination_coverage::Float64
)
    return calculate_endemic_equilibrium_proportions(
        dynamics_params.R_0,
        dynamics_params.beta_mean,
        dynamics_params.sigma,
        dynamics_params.gamma,
        dynamics_params.mu,
        vaccination_coverage
    )
end

function calculate_endemic_equilibrium_proportions(
    R_0::Float64,
    beta_mean::Float64,
    sigma::Float64,
    gamma::Float64,
    mu::Float64,
    vaccination_coverage::Float64
)
    if vaccination_coverage < 0.0 || vaccination_coverage > 1.0
        return Try.Err("Vaccination coverage must be in [0, 1]. Got $vaccination_coverage")
    end
    
    R_eff = R_0 * (1.0 - vaccination_coverage)
    
    if R_eff <= 1.0
        return Try.Err(
            "Effective R₀ must be > 1 for endemic equilibrium. " *
            "Got R_eff = R₀(1 - ρ) = $(round(R_eff; digits = 2))"
        )
    end
    
    s_prop = 1.0 / R_0
    i_prop = mu * (R_eff - 1.0) / beta_mean
    e_prop = i_prop * (gamma + mu) / sigma
    r_prop = 1.0 - s_prop - e_prop - i_prop
    
    if r_prop < 0.0
        return Try.Err("Calculated negative recovered proportion. Check parameter consistency.")
    end
    
    return Try.Ok((
        s_prop = s_prop,
        e_prop = e_prop,
        i_prop = i_prop,
        r_prop = r_prop
    ))
end
```

**Testing Requirements:**
- Test with valid parameters
- Test error cases (R_eff ≤ 1, invalid vaccination coverage)
- Test Try.jl error handling

---

## Phase 2: Vaccination and Noise Optimization

**Key Difference from CSDNoise:**  
OutbreakDetection does not use burnin periods or target Reff calculations. The vaccination coverage for dynamical noise is determined solely by optimization to match a target noise level (Section 2.3). This is simpler than CSDNoise's approach.

**Workflow:**
1. Run target disease simulations → calculate mean incidence
2. Choose target_scaling (e.g., 1.0 for "low noise", 7.0 for "high noise")
3. Optimize to find vaccination coverage where: `mean_noise = target_scaling × mean_target_incidence`
4. Use that vaccination coverage in noise simulations

### 2.1 Simplify Vaccination Utilities (Optional)

**Priority: LOW**  
**Note:** Since OutbreakDetection doesn't use burnin periods or target Reff calculations, the vaccination utilities from CSDNoise are not needed. The vaccination coverage for dynamical noise is determined entirely by the optimization process in Section 2.3.

**If you want basic vaccination sampling utilities for other purposes:**

#### Step 2.1.1: Create Vaccination Directory (Optional)
```bash
mkdir -p OutbreakDetectionCore/src/vaccination
```

#### Step 2.1.2: Implement Vaccination Sampling (Optional)
Only needed if you want to sample vaccination coverage for non-noise purposes:

```julia
# OutbreakDetectionCore/src/vaccination/vaccination-sampling.jl
export sample_vaccination_coverage

"""
Sample vaccination coverage from uniform distribution.

# Arguments
- `min_coverage`: Minimum coverage (0-1)
- `max_coverage`: Maximum coverage (0-1)
- `digits`: Rounding precision (default: 4)

# Returns
- `Float64`: Sampled vaccination coverage
"""
function sample_vaccination_coverage(
    min_coverage::Float64,
    max_coverage::Float64,
    digits::Int = 4
)
    @assert 0.0 <= min_coverage <= 1.0 "Min coverage must be in [0, 1]"
    @assert 0.0 <= max_coverage <= 1.0 "Max coverage must be in [0, 1]"
    @assert min_coverage < max_coverage "Min must be less than max"
    
    return round(
        rand(Distributions.Uniform(min_coverage, max_coverage));
        digits = digits
    )
end
```

**Skip this section if not needed** - The optimization in Section 2.3 will handle finding the correct vaccination coverage automatically.

---

### 2.2 Add Noise Mean Calculation

**Priority: HIGH**
**New files to create:**
- `OutbreakDetectionCore/src/noise/noise-mean-incidence.jl`
- `OutbreakDetectionCore/src/noise/noise-recreation.jl`
- `OutbreakDetectionCore/src/noise/noise-dynamics-parameters.jl`

**Changes:**

#### Step 2.2.1: Implement Mean Incidence Calculation
```julia
# OutbreakDetectionCore/src/noise/noise-mean-incidence.jl
export calculate_mean_incidence

"""
Calculate mean incidence from SEIR results.

Works with both Vector and StructVector of SEIRRun results.

# Arguments
- `seir_results`: StructVector{SEIRRun} or Vector{SEIRRun}

# Returns
- `Float64`: Mean daily incidence across all simulations and time points
"""
function calculate_mean_incidence(seir_results::StructVector{SEIRRun})
    # StructVector provides efficient column access
    # seir_results.incidence is a Vector{Vector{Int64}}
    all_incidence = vcat(seir_results.incidence...)
    return mean(all_incidence)
end

function calculate_mean_incidence(seir_results::Vector{SEIRRun})
    # Fallback for regular Vector
    all_incidence = vcat([result.incidence for result in seir_results]...)
    return mean(all_incidence)
end
```

#### Step 2.2.2: Implement Noise Dynamics Parameter Creation
```julia
# OutbreakDetectionCore/src/noise/noise-dynamics-parameters.jl
export create_noise_dynamics_parameters

"""
Create dynamics parameters for noise simulation.

Adjusts seasonality and transmission parameters based on
correlation type specified in noise specification.

# Arguments
- `noise_spec::DynamicalNoise`: Noise specification with vaccination coverage
- `base_dynamics::DynamicsParameterSpecification`: Base disease dynamics
- `population_N::Int64`: Population size

# Returns
- `DynamicsParameters`: Dynamics parameters for noise simulation
"""
function create_noise_dynamics_parameters(
    noise_spec::DynamicalNoise,
    base_dynamics::DynamicsParameterSpecification,
    population_N::Int64
)
    # Adjust beta_force based on correlation
    noise_beta_force = if noise_spec.correlation == "none"
        0.0
    else
        base_dynamics.beta_force
    end

    # Adjust seasonality based on correlation (using sum types)
    noise_seasonality = if noise_spec.correlation == "out-of-phase"
        # Flip seasonality for out-of-phase correlation
        if LightSumTypes.variant(base_dynamics.seasonality) isa CosineSeasonality
            SeasonalityFunction(SineSeasonality())
        elseif LightSumTypes.variant(base_dynamics.seasonality) isa SineSeasonality
            SeasonalityFunction(CosineSeasonality())
        else
            base_dynamics.seasonality
        end
    else
        base_dynamics.seasonality
    end

    # Calculate noise-specific parameters
    noise_gamma = 1.0 / noise_spec.duration_infection
    noise_sigma = 1.0 / noise_spec.latent_period
    noise_beta_mean = calculate_beta(
        noise_spec.R_0,
        noise_gamma,
        base_dynamics.mu,
        1,
        population_N
    )
    noise_epsilon = calculate_import_rate(
        base_dynamics.mu,
        noise_spec.R_0,
        population_N
    )

    return DynamicsParameters(
        beta_mean = noise_beta_mean,
        beta_force = noise_beta_force,
        seasonality = noise_seasonality,
        sigma = noise_sigma,
        gamma = noise_gamma,
        mu = base_dynamics.mu,
        annual_births_per_k = base_dynamics.annual_births_per_k,
        epsilon = noise_epsilon,
        R_0 = noise_spec.R_0,
        vaccination_coverage = noise_spec.vaccination_coverage
    )
end
```

#### Step 2.2.3: Implement Noise Recreation
```julia
# OutbreakDetectionCore/src/noise/noise-recreation.jl
export recreate_noise_vecs, calculate_mean_dynamical_noise

"""
Recreate noise vectors with updated vaccination coverage.

Uses endemic equilibrium to set initial states for noise simulations.

# Arguments
- `noise_spec::DynamicalNoise`: Noise specification with vaccination coverage
- `ensemble_spec::EnsembleSpecification`: Ensemble parameters
- `base_dynamics::DynamicsParameterSpecification`: Base dynamics

# Keyword Arguments
- `verbose::Bool`: Print warnings (default: false)
- `seed::Int`: Random seed (default: 1234)

# Returns
- `DynamicalNoiseRun` or `PoissonNoiseRun` with noise simulation results
"""
function recreate_noise_vecs(
    noise_spec::DynamicalNoise,
    ensemble_spec::EnsembleSpecification,
    base_dynamics::DynamicsParameterSpecification;
    verbose::Bool = false,
    seed::Int = 1234
)
    @unpack state_parameters, time_parameters, nsims = ensemble_spec
    @unpack N = state_parameters.init_states
    
    # Create noise-specific dynamics
    noise_dynamics = create_noise_dynamics_parameters(
        noise_spec,
        base_dynamics,
        N
    )
    
    # Calculate endemic equilibrium for initial states
    endemic_result = calculate_endemic_equilibrium_proportions(
        noise_dynamics,
        noise_spec.vaccination_coverage
    )
    
    # Handle endemic equilibrium calculation
    endemic_props = if Try.isok(endemic_result)
        Try.unwrap(endemic_result)
    else
        if verbose
            @warn Try.unwrap_err(endemic_result) *
                "\nDefaulting to no initial infections, s_prop = 1 - vaccination_coverage"
        end
        (
            s_prop = 1.0 - noise_spec.vaccination_coverage,
            e_prop = 0.0,
            i_prop = 0.0,
            r_prop = noise_spec.vaccination_coverage
        )
    end
    
    # Create updated state parameters with endemic equilibrium
    noise_state_params = StateParameters(
        N = N,
        s_prop = endemic_props.s_prop,
        e_prop = endemic_props.e_prop,
        i_prop = endemic_props.i_prop
    )
    
    # Run ensemble simulation with noise dynamics
    # Returns StructVector{SEIRRun}
    noise_seir_results = run_ensemble_simulation(
        noise_state_params,
        noise_dynamics,
        time_parameters,
        nsims;
        seed = seed
    )
    
    # Calculate mean dynamical noise from SEIR results
    mean_dynamic_noise = calculate_mean_incidence(noise_seir_results)
    
    # Add Poisson component
    mean_poisson_noise = noise_spec.poisson_component * mean_dynamic_noise
    total_mean_noise = mean_dynamic_noise + mean_poisson_noise
    
    # Extract incidence and Reff vectors
    incidence_vecs = [run.incidence for run in noise_seir_results]
    reff_vecs = [run.Reff for run in noise_seir_results]
    
    return DynamicalNoiseRun(
        incidence = incidence_vecs,
        Reff = reff_vecs,
        mean_noise = total_mean_noise,
        mean_poisson_noise = mean_poisson_noise,
        mean_dynamic_noise = mean_dynamic_noise
    )
end

"""
Calculate mean dynamical noise for given vaccination coverage.

Wrapper function for optimization objective.

# Arguments
- `noise_spec::DynamicalNoise`: Noise specification
- `ensemble_spec::EnsembleSpecification`: Ensemble parameters
- `base_dynamics::DynamicsParameterSpecification`: Base dynamics

# Keyword Arguments
- `verbose::Bool`: Print warnings (default: false)
- `seed::Int`: Random seed (default: 1234)

# Returns
- `Float64`: Mean noise level
"""
function calculate_mean_dynamical_noise(
    noise_spec::DynamicalNoise,
    ensemble_spec::EnsembleSpecification,
    base_dynamics::DynamicsParameterSpecification;
    verbose::Bool = false,
    seed::Int = 1234
)
    result = recreate_noise_vecs(
        noise_spec,
        ensemble_spec,
        base_dynamics;
        verbose = verbose,
        seed = seed
    )
    return result.mean_noise
end
```

**Testing Requirements:**
- Test mean incidence calculation
- Test noise dynamics parameter creation for all correlation types
- Test noise recreation produces expected distributions

---

### 2.3 Implement Vaccination Optimization

**Priority: HIGH**  
**New files to create:**
- `OutbreakDetectionCore/src/noise/noise-parameters-optimization.jl`
- `OutbreakDetectionCore/src/types/optimization-parameters.jl`

**Purpose:**  
This is the core mechanism for determining the vaccination coverage in dynamical noise simulations. The optimization finds the vaccination coverage that produces a mean noise level equal to a specified multiple of the target disease mean incidence. Unlike CSDNoise, there is no burnin period or target Reff - the optimization simply matches noise levels.

**Changes:**

#### Step 2.3.1: Create Optimization Parameters Type
```julia
# OutbreakDetectionCore/src/types/optimization-parameters.jl
export NoiseVaccinationOptimizationParameters

"""
Parameters for optimizing vaccination coverage in dynamical noise.

Uses multistart optimization with Sobol sequences for global search
and NLopt for local optimization.

# Fields
- `n_sobol_points::Int64`: Number of Sobol sequence points for initialization
- `local_algorithm`: NLopt algorithm for local optimization
- `maxeval::Int64`: Maximum function evaluations
- `xtol_rel::Float64`: Relative tolerance on parameters
- `xtol_abs::Float64`: Absolute tolerance on parameters
- `atol::Float64`: Absolute difference tolerance threshold
"""
Base.@kwdef struct NoiseVaccinationOptimizationParameters
    n_sobol_points::Int64 = 100
    local_algorithm = NLopt.LN_BOBYQA
    maxeval::Int64 = 1000
    xtol_rel::Float64 = 1.0e-3
    xtol_abs::Float64 = 1.0e-3
    atol::Float64 = 1.0e-2
end
```

#### Step 2.3.2: Implement Optimization Function
```julia
# OutbreakDetectionCore/src/noise/noise-parameters-optimization.jl
export optimize_dynamic_noise_params

"""
Optimize vaccination coverage to achieve target noise level.

This is the primary function for determining vaccination coverage in dynamical noise
simulations. It finds the vaccination coverage that produces a mean noise level equal
to `target_scaling * mean_target_incidence`.

**Workflow:**
1. User specifies target_scaling (e.g., 7.0 for "high noise")
2. User provides mean_target_incidence from target disease simulations
3. Optimization finds vaccination coverage where:
   mean_noise_incidence = target_scaling × mean_target_incidence

**No burnin or Reff targets needed** - this is purely a noise level matching problem.

Uses multistart optimization with Sobol sequences for global search and local NLopt
optimization to find the optimal vaccination coverage within specified bounds.

# Arguments
- `target_scaling::Float64`: Multiplicative factor for target noise (e.g., 1.0 for "low", 7.0 for "high")
- `mean_target_incidence::Float64`: Mean daily incidence from target disease simulations
- `noise_spec::DynamicalNoiseSpecification`: Noise specification with vaccination_bounds
- `ensemble_spec::EnsembleSpecification`: Ensemble parameters for noise simulations
- `base_dynamics::DynamicsParameterSpecification`: Base disease dynamics (not used for Reff)
- `opt_params::NoiseVaccinationOptimizationParameters`: Optimization settings (optional)

# Keyword Arguments
- `verbose::Bool`: Print optimization progress (default: false)
- `seed::Int`: Random seed for reproducibility (default: 1234)

# Returns
Named tuple with:
- `optimal_vaccination::Float64`: Optimal vaccination coverage that achieves target noise
- `mean_noise::Float64`: Achieved mean noise level
- `target_noise::Float64`: Target noise level (target_scaling × mean_target_incidence)
- `difference::Float64`: Absolute difference from target
- `optimization_result`: Full optimization result object

# Example
```julia
# Define noise specification with vaccination bounds
noise_spec = DynamicalNoiseSpecification(
    R_0 = 5.0,
    latent_period = 7.0,
    duration_infection = 14.0,
    correlation = "in-phase",
    poisson_component = 1.0,
    vaccination_bounds = [0.0, 0.8]  # Search range for optimization
)

# Optimize to find vaccination coverage for "high noise" (7x target incidence)
result = optimize_dynamic_noise_params(
    7.0,   # target_scaling: want noise 7x higher than target disease
    3.0,   # mean_target_incidence: from target disease simulations
    noise_spec,
    ensemble_spec,
    base_dynamics
)

# Result: vaccination coverage that produces mean noise ≈ 21.0 (7.0 × 3.0)
println("Use vaccination coverage: \$(result.optimal_vaccination)")
println("This produces mean noise: \$(result.mean_noise)")
```
"""
function optimize_dynamic_noise_params(
    target_scaling::Float64,
    mean_target_incidence::Float64,
    noise_spec::DynamicalNoiseSpecification,
    ensemble_spec::EnsembleSpecification,
    base_dynamics::DynamicsParameterSpecification,
    opt_params::NoiseVaccinationOptimizationParameters = NoiseVaccinationOptimizationParameters();
    verbose::Bool = false,
    seed::Int = 1234
)
    @unpack vaccination_bounds = noise_spec
    @unpack n_sobol_points, local_algorithm, xtol_rel, xtol_abs, maxeval, atol = opt_params
    
    target_noise = target_scaling * mean_target_incidence

    if verbose
        println("Target noise level: $target_noise")
        println("Vaccination bounds: $vaccination_bounds")
        println("Starting multistart optimization with $n_sobol_points points...")
    end

    # Define objective function
    objective = function(params)
        vaccination_coverage = params[1]

        # Create noise instance with this vaccination coverage
        noise_instance = DynamicalNoise(noise_spec, vaccination_coverage)

        # Calculate mean noise
        noise_level = calculate_mean_dynamical_noise(
            noise_instance,
            ensemble_spec,
            base_dynamics;
            seed = seed
        )

        if verbose
            println("  Vax: $(round(vaccination_coverage; digits=4)), " *
                   "Noise: $(round(noise_level; digits=2)), " *
                   "Target: $(round(target_noise; digits=2))")
        end

        # Return squared error
        return (noise_level - target_noise)^2
    end

    # Setup multistart optimization
    problem = MultistartOptimization.MinimizationProblem(
        objective,
        [vaccination_bounds[1]],  # lower bounds
        [vaccination_bounds[2]]   # upper bounds
    )

    # Configure local method
    local_method = MultistartOptimization.NLopt_local_method(
        local_algorithm;
        xtol_rel = xtol_rel,
        xtol_abs = xtol_abs,
        maxeval = maxeval
    )

    # Configure multistart method
    multistart_method = MultistartOptimization.TikTak(n_sobol_points)

    if verbose
        println("Running multistart optimization...")
    end

    # Run optimization
    result = MultistartOptimization.multistart_minimization(
        multistart_method,
        local_method,
        problem
    )

    # Extract results
    optimal_vaccination = result.location[1]
    achieved_noise = calculate_mean_dynamical_noise(
        DynamicalNoise(noise_spec, optimal_vaccination),
        ensemble_spec,
        base_dynamics;
        seed = seed
    )
    difference = abs(achieved_noise - target_noise)

    # Check convergence
    if !in(result.ret, [:SUCCESS, :XTOL_REACHED, :FTOL_REACHED, :STOPVAL_REACHED]) ||
       difference > atol
        error(
            "Optimization failed or did not converge.\n" *
            "Return code: $(result.ret)\n" *
            "Optimal vaccination: $(optimal_vaccination)\n" *
            "Target noise: $(target_noise)\n" *
            "Achieved noise: $(achieved_noise)\n" *
            "Difference: $(difference)\n" *
            "Tolerance: $(atol)"
        )
    end

    if verbose
        println("\nOptimization successful!")
        println("Optimal vaccination coverage: $(round(optimal_vaccination; digits=4))")
        println("Target noise: $(round(target_noise; digits=2))")
        println("Achieved noise: $(round(achieved_noise; digits=2))")
        println("Difference: $(round(difference; digits=4))")
    end

    return (
        optimal_vaccination = optimal_vaccination,
        mean_noise = achieved_noise,
        target_noise = target_noise,
        difference = difference,
        optimization_result = result
    )
end
```

**Dependencies to add:**
```toml
MultistartOptimization = "0.2"
NLopt = "1.0"
```

**Testing Requirements:**
- Test optimization converges for known cases
- Test with different target scalings
- Test error handling for impossible targets
- Test verbose output
- Benchmark optimization performance

---

### 2.4 Update Noise Generation

**Priority: MEDIUM**
**Files to modify:**
- `OutbreakDetectionCore/src/noise/noise-generation.jl`

**Changes:**

#### Step 2.4.1: Add Dynamical Noise Support
Update existing `create_noise_arr` to handle both Poisson and Dynamical noise:

```julia
# Dispatch on noise type using sum types
function create_noise_arr(
    noise_spec::NoiseSpecification,
    incarr;
    ensemble_spec::EnsembleSpecification,
    base_dynamics::DynamicsParameterSpecification,
    seed::Int = 1234
)
    return create_noise_arr(
        LightSumTypes.variant(noise_spec),
        incarr;
        ensemble_spec = ensemble_spec,
        base_dynamics = base_dynamics,
        seed = seed
    )
end

# Poisson noise (existing implementation)
function create_noise_arr(
    noise_spec::PoissonNoise,
    incarr;
    seed::Int = 1234,
    kwargs...
)
    noise_arr = zeros(Int64, size(incarr))
    add_poisson_noise_arr!(
        noise_arr,
        incarr,
        noise_spec.noise_mean_scaling;
        seed = seed
    )

    mean_poisson_noise = mean(noise_arr)

    return noise_arr, (
        mean_noise = mean_poisson_noise,
        mean_poisson_noise = mean_poisson_noise,
        poisson_noise_prop = 1.0
    )
end

# Dynamical noise (new implementation)
function create_noise_arr(
    noise_spec::DynamicalNoise,
    incarr;
    ensemble_spec::EnsembleSpecification,
    base_dynamics::DynamicsParameterSpecification,
    seed::Int = 1234
)
    # Create noise dynamics
    noise_dynamics = create_noise_dynamics_parameters(
        noise_spec,
        base_dynamics,
        ensemble_spec.state_parameters.init_states.N
    )

    # Run noise simulation
    noise_result = recreate_noise_vecs(
        noise_spec,
        ensemble_spec,
        base_dynamics;
        seed = seed
    )

    return noise_result.noise_vecs, (
        mean_noise = noise_result.mean_noise,
        mean_dynamical_noise = noise_result.mean_dynamical_noise,
        mean_poisson_noise = noise_result.mean_poisson_noise,
        poisson_noise_prop = noise_result.mean_poisson_noise / noise_result.mean_noise
    )
end
```

**Testing Requirements:**
- Test both Poisson and Dynamical noise generation
- Test noise statistics match expected values
- Test correlation types produce expected patterns

---
---

## Phase 3: Integration and Migration

### 3.1 Update Existing Code

**Priority: HIGH**
**Files to update:** All files that use old type constructors

**Changes:**

#### Step 3.1.1: Create Migration Script
Create `scripts/migrate-to-new-types.jl` to help identify usage:

```julia
# Find all uses of old constructors
using Grep

println("Finding DynamicsParameters constructors...")
files = grep("DynamicsParameters(", "OutbreakDetectionCore/src", recursive=true)
println(files)

println("\nFinding constant usage...")
files = grep("BETA_MEAN|MU|EPSILON|SIGMA|GAMMA", "OutbreakDetectionCore/src", recursive=true)
println(files)
```

#### Step 3.1.2: Update Simulation Functions
Update functions that create dynamics parameters:

```julia
# Old way (using constants)
function old_simulation()
    dynamics = DynamicsParameters(SIGMA, GAMMA, R0)
    # ...
end

# New way (explicit parameters)
function new_simulation(;
    R_0 = 16.0,
    latent_period_days = 10.0,
    infectious_duration_days = 8.0,
    beta_force = 0.2
)
    target = TargetDiseaseDynamicsParameters(
        R_0 = R_0,
        latent_period_days = latent_period_days,
        infectious_duration_days = infectious_duration_days,
        beta_force = beta_force
    )
    spec = DynamicsParameterSpecification(target)
    dynamics = DynamicsParameters(spec; vaccination_coverage = 0.8)
    # ...
end
```

#### Step 3.1.3: Add Deprecation Warnings
For backward compatibility, add deprecated constructors:

```julia
# In dynamics-parameters.jl
function DynamicsParameters(sigma::Float64, gamma::Float64, R_0::Float64; kwargs...)
    @warn "DynamicsParameters(sigma, gamma, R_0) is deprecated. " *
          "Use TargetDiseaseDynamicsParameters → DynamicsParameterSpecification → DynamicsParameters"

    # Provide default values for migration
    target = TargetDiseaseDynamicsParameters(
        R_0 = R_0,
        latent_period_days = 1.0 / sigma,
        infectious_duration_days = 1.0 / gamma,
        beta_force = 0.2  # default
    )
    spec = DynamicsParameterSpecification(target)
    return DynamicsParameters(spec; kwargs...)
end
```

### 3.2 Update Scripts

**Priority: MEDIUM**
**Files to update:** All scripts in `scripts/`

**Changes:**

Update each script to use new type system:

```julia
# scripts/optimal-thresholds.jl (example)
using DrWatson
@quickactivate "OutbreakDetection"

# Old way
# dynamics = DynamicsParameters(SIGMA, GAMMA, R0)

# New way
target_dynamics = TargetDiseaseDynamicsParameters(
    R_0 = 16.0,
    latent_period_days = 10.0,
    infectious_duration_days = 8.0,
    beta_force = 0.2,
    life_expectancy_years = 62.5,
    population_N = 500_000
)

dynamics_spec = DynamicsParameterSpecification(target_dynamics)
dynamics = DynamicsParameters(dynamics_spec; vaccination_coverage = 0.8)

# Rest of script...
```

### 3.3 Update Module Exports

**Priority: HIGH**
**Files to modify:**
- `OutbreakDetectionCore/src/OutbreakDetectionCore.jl`

**Changes:**

```julia
# Add new exports
export TargetDiseaseDynamicsParameters,
       DynamicsParameterSpecification,
       DynamicsParameters,
       SeasonalityFunction,
       CosineSeasonality,
       SineSeasonality,
       DynamicalNoiseSpecification,
       DynamicalNoise,
       PoissonNoise,
       NoiseSpecification,
       NoiseVaccinationOptimizationParameters,
       SEIRRun,
       DynamicalNoiseRun,
       PoissonNoiseRun,
       NoiseRun,
       optimize_dynamic_noise_params,
       calculate_mean_dynamical_noise,
       recreate_noise_vecs,
       create_noise_dynamics_parameters,
       calculate_mean_incidence,
       calculate_endemic_equilibrium_proportions

# Include new files
include("types/simulation-results.jl")
include("simulation/endemic-equilibrium.jl")
include("noise/noise-mean-incidence.jl")
include("noise/noise-recreation.jl")
include("noise/noise-dynamics-parameters.jl")
include("noise/noise-parameters-optimization.jl")
include("types/optimization-parameters.jl")

# Optional: Only include if you need vaccination sampling for other purposes
# include("vaccination/vaccination-sampling.jl")
# export sample_vaccination_coverage
```

---

## Phase 4: Dependencies and Project Setup

### 4.1 Update Project.toml

**Priority: CRITICAL**
**File to modify:** `OutbreakDetectionCore/Project.toml`

**Changes:**

```toml
[deps]
# Existing dependencies...
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
LightSumTypes = "0.1"
MultistartOptimization = "0.2"
NLopt = "1.0"
Try = "0.1"  # For error handling
StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"  # For StructVector support
StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"  # For efficient state vectors

[compat]
julia = "1.9"
Distributions = "0.25"
LightSumTypes = "0.1"
MultistartOptimization = "0.2"
NLopt = "1.0"
Try = "0.1"
StructArrays = "0.6"
StaticArrays = "1.9"
```

### 4.2 Update Manifest

**Priority: CRITICAL**

Run after updating Project.toml:
```bash
cd OutbreakDetectionCore
julia --project=. -e 'using Pkg; Pkg.resolve(); Pkg.instantiate()'
```

---

## Phase 5: Documentation and Testing (Optional - Can Be Done Later)
## Implementation Timeline

### Week 1: Core Type System
- [ ] Day 1-2: Implement dynamics parameters hierarchy with seasonality sum types (1.1)
- [ ] Day 2-3: Implement time parameters simplification (1.2)
- [ ] Day 3-4: Implement noise specifications refactoring (1.3)
- [ ] Day 4: Implement test specifications enhancement (1.4)
- [ ] Day 5: Implement ensemble specifications simplification (1.5)

### Week 2: StructVectors and Simulation Infrastructure
- [ ] Day 1-2: Add simulation result types with StructVector support (1.6)
- [ ] Day 2-3: Add endemic equilibrium calculation (1.7)
- [ ] Day 3-4: Implement noise mean calculation with StructVectors (2.2.1)
- [ ] Day 4-5: Implement noise dynamics parameter creation (2.2.2)

### Week 3: Noise Optimization
- [ ] Day 1-2: Implement noise recreation with endemic equilibrium (2.2.3)
- [ ] Day 3-4: Implement vaccination optimization (2.3)
- [ ] Day 5: Update noise generation to use sum types (2.4)

### Week 4: Integration and Migration
- [ ] Day 1-2: Update existing code to use new types (3.1)
- [ ] Day 3: Update scripts (3.2)
- [ ] Day 4: Update module exports (3.3)
- [ ] Day 5: Update dependencies and resolve (4.1, 4.2)

### Week 5: Testing and Validation
- [ ] Day 1-2: Integration testing with StructVectors
- [ ] Day 3-4: Test optimization convergence
- [ ] Day 5: Validate endemic equilibrium calculations

### Week 6: Documentation (Optional - can be done later)
- [ ] Day 1-2: Comprehensive testing (5.2)
- [ ] Day 3-4: Documentation (5.1)
- [ ] Day 5: Final polish and examples

---

## Success Criteria

### Functionality
- [ ] All new types construct correctly with validation
- [ ] Vaccination optimization converges for test cases
- [ ] Noise generation produces expected distributions
- [ ] All correlation types work correctly

### Code Quality
- [ ] All tests pass (including new tests)
- [ ] Code coverage > 80% for new code
- [ ] JET.jl reports no type instabilities
- [ ] Aqua.jl quality checks pass
- [ ] All code formatted with Runic.jl

### Documentation
- [ ] All public functions have docstrings with examples
- [ ] Usage guide complete
- [ ] Migration guide for existing code
- [ ] Examples in documentation run without errors

### Performance
- [ ] Optimization completes in reasonable time (< 5 min for typical case)
- [ ] No performance regressions in existing code
- [ ] Memory usage reasonable for large ensembles

---

## Rollback Plan

If issues arise:

1. **Immediate rollback**: Revert to previous commit
   ```bash
   jj undo
   ```

2. **Partial rollback**: Keep completed phases, revert problematic ones
   - Each phase is independent
   - Can keep Phase 1 (types) while reverting Phase 2 (optimization)

3. **Deprecation period**: Keep old constructors with warnings for 1-2 releases

---

## Notes

- **Backward Compatibility**: Deprecated constructors provide migration path
- **Testing Strategy**: Test each phase independently before moving to next
- **Documentation**: Write docs alongside code, not after
- **Performance**: Profile optimization functions to ensure reasonable runtime
- **Code Review**: Review each phase before proceeding to next

---

## Design Decisions (Resolved)

Based on the CSDNoise implementation, the following design decisions have been made:

1. **SEIR Result Structure**: Use `StructVector{SEIRRun}` for efficient ensemble storage
   - Implemented in Section 1.6

2. **Seasonality Functions**: Use `LightSumTypes.@sumtype` for type-safe seasonality
   - `SeasonalityFunction(CosineSeasonality())` or `SeasonalityFunction(SineSeasonality())`
   - Implemented in Section 1.1

3. **StructVectors**: Fully integrated throughout the codebase
   - All result types designed for StructVector compatibility
   - No `AutoHashEquals` (not compatible with StructVectors)

4. **Initial States for Noise**: Use endemic equilibrium calculation
   - Implemented in Section 1.7
   - Falls back to simple proportions if equilibrium doesn't exist (R_eff ≤ 1)

5. **Ensemble Simulation**: Copy interface from CSDNoise
   - Returns `StructVector{SEIRRun}`
   - Uses same parameter structure

---

## References

- CSDNoise repository structure
- Kent Beck's "Tidy First" principles
- Julia Performance Tips: https://docs.julialang.org/en/v1/manual/performance-tips/
- LightSumTypes.jl: https://github.com/JuliaDiff/LightSumTypes.jl
- MultistartOptimization.jl: https://github.com/tpapp/MultistartOptimization.jl

### 5.1 Comprehensive Documentation

**Priority: MEDIUM**
**Files to update:** All modified type files

**Changes:**

#### Step 5.1.1: Add Detailed Docstrings
For each struct and function:
- Clear description of purpose
- Field/parameter descriptions with types and constraints
- Return value description
- At least one usage example
- Cross-references to related types/functions

Example template:
```julia
"""
    StructName

Brief one-line description.

Detailed description explaining purpose, use cases, and any important
behavioral notes.

# Fields
- `field1::Type`: Description including valid range/constraints
- `field2::Type`: Description

# Constructors
    StructName(; field1, field2)
    StructName(other_struct::OtherType)

# Examples
```julia
# Example 1: Basic usage
obj = StructName(field1 = value1, field2 = value2)

# Example 2: From other struct
obj = StructName(other_obj)
```

# See Also
- [`RelatedType`](@ref): Brief description of relationship
- [`related_function`](@ref): Brief description
"""
```

#### Step 5.1.2: Create Usage Guide
Create `OutbreakDetectionCore/docs/noise-optimization-guide.md`:
- Overview of dynamical noise
- Step-by-step workflow
- Parameter selection guidance
- Troubleshooting common issues
- Performance tips

### 5.2 Comprehensive Testing

**Priority: MEDIUM**  
**New test files to create:**
- `OutbreakDetectionCore/test/type-validation.jl`
- `OutbreakDetectionCore/test/noise-mean-incidence.jl`
- `OutbreakDetectionCore/test/noise-optimization.jl`

**Changes:**

#### Step 5.2.1: Type Validation Tests
```julia
# OutbreakDetectionCore/test/type-validation.jl
@testset "Type Validation" begin
    @testset "DynamicsParameterSpecification" begin
        # Valid construction
        target = TargetDiseaseDynamicsParameters(
            R_0 = 16.0,
            latent_period_days = 10.0,
            infectious_duration_days = 8.0,
            beta_force = 0.2,
            life_expectancy_years = 62.5,
            population_N = 500_000
        )
        spec = DynamicsParameterSpecification(target)
        @test spec.R_0 == 16.0
        @test spec.sigma ≈ 0.1
        @test spec.gamma ≈ 0.125

        # Invalid R_0
        @test_throws AssertionError TargetDiseaseDynamicsParameters(
            R_0 = -1.0,
            latent_period_days = 10.0,
            infectious_duration_days = 8.0,
            beta_force = 0.2
        )

        # Invalid beta_force
        @test_throws AssertionError TargetDiseaseDynamicsParameters(
            R_0 = 16.0,
            latent_period_days = 10.0,
            infectious_duration_days = 8.0,
            beta_force = 1.5
        )
    end

    @testset "NoiseSpecification" begin
        # Valid dynamical noise
        spec = DynamicalNoiseSpecification(
            R_0 = 5.0,
            latent_period = 7.0,
            duration_infection = 14.0,
            correlation = "in-phase",
            poisson_component = 1.0,
            vaccination_bounds = [0.0, 0.8]
        )
        @test spec.R_0 == 5.0

        # Invalid correlation
        @test_throws AssertionError DynamicalNoiseSpecification(
            R_0 = 5.0,
            latent_period = 7.0,
            duration_infection = 14.0,
            correlation = "invalid",
            poisson_component = 1.0
        )

        # Invalid vaccination bounds
        @test_throws AssertionError DynamicalNoiseSpecification(
            R_0 = 5.0,
            latent_period = 7.0,
            duration_infection = 14.0,
            correlation = "in-phase",
            poisson_component = 1.0,
            vaccination_bounds = [0.8, 0.2]  # reversed
        )
    end
end
```

#### Step 5.2.2: Mean Incidence Calculation Tests
```julia
# OutbreakDetectionCore/test/noise-mean-incidence.jl
@testset "Mean Incidence Calculation" begin
    @testset "calculate_mean_incidence" begin
        # Test with mock SEIR results
        # This will depend on your actual SEIR result structure
        # Placeholder test:
        mock_results = [
            (incidence = [1.0, 2.0, 3.0],),
            (incidence = [2.0, 3.0, 4.0],),
        ]
        
        mean_inc = calculate_mean_incidence(mock_results)
        @test mean_inc isa Float64
        @test mean_inc > 0.0
    end
end
```

#### Step 5.2.3: Noise Optimization Tests
```julia
# OutbreakDetectionCore/test/noise-optimization.jl
@testset "Noise Optimization" begin
    @testset "optimize_dynamic_noise_params" begin
        # Setup test parameters
        target_params = TargetDiseaseDynamicsParameters(
            R_0 = 16.0,
            latent_period_days = 10.0,
            infectious_duration_days = 8.0,
            beta_force = 0.2
        )
        base_dynamics = DynamicsParameterSpecification(target_params)

        noise_spec = DynamicalNoiseSpecification(
            R_0 = 5.0,
            latent_period = 7.0,
            duration_infection = 14.0,
            correlation = "in-phase",
            poisson_component = 1.0,
            vaccination_bounds = [0.0, 0.8]
        )

        # Create minimal ensemble spec for testing
        state_params = StateParameters(500_000, Dict(:s_prop => 0.05, :r_prop => 0.95))
        time_params = SimTimeParameters(tmax = 365.0 * 20)
        dynamics_params = DynamicsParameters(base_dynamics; vaccination_coverage = 0.8)
        ensemble_spec = EnsembleSpecification(
            (:test,),
            state_params,
            dynamics_params,
            time_params,
            100  # small nsims for testing
        )

        # Run optimization (with small tolerance for speed)
        opt_params = NoiseVaccinationOptimizationParameters(
            n_sobol_points = 10,
            maxeval = 100,
            atol = 0.5  # relaxed for testing
        )

        result = optimize_dynamic_noise_params(
            7.0,   # target_scaling
            3.0,   # mean_incidence
            noise_spec,
            ensemble_spec,
            base_dynamics,
            opt_params;
            verbose = false
        )

        # Check result structure
        @test haskey(result, :optimal_vaccination)
        @test haskey(result, :mean_noise)
        @test haskey(result, :target_noise)
        @test haskey(result, :difference)

        # Check values are reasonable
        @test 0.0 <= result.optimal_vaccination <= 0.8
        @test result.target_noise ≈ 21.0
        @test result.difference < opt_params.atol
    end
end
```

#### Step 5.2.4: Type Stability Tests
```julia
# Add to OutbreakDetectionCore/test/runtests.jl
using JET

@testset "Type Stability" begin
    # Test key functions for type stability
    target = TargetDiseaseDynamicsParameters(
        R_0 = 16.0,
        latent_period_days = 10.0,
        infectious_duration_days = 8.0,
        beta_force = 0.2
    )

    @test_opt DynamicsParameterSpecification(target)
    @test_call DynamicsParameterSpecification(target)

    spec = DynamicsParameterSpecification(target)
    @test_opt DynamicsParameters(spec; vaccination_coverage = 0.8)
    @test_call DynamicsParameters(spec; vaccination_coverage = 0.8)
end
```

**Testing Requirements:**
- All tests must pass before merging
- Code coverage > 80% for new code
- JET.jl reports no type instabilities
- Aqua.jl quality checks pass

---
