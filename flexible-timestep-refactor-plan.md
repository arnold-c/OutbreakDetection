# Flexible Timestep Refactor Plan

**Goal:** Enable the tau-leaping simulation to support arbitrary timesteps (tstep ≠ 1 day) while maintaining correctness and clarity.

**Current State:** The simulation assumes `tstep = 1.0` (daily timesteps). While `tstep` is configurable in `SimTimeParameters`, many parts of the codebase have hardcoded assumptions about daily timesteps, ambiguous variable names, and parameters that don't scale properly.

**Target State:** Support any positive `tstep` value (e.g., 0.5 for 12-hour, 1/24 for hourly, 7.0 for weekly) with all time-related parameters specified in clear units (using `Dates.Day`) that automatically convert to the appropriate number of timesteps.

---

## Design Decisions

### 1. Seasonality
**Decision:** Pass `tstep` to seasonality functions and use `2π * t * tstep / 365`
- Modify `_calculate_beta_amp` to accept `tstep` parameter
- `t` represents timestep index (1, 2, 3, ...), `t * tstep` gives time in days
- Ensures seasonal peaks occur at correct calendar times regardless of timestep

### 2. Test Result Lag
**Decision:** Use `Dates.Day` period type, convert to timesteps internally
- Store as `test_result_lag::Dates.Day` (e.g., `Dates.Day(2)`)
- Convert to timesteps: `round(Int, Dates.value(lag) / tstep)`
- Clear semantics: 2-day lag is always 2 days, regardless of timestep

### 3. Detection Parameters (Windows, Durations)
**Decision:** Use `Dates.Day` period types
- `moving_average_lag::Dates.Day` instead of `::Int64`
- `min_outbreak_duration::Dates.Day` instead of `::Int64`
- Convert to timesteps when needed using helper function

### 4. Backward Compatibility
**Decision:** Major refactor acceptable for cleaner design
- Breaking changes to type signatures are acceptable
- Will update all scripts and tests
- Provide clear migration guide

---

## Time Conversion Utilities

### Key Functions to Add

```julia
# In OutbreakDetectionUtils/src/time-utils.jl (new file)

"""
    to_days(period::Dates.Period)::Float64

Convert a Dates period to days as Float64.
"""
to_days(period::Dates.Day) = Float64(Dates.value(period))
to_days(period::Dates.Hour) = Dates.value(period) / 24.0
to_days(period::Dates.Minute) = Dates.value(period) / (24.0 * 60.0)
to_days(period::Dates.Second) = Dates.value(period) / (24.0 * 60.0 * 60.0)

"""
    to_timesteps(period::Dates.Period, tstep::Float64)::Int64

Convert a time period to number of simulation timesteps.
"""
function to_timesteps(period::Dates.Period, tstep::Float64)::Int64
    period_days = to_days(period)
    return round(Int, period_days / tstep)
end
```

### Usage Pattern

```julia
# Store parameters in Dates.Day
test_result_lag::Dates.Day = Dates.Day(2)
moving_average_window::Dates.Day = Dates.Day(7)

# Convert to timesteps when needed
lag_timesteps = to_timesteps(test_result_lag, time_params.tstep)
window_timesteps = to_timesteps(moving_average_window, time_params.tstep)
```

---

## Implementation Phases

### Phase 1: Foundation - Time Utilities and Type Updates

#### 1.1 Create Time Conversion Utilities
**File:** `OutbreakDetectionCore/src/utils/time-utils.jl` (new)

```julia
export to_days, to_timesteps

using Dates

# Add the conversion functions shown above
# Include comprehensive docstrings with examples
```

**File:** `OutbreakDetectionCore/src/OutbreakDetectionCore.jl`
- Add `include("utils/time-utils.jl")`

#### 1.2 Update `SimTimeParameters`
**File:** `OutbreakDetectionCore/src/types/time-parameters.jl`

**Changes:**
- Update documentation to clarify `tstep` can be any positive value
- Add examples showing different timesteps (0.5, 1/24, 7.0)
- Document relationship between `tstep` and calendar time

**Example documentation:**
```julia
# Examples
```julia
# Daily timesteps
time_params = SimTimeParameters(tmax = 365.0, tstep = 1.0)

# 12-hour timesteps (twice daily)
time_params = SimTimeParameters(tmax = 365.0, tstep = 0.5)

# Hourly timesteps
time_params = SimTimeParameters(tmax = 365.0, tstep = 1/24)

# Weekly timesteps
time_params = SimTimeParameters(tmax = 365.0, tstep = 7.0)
```
```

#### 1.3 Update `IndividualTestSpecification`
**File:** `OutbreakDetectionCore/src/types/test-specifications.jl`

**Current:**
```julia
struct IndividualTestSpecification
    sensitivity::Float64
    specificity::Float64
    test_result_lag::Integer  # Ambiguous units!
end
```

**Change to:**
```julia
struct IndividualTestSpecification
    sensitivity::Float64
    specificity::Float64
    test_result_lag::Dates.Day  # Clear: always in days
end

# Backward-compatible constructor
IndividualTestSpecification(sens::Float64, spec::Float64, lag::Integer) = 
    IndividualTestSpecification(sens, spec, Dates.Day(lag))
```

**Update all related functions:**
- Update docstrings to reflect `Dates.Day` usage
- Update string representations (e.g., `"$(lag) day lag"` → `"$(Dates.value(lag)) day lag"`)

#### 1.4 Update `DetectionSpecification` Types
**File:** `OutbreakDetectionCore/src/types/detection-specifications.jl`

**Changes needed:**

1. **`DailyThreshold` struct:**
```julia
# Current
struct DailyThreshold <: DetectionMethod
    threshold::Float64
end

# Consider renaming to InstantaneousThreshold or PerTimestepThreshold
# But "Daily" might be okay if it refers to the threshold being applied
# to daily-aggregated data. Review usage and decide.
```

2. **`DailyMovingAverageThreshold` struct:**
```julia
# Current
struct DailyMovingAverageThreshold <: DetectionMethod
    threshold::Float64
    moving_average_lag::Int64  # Default: 7
end

# Change to
struct DailyMovingAverageThreshold <: DetectionMethod
    threshold::Float64
    moving_average_window::Dates.Day  # Renamed for clarity
end

# Constructor with default
DailyMovingAverageThreshold(threshold::Float64) = 
    DailyMovingAverageThreshold(threshold, Dates.Day(7))

# Backward-compatible constructor
DailyMovingAverageThreshold(threshold::Float64, window::Integer) = 
    DailyMovingAverageThreshold(threshold, Dates.Day(window))
```

3. **Update all detection specification constructors:**
- Change all `::Int64` duration/window parameters to `::Dates.Day`
- Add backward-compatible constructors accepting `Integer`
- Update default values from `7` to `Dates.Day(7)`

#### 1.5 Update `NoiseSpecification`
**File:** `OutbreakDetectionCore/src/types/noise-specifications.jl`

**Changes:**
```julia
# Current
struct NoiseSpecification
    # ...
    latent_period::Float64  # Default: 7.0
    # ...
end

# Change to
struct NoiseSpecification
    # ...
    latent_period::Dates.Day  # Default: Dates.Day(7)
    # ...
end
```

Update all constructors and default values.

---

### Phase 2: Seasonality Functions

#### 2.1 Update `_calculate_beta_amp`
**File:** `OutbreakDetectionCore/src/simulation/transmission-functions.jl`

**Current:**
```julia
_calculate_beta_amp(beta_mean, beta_force, t, ::CosineSeasonality) =
    beta_mean * (1 + beta_force * cos(2π * t / 365))

_calculate_beta_amp(beta_mean, beta_force, t, ::SineSeasonality) =
    beta_mean * (1 + beta_force * sin(2π * t / 365))
```

**Change to:**
```julia
_calculate_beta_amp(beta_mean, beta_force, t, tstep, ::CosineSeasonality) =
    beta_mean * (1 + beta_force * cos(2π * t * tstep / 365))

_calculate_beta_amp(beta_mean, beta_force, t, tstep, ::SineSeasonality) =
    beta_mean * (1 + beta_force * sin(2π * t * tstep / 365))
```

**Rationale:**
- `t` is the timestep index (1, 2, 3, ...)
- `t * tstep` gives the time in days
- Seasonal period remains 365 days regardless of timestep

**Update docstrings:**
```julia
"""
    _calculate_beta_amp(beta_mean, beta_force, t, tstep, ::CosineSeasonality)

Calculate beta with cosine seasonality.

# Arguments
- `beta_mean`: Mean transmission rate
- `beta_force`: Amplitude of seasonal forcing (0 to 1)
- `t`: Timestep index (1, 2, 3, ...)
- `tstep`: Timestep size in days
- Seasonality type (CosineSeasonality)

Returns `beta_mean * (1 + beta_force * cos(2π * t * tstep / 365))`.
"""
```

#### 2.2 Update `calculate_beta_vec!`
**File:** `OutbreakDetectionCore/src/simulation/transmission-functions.jl`

**Add `tstep` parameter and pass to `_calculate_beta_amp`:**

```julia
function calculate_beta_vec!(
        beta_vec::AbstractVector{Float64},
        dynamics_params::DynamicsParameters,
        time_params::SimTimeParameters  # Already has tstep!
    )::Nothing
    
    beta_mean = dynamics_params.beta_mean
    beta_force = dynamics_params.beta_force
    seasonality = dynamics_params.seasonality
    tstep = time_params.tstep  # Extract tstep
    
    @inbounds for i in eachindex(beta_vec)
        beta_vec[i] = _calculate_beta_amp(
            beta_mean, beta_force, i, tstep, seasonality  # Pass tstep
        )
    end
    
    return nothing
end
```

**Note:** Check all call sites - should already be passing `time_params`, so no changes needed to callers.

---

### Phase 3: Diagnostic Testing Functions

#### 3.1 Update `calculate_positives_vec!`
**File:** `OutbreakDetectionCore/src/diagnostic-testing/calculate-num-positive.jl`

**Current:**
```julia
function calculate_positives_vec!(
        npos_vec::AbstractVector{Int64},
        tested_vec::AbstractVector{Int64},
        sim_length::Int64,
        lag::Int64,  # Ambiguous: days or timesteps?
        tested_multiplier::Float64
    )
    @inbounds for day in eachindex(npos_vec)
        result_day = day + lag
        if result_day <= sim_length
            npos_vec[result_day] = rand(Distributions.Binomial(tested_vec[day], tested_multiplier))
        end
    end
    return nothing
end
```

**Change to:**
```julia
function calculate_positives_vec!(
        npos_vec::AbstractVector{Int64},
        tested_vec::AbstractVector{Int64},
        sim_length::Int64,
        lag_timesteps::Int64,  # Clear: in timestep units
        tested_multiplier::Float64
    )
    @inbounds for timestep_idx in eachindex(npos_vec)
        result_timestep = timestep_idx + lag_timesteps
        if result_timestep <= sim_length
            npos_vec[result_timestep] = rand(Distributions.Binomial(tested_vec[timestep_idx], tested_multiplier))
        end
    end
    return nothing
end
```

**Add wrapper accepting `Dates.Day`:**
```julia
function calculate_positives_vec!(
        npos_vec::AbstractVector{Int64},
        tested_vec::AbstractVector{Int64},
        sim_length::Int64,
        lag::Dates.Day,
        tstep::Float64,
        tested_multiplier::Float64
    )
    lag_timesteps = to_timesteps(lag, tstep)
    calculate_positives_vec!(npos_vec, tested_vec, sim_length, lag_timesteps, tested_multiplier)
    return nothing
end
```

**Update docstrings:**
```julia
"""
    calculate_positives_vec!(npos_vec, tested_vec, sim_length, lag_timesteps, tested_multiplier)

Calculate test-positive individuals accounting for test result lag.

# Arguments
- `npos_vec`: Output vector for positive test results (modified in-place)
- `tested_vec`: Input vector of individuals tested
- `sim_length`: Length of simulation
- `lag_timesteps`: Test result lag in timesteps
- `tested_multiplier`: Test performance multiplier (sensitivity or 1-specificity)

# Note
For the version accepting `Dates.Day`, use:
`calculate_positives_vec!(npos_vec, tested_vec, sim_length, Dates.Day(2), tstep, tested_multiplier)`
"""
```

#### 3.2 Update Related Functions
**Same file:** `OutbreakDetectionCore/src/diagnostic-testing/calculate-num-positive.jl`

Apply the same pattern to:
- `calculate_positives!` - rename `day` → `timestep_idx`
- `calculate_true_positives!` - update to accept `Dates.Day`
- `calculate_noise_positives!` - update to accept `Dates.Day`

#### 3.3 Update `calculate_num_tested`
**File:** `OutbreakDetectionCore/src/diagnostic-testing/calculate-num-tested.jl`

**Changes:**
- Line 12: Change "each day's incidence" → "each timestep's incidence"
- Update all docstrings to use "timestep" instead of "day"

#### 3.4 Update Test Vector Creation
**File:** `OutbreakDetectionCore/src/diagnostic-testing/create-test-positive-vectors.jl`

**Changes:**
- Update all documentation references to "day lag" → "lag (in days)"
- Update references to "daily incidence" → "per-timestep incidence"
- Update function calls to pass `tstep` parameter where needed
- Lines 41, 124, 183, 185: Update examples and docstrings

---

### Phase 4: Detection Functions

#### 4.1 Update `calculate_moving_average`
**File:** `OutbreakDetectionCore/src/detection/calculate-moving-average.jl`

**Current:**
```julia
function _calculate_float_daily_movingavg(invec, day, avglag)
    @inline moveavg_daystart = _calculate_daily_movingavg_startday(day, avglag)
    return StatsBase.mean(@view(invec[moveavg_daystart:day]))
end

function _calculate_daily_movingavg_startday(day, avglag)
    if day < avglag
        moveavg_daystart = 1
    else
        moveavg_daystart = day - avglag + 1
    end
    return moveavg_daystart
end
```

**Change to:**
```julia
function _calculate_float_movingavg(invec, timestep_idx, window_timesteps)
    @inline window_start = _calculate_movingavg_start_idx(timestep_idx, window_timesteps)
    return StatsBase.mean(@view(invec[window_start:timestep_idx]))
end

function _calculate_movingavg_start_idx(timestep_idx, window_timesteps)
    if timestep_idx < window_timesteps
        window_start = 1
    else
        window_start = timestep_idx - window_timesteps + 1
    end
    return window_start
end
```

**Update public function:**
```julia
function calculate_moving_average(
        invec::AbstractVector,
        window::Dates.Day,
        tstep::Float64
    )
    window_timesteps = to_timesteps(window, tstep)
    outvec = similar(invec, Float64)
    
    @inbounds for timestep_idx in eachindex(invec)
        outvec[timestep_idx] = _calculate_float_movingavg(invec, timestep_idx, window_timesteps)
    end
    
    return outvec
end
```

#### 4.2 Update `classify_outbreaks`
**File:** `OutbreakDetectionCore/src/detection/classify-outbreaks.jl`

**Current:**
```julia
function classify_outbreaks(
        incidence_vec::AbstractVector,
        minoutbreakdur::Int64
    )
    # ...
end
```

**Change to:**
```julia
function classify_outbreaks(
        incidence_vec::AbstractVector,
        min_outbreak_duration::Dates.Day,
        tstep::Float64
    )
    min_duration_timesteps = to_timesteps(min_outbreak_duration, tstep)
    # Use min_duration_timesteps in logic
    # ...
end
```

**Update docstrings:**
- Line 26: "daily incidence" → "per-timestep incidence"
- Line 28: "consecutive days" → "minimum duration (in days, converted to timesteps)"

#### 4.3 Update `calculate_outbreak_threshold`
**File:** `OutbreakDetectionCore/src/detection/outbreak-thresholds-calculation.jl`

**Changes:**
- Line 77: `minimum_outbreak_duration = 7` → `minimum_outbreak_duration = Dates.Day(7)`
- Line 60: Update docstring "Minimum daily incidence" → "Minimum per-timestep incidence"
- Line 141: Update docstring similarly
- Update all function signatures to accept `Dates.Day` for duration parameters

#### 4.4 Update Detection Metric Functions
**File:** `OutbreakDetectionCore/src/detection/detection-metric-functions.jl`

**Changes:**
- Line 34: Update docstring to clarify delays are in timesteps
- Lines 180-181: Update examples to be timestep-agnostic or note they assume tstep=1

---

### Phase 5: Global Renaming and Documentation

#### 5.1 Variable Renaming

**Search and replace throughout codebase:**

| Old Name | New Name | Context |
|----------|----------|---------|
| `day` | `timestep_idx` | When used as loop index over simulation time |
| `result_day` | `result_timestep` | Timestep index after applying lag |
| `moveavg_daystart` | `moveavg_start_idx` | Start index for moving average window |
| `avglag` | `window_timesteps` | Moving average window size in timesteps |
| `minoutbreakdur` | `min_outbreak_duration` | Minimum outbreak duration parameter |

**Important:** Keep "day" when actually referring to calendar days (e.g., `Dates.Day`, documentation about real-world time).

#### 5.2 Function Renaming

| Old Name | New Name | Reason |
|----------|----------|--------|
| `_calculate_float_daily_movingavg` | `_calculate_float_movingavg` | "daily" is misleading |
| `_calculate_daily_movingavg_startday` | `_calculate_movingavg_start_idx` | More accurate |

**Note:** Consider whether `DailyThreshold` and `DailyMovingAverageThreshold` should be renamed. If "Daily" refers to the data being daily-aggregated (as opposed to the threshold being applied daily), the names might be okay. Review usage context and decide.

#### 5.3 Documentation Updates

**Update all docstrings containing:**
- "daily incidence" → "per-timestep incidence"
- "day lag" → "lag in days (converted to timesteps internally)"
- "consecutive days" → "consecutive timesteps (specified in days)"
- "7-day moving average" → "moving average window (specified in days, e.g., Dates.Day(7))"
- "each day" → "each timestep"

**Files to review:**
- All files in `OutbreakDetectionCore/src/diagnostic-testing/`
- All files in `OutbreakDetectionCore/src/detection/`
- All files in `OutbreakDetectionCore/src/simulation/`
- `OutbreakDetectionCore/src/types/*.jl`

---

### Phase 6: Update Call Sites and Integration

#### 6.1 Update `seir_mod!`
**File:** `OutbreakDetectionCore/src/simulation/seir-model.jl`

**Verify:**
- `calculate_beta_vec!` call passes `time_params` (should already be correct)
- Rate-to-probability conversions work for larger timesteps
- Consider adding validation warning

**Add validation (optional but recommended):**
```julia
function seir_mod!(...)
    # ... existing code ...
    
    # Validate that rate * tstep is small enough for Euler approximation
    max_rate = max(mu_timestep, sigma_timestep, gamma_timestep, epsilon_timestep)
    if max_rate > 0.5
        @warn "Large timestep detected: max(rate * tstep) = $max_rate > 0.5. " *
              "Euler approximation may be inaccurate. Consider reducing tstep."
    end
    
    # ... rest of function ...
end
```

#### 6.2 Update Test Specification Constructors
**File:** `OutbreakDetectionCore/src/types/test-specifications.jl`

**Ensure all constructors accept both `Integer` and `Dates.Day`:**
```julia
# Example for IndividualTestSpecification
function IndividualTestSpecification(
        sensitivity::Float64,
        specificity::Float64,
        lag::Integer  # Backward compatibility
    )
    return IndividualTestSpecification(sensitivity, specificity, Dates.Day(lag))
end
```

#### 6.3 Update Detection Specification Constructors
**File:** `OutbreakDetectionCore/src/types/detection-specifications.jl`

**Update all constructors to:**
1. Accept `Dates.Day` for duration/window parameters
2. Provide backward-compatible constructors accepting `Integer`
3. Update default values from integers to `Dates.Day`

**Example:**
```julia
function DailyMovingAverageThreshold(
        threshold::Float64,
        window::Integer  # Backward compatibility
    )
    return DailyMovingAverageThreshold(threshold, Dates.Day(window))
end
```

#### 6.4 Update Ensemble Specifications
**File:** `OutbreakDetectionCore/src/types/ensemble-specifications.jl`

**Verify:**
- Line 155: Directory path generation handles different tsteps correctly (should already be fine)
- Any other references to time parameters are correct

---

### Phase 7: Testing and Validation

#### 7.1 Create New Test File
**File:** `OutbreakDetectionCore/test/flexible-timesteps.jl` (new)

**Test cases:**

```julia
using Test
using OutbreakDetectionCore
using Dates

@testset "Flexible Timesteps" begin
    @testset "Time Conversion Utilities" begin
        @test to_days(Dates.Day(2)) == 2.0
        @test to_days(Dates.Hour(12)) == 0.5
        @test to_days(Dates.Hour(6)) == 0.25
        @test to_days(Dates.Minute(720)) == 0.5
        
        @test to_timesteps(Dates.Day(2), 1.0) == 2
        @test to_timesteps(Dates.Day(2), 0.5) == 4
        @test to_timesteps(Dates.Day(7), 7.0) == 1
        @test to_timesteps(Dates.Hour(12), 0.5) == 1
    end
    
    @testset "Seasonality with Different Timesteps" begin
        # Test that seasonal peak occurs at same calendar time
        # regardless of timestep
        
        # Daily timesteps
        time_params_daily = SimTimeParameters(tmax = 365.0, tstep = 1.0)
        # ... create dynamics params with seasonality ...
        # ... run simulation ...
        # ... find peak time ...
        
        # Weekly timesteps
        time_params_weekly = SimTimeParameters(tmax = 365.0, tstep = 7.0)
        # ... run same simulation ...
        # ... verify peak occurs at same calendar day ...
    end
    
    @testset "Test Lag with Different Timesteps" begin
        # Verify 2-day lag is always 2 days regardless of timestep
        
        test_spec = IndividualTestSpecification(0.9, 0.95, Dates.Day(2))
        
        # Daily
        time_params = SimTimeParameters(tmax = 100.0, tstep = 1.0)
        lag_ts = to_timesteps(test_spec.test_result_lag, time_params.tstep)
        @test lag_ts == 2
        
        # 12-hour
        time_params = SimTimeParameters(tmax = 100.0, tstep = 0.5)
        lag_ts = to_timesteps(test_spec.test_result_lag, time_params.tstep)
        @test lag_ts == 4
        
        # Hourly
        time_params = SimTimeParameters(tmax = 100.0, tstep = 1/24)
        lag_ts = to_timesteps(test_spec.test_result_lag, time_params.tstep)
        @test lag_ts == 48
    end
    
    @testset "Moving Average Window" begin
        # Verify 7-day window covers same time period regardless of tstep
        
        window = Dates.Day(7)
        
        @test to_timesteps(window, 1.0) == 7
        @test to_timesteps(window, 0.5) == 14
        @test to_timesteps(window, 7.0) == 1
    end
    
    @testset "Full Simulation with Different Timesteps" begin
        # Run identical simulation with different timesteps
        # Verify results are comparable (accounting for discretization)
        
        # Setup common parameters
        # ...
        
        # Run with tstep = 1.0
        # ...
        
        # Run with tstep = 0.5
        # ...
        
        # Compare results (may need tolerance for discretization differences)
    end
end
```

#### 7.2 Update Existing Tests
**Files to update:**
- `OutbreakDetectionCore/test/SEIR-model.jl`
- All test files that create `IndividualTestSpecification`
- All test files that create detection specifications

**Changes:**
- Use `Dates.Day()` constructors instead of raw integers
- Add test cases with different timesteps where appropriate
- Verify backward compatibility (integer constructors still work)

#### 7.3 Add Validation Tests
**File:** `OutbreakDetectionCore/test/validation.jl` (new or add to existing)

**Test edge cases:**
```julia
@testset "Timestep Validation" begin
    @testset "Non-divisible periods" begin
        # Test warning when period not evenly divisible by tstep
        # e.g., Dates.Day(7) with tstep = 2.0 → 3.5 timesteps → rounds to 4
        
        window = Dates.Day(7)
        tstep = 2.0
        
        # Should warn but not error
        @test_logs (:warn,) to_timesteps(window, tstep)
        @test to_timesteps(window, tstep) == 4  # Rounds 3.5 → 4
    end
    
    @testset "Large timestep approximation" begin
        # Test that large tstep triggers warning in seir_mod!
        # when rate * tstep > 0.5
        
        # ... setup parameters with large tstep ...
        # ... verify warning is issued ...
    end
end
```

#### 7.4 Integration Tests
**Verify:**
- Scripts in `scripts/` run without errors
- Results are reasonable for different timesteps
- Performance is acceptable (smaller tstep = more iterations)

---

### Phase 8: Documentation and Migration Guide

#### 8.1 Update README
**File:** `README.md` or `OutbreakDetectionCore/README.md`

**Add section:**
```markdown
## Flexible Timesteps

The simulation supports arbitrary timestep sizes via the `tstep` parameter in `SimTimeParameters`:

- `tstep = 1.0` - Daily timesteps (default)
- `tstep = 0.5` - 12-hour timesteps
- `tstep = 0.25` - 6-hour timesteps
- `tstep = 1/24` - Hourly timesteps
- `tstep = 7.0` - Weekly timesteps

All time-related parameters (test lag, moving average windows, outbreak durations) are specified using `Dates.Day` and automatically converted to the appropriate number of timesteps.

See `docs/flexible-timesteps.md` for detailed usage guide.
```

#### 8.2 Create Migration Guide
**File:** `docs/flexible-timesteps.md` (new)

**Contents:**

```markdown
# Flexible Timesteps Guide

## Overview

This guide explains how to use flexible timesteps in the OutbreakDetection simulation framework.

## Specifying Timestep Size

The timestep size is controlled by the `tstep` parameter in `SimTimeParameters`:

```julia
# Daily timesteps (default)
time_params = SimTimeParameters(tmax = 365.0 * 10, tstep = 1.0)

# 12-hour timesteps
time_params = SimTimeParameters(tmax = 365.0 * 10, tstep = 0.5)

# Hourly timesteps
time_params = SimTimeParameters(tmax = 365.0 * 10, tstep = 1/24)

# Weekly timesteps
time_params = SimTimeParameters(tmax = 365.0 * 10, tstep = 7.0)
```

## Time-Related Parameters

All time-related parameters use `Dates.Day` for clarity:

```julia
# Test specifications
test_spec = IndividualTestSpecification(
    0.9,              # sensitivity
    0.95,             # specificity
    Dates.Day(2)      # 2-day lag
)

# Detection specifications
detection_spec = DailyMovingAverageThreshold(
    5.0,              # threshold
    Dates.Day(7)      # 7-day moving average window
)
```

These parameters automatically convert to the appropriate number of timesteps based on your `tstep` setting.

## How It Works

Internally, the simulation converts `Dates.Day` parameters to timesteps:

```julia
lag_timesteps = round(Int, Dates.value(lag) / tstep)
```

Examples:
- 2-day lag with `tstep=1.0` → 2 timesteps
- 2-day lag with `tstep=0.5` → 4 timesteps
- 2-day lag with `tstep=1/24` → 48 timesteps

## Choosing an Appropriate Timestep

### Considerations

1. **Accuracy**: Smaller timesteps provide more accurate approximations but require more computation
2. **Performance**: Timestep size directly affects simulation time (tstep=0.5 takes ~2x longer than tstep=1.0)
3. **Seasonality**: For best results, choose tstep that divides evenly into 365 days
4. **Data granularity**: Match your timestep to your data collection frequency

### Recommendations

- **Daily data**: Use `tstep = 1.0`
- **Sub-daily monitoring**: Use `tstep = 0.5` (12-hour) or `tstep = 0.25` (6-hour)
- **Weekly aggregation**: Use `tstep = 7.0`
- **High-resolution studies**: Use `tstep = 1/24` (hourly) but be aware of performance cost

### Limitations

- Very large timesteps (e.g., `tstep > 7.0`) may cause inaccuracies in the Euler approximation
- The simulation will warn if `rate * tstep > 0.5` (approximation breaks down)
- Non-divisor timesteps (e.g., `tstep = 3.5`) may cause seasonal peaks to drift slightly

## Migration from Previous Versions

If you have existing code using integer lag/window parameters:

### Old Code
```julia
test_spec = IndividualTestSpecification(0.9, 0.95, 2)  # Integer lag
detection_spec = DailyMovingAverageThreshold(5.0, 7)   # Integer window
```

### New Code (Recommended)
```julia
test_spec = IndividualTestSpecification(0.9, 0.95, Dates.Day(2))
detection_spec = DailyMovingAverageThreshold(5.0, Dates.Day(7))
```

**Note:** Backward-compatible constructors are provided, so old code will still work, but using `Dates.Day` is recommended for clarity.

## Examples

### Example 1: High-Resolution Outbreak Detection

```julia
# Hourly timesteps for detailed outbreak dynamics
time_params = SimTimeParameters(
    tmax = 365.0,
    tstep = 1/24  # Hourly
)

# 2-day test lag (48 timesteps)
test_spec = IndividualTestSpecification(0.9, 0.95, Dates.Day(2))

# 7-day moving average (168 timesteps)
detection_spec = DailyMovingAverageThreshold(5.0, Dates.Day(7))
```

### Example 2: Weekly Surveillance

```julia
# Weekly timesteps for long-term trends
time_params = SimTimeParameters(
    tmax = 365.0 * 20,  # 20 years
    tstep = 7.0         # Weekly
)

# 2-week test lag (2 timesteps)
test_spec = IndividualTestSpecification(0.9, 0.95, Dates.Day(14))

# 4-week moving average (4 timesteps)
detection_spec = DailyMovingAverageThreshold(10.0, Dates.Day(28))
```

## Troubleshooting

### Warning: "Large timestep detected"

If you see this warning, your timestep may be too large for accurate Euler approximation. Consider:
1. Reducing `tstep`
2. Adjusting disease dynamics parameters
3. Accepting the approximation error if it's acceptable for your use case

### Non-integer timestep conversions

If a period doesn't divide evenly into your timestep (e.g., 7-day window with tstep=2.0), the conversion will round to the nearest integer. A warning will be issued.

## Performance Tips

- Use the largest timestep that provides acceptable accuracy for your application
- Profile your code to identify bottlenecks before reducing timestep
- Consider using weekly timesteps for long-term (multi-year) simulations
- Use hourly/sub-daily timesteps only when necessary for your research question
```

#### 8.3 Update Example Scripts
**File:** `scripts/optimal-thresholds_optims.jl`

**Add examples with different timesteps:**
```julia
# Example 1: Daily timesteps (original)
time_params_daily = SimTimeParameters(
    tmin = 0.0,
    tmax = 365.0 * 10,
    tstep = 1.0
)

# Example 2: 12-hour timesteps
time_params_12h = SimTimeParameters(
    tmin = 0.0,
    tmax = 365.0 * 10,
    tstep = 0.5
)

# Example 3: Weekly timesteps
time_params_weekly = SimTimeParameters(
    tmin = 0.0,
    tmax = 365.0 * 10,
    tstep = 7.0
)
```

#### 8.4 Add Docstring Examples
**Update key functions with timestep examples:**

```julia
"""
    SimTimeParameters(; tmin = 0.0, tmax = 365.0 * 100.0, tstep = 1.0)

# Examples with Different Timesteps
```julia
# Daily timesteps (default)
params = SimTimeParameters(tmax = 365.0 * 10, tstep = 1.0)

# 12-hour timesteps (twice daily)
params = SimTimeParameters(tmax = 365.0 * 10, tstep = 0.5)

# 6-hour timesteps
params = SimTimeParameters(tmax = 365.0 * 10, tstep = 0.25)

# Hourly timesteps
params = SimTimeParameters(tmax = 365.0 * 10, tstep = 1/24)

# Weekly timesteps
params = SimTimeParameters(tmax = 365.0 * 10, tstep = 7.0)
```
"""
```

---

## Implementation Checklist

### Phase 1: Foundation
- [ ] Create `OutbreakDetectionCore/src/utils/time-utils.jl`
- [ ] Add `to_days()` function with multiple dispatch
- [ ] Add `to_timesteps()` function
- [ ] Update `SimTimeParameters` documentation
- [ ] Update `IndividualTestSpecification` to use `Dates.Day`
- [ ] Update `DetectionSpecification` types to use `Dates.Day`
- [ ] Update `NoiseSpecification` to use `Dates.Day`
- [ ] Add backward-compatible constructors for all types

### Phase 2: Seasonality
- [ ] Update `_calculate_beta_amp` (cosine) to accept `tstep`
- [ ] Update `_calculate_beta_amp` (sine) to accept `tstep`
- [ ] Update `calculate_beta_vec!` to pass `tstep`
- [ ] Update docstrings
- [ ] Verify all call sites

### Phase 3: Diagnostic Testing
- [ ] Update `calculate_positives_vec!` - rename variables
- [ ] Add `Dates.Day` wrapper for `calculate_positives_vec!`
- [ ] Update `calculate_positives!`
- [ ] Update `calculate_true_positives!`
- [ ] Update `calculate_noise_positives!`
- [ ] Update `calculate_num_tested` documentation
- [ ] Update `create-test-positive-vectors.jl` documentation

### Phase 4: Detection
- [ ] Rename `_calculate_float_daily_movingavg` → `_calculate_float_movingavg`
- [ ] Rename `_calculate_daily_movingavg_startday` → `_calculate_movingavg_start_idx`
- [ ] Update variable names in moving average functions
- [ ] Add `Dates.Day` wrapper for `calculate_moving_average`
- [ ] Update `classify_outbreaks` to accept `Dates.Day`
- [ ] Update `calculate_outbreak_threshold` defaults
- [ ] Update detection metric function documentation

### Phase 5: Documentation
- [ ] Global search/replace: `day` → `timestep_idx` (in loop contexts)
- [ ] Update all "daily incidence" → "per-timestep incidence"
- [ ] Update all "day lag" → "lag in days"
- [ ] Update all "consecutive days" → "consecutive timesteps"
- [ ] Review and update function names if needed

### Phase 6: Integration
- [ ] Verify `seir_mod!` works with updated `calculate_beta_vec!`
- [ ] Add validation warning for large timesteps (optional)
- [ ] Update all test specification constructors
- [ ] Update all detection specification constructors
- [ ] Verify ensemble specifications handle different tsteps

### Phase 7: Testing
- [ ] Create `test/flexible-timesteps.jl`
- [ ] Add time conversion utility tests
- [ ] Add seasonality tests with different timesteps
- [ ] Add test lag tests with different timesteps
- [ ] Add moving average window tests
- [ ] Add full simulation comparison tests
- [ ] Update existing test files to use `Dates.Day`
- [ ] Add validation tests for edge cases
- [ ] Run full test suite: `just tests`
- [ ] Verify all tests pass

### Phase 8: Documentation
- [ ] Update README with flexible timestep section
- [ ] Create `docs/flexible-timesteps.md` migration guide
- [ ] Update example scripts with different timesteps
- [ ] Add timestep examples to key docstrings
- [ ] Update AGENTS.md if needed

### Final Steps
- [ ] Run `runic -i` on all modified files
- [ ] Run full test suite one more time
- [ ] Build manuscript to verify no breakage: `just manuscript`
- [ ] Create comprehensive commit message
- [ ] Update CHANGELOG if applicable

---

## Potential Issues and Mitigations

### Issue 1: Fractional Timesteps in Lag Conversion
**Problem:** `lag_days / tstep` may not be integer (e.g., 7 days / 2.0 tstep = 3.5 timesteps)

**Mitigation:**
- Use `round(Int, ...)` for conversion
- Document rounding behavior in docstrings
- Optionally add validation warning when rounding occurs
- Consider adding `validate_period_divisibility()` helper function

**Example validation function:**
```julia
function validate_period_divisibility(period::Dates.Day, tstep::Float64)
    period_days = Dates.value(period)
    exact_timesteps = period_days / tstep
    rounded_timesteps = round(Int, exact_timesteps)
    
    if !isapprox(exact_timesteps, rounded_timesteps, atol=1e-10)
        @warn "Period $period not evenly divisible by tstep=$tstep. " *
              "Rounding $exact_timesteps to $rounded_timesteps timesteps."
    end
    
    return rounded_timesteps
end
```

### Issue 2: Performance with Small Timestep
**Problem:** `tstep = 0.1` means 10x more iterations, significantly slower

**Mitigation:**
- Document performance implications clearly
- Suggest profiling before choosing very small timesteps
- Consider adding performance benchmarks for different timesteps
- Recommend starting with larger timesteps and reducing only if needed

### Issue 3: Probability Approximation Breaks Down
**Problem:** `convert_rate_to_prob(rate) = min(rate, 1.0)` assumes small `rate * tstep`

**Current implementation** (seir-model.jl:249-250):
```julia
convert_rate_to_prob(rate::Float64)::Float64 = min(rate, 1.0)
```

**Mitigation:**
- Add validation in `seir_mod!` to warn if `rate * tstep > 0.5`
- Document limitation in `SimTimeParameters` docstring
- Future improvement: use exact Poisson/binomial conversions

**Validation code:**
```julia
function seir_mod!(...)
    # ... existing setup ...
    
    # Validate Euler approximation
    max_rate_timestep = max(mu_timestep, sigma_timestep, gamma_timestep, epsilon_timestep)
    if max_rate_timestep > 0.5
        @warn "Large rate*tstep detected: max = $max_rate_timestep > 0.5. " *
              "Euler approximation may be inaccurate. Consider reducing tstep or adjusting parameters."
    end
    
    # ... rest of function ...
end
```

### Issue 4: Seasonality with Non-Divisor Timesteps
**Problem:** If `tstep = 3.5`, seasonality period won't align perfectly with years

**Mitigation:**
- Document that tstep should ideally divide 365 evenly
- List recommended timesteps: 1.0, 0.5, 0.25, 1/24, 7.0, etc.
- Explain that non-divisor timesteps may cause slight drift
- For most applications, the drift will be negligible

**Recommended timesteps:**
- 1.0 (daily)
- 0.5 (12-hour)
- 0.25 (6-hour)
- 1/24 ≈ 0.0417 (hourly)
- 1/48 ≈ 0.0208 (30-minute)
- 7.0 (weekly)
- 365/12 ≈ 30.42 (monthly - but note this is approximate)

### Issue 5: Backward Compatibility
**Problem:** Existing code uses integer lag/window parameters

**Mitigation:**
- Provide backward-compatible constructors accepting `Integer`
- Auto-convert `Integer` → `Dates.Day` in constructors
- Document migration path clearly
- Add deprecation warnings if desired (optional)

**Example:**
```julia
# Backward-compatible constructor
function IndividualTestSpecification(
        sensitivity::Float64,
        specificity::Float64,
        lag::Integer
    )
    @warn "Using Integer for lag is deprecated. Use Dates.Day($lag) instead." maxlog=1
    return IndividualTestSpecification(sensitivity, specificity, Dates.Day(lag))
end
```

---

## Testing Strategy

### Unit Tests
- Test time conversion utilities in isolation
- Test each updated function with different timesteps
- Test backward compatibility (integer constructors)
- Test edge cases (non-divisible periods, very large/small timesteps)

### Integration Tests
- Run full simulations with different timesteps
- Compare results across timesteps (accounting for discretization)
- Verify seasonality peaks occur at correct calendar times
- Verify test lags represent correct time periods

### Regression Tests
- Ensure existing tests still pass
- Verify default behavior unchanged (tstep=1.0)
- Check that results match previous version for tstep=1.0

### Performance Tests
- Benchmark simulation time for different timesteps
- Document performance scaling (should be approximately linear with 1/tstep)
- Identify any unexpected performance bottlenecks

---

## Success Criteria

### Functional Requirements
- ✅ Simulation runs correctly with any positive tstep value
- ✅ Seasonality peaks occur at correct calendar times regardless of tstep
- ✅ Test lags represent correct time periods (e.g., 2 days is always 2 days)
- ✅ Moving average windows cover correct time periods
- ✅ Detection thresholds work correctly with different timesteps
- ✅ All existing tests pass
- ✅ New tests cover flexible timestep functionality

### Code Quality Requirements
- ✅ All variable names clearly indicate units (timestep_idx vs. days)
- ✅ All functions have clear docstrings with examples
- ✅ No ambiguous parameters (all use Dates.Day for time periods)
- ✅ Code formatted with Runic.jl
- ✅ No type instabilities (JET.jl tests pass)
- ✅ All Aqua.jl quality checks pass

### Documentation Requirements
- ✅ Migration guide explains how to use flexible timesteps
- ✅ Examples demonstrate different timestep scenarios
- ✅ Limitations and performance implications documented
- ✅ Troubleshooting guide for common issues
- ✅ README updated with flexible timestep section

---

## Timeline Estimate

**Phase 1 (Foundation):** 2-3 hours
- Type updates and utility functions
- Straightforward changes

**Phase 2 (Seasonality):** 1 hour
- Simple function signature changes
- Few call sites to update

**Phase 3 (Diagnostic Testing):** 2-3 hours
- Multiple functions to update
- Need to add wrapper functions

**Phase 4 (Detection):** 2-3 hours
- Function renaming and updates
- Multiple files affected

**Phase 5 (Documentation):** 1-2 hours
- Search/replace and manual review
- Careful to not change unintended instances

**Phase 6 (Integration):** 1-2 hours
- Verification and validation
- Should be mostly working already

**Phase 7 (Testing):** 3-4 hours
- Writing comprehensive tests
- Debugging any issues found

**Phase 8 (Documentation):** 2-3 hours
- Writing migration guide
- Updating examples and docstrings

**Total Estimated Time:** 14-21 hours

---

## Notes

### Key Design Principles
1. **Clarity over brevity**: Use `Dates.Day` even though it's more verbose
2. **Type safety**: Let the compiler catch unit mismatches
3. **User-friendly**: Parameters in intuitive units (days), conversion is internal
4. **Backward compatible**: Provide migration path for existing code
5. **Well-documented**: Clear examples and migration guide

### Alternative Approaches Considered

**Alternative 1: Store everything in timesteps**
- Rejected: Less intuitive for users (lag=2 means different things for different tsteps)

**Alternative 2: Use separate time-in-days tracking**
- Rejected: More complex, harder to maintain, not necessary

**Alternative 3: Keep current design, just document better**
- Rejected: Doesn't solve the core ambiguity problem

### Future Enhancements
- Consider adaptive timestep methods for better accuracy/performance tradeoff
- Implement exact Poisson/binomial conversions instead of Euler approximation
- Add automatic timestep selection based on disease dynamics parameters
- Support for time-varying timesteps (smaller during outbreaks, larger during endemic periods)

---

## Questions to Resolve Before Implementation

1. **Deprecation warnings:** Should we add deprecation warnings for integer constructors, or just silently convert?
   - Recommendation: Silent conversion for now, add warnings in future version if needed

2. **Validation strictness:** Should non-divisible periods error, warn, or silently round?
   - Recommendation: Silent round, optional warning with `@warn ... maxlog=1`

3. **Function naming:** Should `DailyThreshold` be renamed to `InstantaneousThreshold`?
   - Recommendation: Review usage context first, rename only if "Daily" is truly misleading

4. **Testing scope:** How comprehensive should timestep comparison tests be?
   - Recommendation: Test key scenarios (0.5, 1.0, 7.0), document known discretization differences

5. **Performance benchmarks:** Should we add formal performance benchmarks?
   - Recommendation: Add simple benchmarks in tests, document scaling behavior

---

## References

- Julia Dates documentation: https://docs.julialang.org/en/v1/stdlib/Dates/
- Current codebase analysis (from exploration agent)
- Kent Beck's "Tidy First" principles (from AGENTS.md)
- DrWatson integration patterns (from project structure)

---

**Document Version:** 1.0  
**Last Updated:** 2026-01-22  
**Author:** AI Assistant (based on codebase analysis and user requirements)
