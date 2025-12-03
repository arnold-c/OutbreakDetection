export SEIRRun,
    NoiseRun,
    DynamicalNoiseRun,
    PoissonNoiseRun

"""
    SEIRRun

Results from a single SEIR epidemic simulation run.

This type is designed to work efficiently with StructVectors for ensemble storage.
Using StructVector{SEIRRun} provides column-wise access to simulation results,
which is more memory-efficient and faster for analysis than row-wise storage.

# Fields
- `states::Vector{StaticArrays.SVector{5, Int64}}`: Time series of [S, E, I, R, N]
- `incidence::Vector{Int64}`: Time series of new infections
- `Reff::Vector{Float64}`: Time series of effective reproduction number

# Examples
```julia
# Create a single run result
run = SEIRRun(
    states = [SVector{5, Int64}(450000, 100, 50, 49850, 500000), ...],
    incidence = [10, 12, 15, ...],
    Reff = [1.2, 1.3, 1.1, ...]
)

# Store ensemble results in StructVector for efficiency
ensemble_results = StructVector{SEIRRun}([run1, run2, run3, ...])

# Access all incidence time series efficiently
all_incidence = ensemble_results.incidence  # Vector{Vector{Int64}}
```

# See Also
- [`DynamicalNoiseRun`](@ref): Results from dynamical noise simulations
- [`PoissonNoiseRun`](@ref): Results from Poisson noise simulations
"""
Base.@kwdef struct SEIRRun
    states::Vector{StaticArrays.SVector{5, Int64}}
    incidence::Vector{Int64}
    Reff::Vector{Float64}
end

"""
    DynamicalNoiseRun

Results from noise simulations with both Poisson and dynamical noise components.

This type stores the results of dynamical noise simulations, which combine
a dynamical SEIR process with optional Poisson noise. It includes both the
noise time series and summary statistics.

# Fields
- `incidence::Vector{Vector{Int64}}`: Collection of noisy incidence time series
- `Reff::Vector{Vector{Float64}}`: Collection of noisy Reff time series
- `mean_noise::Float64`: Overall mean noise level (dynamical + Poisson)
- `mean_poisson_noise::Float64`: Mean Poisson noise component
- `mean_dynamic_noise::Float64`: Mean dynamical noise component

# Examples
```julia
# Create dynamical noise run results
noise_run = DynamicalNoiseRun(
    incidence = [[10, 12, 15, ...], [11, 13, 14, ...], ...],
    Reff = [[1.2, 1.3, 1.1, ...], [1.1, 1.2, 1.0, ...], ...],
    mean_noise = 15.5,
    mean_poisson_noise = 5.5,
    mean_dynamic_noise = 10.0
)

# Store in StructVector for ensemble analysis
ensemble_noise = StructVector{DynamicalNoiseRun}([noise_run1, noise_run2, ...])
```

# See Also
- [`SEIRRun`](@ref): Single SEIR simulation results
- [`PoissonNoiseRun`](@ref): Simpler Poisson-only noise results
- [`DynamicalNoise`](@ref): Noise specification type
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

This type stores the results of simple Poisson noise simulations, which
add random Poisson-distributed noise to a baseline incidence pattern.

# Fields
- `incidence::Vector{Vector{Int64}}`: Collection of noisy incidence time series
- `mean_noise::Float64`: Mean Poisson noise level

# Examples
```julia
# Create Poisson noise run results
noise_run = PoissonNoiseRun(
    incidence = [[10, 12, 15, ...], [11, 13, 14, ...], ...],
    mean_noise = 12.5
)

# Store in StructVector for ensemble analysis
ensemble_noise = StructVector{PoissonNoiseRun}([noise_run1, noise_run2, ...])
```

# See Also
- [`DynamicalNoiseRun`](@ref): More complex dynamical noise results
- [`PoissonNoise`](@ref): Noise specification type
"""
Base.@kwdef struct PoissonNoiseRun
    incidence::Vector{Vector{Int64}}
    mean_noise::Float64
end

# Sum type for noise runs
abstract type AbstractNoiseRun end

LightSumTypes.@sumtype NoiseRun(DynamicalNoiseRun, PoissonNoiseRun) <: AbstractNoiseRun
