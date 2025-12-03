export SimTimeParameters

"""
    SimTimeParameters

Parameters defining the temporal structure of a simulation.

This struct encapsulates all time-related parameters for SEIR simulations,
including the time range, step size, and derived quantities.

# Fields
- `tmin::Float64`: Start time of the simulation (default: 0.0)
- `tmax::Float64`: End time of the simulation in days (default: 36,500.0 = 100 years)
- `tstep::Float64`: Time step size in days (default: 1.0)
- `trange::StepRangeLen{Float64, Float64, Float64, Int64}`: Range from tmin to tmax by tstep
- `tspan::Tuple{Float64, Float64}`: Tuple of (tmin, tmax) for ODE solvers
- `tlength::Int64`: Number of time steps in trange

# Constructors
    SimTimeParameters(; tmin = 0.0, tmax = 365.0 * 100.0, tstep = 1.0)

Create a SimTimeParameters instance with specified time bounds and step size.
Automatically calculates trange, tspan, and tlength.

# Validation
- tmax must be greater than tmin + tstep
- tstep must be positive

# Examples
```julia
# Default: 100-year simulation with daily steps
time_params = SimTimeParameters()

# Custom: 20-year simulation with daily steps
time_params = SimTimeParameters(tmax = 365.0 * 20)

# Custom: 10-year simulation with weekly steps
time_params = SimTimeParameters(tmax = 365.0 * 10, tstep = 7.0)
```

# See Also
- [`DynamicsParameters`](@ref): Disease dynamics parameters
- [`StateParameters`](@ref): Initial state parameters
"""
struct SimTimeParameters
    tmin::Float64
    tmax::Float64
    tstep::Float64
    trange::StepRangeLen{Float64, Float64, Float64, Int64}
    tspan::Tuple{Float64, Float64}
    tlength::Int64

    function SimTimeParameters(tmin, tmax, tstep, trange, tspan, tlength)
        @assert tmax > tmin + tstep "tmax must be greater than tmin + tstep"
        @assert tstep > 0.0 "tstep must be positive"
        return new(tmin, tmax, tstep, trange, tspan, tlength)
    end
end

"""
    SimTimeParameters(; tmin = 0.0, tmax = 365.0 * 100.0, tstep = 1.0)

Keyword constructor for SimTimeParameters.

Automatically calculates derived fields (trange, tspan, tlength) from the
provided time bounds and step size.

# Arguments
- `tmin::Float64`: Start time (default: 0.0)
- `tmax::Float64`: End time in days (default: 36,500.0 = 100 years)
- `tstep::Float64`: Time step in days (default: 1.0)

# Returns
- `SimTimeParameters`: Fully initialized time parameters

# Examples
```julia
# Default 100-year simulation
params = SimTimeParameters()

# 20-year simulation
params = SimTimeParameters(tmax = 365.0 * 20)

# 5-year simulation starting at year 10
params = SimTimeParameters(tmin = 365.0 * 10, tmax = 365.0 * 15)
```
"""
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
