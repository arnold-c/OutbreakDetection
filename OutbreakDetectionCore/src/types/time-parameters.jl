export SimTimeParameters

"""
    SimTimeParameters

Parameters defining the temporal structure of a simulation.

# Fields
- `tmin::AbstractFloat`: Start time of the simulation
- `tmax::AbstractFloat`: End time of the simulation in days
- `tstep::AbstractFloat`: Time step size in days
- `trange::StepRangeLen`: Range from tmin to tmax by tstep
- `tspan::Tuple{AbstractFloat, AbstractFloat}`: Tuple of (tmin, tmax)
- `tlength::Int`: Number of time steps in trange

# Constructor
    SimTimeParameters(; tmin=0.0, tmax=365.0*100.0, tstep=1.0)

Creates a `SimTimeParameters` instance with default values suitable for a 100-year
simulation with daily time steps.
"""
struct SimTimeParameters{
        T1 <: AbstractFloat,
        T2 <: StepRangeLen,
        T3 <: Tuple{T1, T1},
        T4 <: Int,
    }
    tmin::T1
    tmax::T1
    tstep::T1
    trange::T2
    tspan::T3
    tlength::T4
end

function SimTimeParameters(; tmin = 0.0, tmax = 365.0 * 100.0, tstep = 1.0)
    return SimTimeParameters(
        tmin, tmax, tstep, tmin:tstep:tmax, (tmin, tmax),
        length(tmin:tstep:tmax),
    )
end
