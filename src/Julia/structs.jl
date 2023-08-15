struct EnsembleTransmissionParameters
    R₀::Float64
    σ::Float64
    γ::Float64
end

struct EnsembleTimeParameters
    tlower
    tmax
    tstep
    trange
    tlength
    # Use an inner constructor to calculate the trange and tlength, and ensure
    # that the struct pointer can be redefined below without rewriting a const method
    function EnsembleTimeParameters(tmin, tmax, tstep)
        return new(tmin, tmax, tstep, tmin:tstep:tmax, length(tmin:tstep:tmax))
    end
end

struct EnsembleSpecification
    modeltypes::Tuple
    N::Int64
    Rinit_prop::Float64
    nsims::Int64
end
