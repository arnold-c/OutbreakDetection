struct EnsembleTransmissionParameters
    R₀::Float64
    sigma::Float64
    γ::Float64
end

struct SimTimeParameters
    tlower
    tmax
    tstep
    trange
    tlength
    # Use an inner constructor to calculate the trange and tlength, and ensure
    # that the struct pointer can be redefined below without rewriting a const method
    function SimTimeParameters(tmin, tmax, tstep)
        return new(tmin, tmax, tstep, tmin:tstep:tmax, length(tmin:tstep:tmax))
    end
end

struct EnsembleSpecification
    modeltypes::Tuple
    N::Int64
    Rinit_prop::Float64
    nsims::Int64
    births_per_k::Int64
    beta_force::Float64
    time_parameters::SimTimeParameters
end

