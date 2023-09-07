# module ODStructs
#
# export SimTimeParameters, EnsembleSpecification, DynamicsParameters,
#     StateParameters, OutbreakThresholdChars, OutbreakDetectionSpecification,
#     OutbreakSpecification, IndividualTestSpecification, NoiseSpecification

using LabelledArrays

# include("transmission-functions.jl")
# using .TransmissionFunctions

struct SimTimeParameters
    tmin
    tmax
    tstep
    trange
    tspan
    tlength
end

function SimTimeParameters(; tmin = 0.0, tmax = 365.0 * 100.0, tstep = 1.0)
    return SimTimeParameters(
        tmin, tmax, tstep, tmin:tstep:tmax, (tmin, tmax),
        length(tmin:tstep:tmax),
    )
end

const POPULATION_N = 500_000
const LATENT_PER_DAYS = 8
const DUR_INF_DAYS = 5
const R0 = 10.0
const SIGMA = 1 / LATENT_PER_DAYS
const GAMMA = 1 / DUR_INF_DAYS
const LIFE_EXPECTANCY_YEARS = 62.5
const MU = 1 / (LIFE_EXPECTANCY_YEARS * 365)
const BETA_MEAN = calculate_beta(R0, GAMMA, MU, 1, POPULATION_N)
const BETA_FORCE = 0.2
const EPSILON = calculate_import_rate(MU, R0, POPULATION_N)

struct DynamicsParameters
    beta_mean::Float64
    beta_force::Float64
    sigma::Float64
    gamma::Float64
    mu::Float64
    annual_births_per_k::Union{Int64,Float64}
    epsilon::Float64
    R_0::Float64
end

function DynamicsParameters(sigma::Float64, gamma::Float64, R_0::Float64)
    annual_births_per_k = 1000 / LIFE_EXPECTANCY_YEARS

    return DynamicsParameters(
        BETA_MEAN,
        BETA_FORCE,
        sigma,
        gamma,
        MU,
        annual_births_per_k,
        EPSILON,
        R_0,
    )
end

function DynamicsParameters(
    N::Int64, annual_births_per_k::Int64, beta_force::Float64
)
    mu = calculate_mu(annual_births_per_k)
    beta_mean = calculate_beta(R0, GAMMA, mu, 1, N)
    epsilon = calculate_import_rate(mu, R0, N)

    return DynamicsParameters(
        beta_mean,
        beta_force,
        SIGMA,
        GAMMA,
        mu,
        annual_births_per_k,
        epsilon,
        R0,
    )
end

struct StateParameters
    init_states
    init_state_props
end

function StateParameters(N::Int64, init_state_props::Dict)
    return StateParameters(;
        N = N,
        s_prop = init_state_props[:s_prop],
        e_prop = init_state_props[:e_prop],
        i_prop = init_state_props[:i_prop]
    )
end

function StateParameters(;
    N = 500_00, s_prop = 0.1, e_prop = 0.01, i_prop = 0.01
)
    r_prop = 1 - (s_prop + e_prop + i_prop)
    states = @LArray [
        map(x -> Int64(round(x * N)), [s_prop, e_prop, i_prop, r_prop])..., N
    ] (:S, :E, :I, :R, :N)
    state_props = @LArray [s_prop, e_prop, i_prop, r_prop] (
        :s_prop, :i_prop, :e_prop, :r_prop
    )

    return StateParameters(
        states, state_props
    )
end

struct EnsembleSpecification
    modeltypes::Tuple
    state_parameters::StateParameters
    dynamics_parameters::DynamicsParameters
    time_parameters::SimTimeParameters
    nsims::Int64
    dirpath::String
end

function EnsembleSpecification(
    modeltypes::Tuple,
    state_parameters::StateParameters,
    dynamics_parameters::DynamicsParameters,
    time_parameters::SimTimeParameters,
    nsims::Int64,
)
    dirpath = datadir(
        modeltypes...,
        "N_$(state_parameters.init_states.N)",
        "r_$(state_parameters.init_state_props.r_prop)",
        "nsims_$(nsims)",
        "births_per_k_$(dynamics_parameters.annual_births_per_k)",
        "beta_force_$(dynamics_parameters.beta_force)",
        "tmax_$(time_parameters.tmax)",
        "tstep_$(time_parameters.tstep)",
    )

    return EnsembleSpecification(
        modeltypes,
        state_parameters,
        dynamics_parameters,
        time_parameters,
        nsims,
        dirpath,
    )
end

struct OutbreakThresholdChars{A,B,C,D}
    crosstab::A
    tp::B
    tn::B
    fp::B
    fn::B
    sensitivity::C
    specificity::C
    ppv::C
    npv::C
    noutbreaks::B
    ndetectoutbreaks::B
    outbreakbounds::D
    detectoutbreakbounds::D
end

struct OutbreakSpecification
    outbreak_threshold
    minimum_outbreak_duration
    minimum_outbreak_size
end

struct OutbreakDetectionSpecification
    detection_threshold
    moving_average_lag
    percent_tested
    test_result_lag
end

function OutbreakDetectionSpecification(
    detection_threshold,
    moving_average_lag,
    percent_clinic,
    percent_clinic_tested,
    test_result_lag,
)
    return OutbreakDetectionSpecification(
        detection_threshold,
        moving_average_lag,
        percent_clinic * percent_clinic_tested,
        test_result_lag,
    )
end

struct IndividualTestSpecification
    sensitivity
    specificity
end

struct NoiseSpecification
    noise_type
    noise_array
end

struct ScenarioSpecification
    ensemble_specification::EnsembleSpecification
    outbreak_specification::OutbreakSpecification
    noise_specification::NoiseSpecification
    outbreak_detection_specification::OutbreakDetectionSpecification
    individual_test_specification::IndividualTestSpecification
    dirpath::String
end

function ScenarioSpecification(
    ensemble_specification::EnsembleSpecification,
    outbreak_specification::OutbreakSpecification,
    noise_specification::NoiseSpecification,
    outbreak_detection_specification::OutbreakDetectionSpecification,
    individual_test_specification::IndividualTestSpecification,
)
    dirpath = joinpath(
        ensemble_specification.dirpath,
        "noise_$(noise_specification.noise_type)",
        "min_outbreak_dur_$(outbreak_specification.minimum_outbreak_duration)",
        "min_outbreak_size_$(outbreak_specification.minimum_outbreak_size)",
        "outbreak_threshold_$(outbreak_specification.outbreak_threshold)",
        "detectthreshold_$(outbreak_detection_specification.detection_threshold)",
        "testlag_$(outbreak_detection_specification.test_result_lag)",
        "moveavglag_$(outbreak_detection_specification.moving_average_lag)",
        "perc_tested_$(outbreak_detection_specification.percent_tested)",
        "testsens_$(individual_test_specification.sensitivity)",
        "testspec_$(individual_test_specification.specificity)",
    )

    return ScenarioSpecification(
        ensemble_specification,
        outbreak_specification,
        noise_specification,
        outbreak_detection_specification,
        individual_test_specification,
        dirpath,
    )
end

# end
