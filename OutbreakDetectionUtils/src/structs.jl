using LabelledArrays: SLVector, SLArray
using StructArrays: StructVector
using Match: Match

struct SimTimeParameters{
    T1<:AbstractFloat,
    T2<:StepRangeLen,
    T3<:Tuple{T1,T1},
    T4<:Int,
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

struct DynamicsParameters{
    T1<:AbstractFloat,T2<:Union{<:Integer,T1},T3<:Function
}
    beta_mean::T1
    beta_force::T1
    seasonality::T3
    sigma::T1
    gamma::T1
    mu::T1
    annual_births_per_k::T2
    epsilon::T1
    R_0::T1
    vaccination_coverage::T1
end

function DynamicsParameters(
    sigma::Float64,
    gamma::Float64,
    R_0::Float64;
    vaccination_coverage::Float64 = 0.8,
)
    return DynamicsParameters(
        BETA_MEAN,
        BETA_FORCE,
        cos,
        sigma,
        gamma,
        MU,
        ANNUAL_BIRTHS_PER_K,
        EPSILON,
        R_0,
        vaccination_coverage,
    )
end

function DynamicsParameters(
    N::Int64,
    annual_births_per_k::Int64,
    beta_force::Float64,
    sigma::Float64,
    gamma::Float64,
    R_0::Float64,
    vaccination_coverage::Float64;
    seasonality::Function = cos,
)
    mu = calculate_mu(annual_births_per_k)
    beta_mean = calculate_beta(R_0, gamma, mu, 1, N)
    epsilon = calculate_import_rate(mu, R_0, N)

    return DynamicsParameters(
        beta_mean,
        beta_force,
        seasonality,
        sigma,
        gamma,
        mu,
        annual_births_per_k,
        epsilon,
        R_0,
        vaccination_coverage,
    )
end

function DynamicsParameters(
    N::Int64, annual_births_per_k::Int64, beta_force::Float64;
    vaccination_coverage::Float64 = 0.8,
)
    mu = calculate_mu(annual_births_per_k)
    beta_mean = calculate_beta(R0, GAMMA, mu, 1, N)
    epsilon = calculate_import_rate(mu, R0, N)

    return DynamicsParameters(
        beta_mean,
        beta_force,
        cos,
        SIGMA,
        GAMMA,
        mu,
        annual_births_per_k,
        epsilon,
        R0,
        vaccination_coverage,
    )
end

struct StateParameters{T1<:SLArray,T2<:SLArray}
    init_states::T1
    init_state_props::T2
end

function StateParameters(N::Int64, init_state_props::Dict)
    return StateParameters(;
        N = N,
        s_prop = init_state_props[:s_prop],
        e_prop = init_state_props[:e_prop],
        i_prop = init_state_props[:i_prop],
    )
end

function StateParameters(;
    N = 500_00, s_prop = 0.1, e_prop = 0.01, i_prop = 0.01
)
    r_prop = 1 - (s_prop + e_prop + i_prop)

    states = SLVector(;
        S = Int64(round(s_prop * N)),
        E = Int64(round(e_prop * N)),
        I = Int64(round(i_prop * N)),
        R = Int64(round(r_prop * N)),
        N = N,
    )
    state_props = SLVector(;
        s_prop = s_prop,
        e_prop = e_prop,
        i_prop = i_prop,
        r_prop = r_prop,
    )

    return StateParameters(
        states, state_props
    )
end

struct EnsembleSpecification{
    T1<:Tuple,
    T2<:StateParameters,
    T3<:DynamicsParameters,
    T4<:SimTimeParameters,
    T5<:Integer,
    T6<:AbstractString,
}
    modeltypes::T1
    state_parameters::T2
    dynamics_parameters::T3
    time_parameters::T4
    nsims::T5
    dirpath::T6
end

function EnsembleSpecification(
    modeltypes::Tuple,
    state_parameters::StateParameters,
    dynamics_parameters::DynamicsParameters,
    time_parameters::SimTimeParameters,
    nsims::Int64,
)
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

struct OutbreakThresholdChars{
    T1<:AbstractFloat,
    T2<:Integer,
    T3<:Vector{<:AbstractFloat},
    T4<:Vector{<:Integer},
    T5<:AbstractMatrix{<:Integer},
}
    daily_sensitivity::T1
    daily_specificity::T1
    daily_ppv::T1
    daily_npv::T1
    accuracy::T1
    f1_score::T1
    matchedoutbreakbounds::T5
    noutbreaks::T2
    nalerts::T2
    outbreak_duration_vec::T4
    alert_duration_vec::T4
    detected_outbreak_size::T4
    missed_outbreak_size::T4
    n_true_outbreaks_detected::T2
    n_missed_outbreaks::T2
    n_correct_alerts::T2
    n_false_alerts::T2
    n_alerts_per_outbreak::T4
    period_sum_per_outbreak::T4
    perc_true_outbreaks_detected::T1
    perc_true_outbreaks_missed::T1
    falsealert_trueoutbreak_prop::T1
    correctalert_trueoutbreak_prop::T1
    trueoutbreak_alerts_prop::T1
    outbreaksmissed_alerts_prop::T1
    perc_alerts_false::T1
    perc_alerts_correct::T1
    detectiondelays::T4
    cases_before_alerts::T4
    cases_perc_before_alerts::T3
    cases_after_alerts::T4
    cases_perc_after_alerts::T3
    unavoidable_cases::T2
    avoidable_cases::T2
    n_outbreak_cases::T2
    n_tests::T2
    n_outbreak_tests::T2
    mean_noise_incidence_ratio::T1
    mean_poisson_noise::T1
    poisson_noise_prop::T1
    proportion_timeseries_in_outbreak::T1
    proportion_timeseries_in_alert::T1
    alert_outbreak_timeseries_prop_diff::T1
end

struct OutbreakSpecification{T1<:Integer,T2<:AbstractString}
    outbreak_threshold::T1
    minimum_outbreak_duration::T1
    minimum_outbreak_size::T1
    dirpath::T2
end

function OutbreakSpecification(
    outbreak_threshold, minimum_outbreak_duration, minimum_outbreak_size
)
    dirpath = joinpath(
        "min_outbreak_dur_$(minimum_outbreak_duration)",
        "min_outbreak_size_$(minimum_outbreak_size)",
        "outbreak_threshold_$(outbreak_threshold)",
    )

    return OutbreakSpecification(
        outbreak_threshold,
        minimum_outbreak_duration,
        minimum_outbreak_size,
        dirpath,
    )
end

struct AlertMethod{T1<:AbstractString}
    method_name::T1
    function AlertMethod(method_name::T1) where {T1<:AbstractString}
        available_test_methods = [
            "dailythreshold", "movingavg", "dailythreshold_movingavg"
        ]
        if !in(method_name, available_test_methods)
            error(
                "$(method_name) is not a valid test method. It must be one of $(available_test_methods)"
            )
        end
        return new{T1}(method_name)
    end
end

struct OutbreakDetectionSpecification{
    TReal<:Real,T1<:Integer,T2<:AbstractFloat,T3<:AlertMethod,T4<:AbstractString
}
    alert_threshold::TReal
    moving_average_lag::T1
    percent_visit_clinic::T2
    percent_clinic_tested::T2
    percent_tested::T2
    alert_method::T3
    dirpath::T4
end

function OutbreakDetectionSpecification(
    alert_threshold,
    moving_average_lag,
    percent_visit_clinic,
    percent_clinic_tested,
    alert_method,
)
    alertdirpath = joinpath(
        "alertmethod_$(alert_method)", "alertthreshold_$(alert_threshold)"
    )
    testingdirpath = joinpath(
        "perc_visit_clinic_$(percent_visit_clinic)",
        "perc_clinic_tested_$(percent_clinic_tested)",
    )

    dirpath = Match.@match alert_method begin
        "dailythreshold" => joinpath(
            alertdirpath,
            testingdirpath,
        )
        _ => joinpath(
            alertdirpath,
            "moveavglag_$(moving_average_lag)",
            testingdirpath,
        )
    end

    return OutbreakDetectionSpecification(
        alert_threshold,
        moving_average_lag,
        percent_visit_clinic,
        percent_clinic_tested,
        percent_visit_clinic * percent_clinic_tested,
        AlertMethod(alert_method),
        dirpath,
    )
end

struct IndividualTestSpecification{T1<:AbstractFloat,T2<:Integer}
    sensitivity::T1
    specificity::T1
    test_result_lag::T2
end

function get_test_description(test_specification::IndividualTestSpecification)
    description = Match.@match test_specification begin
        IndividualTestSpecification(1.0, 0.0, 0) => "Clinical case definition"
        IndividualTestSpecification(0.98, 0.98, x::Int) || IndividualTestSpecification(0.95, 0.98, x::Int) => "ELISA-like ($(test_specification.sensitivity * 100)% sens, $(test_specification.specificity * 100)% spec, $(x) day lag)"
        IndividualTestSpecification(1.0, 1.0, x::Int) => "Perfect test ($(x) day lag)"
        IndividualTestSpecification(x::AbstractFloat, x::AbstractFloat, 0) where {x<1.0} => "RDT-like ($(test_specification.sensitivity * 100)% sens/spec)"
    end
    return description
end

function table_test_type(sensitivity, specificity, test_lag)
    return Match.@match (sensitivity, specificity, test_lag) begin
        (1.0, 0.0, 0) => "Clinical Case Definition"
        (x::AbstractFloat, x::AbstractFloat, 0) where {x<1.0} => "Imperfect Test ($(Int64(round(sensitivity * 100; digits = 0)))%)"
        (1.0, 1.0, x::Int) => "Perfect Test"
    end
end

function plot_test_description(test_specification::IndividualTestSpecification)
    return Match.@match test_specification begin
        IndividualTestSpecification(1.0, 0.0, 0) => "Clinical Case Definition"
        IndividualTestSpecification(x::AbstractFloat, x::AbstractFloat, 0) where {x<1.0} => "Imperfect Test ($(Int64(round(x * 100; digits = 0)))%)"
        IndividualTestSpecification(1.0, 1.0, x::Int) => "Perfect Test ($x-day lag)"
    end
end

abstract type NoiseSpecification end

struct PoissonNoiseSpecification{
    T1<:AbstractString,T2<:AbstractFloat
} <: NoiseSpecification
    noise_type::T1
    noise_mean_scaling::T2
end

function PoissonNoiseSpecification(noise_mean_scaling)
    return PoissonNoiseSpecification("Poisson", noise_mean_scaling)
end

struct DynamicalNoiseSpecification{
    T1<:AbstractString,T2<:AbstractFloat,T3<:Integer
} <: NoiseSpecification
    noise_type::T1
    R_0::T2
    latent_period::T3
    duration_infection::T3
    correlation::T1
    noise_mean_scaling::T2
    vaccination_coverage::T2
end

function DynamicalNoiseSpecification()
end

function get_noise_description(
    noise_specification::T
) where {T<:NoiseSpecification}
    return noise_specification.noise_type
end

function get_noise_description(noise_specification::DynamicalNoiseSpecification)
    return string(
        noise_specification.noise_type, ", ", noise_specification.correlation
    )
end

function get_noise_magnitude(
    noise_specification::T
) where {T<:NoiseSpecification}
    return string("Poisson scaling: ", noise_specification.noise_mean_scaling)
end

function get_noise_magnitude(
    noise_specification::DynamicalNoiseSpecification
)
    return string("Rubella vax: ", noise_specification.vaccination_coverage)
end

function getdirpath(spec::NoiseSpecification)
    return reduce(
        joinpath,
        map(
            p -> "$(p)_$(getproperty(spec, p))",
            propertynames(spec),
        ),
    )
end

struct ScenarioSpecification{
    T1<:EnsembleSpecification,
    T2<:OutbreakSpecification,
    T3<:NoiseSpecification,
    T4<:OutbreakDetectionSpecification,
    T5<:IndividualTestSpecification,
    T6<:AbstractString,
}
    ensemble_specification::T1
    outbreak_specification::T2
    noise_specification::T3
    outbreak_detection_specification::T4
    individual_test_specification::T5
    dirpath::T6
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
        outbreak_specification.dirpath,
        getdirpath(noise_specification),
        outbreak_detection_specification.dirpath,
        "testsens_$(individual_test_specification.sensitivity)",
        "testspec_$(individual_test_specification.specificity)",
        "testlag_$(individual_test_specification.test_result_lag)",
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

struct TestPositivity{T1<:AbstractArray{<:AbstractFloat}}
    one_day::T1
    seven_day::T1
    fourteen_day::T1
    thirty_day::T1
end

function TestPositivity(true_positive_vec, total_test_vec, alert_vec)
    return TestPositivity(
        calculate_test_positivity(
            true_positive_vec, total_test_vec, alert_vec, 1
        ),
        calculate_test_positivity(
            true_positive_vec, total_test_vec, alert_vec, 7
        ),
        calculate_test_positivity(
            true_positive_vec, total_test_vec, alert_vec, 14
        ),
        calculate_test_positivity(
            true_positive_vec, total_test_vec, alert_vec, 30
        ),
    )
end

struct OptimalThresholdCharacteristics{
    T1<:StructVector{<:OutbreakThresholdChars},
    T2<:IndividualTestSpecification,
    T3<:NoiseSpecification,
    T4<:AbstractFloat,
    TReal<:Real,
}
    outbreak_threshold_chars::T1
    individual_test_specification::T2
    noise_specification::T3
    percent_clinic_tested::T4
    alert_threshold::TReal
    accuracy::T4
end

abstract type OptimizationMethods end

# QuadDIRECT optimization
struct QD <: OptimizationMethods end
# MultistartOptimization
struct MSO <: OptimizationMethods end
