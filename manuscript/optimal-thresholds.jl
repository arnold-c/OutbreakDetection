#%%
using DrWatson
@quickactivate "OutbreakDetection"

using ProgressMeter
using FLoops
using NaNMath: NaNMath
using CSV: CSV
using DataFrames
using CategoricalArrays
using DataFramesMeta
using Statistics
using Match

using OutbreakDetectionUtils

includet(srcdir("makie-plotting-setup.jl"))
include(srcdir("ensemble-parameters.jl"))

#%%
rdt_spec_vec = [
    IndividualTestSpecification(0.85, 0.85, 0),
    IndividualTestSpecification(0.9, 0.9, 0),
    IndividualTestSpecification(1.0, 0.0, 0),
]

elisa_test_spec_vec = [
    IndividualTestSpecification(1.0, 1.0, 0),
    IndividualTestSpecification(1.0, 1.0, 3),
    IndividualTestSpecification(1.0, 1.0, 7),
    IndividualTestSpecification(1.0, 1.0, 14),
]

optimal_threshold_alertthreshold_vec = collect(1:1:15)

R_0_vec = collect(16.0)

ensemble_dynamics_spec_vec = create_combinations_vec(
    DynamicsParameters,
    (
        [ensemble_state_specification.init_states.N],
        [27],
        [0.2],
        [SIGMA],
        [GAMMA],
        R_0_vec,
        [0.8],
    ),
)

ensemble_spec_vec = create_combinations_vec(
    EnsembleSpecification,
    (
        [ensemble_model_type],
        [ensemble_state_specification],
        ensemble_dynamics_spec_vec,
        [ensemble_time_specification],
        [ensemble_nsims],
    ),
)

# alert_method_vec = ["movingavg", "dailythreshold_movingavg"]
alert_method_vec = ["movingavg"]

#%%
cfr_df = CSV.read(
    datadir("CFR_2022.csv"),
    DataFrame; delim = ',',
    header = true,
    types = [String, Int64, Float64],
    silencewarnings = true,
)
dropmissing!(cfr_df)

gha_cfr = only(cfr_df[cfr_df.country .== "GHA", :CFR])

population_df = CSV.read(
    datadir("input-populations.csv"),
    DataFrame; delim = ',',
    header = true,
)

gha_2022_pop = only(
    population_df[population_df.ISO3_code .== "GHA", "2022"]
)
gha_2022_scale_population =
    gha_2022_pop / ensemble_state_specification.init_states.N

countries = [
    (;
        name = "Ghana",
        code = "GHA",
        cfr = gha_cfr,
        year = "2022",
        population_size = gha_2022_pop,
        scale_population = gha_2022_scale_population,
    ),
]

#%%
noise_spec_vec = [
    PoissonNoiseSpecification(8.0),
    DynamicalNoiseSpecification(
        "dynamical", 5.0, 7, 14, "in-phase", 0.15, 0.050
    ),
]
ensemble_specification = ensemble_spec_vec[1]
alertmethod = alert_method_vec[1]

# #%%
# thresholds_df = DataFrame(
#     "percent_clinic_tested" => Float64[],
#     "sensitivity" => Float64[],
#     "specificity" => Float64[],
#     "test_lag" => Int64[],
#     "alert_threshold" => Int64[],
#     "accuracy" => Float64[],
#     "noise_spec" => String[],
# )

#%%
thresholds_df = mapreduce(
    vcat,
    enumerate(
        Iterators.product(
            noise_spec_vec, ensemble_spec_vec, alert_method_vec
        ),
    ),
) do (
    i, (ensemble_noise_specification, ensemble_specification, alertmethod)
)
    optimal_threshold_comparison_params = (
        alertthreshold_vec = optimal_threshold_alertthreshold_vec,
        ensemble_specification = ensemble_specification,
        noise_specification = ensemble_noise_specification,
        outbreak_specification = ensemble_outbreak_specification,
        moving_avg_detection_lag = ensemble_moving_avg_detection_lag,
        percent_visit_clinic = ensemble_percent_visit_clinic,
        alertmethod = alertmethod,
    )

    rdt_thresholds_vec = calculate_OptimalThresholdCharacteristics(
        ensemble_percent_clinic_tested_vec,
        rdt_spec_vec,
        optimal_threshold_comparison_params,
    )

    optimal_thresholds_df = create_optimal_thresholds_df(
        rdt_thresholds_vec
    )
    optimal_thresholds_df[!, :noise_spec] .= get_noise_description(
        ensemble_noise_specification
    )

    if i == 1
        elisa_thresholds_vec = calculate_OptimalThresholdCharacteristics(
            ensemble_percent_clinic_tested_vec,
            elisa_test_spec_vec,
            optimal_threshold_comparison_params,
        )

        elisa_thresholds_df = create_optimal_thresholds_df(
            elisa_thresholds_vec
        )
        elisa_thresholds_df[!, :noise_spec] .= "Common"
        optimal_thresholds_df = vcat(
            optimal_thresholds_df, elisa_thresholds_df
        )
    end

    return optimal_thresholds_df
end

replace!(
    thresholds_df[!, :noise_spec],
    "dynamical, in-phase" => "Dynamical noise: in-phase",
    "Poisson" => "Poisson noise",
)

# #%%
# for (i, (ensemble_noise_specification, ensemble_specification, alertmethod)) in
#     enumerate(
#     Iterators.product(
#         noise_spec_vec, ensemble_spec_vec, alert_method_vec
#     ),
# )
#     @info "Creating plots and tables for R0: $(ensemble_specification.dynamics_parameters.R_0), $(getdirpath(ensemble_noise_specification)), $(alertmethod)"
#     println("==============================================")
#
#     optimal_threshold_comparison_params = (
#         alertthreshold_vec = optimal_threshold_alertthreshold_vec,
#         ensemble_specification = ensemble_specification,
#         noise_specification = ensemble_noise_specification,
#         outbreak_specification = ensemble_outbreak_specification,
#         moving_avg_detection_lag = ensemble_moving_avg_detection_lag,
#         percent_visit_clinic = ensemble_percent_visit_clinic,
#         alertmethod = alertmethod,
#     )
#
#     rdt_thresholds_vec = calculate_OptimalThresholdCharacteristics(
#         ensemble_percent_clinic_tested_vec,
#         rdt_spec_vec,
#         optimal_threshold_comparison_params,
#     )
#
#     optimal_thresholds_df = create_optimal_thresholds_df(
#         rdt_thresholds_vec
#     )
#     optimal_thresholds_df[!, :noise_spec] .= get_noise_description(
#         ensemble_noise_specification
#     )
#
#     if i == 1
#         elisa_thresholds_vec = calculate_OptimalThresholdCharacteristics(
#             ensemble_percent_clinic_tested_vec,
#             elisa_test_spec_vec,
#             optimal_threshold_comparison_params,
#         )
#
#         elisa_thresholds_df = create_optimal_thresholds_df(
#             elisa_thresholds_vec
#         )
#         elisa_thresholds_df[!, :noise_spec] .= "Common"
#         optimal_thresholds_df = vcat(
#             optimal_thresholds_df, elisa_thresholds_df
#         )
#     end
# end
# thresholds_df = vcat(thresholds_df, optimal_thresholds_df)

#%%
function create_wide_df(
    long_df,
    outcome::Symbol;
    noise_order = ["Dynamical noise: in-phase", "Poisson noise", "Common"],
    digits = 3,
)
    if digits == 0
        long_df[!, outcome] .= Int64.(long_df[!, outcome])
    else
        long_df[!, outcome] = round.(long_df[!, outcome]; digits = digits)
    end

    wide_df =
        map(
            collect(groupby(long_df, :noise_spec))) do df
            unstack(
                df,
                [:noise_spec, :sensitivity, :specificity, :test_lag],
                :percent_clinic_tested,
                outcome,
            )
        end |>
        x -> vcat(x...; cols = :union)

    wide_df[!, :noise_spec] = categorical(
        wide_df.noise_spec
    )

    levels!(
        wide_df.noise_spec,
        noise_order,
    )

    sort!(wide_df, [:noise_spec, order(:specificity; rev = false)])

    wide_df[!, :test_type] =
        get_test_type.(
            wide_df.sensitivity,
            wide_df.specificity,
            wide_df.test_lag,
        )

    select!(
        wide_df,
        :noise_spec,
        :test_type,
        :test_lag,
        Not(
            :sensitivity,
            :specificity,
            :test_lag,
        ),
    )
    rename!(
        wide_df,
        Dict(
            :noise_spec => "Noise Type",
            :test_type => "Test Type",
            :test_lag => "Test Lag",
        ),
    )
    rename!(
        x -> "$(Int64(parse(Float64, x)*100))%",
        wide_df;
        cols = contains('0'),
    )
    return wide_df
end

function get_test_type(sensitivity, specificity, test_lag)
    return Match.@match (sensitivity, specificity, test_lag) begin
        (1.0, 0.0, 0) => "Clinical Case Definition"
        (x::AbstractFloat, x::AbstractFloat, 0) where {x<1.0} => "RDT Equivalent ($(sensitivity * 100)%)"
        (1.0, 1.0, x::Int) => "ELISA Equivalent"
    end
end

#%%
wide_thresholds_df = create_wide_df(thresholds_df, :alert_threshold; digits = 0)
CSV.write(projectdir("manuscript/optimal-thresholds.csv"), wide_thresholds_df)

#%%
wide_accuracy_df = create_wide_df(thresholds_df, :accuracy; digits = 2)
CSV.write(
    projectdir("manuscript/optimal-thresholds_accuracy.csv"), wide_accuracy_df
)
