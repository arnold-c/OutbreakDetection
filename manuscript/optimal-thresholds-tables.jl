#%%
using DrWatson
@quickactivate "OutbreakDetection"

using DataFrames
using CategoricalArrays
using Match: Match
using CSV: CSV

using OutbreakDetectionUtils: create_optimal_thresholds_df

include(projectdir("manuscript", "optimal-thresholds-loading.jl"));

if false
    include("optimal-thresholds-loading.jl")
end

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
dynamical_noise_rdt_optimal_solutions = filter(
    chars -> in(chars.individual_test_specification, rdt_test_spec_vec),
    dynamical_noise_optimal_solutions,
);

poisson_noise_rdt_optimal_solutions = filter(
    chars -> in(chars.individual_test_specification, rdt_test_spec_vec),
    poisson_noise_optimal_solutions,
);

elisa_optimal_solutions = filter(
    chars -> in(chars.individual_test_specification, elisa_test_spec_vec),
    dynamical_noise_optimal_solutions,
);

#%%
elisa_df = create_optimal_thresholds_df(elisa_optimal_solutions)
elisa_df[!, :noise_spec] .= "Common"

rdt_dynamical_df = create_optimal_thresholds_df(
    dynamical_noise_rdt_optimal_solutions
)
rdt_dynamical_df[!, :noise_spec] .= "Dynamical noise: in-phase"

rdt_poisson_df = create_optimal_thresholds_df(
    poisson_noise_rdt_optimal_solutions
)
rdt_poisson_df[!, :noise_spec] .= "Poisson noise"

thresholds_df = vcat(
    rdt_dynamical_df, rdt_poisson_df, elisa_df
)

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
CSV.write(projectdir("manuscript/optimal-thresholds.csv"), wide_thresholds_df);

#%%
wide_accuracy_df = create_wide_df(thresholds_df, :accuracy; digits = 2)
CSV.write(
    projectdir("manuscript/optimal-thresholds_accuracy.csv"), wide_accuracy_df
);
