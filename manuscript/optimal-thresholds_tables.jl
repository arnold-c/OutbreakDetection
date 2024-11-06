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

gha_2022_scale_population_per_annum = gha_2022_scale_population / 100

#%%
dynamical_noise_rdt_optimal_solutions = filter(
    chars -> in(chars.individual_test_specification, rdt_test_spec_vec),
    dynamical_noise_optimal_solutions,
);

poisson_noise_rdt_optimal_solutions = filter(
    chars -> in(chars.individual_test_specification, rdt_test_spec_vec),
    poisson_noise_optimal_solutions,
);

# dynamical_noise_elisa_optimal_solutions = filter(
#     chars -> in(chars.individual_test_specification, elisa_test_spec_vec),
#     dynamical_noise_optimal_solutions,
# );
#
# poisson_noise_elisa_optimal_solutions = filter(
#     chars -> in(chars.individual_test_specification, elisa_test_spec_vec),
#     poisson_noise_optimal_solutions,
# );

perfect_test_optimal_solutions = filter(
    chars -> in(chars.individual_test_specification, perfect_test_spec_vec),
    poisson_noise_optimal_solutions,
);

#%%
perfect_test_df = create_optimal_threshold_summary_df(
    perfect_test_optimal_solutions,
    [:unavoidable_cases, :detectiondelays, :alert_duration_vec];
    percentiles = nothing,
    nboots = nothing,
)
perfect_test_df[!, :noise_spec] .= "All noise structures"
perfect_test_df

rdt_dynamical_df = create_optimal_threshold_summary_df(
    dynamical_noise_rdt_optimal_solutions,
    [:unavoidable_cases, :detectiondelays, :alert_duration_vec];
    percentiles = nothing,
    nboots = nothing,
)
rdt_dynamical_df[!, :noise_spec] .= "Dynamical noise: in-phase"

rdt_poisson_df = create_optimal_threshold_summary_df(
    poisson_noise_rdt_optimal_solutions,
    [:unavoidable_cases, :detectiondelays, :alert_duration_vec];
    percentiles = nothing,
    nboots = nothing,
)
rdt_poisson_df[!, :noise_spec] .= "Poisson noise"

# elisa_dynamical_df = create_optimal_threshold_summary_df(
#     dynamical_noise_elisa_optimal_solutions,
#     [:unavoidable_cases, :detectiondelays, :alert_duration_vec];
#     percentiles = nothing,
#     nboots = nothing,
# )
# elisa_dynamical_df[!, :noise_spec] .= "Dynamical noise: in-phase"
#
# elisa_poisson_df = create_optimal_threshold_summary_df(
#     poisson_noise_elisa_optimal_solutions,
#     [:unavoidable_cases, :detectiondelays, :alert_duration_vec];
#     percentiles = nothing,
#     nboots = nothing,
# )
# elisa_poisson_df[!, :noise_spec] .= "Poisson noise"

thresholds_df = vcat(
    rdt_dynamical_df,
    rdt_poisson_df,
    # elisa_dynamical_df,
    # elisa_poisson_df,
    perfect_test_df,
)

thresholds_df[!, :unavoidable_cases] =
    Int64.(
        round.(
            thresholds_df[!, :unavoidable_cases] *
            gha_2022_scale_population_per_annum;
            digits = 0,
        )
    )

#%%
function create_wide_df(
    long_df,
    outcome::Symbol;
    noise_order = [
        "Dynamical noise: in-phase", "Poisson noise", "All noise structures"
    ],
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
        (0.98, 0.98, x::Int) => "ELISA Equivalent ($(sensitivity * 100)%)"
        (x::AbstractFloat, x::AbstractFloat, 0) where {x<1.0} => "RDT Equivalent ($(sensitivity * 100)%)"
        (1.0, 1.0, x::Int) => "Perfect Test"
    end
end

#%%
wide_thresholds_df = create_wide_df(thresholds_df, :alert_threshold; digits = 0)
CSV.write(
    manuscript_tabledir("optimal-thresholds.csv"),
    wide_thresholds_df,
);

#%%
wide_accuracy_df = create_wide_df(thresholds_df, :accuracy; digits = 2)
CSV.write(
    manuscript_tabledir("optimal-thresholds_accuracy.csv"), wide_accuracy_df
);

#%%
wide_unavoidable_df = create_wide_df(
    thresholds_df, :unavoidable_cases; digits = 2
)
transform!(
    wide_unavoidable_df,
    Not(r".*Type") .=> (x -> Int64.(x));
    renamecols = false,
)
CSV.write(
    manuscript_tabledir("optimal-thresholds_unavoidable-cases.csv"),
    wide_unavoidable_df,
);

#%%
wide_delays_df = create_wide_df(thresholds_df, :detectiondelays; digits = 2)
CSV.write(
    manuscript_tabledir("optimal-thresholds_detection-delays.csv"),
    wide_delays_df,
);
