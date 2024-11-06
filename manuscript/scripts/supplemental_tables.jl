#%%
wide_accuracy_df = create_wide_df(thresholds_df, :accuracy; digits = 2)
CSV.write(
    appendix_tabledir("optimal-thresholds_accuracy.csv"), wide_accuracy_df
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
    appendix_tabledir("optimal-thresholds_unavoidable-cases.csv"),
    wide_unavoidable_df,
);

#%%
wide_delays_df = create_wide_df(thresholds_df, :detectiondelays; digits = 2)
CSV.write(
    appendix_tabledir("optimal-thresholds_detection-delays.csv"),
    wide_delays_df,
);
