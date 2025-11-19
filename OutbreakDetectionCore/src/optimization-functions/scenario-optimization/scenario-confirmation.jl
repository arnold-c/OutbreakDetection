export check_missing_scenario_optimizations

"""
    check_missing_scenario_optimizations(optim_df, combinations_to_run; 
                                        disable_time_check=false, 
                                        time_per_run_s=45, 
                                        scenario_parameter_symbols=...)

Check for missing scenario optimizations and prompt user to continue.

Returns Try.Ok(missing_combinations) or Try.Err(message).
"""
function check_missing_scenario_optimizations(
        optim_df,
        combinations_to_run;
        disable_time_check = false,
        time_per_run_s = 45,
        scenario_parameter_symbols = [
            :ensemble_spec,
            :outbreak_spec,
            :noise_spec,
            :outbreak_detection_spec,
            :test_spec,
            :optimization_method,
            :accuracy_function,
        ],
    )
    missing_combinations = DataFrames.antijoin(
        combinations_to_run,
        optim_df;
        on = scenario_parameter_symbols,
    )

    missing_runs = DataFrames.nrow(missing_combinations)
    if missing_runs == 0
        return Try.Err("No missing simulations")
    end

    if missing_runs > 0 && !disable_time_check
        nrun_time_s = missing_runs * time_per_run_s
        nrun_time_minutes = round(nrun_time_s / 60; digits = 2)
        nrun_time_hours = divrem(nrun_time_minutes, 60)
        nrun_time_message = if nrun_time_s < 10
            "less than 10 seconds"
        elseif nrun_time_s < 60
            "approximately $(round(nrun_time_s; digits = 0)) seconds"
        elseif nrun_time_minutes < 120
            "approximately $(nrun_time_minutes) minutes"
        else
            "approximately $(nrun_time_hours[1]) hours and $(nrun_time_hours[2]) minutes"
        end
        choice = request(
            "There are $(missing_runs) missing simulations. This is estimated to take $(nrun_time_message). Do you want to continue?",
            RadioMenu(["No", "Yes"]; ctrl_c_interrupt = false),
        )

        if choice != 2
            return Try.Err("User aborted")
        end

        println("Continuing ...")
    end

    return Try.Ok(missing_combinations)
end
