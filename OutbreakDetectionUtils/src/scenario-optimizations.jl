using DataFrames: DataFrames
using DrWatson: @dict, @tagsave
using DataFramesMeta: exec
using StructArrays: StructVector
using ProgressMeter: Progress, next!
using FLoops: FLoops
using Try: Try
using TryExperimental: TryExperimental
using REPL: REPL
using REPL.TerminalMenus: RadioMenu, request
using StyledStrings
using Dates: Dates
using Match: Match

function run_scenario_optimizations(
    ensemble_specifications,
    outbreak_specifications,
    noise_specifications,
    outbreak_detection_specifications,
    individual_test_specifications,
    optim_method::TMethod = MSO,
    accuracy_functions = [arithmetic_mean, calculate_f_beta_score];
    seed = 1234,
    executor = FLoops.SequentialEx(),
    filedir = outdir("ensemble", "scenario-optimizations"),
    optimization_filename_base = "alert-threshold-optimization.jld2",
    optimization_output_filepath = joinpath(
        filedir,
        string(Dates.now()) * "_" * optimization_filename_base,
    ),
    force = false,
    return_df = true,
    disable_time_check = false,
    time_per_run_s = 98,
    kwargs...,
) where {TMethod<:Type{<:OptimizationMethods}}
    if !isdir(filedir)
        mkpath(filedir)
    end

    load_filepath = get_most_recent_optimization_filepath(
        optimization_filename_base,
        filedir,
    )

    optim_df = nothing

    if Try.isok(load_filepath) && !force
        optim_df = JLD2.load(Try.unwrap(load_filepath))["threshold_optim_df"]
    else
        optim_df = DataFrames.DataFrame((
            ensemble_spec = typeof(ensemble_specifications[1])[],
            outbreak_spec = typeof(outbreak_specifications[1])[],
            noise_spec = Union{
                PoissonNoiseSpecification,DynamicalNoiseSpecification
            }[],
            outbreak_detection_spec = typeof(
                outbreak_detection_specifications[1]
            )[],
            test_spec = typeof(individual_test_specifications[1])[],
            optimal_threshold = Float64[],
            optimal_accuracy = Float64[],
            optimization_method = Union{Type{QD},Type{MSO}}[],
            accuracy_function = Symbol[],
            OT_chars = StructVector{<:OutbreakThresholdChars}[],
        ))
    end

    missing_optimizations = check_missing_scenario_optimizations(
        optim_df,
        ensemble_specifications,
        outbreak_specifications,
        noise_specifications,
        outbreak_detection_specifications,
        individual_test_specifications,
        optim_method,
        accuracy_functions;
        disable_time_check = disable_time_check,
        time_per_run_s = time_per_run_s,
    )

    if Try.iserr(missing_optimizations)
        println(Try.unwrap_err(missing_optimizations))
        if return_df
            println("Returning previously computed threshold_optim_df")
            return optim_df
        end
        return nothing
    end

    missing_optims_df = Try.unwrap(missing_optimizations)

    if DataFrames.nrow(missing_optims_df) == 0
        @info "🟨 Previous optimal ews_df is the same as the current. No need to save a new version. 🟨"
    end

    run_missing_scenario_optimizations!(
        optim_df,
        missing_optims_df;
        seed = seed,
        executor = executor,
        kwargs...,
    )

    @tagsave(
        optimization_output_filepath,
        Dict("threshold_optim_df" => optim_df)
    )
    @info "🟢 Saved optimal ews_df to $(optimization_output_filepath) 🟢"

    if return_df
        return optim_df
    end

    return nothing
end

function check_missing_scenario_optimizations(
    optim_df,
    ensemble_specifications,
    outbreak_specifications,
    noise_specifications,
    outbreak_detection_specifications,
    individual_test_specifications,
    optim_method::TMethod = MSO,
    accuracy_functions = [arithmetic_mean, calculate_f_beta_score];
    disable_time_check = false,
    time_per_run_s = 45,
) where {TMethod<:Type{<:OptimizationMethods}}
    scenario_parameter_symbols = [
        :ensemble_spec,
        :outbreak_spec,
        :noise_spec,
        :outbreak_detection_spec,
        :test_spec,
        :optimization_method,
        :accuracy_function,
    ]

    combinations_to_run = DataFrames.DataFrame(
        Iterators.product(
            ensemble_specifications,
            outbreak_specifications,
            noise_specifications,
            outbreak_detection_specifications,
            individual_test_specifications,
            [optim_method],
            nameof.(accuracy_functions),
        ),
        scenario_parameter_symbols,
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

function run_missing_scenario_optimizations!(
    optim_df,
    missing_optims_df;
    seed = 1234,
    executor = FLoops.SequentialEx(),
    kwargs...,
)
    prog = Progress(DataFrames.nrow(missing_optims_df))

    incidence_sim_grouped_df = DataFrames.groupby(
        missing_optims_df,
        [:ensemble_spec, :outbreak_spec],
    )

    lk = ReentrantLock()
    for inc_gp in incidence_sim_grouped_df
        ensemble_spec = inc_gp[1, :ensemble_spec]
        outbreak_spec = inc_gp[1, :outbreak_spec]

        base_param_dict = @dict(
            ensemble_spec = ensemble_spec,
            outbreak_spec = outbreak_spec,
            seed = seed,
        )

        ensemble_inc_arr, thresholds_vec = setup_optimization(
            base_param_dict
        )

        FLoops.@floop executor for noise_gp in
                                   DataFrames.groupby(inc_gp, [:noise_spec])
            noise_spec = noise_gp[1, :noise_spec]

            println(
                styled"{green:\n=================================================================}"
            )
            println(
                styled"Noise type: {green,inverse: $(getdirpath(noise_spec))}"
            )

            noise_array, noise_means = create_noise_arr(
                noise_spec,
                ensemble_inc_arr;
                ensemble_specification = ensemble_spec,
                seed = seed,
            )

            for detect_test_gp in DataFrames.groupby(
                noise_gp, [:outbreak_detection_spec, :test_spec]
            )
                outbreak_detection_spec = detect_test_gp[
                    1, :outbreak_detection_spec
                ]
                individual_test_spec = detect_test_gp[1, :test_spec]

                println(
                    styled"\t\t-> Test specification: {blue: $(get_test_description(individual_test_spec))}, Percent tested: {red,inverse: $(outbreak_detection_spec.percent_tested)}"
                )

                for accuracy_function_df in DataFrames.groupby(
                    detect_test_gp, [:accuracy_function]
                )
                    accuracy_function_symbol = accuracy_function_df[
                        1, :accuracy_function
                    ]
                    accuracy_function = Match.@match accuracy_function_symbol begin
                        :arithmetic_mean => arithmetic_mean
                        :calculate_f_beta_score => calculate_f_beta_score
                    end

                    println(
                        styled"\t\t\t\t-> Accuracy Function: {yellow,inverse: $(accuracy_function)}"
                    )

                    obj_inputs = (;
                        ensemble_inc_arr,
                        noise_array,
                        outbreak_detection_specification = outbreak_detection_spec,
                        individual_test_specification = individual_test_spec,
                        thresholds_vec,
                        accuracy_function,
                    )

                    objective_function_closure =
                        x -> objective_function(x, obj_inputs)

                    @assert length(
                        unique(detect_test_gp[:, :optimization_method])
                    ) == 1

                    optim_method = detect_test_gp[1, :optimization_method]

                    optim_minimizer, optim_minimum = optimization_wrapper(
                        objective_function_closure,
                        optim_method;
                        kwargs...,
                    )

                    optimal_outbreak_detection_spec = OutbreakDetectionSpecification(
                        optim_minimizer,
                        outbreak_detection_spec.moving_average_lag,
                        outbreak_detection_spec.percent_visit_clinic,
                        outbreak_detection_spec.percent_clinic_tested,
                        outbreak_detection_spec.alert_method.method_name,
                    )

                    testarr = create_testing_arrs(
                        ensemble_inc_arr,
                        noise_array,
                        optimal_outbreak_detection_spec,
                        individual_test_spec,
                    )[1]

                    OT_chars = calculate_OutbreakThresholdChars(
                        testarr, ensemble_inc_arr, thresholds_vec, noise_means
                    )

                    lock(lk)
                    push!(
                        optim_df,
                        (
                            ensemble_spec,
                            outbreak_spec,
                            noise_spec,
                            outbreak_detection_spec,
                            individual_test_spec,
                            optim_minimizer,
                            1 - optim_minimum,
                            optim_method,
                            accuracy_function_symbol,
                            OT_chars,
                        ),
                    )
                    unlock(lk)

                    next!(prog)
                    println("")
                end
            end
        end
    end

    return nothing
end

function get_most_recent_optimization_filepath(
    filename_base,
    filedir,
)
    @assert isdir(filedir)
    optimization_files = readdir(filedir)

    if length(optimization_files) == 0
        return Try.Err("No optimization files found.")
    end

    filter_regex = Regex("(.*)$(filename_base)\$")

    filtered_optimization_files = filter(
        f -> contains(f, filter_regex),
        optimization_files,
    )

    if length(filtered_optimization_files) == 0
        return Try.Err("No optimization files found.")
    end

    filtered_optimization_datetimes = Vector{Union{Try.Ok,Try.Err}}(
        undef, length(filtered_optimization_files)
    )

    for (i, f) in pairs(filtered_optimization_files)
        matches = match(filter_regex, f)
        if isnothing(matches)
            filtered_optimization_datetimes[i] = Try.Err(
                "No matches for filename $(f)"
            )
            continue
        end

        filtered_optimization_datetimes[i] = Try.Ok(
            tryparse(
                Dates.DateTime,
                strip(
                    matches[1],
                    '_',
                ),
            ),
        )
    end

    filtered_optimization_datetimes = filter(
        Try.isok, filtered_optimization_datetimes
    )

    if length(filtered_optimization_datetimes) == 0
        return Try.Err("No optimization files found.")
    end

    most_recent_optimization_datetime = sort(
        Try.unwrap.(filtered_optimization_datetimes)
    )[end]

    most_recent_filepath = joinpath(
        filedir,
        string(most_recent_optimization_datetime) *
        "_$(filename_base)",
    )
    return Try.Ok(most_recent_filepath)
end
