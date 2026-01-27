#%%
using DrWatson
@quickactivate "OutbreakDetection"

using OutbreakDetection: visualize_timeseries
using OutbreakDetectionCore
using Try
using CairoMakie

#%%
# Load optimized results with default action to continue despite struct changes
optimized_filedir = OutbreakDetectionCore.outdir("ensemble", "threshold-optimization")
optimized_threshold_results = Try.unwrap(
    OutbreakDetectionCore.load_previous_optimization_results_structvector(
        optimized_filedir,
        "threshold-optimization.jld2",
        joinpath(optimized_filedir, "checkpoints");
        default_action = :continue
    )
);

println("Loaded $(length(optimized_threshold_results)) optimization results")
println("First result has $(optimized_threshold_results[1].ensemble_specification.nsims) simulations")

#%%
dynamic_optimized_results = filter(res -> res.noise_type_description == :dynamic, optimized_threshold_results)
for noise_level in unique(optimized_threshold_results.noise_level)
    noise_filtered_results = filter(res -> res.noise_level == noise_level, dynamic_optimized_results)
    noise_plotpath = DrWatson.plotsdir("noise-type_$(string(noise_type))", "noise-level_$noise_level")
    mkpath(noise_plotpath)
    for test in unique(noise_filtered_results.test_specification), percent_tested in unique(noise_filtered_results.percent_tested)
        test_filtered_results = filter(
            res -> res.test_specification == test && res.percent_tested == percent_tested,
            noise_filtered_results
        )
        test_type_description = replace(
            get_test_description(test),
            r"[\s()]+" => "-"
        ) |> r -> rstrip(r, '-')
        test_plotpath = joinpath(noise_plotpath, "test-type_$test_type_description", "percent-tested_$percent_tested")
        mkpath(test_plotpath)
        for accuracy_metric in unique(test_filtered_results.accuracy_metric)
            accuracy_filtered_results = filter(
                res -> res.accuracy_metric == accuracy_metric,
                test_filtered_results
            )
            accuracy_metric_description = replace(
                string(accuracy_metric),
                r"AccuracyMetric\((.*?)\(\)\)" => s"\1"
            )
            accuracy_plotpath = joinpath(test_plotpath, "accuracy_metric_$accuracy_metric_description")
            mkpath(accuracy_plotpath)
            @assert length(accuracy_filtered_results) == 1
            for sim in [1, 10, 49, 100]
                println("Creating visualization for simulation $sim...")
                sim_fig = visualize_timeseries(
                    accuracy_filtered_results[1];
                    sim_number = sim,
                )

                save(
                    joinpath(accuracy_plotpath, "sim-$(sim)_timeseries.svg"),
                    sim_fig
                )

                println("Visualization created successfully!")
            end
        end
    end
end
