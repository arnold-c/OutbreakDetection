#%%
using DrWatson: DrWatson
using OutbreakDetectionCore: OutbreakDetectionCore
using OutbreakDetection
using Try: Try
using StructArrays: StructVector


# Make sure these values are present in the optimization script
alert_method = OutbreakDetectionCore.AlertMethod(OutbreakDetectionCore.MovingAverage(7))
accuracy_metric = OutbreakDetectionCore.AccuracyMetric(OutbreakDetectionCore.BalancedAccuracy())
threshold_bounds = (; lower = 0.0, upper = 20.0)
alert_filtering_strategy = OutbreakDetectionCore.AlertFilteringStrategy(OutbreakDetectionCore.AllAlerts())
plotdirpath = DrWatson.plotsdir()

#%%
optimized_filedir = OutbreakDetectionCore.outdir("ensemble", "threshold-optimization")
optimized_threshold_results = Try.unwrap(
    OutbreakDetectionCore.load_previous_optimization_results_structvector(
        optimized_filedir,
        "threshold-optimization.jld2",
        joinpath(optimized_filedir, "checkpoints");
        default_action = :continue
    )
)

#%%
filtered_results = filter(
    r -> r.alert_method == alert_method &&
        r.accuracy_metric == accuracy_metric &&
        r.threshold_bounds == threshold_bounds &&
        r.alert_filtering_strategy == alert_filtering_strategy,
    optimized_threshold_results
)

perfect_test_specs = [
    OutbreakDetectionCore.IndividualTestSpecification(1.0, 1.0, 14),
    OutbreakDetectionCore.IndividualTestSpecification(1.0, 1.0, 0),
]

rdt_test_specs = [
    OutbreakDetectionCore.IndividualTestSpecification(0.9, 0.9, 0),
    OutbreakDetectionCore.IndividualTestSpecification(0.85, 0.85, 0),
]

function filter_by_noise_and_tests(results, noise_type, test_specs)
    return filter(
        r -> r.noise_type_description == noise_type &&
            any(==(r.test_specification), test_specs),
        results,
    )
end

function group_results_by_noise_level(results)
    noise_levels = sort(unique(results.noise_level))
    grouped_results = map(noise_levels) do noise_level
        StructVector(filter(r -> r.noise_level == noise_level, results))
    end
    return grouped_results, noise_levels
end

all_perfect_test_optimal_solutions, static_noise_levels =
    filter_by_noise_and_tests(filtered_results, :static, perfect_test_specs) |>
    group_results_by_noise_level

all_static_noise_rdt_test_optimal_solutions, static_rdt_noise_levels =
    filter_by_noise_and_tests(filtered_results, :static, rdt_test_specs) |>
    group_results_by_noise_level

all_dynamical_noise_rdt_test_optimal_solutions, dynamic_noise_levels =
    filter_by_noise_and_tests(filtered_results, :dynamic, rdt_test_specs) |>
    group_results_by_noise_level

@assert length(all_perfect_test_optimal_solutions) ==
    length(all_static_noise_rdt_test_optimal_solutions) ==
    length(all_dynamical_noise_rdt_test_optimal_solutions)
@assert static_noise_levels == static_rdt_noise_levels
@assert static_noise_levels == dynamic_noise_levels

#%%
colors = [:red, :blue, :green, :magenta, :cyan]
for i in eachindex(all_perfect_test_optimal_solutions)
    noise_level = static_noise_levels[i]
    color = colors[mod1(i, length(colors))]

    println("\n", "="^80)
    println("Comparison: Perfect test + Static-noise RDT")
    println("="^80)
    compare_optimal_solution_mean_extrema(
        [all_perfect_test_optimal_solutions[i]],
        [all_static_noise_rdt_test_optimal_solutions[i]],
        [
            # :detection_delays,
            # :alert_durations,
            :accuracies,
        ],
        [
            OutbreakDetectionCore.IndividualTestSpecification(1.0, 1.0, 14),
            OutbreakDetectionCore.IndividualTestSpecification(1.0, 1.0, 0),
            OutbreakDetectionCore.IndividualTestSpecification(0.9, 0.9, 0),
            OutbreakDetectionCore.IndividualTestSpecification(0.85, 0.85, 0),
        ];
        digits = 2,
        color = color,
    )

    println("\n", "="^80)
    println("Comparison: Perfect test + Dynamical-noise RDT")
    println("="^80)
    compare_optimal_solution_mean_extrema(
        [all_perfect_test_optimal_solutions[i]],
        [all_dynamical_noise_rdt_test_optimal_solutions[i]],
        [
            # :detection_delays,
            # :alert_durations,
            :accuracies,
        ],
        [
            OutbreakDetectionCore.IndividualTestSpecification(1.0, 1.0, 14),
            OutbreakDetectionCore.IndividualTestSpecification(1.0, 1.0, 0),
            OutbreakDetectionCore.IndividualTestSpecification(0.9, 0.9, 0),
            OutbreakDetectionCore.IndividualTestSpecification(0.85, 0.85, 0),
        ];
        digits = 2,
        color = color,
    )
end
