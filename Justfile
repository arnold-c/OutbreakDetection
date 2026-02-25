default: optimizations plots_worker supplemental-plots_worker

julia_script := "julia --project=. --startup-file=no --history-file=no"



optimizations:
	{{ julia_script }}  ./scripts/optimal-thresholds_optims.jl

plots: optimizations plots_worker
plots_worker:
	{{ julia_script }} ./scripts/optimal-thresholds_plots.jl

supplemental-plots: optimizations supplemental-plots_worker
supplemental-plots_worker:
	{{ julia_script }} ./scripts/supplemental_plots.jl

tables: optimizations tables_worker
tables_worker:
	{{ julia_script }} ./scripts/optimal-thresholds_tables.jl

plot-timeseries: optimizations plot-timeseries_worker
plot-timeseries_worker:
	{{ julia_script }} ./scripts/visualize-timeseries.jl

optimizations-comparison: optimizations optimization-comparison_worker
optimization-comparison_worker:
	{{ julia_script }} ./scripts/optimal-thresholds_comparisons.jl











manuscript: optimizations plots_worker tables_worker plot-timeseries_worker optimization-comparison_worker
	typst compile ./manuscript/combined-manuscript.typ

tests:
	julia --project=./OutbreakDetectionCore -e "using Pkg; Pkg.test()"
	julia --project=. -e "using Pkg; Pkg.test()"
