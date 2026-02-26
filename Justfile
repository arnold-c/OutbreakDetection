default: optimizations plots_worker supplemental-plots_worker

julia_script := "julia --project=. --startup-file=no --history-file=no"



optimizations:
	{{ julia_script }}  ./scripts/optimal-thresholds_optims.jl

plots: optimizations plots_worker
plots_worker:
	{{ julia_script }} ./scripts/optimal-thresholds_plots.jl

supplemental-plots: optimizations supplemental-plots_worker
supplemental-plots_worker:
	{{ julia_script }} ./scripts/optimal-thresholds_supplement-plots.jl

supplemental-tables: optimizations supplemental-tables_worker
supplemental-tables_worker:
	{{ julia_script }} ./scripts/optimal-thresholds_supplement-tables.jl

plot-timeseries: optimizations plot-timeseries_worker
plot-timeseries_worker:
	{{ julia_script }} ./scripts/visualize-timeseries.jl

plot-schematic: plot-schematic_worker
plot-schematic_worker:
	{{ julia_script }} ./scripts/schematic-plot.jl

optimizations-comparison: optimizations optimization-comparison_worker
optimization-comparison_worker:
	{{ julia_script }} ./scripts/optimal-thresholds_comparisons.jl











manuscript: optimizations \
	plots_worker \
	supplemental-tables_worker \
	supplemental-plots_worker \
	plot-schematic_worker \
	optimization-comparison_worker
	typst compile ./manuscript/combined-manuscript.typ

tests:
	julia --project=./OutbreakDetectionCore -e "using Pkg; Pkg.test()"
	julia --project=. -e "using Pkg; Pkg.test()"
