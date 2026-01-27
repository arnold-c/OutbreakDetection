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



plot-timeseries: optimizations plot-timeseries_worker

plot-timeseries_worker:
	{{ julia_script }} ./scripts/visualize-timeseries.jl



manuscript:
	julia ./scripts/manuscript/optimal-thresholds.jl
	typst compile ./manuscript/combined-manuscript.typ

tests:
	julia --project=./OutbreakDetectionCore -e "using Pkg; Pkg.test()"
	julia --project=. -e "using Pkg; Pkg.test()"
