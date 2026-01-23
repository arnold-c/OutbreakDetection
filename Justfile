default: tests manuscript

julia_script := "julia --project=. --startup-file=no --history-file=no"

optimizations:
	{{ julia_script }}  ./scripts/optimal-thresholds_optims.jl

plots: optimizations
	{{ julia_script }} ./scripts/optimal-thresholds_plots.jl


manuscript:
	julia ./scripts/manuscript/optimal-thresholds.jl
	typst compile ./manuscript/combined-manuscript.typ

tests:
	julia --project=./OutbreakDetectionCore -e "using Pkg; Pkg.test()"
	julia --project=. -e "using Pkg; Pkg.test()"
