default: tests manuscript

optimizations:
	julia --project=. --startup-file=no --history-file=no ./scripts/optimal-thresholds_optims.jl

manuscript:
	julia ./scripts/manuscript/optimal-thresholds.jl
	typst compile ./manuscript/combined-manuscript.typ

tests:
	julia --project=./OutbreakDetectionUtils -e "using Pkg; Pkg.test()"
