default: tests manuscript

manuscript:
	julia ./manuscript/scripts/optimal-thresholds.jl
	typst compile ./manuscript/combined-manuscript.typ

tests:
	julia --project=./OutbreakDetectionUtils -e "using Pkg; Pkg.test()"
