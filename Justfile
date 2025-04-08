default: tests manuscript

manuscript:
	julia ./manuscript/scripts/optimal-thresholds.jl

tests:
	julia ./OutbreakDetectionUtils/test/runtests.jl
