module OutbreakDetection

include("plotting-functions.jl")
export single_seir_plot

@static if false
    include("../scripts/single-sim_plots.jl")
    include("../scripts/ensemble-diag-testing_scenarios_plots.jl")
end

end
