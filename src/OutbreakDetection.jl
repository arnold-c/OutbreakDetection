"""
Note that each file is a module that defines the exported functions.
They are just listed here for convenience of sourcing one file.
"""
module OutbreakDetection


include("dynamics-constants.jl")
export POPULATION_N, LATENT_PER_DAYS, DUR_INF_DAYS, R0, SIGMA, GAMMA,
    LIFE_EXPECTANCY_YEARS, ANNUAL_BIRTHS_PER_K, VACCINATION_COVERAGE,
    MU, BETA_MEAN, BETA_FORCE, EPSILON

include("test-constants.jl")
export CLINICAL_CASE_TEST_SPEC, EPI_LINKED_CASE_TEST_SPEC, CLINICAL_TEST_SPECS

@static if false
    include("../scripts/single-sim.jl")
    include("../scripts/ensemble-sim.jl")
    include("../scripts/ensemble-sim_single-scenario.jl")
    include("../scripts/ensemble-sim_noise-visualizations.jl")
    include("../scripts/ensemble-diag-testing_scenarios_plots.jl")
    include("../scripts/ensemble-diag-testing_optimal-thresholds.jl")
    include("../scripts/ensemble-diag-testing_constant-thresholds.jl")
    include("../scripts/ensemble-diag-testing_optimal-thresholds_single-timeseries.jl")
    include("../scripts/single-sim_plots.jl")
    # include("../scripts/debugging.jl")
    include("../scripts/outbreak-detection-schematic.jl")
end

end
