#%%
using DrWatson
@quickactivate "OutbreakDetection"

using OutbreakDetectionUtils
using OutbreakDetection
using Optimization
using OptimizationBBO

include(srcdir("ensemble-parameters.jl"))

if false
    include("../src/ensemble-parameters.jl")
end

#%%
using DrWatson
@quickactivate "OutbreakDetection"


#%%
model_types_vec = [("seasonal-infectivity-import", "tau-leaping")]
seed = 1234

#%%
ensemble_dynamics_spec = DynamicsParameters(
        ensemble_state_specification.init_states.N,
        27,
        0.2,
        SIGMA,
        GAMMA,
        16.0,
        0.8,
    )


ensemble_spec = EnsembleSpecification(
	ensemble_model_type,
	ensemble_state_specification,
	ensemble_dynamics_spec,
	ensemble_time_specification,
	ensemble_nsims,
)

outbreak_spec = OutbreakSpecification(5, 30, 500)

base_param_dict = @dict(
    ensemble_spec = ensemble_spec,
	outbreak_spec = outbreak_spec,
    seed = seed,
)

#%%
ensemble_inc_arr, thresholds_vec = setup_optimization(base_param_dict)

#%%
noise_spec = PoissonNoiseSpecification(1.0)
individual_test_specification = IndividualTestSpecification(0.85, 0.85, 0)
outbreak_detection_specification = OutbreakDetectionSpecification(
	2,
	7,
	0.6,
	0.6,
	"movingavg",
)


scenario_spec = ScenarioSpecification(
	ensemble_spec,
	outbreak_spec,
	noise_spec,
	outbreak_detection_specification,
	individual_test_specification
)

OT_chars_param_dict = (;
	scenario_spec,
	ensemble_inc_arr,
	thresholds_vec,
	seed = seed
)

#%%
noise_array = create_noise_arr(
	noise_spec,
	ensemble_inc_arr;
	ensemble_specification = scenario_spec.ensemble_specification,
	seed = seed,
)[1]

output_df = DataFrames.DataFrame(
	"threshold" => Float64[],
	"loss" => Float64[]
)

obj_inputs = (;
	ensemble_inc_arr,
	noise_array,
	outbreak_detection_specification,
	individual_test_specification,
	thresholds_vec,
	output_df
)

f = Optimization.OptimizationFunction(objective_function)

prob = Optimization.OptimizationProblem(
	f,
	[2.0],
	obj_inputs;
	lb = [0.5],
	ub = [50.0],
)

function callback!(state, loss)

	return false
end

sol = Optimization.solve(
	prob,
	OptimizationBBO.BBO_adaptive_de_rand_1_bin_radiuslimited();
	maxiters = 1000,
	maxtime = 30.0,
	abstol = 0.01,
)

#%%
# run_optimization(
# 	OT_chars_param_dict;
# 	guess_threshold = 2.0,
# 	threshold_lower_bound = 0.5,
# 	threshold_upper_bound = 50.0,
# 	maxiters = 100000,
# 	maxtime	= 10.0,
# )

