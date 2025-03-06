#%%
using DrWatson
@quickactivate "OutbreakDetection"

using OutbreakDetectionUtils
using OutbreakDetection
using Optimization
using OptimizationBBO
using DataFrames
using GLMakie

include(srcdir("makie-plotting-setup.jl"))

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
# noise_spec = PoissonNoiseSpecification(8.0)
noise_spec = DynamicalNoiseSpecification(
	"dynamical",
	5.0,
	7,
	14,
	"in-phase",
	0.15,
	0.0492 #0.8538, 0.7383, 0.5088, 0.2789, 0.0492
)

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
    individual_test_specification,
)

# OT_chars_param_dict = (;
#     scenario_spec,
#     ensemble_inc_arr,
#     thresholds_vec,
#     seed = seed,
# )

#%%
noise_array = create_noise_arr(
    noise_spec,
    ensemble_inc_arr;
    ensemble_specification = scenario_spec.ensemble_specification,
    seed = seed,
)[1]

#%%
output_df = DataFrames.DataFrame(
    "threshold" => Float64[],
    "loss" => Float64[],
)

obj_inputs = (;
    ensemble_inc_arr,
    noise_array,
    outbreak_detection_specification,
    individual_test_specification,
    thresholds_vec,
    output_df,
)

# f = Optimization.OptimizationFunction(
# 	objective_function,
# 	Optimization.AutoForwardDiff()
# )
#
# prob = Optimization.OptimizationProblem(
#     f,
#     [3.0],
#     obj_inputs;
#     lb = [0.5],
#     ub = [50.0],
# )
#
# function callback!(state, loss; doplot = true)
# 	push!(threshold_trace, deepcopy(state.u[1]))
# 	push!(obj_trace, deepcopy(state.objective))
#     return false
# end
#
# threshold_trace = Float64[]
# obj_trace = Float64[]
# sol = Optimization.solve(
#     prob,
# 	Optim.LBFGS(),
# 	# NLopt.GN_DIRECT(),
# 	# QuadDirect(),
#     # OptimizationBBO.BBO_adaptive_de_rand_1_bin_radiuslimited();
#     maxiters = 10000,
#     maxtime = 240.0,
#     abstol = 0.02,
# 	# splits = ([0.5, 2.0, 5.0]),
# 	callback = callback!
# )
#
#
@showprogress for threshold in collect(0.1:0.1:30.0)
	if length(filter(t -> t.==threshold, output_df.threshold)) > 0
	       continue
	   end

    loss = objective_function(
        [threshold],
        obj_inputs,
    )

    DataFrames.push!(
        output_df,
        (
        	threshold = threshold,
        	loss = loss
        );
        cols = :union,
    )
end

#%%
sort!(output_df, :threshold)
scatter(output_df.threshold, output_df.loss; markersize = 20, alpha = 0.2)
# scatter!(threshold_trace, obj_trace, color = :red)
#
# lines(eachindex(threshold_trace), threshold_trace)
#
# lines(eachindex(obj_trace), obj_trace)

#%%
# run_optimization(
# 	OT_chars_param_dict;
# 	guess_threshold = 2.0,
# 	threshold_lower_bound = 0.5,
# 	threshold_upper_bound = 50.0,
# 	maxiters = 100000,
# 	maxtime	= 10.0,
# )


#%%
optim_sol = Optim.optimize(
	t -> objective_function(t, obj_inputs),
	0.0,
	100.0,
	Brent();
	store_trace = true
)
Optim.minimizer(optim_sol)
Optim.minimum(optim_sol)

#%%
optim_sol2 = Optim.optimize(
	t -> objective_function(t, obj_inputs),
	0.0,
	100.0,
	GoldenSection();
	store_trace = true
)
Optim.minimizer(optim_sol2)
Optim.minimum(optim_sol2)

#%%
optim_sol3 = Optim.optimize(
	t -> objective_function(t, obj_inputs),
	[0.0],
	[100.0],
	[7.0],
	SAMIN(),
	Optim.Options(
		store_trace = true,
		extended_trace = true,
		iterations = 10^6
	)
)

Optim.minimizer(optim_sol3)
Optim.minimum(optim_sol3)

#%%
optim_sol4 = Optim.optimize(
	t -> objective_function(t, obj_inputs),
	[7.0],
	NelderMead(),
	Optim.Options(
		store_trace = true,
		extended_trace = true,
		iterations = 10^6
	)
)

Optim.minimizer(optim_sol4)
Optim.minimum(optim_sol4)

#%%
optim_sol5 = Optim.optimize(
	t -> objective_function(t, obj_inputs),
	[7.0],
	ParticleSwarm(),
	Optim.Options(
		store_trace = true,
		extended_trace = true,
		iterations = 10^6
	)
)

Optim.minimizer(optim_sol5)
Optim.minimum(optim_sol5)

#%%
optim_sol6 = Optim.optimize(
	t -> objective_function(t, obj_inputs),
	[7.0],
	SimulatedAnnealing(),
	Optim.Options(
		store_trace = true,
		extended_trace = true,
		iterations = 10^6
	)
)

Optim.minimizer(optim_sol6)
Optim.minimum(optim_sol6)


#%%
sort!(output_df, :threshold)
scatter(output_df.threshold, output_df.loss; markersize = 20, alpha = 0.2)
scatter!(Optim.x_trace(optim_sol), Optim.f_trace(optim_sol), color = :red)
scatter!(Optim.x_trace(optim_sol2), Optim.f_trace(optim_sol2), color = :purple)
scatter!(reduce(vcat, Optim.x_trace(optim_sol3)), reduce(vcat, Optim.f_trace(optim_sol3)), color = :green)

#%%
lines(eachindex(Optim.x_trace(optim_sol)), Optim.x_trace(optim_sol))
lines!(eachindex(Optim.x_trace(optim_sol2)), Optim.x_trace(optim_sol2))
lines!(eachindex(reduce(vcat, Optim.x_trace(optim_sol3))), reduce(vcat, Optim.x_trace(optim_sol3)), color = :green)

#%%
lines(eachindex(Optim.f_trace(optim_sol)), Optim.f_trace(optim_sol))
lines!(eachindex(Optim.f_trace(optim_sol2)), Optim.f_trace(optim_sol2))
lines!(eachindex(reduce(vcat, Optim.f_trace(optim_sol3))), reduce(vcat, Optim.f_trace(optim_sol3)), color = :green)

