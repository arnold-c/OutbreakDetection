#%%
using DrWatson
@quickactivate "OutbreakDetection"

using ProgressMeter
using DifferentialEquations

include("../src/OutbreakDetection.jl")
using .OutbreakDetection

#%%
init_noise = [10.0]
sde_cb = DiscreteCallback(
    sde_condition, sde_affect!; save_positions = (false, false)
)

ensemble_single_scenario_noise_spec = create_static_NoiseSpecification(
    init_noise,
    SimTimeParameters(; tmin = 0.0, tmax = 365.0 * 100, tstep = 1.0),
    0.0,
    0.1,
    1_000;
    callback = sde_cb
)
