#%%
using DrWatson
@quickactivate "OutbreakDetection"

using ProgressMeter
using FLoops

using OutbreakDetection

include(srcdir("makie-plotting-setup.jl"))

#%%
@unpack singlesim_states_p, singlesim_time_p, singlesim_dynamics_p = load(
    "data/singlesim/single-sim_setup.jld2"
)
@unpack init_states = singlesim_states_p
@unpack tlength = singlesim_time_p

#%%
beta_arr = Vector{Float64}(undef, tlength);
rates = Vector{Float64}(undef, 9);
years = (40 * 365):365:(tlength - 365)
n_transitions = (size(init_states, 1) - 1) * 2 + 1

#%%
annual_birth_rate_per_k_min = 10
annual_birth_rate_per_k_max = 100
annual_birth_rate_per_k_step = 2.0
annual_birth_rate_per_k_vec =
    annual_birth_rate_per_k_min:annual_birth_rate_per_k_step:annual_birth_rate_per_k_max

n_annual_birth_rate_per_k = length(annual_birth_rate_per_k_vec)

#%%
bifurc_mu_seir_arr = zeros(
    Float64, tlength, size(init_states, 1), n_annual_birth_rate_per_k
);
bifurc_mu_change_arr = similar(bifurc_mu_seir_arr);
bifurc_mu_jump_arr = zeros(
    Float64,
    tlength,
    n_transitions,
    n_annual_birth_rate_per_k
);

#%%
birth_rate_bifurcation_simulation!(
    bifurc_mu_seir_arr,
    bifurc_mu_change_arr,
    bifurc_mu_jump_arr,
    beta_arr,
    init_states,
    rates,
    annual_birth_rate_per_k_vec,
    singlesim_dynamics_p,
    singlesim_time_p;
)

#%%
bifurc_mu_annual_summary = bifurcation_summary(
    bifurc_mu_seir_arr,
    annual_birth_rate_per_k_vec,
    years;
    state_labels = seir_state_labels,
)

#%%
birth_rate_bifurcation_plot = bifurcation_plot(
    annual_birth_rate_per_k_vec,
    bifurc_mu_annual_summary;
    years = years,
    xlabel = "Birth rate (per 1_000, per annum)",
    ylabel = "Max. I",
)

save(
    plotsdir("singlesim/bifurcation_birth-rate.png"),
    birth_rate_bifurcation_plot,
)

#%%
beta_force_min = 0.0
beta_force_max = 1.0
beta_force_step = 0.02
n_beta_forces = length(beta_force_min:beta_force_step:beta_force_max)
beta_force_vec = zeros(Float64, n_beta_forces)
beta_force_vec .= collect(beta_force_min:beta_force_step:beta_force_max)

bifurc_beta_force_seir_arr = zeros(
    Float64, tlength, size(init_states, 1), n_beta_forces
);
bifurc_beta_force_change_arr = similar(bifurc_beta_force_seir_arr);
bifurc_beta_force_jump_arr = zeros(
    Float64, tlength, n_transitions, n_beta_forces
);

beta_force_annual_birth_rate_per_k = 50

beta_force_bifurcation_simulation!(
    bifurc_beta_force_seir_arr,
    bifurc_beta_force_change_arr,
    bifurc_beta_force_jump_arr,
    beta_arr,
    init_states,
    rates,
    beta_force_vec,
    singlesim_dynamics_p,
    singlesim_time_p;
    birth_rate = beta_force_annual_birth_rate_per_k,
)

bifurc_beta_force_annual_summary = bifurcation_summary(
    bifurc_beta_force_seir_arr,
    beta_force_vec,
    years;
    state_labels = seir_state_labels,
)

#%%
beta_force_bifurcation_plot = bifurcation_plot(
    beta_force_vec,
    bifurc_beta_force_annual_summary;
    years = years,
    xlabel = "beta_force (seasonality)",
    ylabel = "Max. I",
)

save(
    plotsdir("singlesim/bifurcation_beta-force.png"),
    beta_force_bifurcation_plot,
)

#%%
bifurc_mu_beta_force_seir_arr = zeros(
    Float64, tlength, size(init_states, 1), n_annual_birth_rate_per_k,
    n_beta_forces,
);
bifurc_mu_beta_force_change_arr = similar(bifurc_mu_beta_force_seir_arr);
bifurc_mu_beta_force_jump_arr = zeros(
    Float64, tlength, n_transitions, n_annual_birth_rate_per_k, n_beta_forces
);

#%%
birth_rate_beta_force_bifurcation_simulation!(
    bifurc_mu_beta_force_seir_arr,
    bifurc_mu_beta_force_change_arr,
    bifurc_mu_beta_force_jump_arr,
    beta_arr,
    init_states,
    rates,
    annual_birth_rate_per_k_vec,
    beta_force_vec,
    singlesim_dynamics_p,
    singlesim_time_p,
)

bifurc_mu_beta_force_seir_arr[:, :, 2, 1]
sum(bifurc_mu_beta_force_seir_arr[:, :, 2, 1][end, 1:4])

#%%
bifurc_mu_beta_force_annual_summary = birth_rate_beta_force_bifurcation_annual_summary(
    bifurc_mu_beta_force_seir_arr,
    annual_birth_rate_per_k_vec,
    beta_force_vec,
    years,
)

bifurc_mu_beta_force_cycle_summary = birth_rate_beta_force_bifurcation_cycle_summary(
    bifurc_mu_beta_force_annual_summary,
    annual_birth_rate_per_k_vec,
    beta_force_vec,
)

#%%
test_seir_array, test_change_array, test_jump_array, test_beta_arr = seir_mod(
    init_states,
    DynamicsParameters(
        500_000, Integer(annual_birth_rate_per_k_vec[1]), beta_force_vec[1]
    ), singlesim_time_p;
    type = "det", seed = 1234,
);

test_seir_array == test_seir_array1
test_change_array == test_change_array1
test_jump_array == test_jump_array1

lines(test_beta_arr[1:(365 * 2)])

#%%
comb = Iterators.product(
    pairs(annual_birth_rate_per_k_vec[1]), pairs(beta_force_vec[1])
)

function test_fun(
    beta_force_seir_arr,
    beta_force_change_arr,
    beta_force_jump_arr,
    beta_arr,
    initial_states,
    rates,
    birth_rate_vec,
    beta_force_vec,
    dynamics_parameters,
    time_parameters,
)
    for birth_rate_pair in pairs(birth_rate_vec),
        beta_force_pair in pairs(beta_force_vec)
        k = birth_rate_pair[1]
        birth_rate_run = birth_rate_pair[2]
        mu_run = calculate_mu(birth_rate_run)

        l = beta_force_pair[1]
        beta_force_run = beta_force_pair[2]
        epsilon_run = calculate_import_rate(
            mu_run, dynamics_parameters.R_0, initial_states.N
        )

        run_dynamics_parameters = DynamicsParameters(
            dynamics_parameters.beta_mean,
            beta_force_run,
            dynamics_parameters.sigma,
            dynamics_parameters.gamma,
            mu_run,
            birth_rate_run,
            epsilon_run,
            dynamics_parameters.R_0,
        )

        seir_mod!(
            beta_force_seir_arr,
            beta_force_change_arr,
            beta_force_jump_arr,
            beta_arr,
            initial_states,
            rates,
            run_dynamics_parameters,
            time_parameters;
            type = "det",
        )
    end
    # for (birth_rate_pair, beta_force_pair) in comb
    #     k = birth_rate_pair[1]
    #     birth_rate_run = birth_rate_pair[2]
    #     mu_run = calculate_mu(birth_rate_run)
    #
    #     l = beta_force_pair[1]
    #     beta_force_run = beta_force_pair[2]
    #     epsilon_run = calculate_import_rate(
    #         mu_run, singlesim_dynamics_p.R_0, init_states.N
    #     )
    #
    #     run_dynamics_parameters = DynamicsParameters(
    #         singlesim_dynamics_p.beta_mean,
    #         beta_force_run,
    #         singlesim_dynamics_p.sigma,
    #         singlesim_dynamics_p.gamma,
    #         mu_run,
    #         birth_rate_run,
    #         epsilon_run,
    #         singlesim_dynamics_p.R_0,
    #     )
    #
    #     seir_mod!(
    #         test_seir_array,
    #         test_change_array,
    #         test_jump_array,
    #         test_beta_arr,
    #         init_states,
    #         Vector{Float64}(undef, 9),
    #         run_dynamics_parameters,
    #         singlesim_time_p;
    #         type = "det",
    #     )
    # end
    return nothing
end

test_seir_array .= 0
test_change_array .= 0
test_jump_array .= 0

test_fun(
    test_seir_array,
    test_change_array,
    test_jump_array,
    test_beta_arr,
    init_states,
    rates,
    annual_birth_rate_per_k_vec[1:2],
    beta_force_vec[1:2],
    singlesim_dynamics_p,
    singlesim_time_p,
)

test_seir_array

test_seir_df = create_sir_df(
    test_seir_array, singlesim_time_p.trange, [:S, :E, :I, :R, :N]
)
draw_sir_plot(test_seir_df;
    colors = ["dodgerblue4", "green", "firebrick3", "chocolate2", "purple"],
    labels = ["S", "E", "I", "R", "N"])

bifurc_wide = load("tmp/bifurc_mu_beta_force_arrays.jld2")

bifurc_mu_beta_force_seir_arr[:, :, 1, 1]

sum(bifurc_mu_beta_force_seir_arr[:, :, 1, 1][end, 1:4])

bifurc_mu_beta_force_change_arr[:, :, 1, 1]
sum(bifurc_mu_beta_force_change_arr[:, :, 1, 1][end, 1:4])
beta_arr

bifurc_wide["bifurc_mu_beta_force_seir_arr"][:, :, 1, 1]'
sum(bifurc_wide["bifurc_mu_beta_force_seir_arr"][:, :, 1, 1]'[end, 1:4])

bifurc_wide["bifurc_mu_beta_force_change_arr"][:, :, 1, 1]'

bifurc_mu_beta_force_annual_summary[:, :, 1, 1]

bifurc_wide["bifurc_mu_beta_force_annual_summary"][:, :, 1, 1]

bifurc_mu_beta_force_change_arr[:, :, 1, 1]
bifurc_wide["bifurc_mu_beta_force_change_arr"][:, :, 1, 1]'

#%%
birth_rate_beta_force_bifurcation_heatmap = bifurcation_heatmap(
    annual_birth_rate_per_k_vec,
    beta_force_vec,
    bifurc_mu_beta_force_cycle_summary[:, :, 2],
)

save(
    plotsdir("singlesim/bifurcation_birth-rate_beta-force.png"),
    birth_rate_beta_force_bifurcation_heatmap,
)
