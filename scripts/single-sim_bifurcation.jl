#%%
using DrWatson
@quickactivate "OutbreakDetection"

using ProgressMeter
using FLoops

include(srcdir("makie-plotting-setup.jl"))

include("single-sim.jl")

#%%
@unpack init_states = singlesim_states_p
@unpack tlength = singlesim_time_p

beta_arr = Vector{Float64}(undef, tlength);
rates = Vector{Float64}(undef, 9);
years = (40 * 365):365:(tlength - 365)

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
bifurc_mu_change_arr = zeros(
    Float64, tlength, size(init_states, 1), n_annual_birth_rate_per_k
);
bifurc_mu_jump_arr = zeros(
    Float64,
    tlength,
    (size(init_states, 1) - 1) * 2 + 1,
    n_annual_birth_rate_per_k,
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
bifurc_mu_annual_summary = birth_rate_birucation_summary(
    bifurc_mu_seir_arr,
    annual_birth_rate_per_k_vec,
    years;
    state_labels = seir_state_labels,
)

#%%
birth_rate_bifurcation_plot(
    annual_birth_rate_per_k_vec,
    bifurc_mu_annual_summary;
    years = years,
    xlabel = "Birth rate (per 1_000, per annum)",
    ylabel = "Max. I",
)

#%%
beta_force_min = 0.0
beta_force_max = 1.0
beta_force_step = 0.02
n_beta_forces = length(beta_force_min:beta_force_step:beta_force_max)
beta_force_vec = zeros(Float64, n_beta_forces)
beta_force_vec .= collect(beta_force_min:beta_force_step:beta_force_max)

bifurc_beta_force_seir_arr = zeros(
    Float64, size(init_states, 1), tlength, n_beta_forces
);
bifurc_beta_force_change_arr = zeros(
    Float64, size(init_states, 1), tlength, n_beta_forces
);
bifurc_beta_force_jump_arr = zeros(Float64, 9, tlength, n_beta_forces);

beta_force_annual_birth_rate_per_k = 50
beta_force_mu = calculate_mu(beta_force_annual_birth_rate_per_k)
beta_force_epsilon = calculate_import_rate(beta_force_mu, R_0, init_states.N)

@showprogress for (k, beta_force_run) in pairs(beta_force_vec)
    seir = @view bifurc_beta_force_seir_arr[:, :, k]
    change = @view bifurc_beta_force_change_arr[:, :, k]
    jump = @view bifurc_beta_force_jump_arr[:, :, k]

    bifurc_beta_force_dynamics_p = DynamicsParameters(
        singlesim_dynamics_p.beta_mean,
        beta_force_run,
        singlesim_dynamics_p.sigma,
        singlesim_dynamics_p.gamma,
        beta_force_mu,
        beta_force_annual_birth_rate_per_k,
        beta_force_epsilon,
        singlesim_dynamics_p.R_0,
    )

    seir_mod!(
        seir, change, jump, init_states, bifurc_beta_force_dynamics_p,
        singlesim_time_p; type = "det",
    )
end

years = (40 * 365):365:(tlength - 365)
bifurc_beta_force_annual_summary = zeros(
    Float64, 5, length(years), n_beta_forces
)

@floop for (k, beta_force) in pairs(beta_force_vec),
    (year, day) in pairs(years),
    state in eachindex(seir_state_labels)

    bifurc_beta_force_annual_summary[state, year, k] = maximum(
        bifurc_beta_force_seir_arr[state, day:(day + 364), k]
    )
end

bifurc_beta_force_seir_arr[2, (40 * 365):(40 * 365 + 364), 1] ==
bifurc_beta_force_seir_arr[2, (40 * 365):(40 * 365 + 364), 11]

#%%
bifurc_beta_force_fig = Figure()
bifurc_beta_force_ax = Axis(
    bifurc_beta_force_fig[1, 1]; xlabel = "beta_force (seasonality)",
    ylabel = "Max. I",
)

for year in eachindex(years)
    scatter!(
        bifurc_beta_force_ax,
        beta_force_min:beta_force_step:beta_force_max,
        bifurc_beta_force_annual_summary[2, year, :];
        markersize = 4,
        color = :black,
    )
end

bifurc_beta_force_fig

#%%
bifurc_mu_beta_force_seir_arr = zeros(
    Float64, size(init_states, 1), tlength, n_annual_birth_rate_per_k,
    n_beta_forces,
);
bifurc_mu_beta_force_change_arr = zeros(
    Float64, size(init_states, 1), tlength, n_annual_birth_rate_per_k,
    n_beta_forces,
);
bifurc_mu_beta_force_jump_arr = zeros(
    Float64, 9, tlength, n_annual_birth_rate_per_k, n_beta_forces
);

prog = Progress(length(annual_birth_rate_per_k_vec) * length(beta_force_vec))
@floop for (annual_birth_rate_per_k_pair, beta_force_pair) in
           Iterators.product(
    pairs(annual_birth_rate_per_k_vec), pairs(beta_force_vec)
)
    k = annual_birth_rate_per_k_pair[1]
    annual_birth_rate_per_k_run = annual_birth_rate_per_k_pair[2]
    mu_run = calculate_mu(annual_birth_rate_per_k_run)

    l = beta_force_pair[1]
    beta_force_run = beta_force_pair[2]
    epsilon_run = calculate_import_rate(
        mu_run, singlesim_dynamics_p.R_0, init_states.N
    )

    seir = @view bifurc_mu_beta_force_seir_arr[:, :, k, l]
    change = @view bifurc_mu_beta_force_change_arr[:, :, k, l]
    jump = @view bifurc_mu_beta_force_jump_arr[:, :, k, l]

    bifurc_mu_beta_force_dynamics_p = DynamicsParameters(
        singlesim_dynamics_p.beta_mean,
        beta_force_run,
        singlesim_dynamics_p.sigma,
        singlesim_dynamics_p.gamma,
        mu_run,
        annual_birth_rate_per_k_run,
        beta_force_epsilon,
        singlesim_dynamics_p.R_0,
    )

    seir_mod!(
        seir,
        change,
        jump,
        init_states,
        bifurc_mu_beta_force_dynamics_p,
        singlesim_time_p;
        type = "det",
    )

    next!(prog)
end

#%%
years = (40 * 365):365:(tlength - 365)

bifurc_mu_beta_force_annual_summary = zeros(
    Float64, size(init_states, 1), length(years), n_annual_birth_rate_per_k,
    n_beta_forces,
);
prog = Progress(length(years) * 5 * length(beta_force_vec))
@floop for (beta_force, state, years) in
           Iterators.product(
    eachindex(beta_force_vec), eachindex(init_states), pairs(years)
)
    year = years[1]
    day = years[2]

    bifurc_mu_beta_force_annual_summary[state, year, :, beta_force] = maximum(
        bifurc_mu_beta_force_seir_arr[state, day:(day + 364), :, beta_force];
        dims = 1,
    )
    next!(prog)
end

bifurc_mu_beta_force_cycle_summary = zeros(
    Float64, n_beta_forces, n_annual_birth_rate_per_k, 5
);
prog = Progress(
    5 * length(beta_force_vec) * length(annual_birth_rate_per_k_vec)
)
@floop for (beta_force, annual_birth_rate_per_k, state) in
           Iterators.product(
    eachindex(beta_force_vec), eachindex(annual_birth_rate_per_k_vec), 1:5
)
    bifurc_mu_beta_force_cycle_summary[beta_force, annual_birth_rate_per_k, state] = length(
        Set(
            round.(
                bifurc_mu_beta_force_annual_summary[
                    state, :, annual_birth_rate_per_k, beta_force
                ]
            ),
        ),
    )
    next!(prog)
end

#%%
# Because Julia is column-major, we need to transpose the heatmap as it reads the data one column at a time, and in our original matrix, each column represents a different beta_force value.
bifurc_mu_beta_force_fig, bifurc_mu_beta_force_ax, bifurc_mu_beta_force_hm = heatmap(
    annual_birth_rate_per_k_vec,
    beta_force_vec,
    bifurc_mu_beta_force_cycle_summary[:, :, 2]',
)
Colorbar(
    bifurc_mu_beta_force_fig[:, end + 1],
    bifurc_mu_beta_force_hm;
    label = "Periodicity",
)

bifurc_mu_beta_force_ax.xlabel = "Birth rate (per 1_000, per annum)"
bifurc_mu_beta_force_ax.ylabel = "beta_force (seasonality)"

bifurc_mu_beta_force_fig
