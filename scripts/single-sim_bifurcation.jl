#%%
using DrWatson
@quickactivate "OutbreakDetection"

using ProgressMeter
using IterTools

includet(scriptsdir("single-sim.jl"))

#%%
mu_min = 10
mu_max = 100
mu_step = 2.0
n_mus = length(mu_min:mu_step:mu_max)
mu_vec = zeros(Float64, n_mus)
mu_vec .= collect(mu_min:mu_step:mu_max) ./ (1_000 * 365)

@unpack init_states = singlesim_init_states
@unpack tlength = singlesim_time_p

epsilon_vec = zeros(Float64, n_mus)
for (i, mu) in pairs(mu_vec)
    epsilon_vec[i] = calculate_import_rate(mu, singlesim_dynamics_p.R_0, init_states.N)
end

bifurc_mu_seir_arr = zeros(Float64, size(init_states, 1), tlength, n_mus);
bifurc_mu_change_arr = zeros(Float64, size(init_states, 1), tlength, n_mus);
bifurc_mu_jump_arr = zeros(Float64, 9, tlength, n_mus);

prog = Progress(n_mus)
@floop for (k, mu_run) in pairs(mu_vec)
    epsilon_run = epsilon_vec[k]

    seir = @view bifurc_mu_seir_arr[:, :, k]
    change = @view bifurc_mu_change_arr[:, :, k]
    jump = @view bifurc_mu_jump_arr[:, :, k]

    bifurc_mu_dynamics_p = DynamicsParameters(;
        beta_mean = singlesim_dynamics_p.beta_mean,
        beta_force = singlesim_dynamics_p.beta_force,
        sigma = singlesim_dynamics_p.sigma,
        gamma = singlesim_dynamics_p.gamma,
        mu = mu_run,
        epsilon = epsilon_run,
        R_0 = singlesim_dynamics_p.R_0,
    )

    seir_mod!(seir, change, jump,
        init_states, bifurc_mu_dynamics_p, singlesim_time_p; type = "det",
    )
    next!(prog)
end

years = (40 * 365):365:(tlength - 365)
bifurc_mu_annual_summary = zeros(Float64, 5, length(years), n_mus);

for mu in eachindex(mu_vec), state in eachindex(seir_state_labels),
    (year, day) in pairs(years)

    bifurc_mu_annual_summary[state, year, mu] = maximum(
        bifurc_mu_seir_arr[state, day:(day + 364), mu]
    )
end

bifurc_mu_seir_arr[2, (40 * 365):(40 * 365 + 364), 1] ==
bifurc_mu_seir_arr[2, (40 * 365):(40 * 365 + 364), 10]

#%%
bifurc_mu_fig = Figure()
bifurc_mu_ax = Axis(
    bifurc_mu_fig[1, 1]; xlabel = "mu (per 1_000, per annum)", ylabel = "Max. I"
)

for year in eachindex(years)
    scatter!(
        bifurc_mu_ax,
        mu_min:mu_step:mu_max,
        bifurc_mu_annual_summary[2, year, :];
        markersize = 4,
        color = :black,
    )
end

bifurc_mu_fig

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

beta_force_mu = 50 / (1_000 * 365)
beta_force_epsilon = calculate_import_rate(mu, R_0, init_states.N)

@showprogress for (k, beta_force_run) in pairs(beta_force_vec)
    seir = @view bifurc_beta_force_seir_arr[:, :, k]
    change = @view bifurc_beta_force_change_arr[:, :, k]
    jump = @view bifurc_beta_force_jump_arr[:, :, k]

    bifurc_beta_force_dynamics_p = DynamicsParameters(;
        beta_mean = singlesim_dynamics_p.beta_mean,
        beta_force = beta_force_run,
        sigma = singlesim_dynamics_p.sigma,
        gamma = singlesim_dynamics_p.gamma,
        mu = beta_force_mu,
        epsilon = beta_force_epsilon,
        R_0 = singlesim_dynamics_p.R_0,
    )

    seir_mod!(
        seir, change, jump, init_states, bifurc_beta_force_dynamics_p, singlesim_time_p; type = "det"
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
    Float64, size(init_states, 1), tlength, n_mus, n_beta_forces
);
bifurc_mu_beta_force_change_arr = zeros(
    Float64, size(init_states, 1), tlength, n_mus, n_beta_forces
);
bifurc_mu_beta_force_jump_arr = zeros(Float64, 9, tlength, n_mus, n_beta_forces);

prog = Progress(length(mu_vec) * length(beta_force_vec))
@floop for (beta_force_pair, mu_pair) in
           IterTools.product(pairs(beta_force_vec), pairs(mu_vec))
    k = mu_pair[1]
    mu_run = mu_pair[2]
    l = beta_force_pair[1]
    beta_force_run = beta_force_pair[2]
    epsilon_run = epsilon_vec[k]

    seir = @view bifurc_mu_beta_force_seir_arr[:, :, k, l]
    change = @view bifurc_mu_beta_force_change_arr[:, :, k, l]
    jump = @view bifurc_mu_beta_force_jump_arr[:, :, k, l]

    bifurc_mu_beta_force_dynamics_p = DynamicsParameters(;
        beta_mean = singlesim_dynamics_p.beta_mean,
        beta_force = beta_force_run,
        sigma = singlesim_dynamics_p.sigma,
        gamma = singlesim_dynamics_p.gamma,
        mu = mu_run,
        epsilon = epsilon_run,
        R_0 = singlesim_dynamics_p.R_0,
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
    Float64, size(init_states, 1), length(years), n_mus, n_beta_forces
);
prog = Progress(length(years) * 5 * length(beta_force_vec))
@floop for (beta_force, state, years) in
           IterTools.product(
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

bifurc_mu_beta_force_cycle_summary = zeros(Float64, n_beta_forces, n_mus, 5);
prog = Progress(5 * length(beta_force_vec) * length(mu_vec))
@floop for (beta_force, mu, state) in
           IterTools.product(eachindex(beta_force_vec), eachindex(mu_vec), 1:5)
    bifurc_mu_beta_force_cycle_summary[beta_force, mu, state] = length(
        Set(
            round.(
                bifurc_mu_beta_force_annual_summary[state, :, mu, beta_force]
            ),
        ),
    )
    next!(prog)
end

#%%
# Because Julia is column-major, we need to transpose the heatmap as it reads the data one column at a time, and in our original matrix, each column represents a different beta_force value.
bifurc_mu_beta_force_fig, bifurc_mu_beta_force_ax, bifurc_mu_beta_force_hm = heatmap(
    mu_vec .* (1_000 * 365), beta_force_vec,
    bifurc_mu_beta_force_cycle_summary[:, :, 2]',
)
Colorbar(bifurc_mu_beta_force_fig[:, end + 1], bifurc_mu_beta_force_hm)

bifurc_mu_beta_force_ax.xlabel = "mu (per 1_000, per annum)"
bifurc_mu_beta_force_ax.ylabel = "beta_force (seasonality)"

bifurc_mu_beta_force_fig
