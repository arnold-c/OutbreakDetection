#%%
using DrWatson
@quickactivate "OutbreakDetection"

using ProgressMeter
using IterTools

includet(scriptsdir("single-sim.jl"))

#%%
μ_min = 10
μ_max = 100
μ_step = 2.0
n_μs = length(μ_min:μ_step:μ_max)
μ_vec = zeros(Float64, n_μs)
μ_vec .= collect(μ_min:μ_step:μ_max) ./ (1000 * 365)

ε_vec = (1.06 .* μ_vec .* (R₀ - 1)) ./ sqrt(N)

bifurc_μ_seir_arr = zeros(Float64, size(u₀, 1), tlength, n_μs);
bifurc_μ_change_arr = zeros(Float64, size(u₀, 1), tlength, n_μs);
bifurc_μ_jump_arr = zeros(Float64, 9, tlength, n_μs);

prog = Progress(n_μs)
@floop for (k, μ_run) in pairs(μ_vec)
    ε = ε_vec[k]

    seir = @view bifurc_μ_seir_arr[:, :, k]
    change = @view bifurc_μ_change_arr[:, :, k]
    jump = @view bifurc_μ_jump_arr[:, :, k]

    params = (beta_mean, beta_force, sigma, γ, μ_run, ε, R₀)
    seir_mod!(seir, change, jump,
        u₀, params, trange; dt = τ, type = "det",
    )
    next!(prog)
end

years = (40 * 365):365:(tlength - 365)
bifurc_μ_annual_summary = zeros(Float64, 5, length(years), n_μs);

for μ in eachindex(μ_vec), state in eachindex(state_labels),
    (year, day) in pairs(years)

    bifurc_μ_annual_summary[state, year, μ] = maximum(
        bifurc_μ_seir_arr[state, day:(day + 364), μ]
    )
end

bifurc_μ_seir_arr[2, (40 * 365):(40 * 365 + 364), 1] ==
bifurc_μ_seir_arr[2, (40 * 365):(40 * 365 + 364), 10]

#%%
bifurc_μ_fig = Figure()
bifurc_μ_ax = Axis(
    bifurc_μ_fig[1, 1]; xlabel = "μ (per 1000, per annum)", ylabel = "Max. I"
)

for year in eachindex(years)
    scatter!(
        bifurc_μ_ax,
        μ_min:μ_step:μ_max,
        bifurc_μ_annual_summary[2, year, :];
        markersize = 4,
        color = :black,
    )
end

bifurc_μ_fig

#%%
beta_force_min = 0.0
beta_force_max = 1.0
beta_force_step = 0.02
n_beta_forces = length(beta_force_min:beta_force_step:beta_force_max)
beta_force_vec = zeros(Float64, n_beta_forces)
beta_force_vec .= collect(beta_force_min:beta_force_step:beta_force_max)

bifurc_beta_force_seir_arr = zeros(Float64, size(u₀, 1), tlength, n_beta_forces);
bifurc_beta_force_change_arr = zeros(Float64, size(u₀, 1), tlength, n_beta_forces);
bifurc_beta_force_jump_arr = zeros(Float64, 9, tlength, n_beta_forces);

μ = 50 / (1000 * 365)
ε = 1.06 * μ * (R₀ - 1) / sqrt(N)

@showprogress for (k, beta_force) in pairs(beta_force_vec)
    seir = @view bifurc_beta_force_seir_arr[:, :, k]
    change = @view bifurc_beta_force_change_arr[:, :, k]
    jump = @view bifurc_beta_force_jump_arr[:, :, k]

    local p = (beta_mean, beta_force, sigma, γ, μ, ε, R₀)

    seir_mod!(seir, change, jump, u₀, p, trange; dt = τ, type = "det")
end

years = (40 * 365):365:(tlength - 365)
bifurc_beta_force_annual_summary = zeros(Float64, 5, length(years), n_beta_forces)

@floop for (k, beta_force) in pairs(beta_force_vec), (year, day) in pairs(years),
    state in eachindex(state_labels)

    bifurc_beta_force_annual_summary[state, year, k] = maximum(
        bifurc_beta_force_seir_arr[state, day:(day + 364), k]
    )
end

bifurc_beta_force_seir_arr[2, (40 * 365):(40 * 365 + 364), 1] ==
bifurc_beta_force_seir_arr[2, (40 * 365):(40 * 365 + 364), 11]

#%%
bifurc_beta_force_fig = Figure()
bifurc_beta_force_ax = Axis(
    bifurc_beta_force_fig[1, 1]; xlabel = "beta_force (seasonality)", ylabel = "Max. I"
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
bifurc_μ_beta_force_seir_arr = zeros(Float64, size(u₀, 1), tlength, n_μs, n_beta_forces);
bifurc_μ_beta_force_change_arr = zeros(Float64, size(u₀, 1), tlength, n_μs, n_beta_forces);
bifurc_μ_beta_force_jump_arr = zeros(Float64, 9, tlength, n_μs, n_beta_forces);

prog = Progress(length(μ_vec) * length(beta_force_vec))
@floop for (beta_force_pair, μ_pair) in
           IterTools.product(pairs(beta_force_vec), pairs(μ_vec))
    k = μ_pair[1]
    μ = μ_pair[2]
    l = beta_force_pair[1]
    beta_force = beta_force_pair[2]
    ε = ε_vec[k]

    seir = @view bifurc_μ_beta_force_seir_arr[:, :, k, l]
    change = @view bifurc_μ_beta_force_change_arr[:, :, k, l]
    jump = @view bifurc_μ_beta_force_jump_arr[:, :, k, l]

    params = (beta_mean, beta_force, sigma, γ, μ, ε, R₀)

    seir_mod!(seir, change, jump, u₀, params, trange; dt = τ, type = "det")
    next!(prog)
end

#%%
years = (40 * 365):365:(tlength - 365)

bifurc_μ_beta_force_annual_summary = zeros(
    Float64, size(u₀, 1), length(years), n_μs, n_beta_forces
);
prog = Progress(length(years) * 5 * length(beta_force_vec))
@floop for (beta_force, state, years) in
           IterTools.product(eachindex(beta_force_vec), eachindex(u₀), pairs(years))
    year = years[1]
    day = years[2]

    bifurc_μ_beta_force_annual_summary[state, year, :, beta_force] = maximum(
        bifurc_μ_beta_force_seir_arr[state, day:(day + 364), :, beta_force]; dims = 1
    )
    next!(prog)
end

bifurc_μ_beta_force_cycle_summary = zeros(Float64, n_beta_forces, n_μs, 5);
prog = Progress(5 * length(beta_force_vec) * length(μ_vec))
@floop for (beta_force, μ, state) in
           IterTools.product(eachindex(beta_force_vec), eachindex(μ_vec), 1:5)
    bifurc_μ_beta_force_cycle_summary[beta_force, μ, state] = length(
        Set(round.(bifurc_μ_beta_force_annual_summary[state, :, μ, beta_force]))
    )
    next!(prog)
end

#%%
# Because Julia is column-major, we need to transpose the heatmap as it reads the data one column at a time, and in our original matrix, each column represents a different beta_force value.
bifurc_μ_beta_force_fig, bifurc_μ_beta_force_ax, bifurc_μ_beta_force_hm = heatmap(
    μ_vec .* (1000 * 365), beta_force_vec, bifurc_μ_beta_force_cycle_summary[:, :, 2]'
)
Colorbar(bifurc_μ_beta_force_fig[:, end + 1], bifurc_μ_beta_force_hm)

bifurc_μ_beta_force_ax.xlabel = "μ (per 1000, per annum)"
bifurc_μ_beta_force_ax.ylabel = "beta_force (seasonality)"

bifurc_μ_beta_force_fig

