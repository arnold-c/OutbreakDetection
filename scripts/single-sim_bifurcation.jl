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

    params = (β₀, β₁, σ, γ, μ_run, ε, R₀)
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
β₁_min = 0.0
β₁_max = 1.0
β₁_step = 0.02
n_β₁s = length(β₁_min:β₁_step:β₁_max)
β₁_vec = zeros(Float64, n_β₁s)
β₁_vec .= collect(β₁_min:β₁_step:β₁_max)

bifurc_β₁_seir_arr = zeros(Float64, size(u₀, 1), tlength, n_β₁s);
bifurc_β₁_change_arr = zeros(Float64, size(u₀, 1), tlength, n_β₁s);
bifurc_β₁_jump_arr = zeros(Float64, 9, tlength, n_β₁s);

μ = 50 / (1000 * 365)
ε = 1.06 * μ * (R₀ - 1) / sqrt(N)

@showprogress for (k, β₁) in pairs(β₁_vec)
    seir = @view bifurc_β₁_seir_arr[:, :, k]
    change = @view bifurc_β₁_change_arr[:, :, k]
    jump = @view bifurc_β₁_jump_arr[:, :, k]

    local p = (β₀, β₁, σ, γ, μ, ε, R₀)

    seir_mod!(seir, change, jump, u₀, p, trange; dt = τ, type = "det")
end

years = (40 * 365):365:(tlength - 365)
bifurc_β₁_annual_summary = zeros(Float64, 5, length(years), n_β₁s)

@floop for (k, β₁) in pairs(β₁_vec), (year, day) in pairs(years),
    state in eachindex(state_labels)

    bifurc_β₁_annual_summary[state, year, k] = maximum(
        bifurc_β₁_seir_arr[state, day:(day + 364), k]
    )
end

bifurc_β₁_seir_arr[2, (40 * 365):(40 * 365 + 364), 1] ==
bifurc_β₁_seir_arr[2, (40 * 365):(40 * 365 + 364), 11]

#%%
bifurc_β₁_fig = Figure()
bifurc_β₁_ax = Axis(
    bifurc_β₁_fig[1, 1]; xlabel = "β₁ (seasonality)", ylabel = "Max. I"
)

for year in eachindex(years)
    scatter!(
        bifurc_β₁_ax,
        β₁_min:β₁_step:β₁_max,
        bifurc_β₁_annual_summary[2, year, :];
        markersize = 4,
        color = :black,
    )
end

bifurc_β₁_fig

#%%
bifurc_μ_β₁_seir_arr = zeros(Float64, size(u₀, 1), tlength, n_μs, n_β₁s);
bifurc_μ_β₁_change_arr = zeros(Float64, size(u₀, 1), tlength, n_μs, n_β₁s);
bifurc_μ_β₁_jump_arr = zeros(Float64, 9, tlength, n_μs, n_β₁s);

prog = Progress(length(μ_vec) * length(β₁_vec))
@floop for (β₁_pair, μ_pair) in
           IterTools.product(pairs(β₁_vec), pairs(μ_vec))
    k = μ_pair[1]
    μ = μ_pair[2]
    l = β₁_pair[1]
    β₁ = β₁_pair[2]
    ε = ε_vec[k]

    seir = @view bifurc_μ_β₁_seir_arr[:, :, k, l]
    change = @view bifurc_μ_β₁_change_arr[:, :, k, l]
    jump = @view bifurc_μ_β₁_jump_arr[:, :, k, l]

    params = (β₀, β₁, σ, γ, μ, ε, R₀)

    seir_mod!(seir, change, jump, u₀, params, trange; dt = τ, type = "det")
    next!(prog)
end

#%%
years = (40 * 365):365:(tlength - 365)

bifurc_μ_β₁_annual_summary = zeros(
    Float64, size(u₀, 1), length(years), n_μs, n_β₁s
);
prog = Progress(length(years) * 5 * length(β₁_vec))
@floop for (β₁, state, years) in
           IterTools.product(eachindex(β₁_vec), eachindex(u₀), pairs(years))
    year = years[1]
    day = years[2]

    bifurc_μ_β₁_annual_summary[state, year, :, β₁] = maximum(
        bifurc_μ_β₁_seir_arr[state, day:(day + 364), :, β₁]; dims = 1
    )
    next!(prog)
end

bifurc_μ_β₁_cycle_summary = zeros(Float64, n_β₁s, n_μs, 5);
prog = Progress(5 * length(β₁_vec) * length(μ_vec))
@floop for (β₁, μ, state) in
           IterTools.product(eachindex(β₁_vec), eachindex(μ_vec), 1:5)
    bifurc_μ_β₁_cycle_summary[β₁, μ, state] = length(
        Set(round.(bifurc_μ_β₁_annual_summary[state, :, μ, β₁]))
    )
    next!(prog)
end

#%%
# Because Julia is column-major, we need to transpose the heatmap as it reads the data one column at a time, and in our original matrix, each column represents a different β₁ value.
bifurc_μ_β₁_fig, bifurc_μ_β₁_ax, bifurc_μ_β₁_hm = heatmap(
    μ_vec .* (1000 * 365), β₁_vec, bifurc_μ_β₁_cycle_summary[:, :, 2]'
)
Colorbar(bifurc_μ_β₁_fig[:, end + 1], bifurc_μ_β₁_hm)

bifurc_μ_β₁_ax.xlabel = "μ (per 1000, per annum)"
bifurc_μ_β₁_ax.ylabel = "β₁ (seasonality)"

bifurc_μ_β₁_fig

