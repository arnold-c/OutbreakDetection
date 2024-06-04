using GLMakie

function visualize_ensemble_noise(
    ensemble_noise_arr, poisson_noise_prop_vec, timespecification, noisedir;
    xlabel="Time (years)", ylabel="Noise Incidence",
)
    times = collect(timespecification.trange) ./ 365
    meanline = vec(mean(ensemble_noise_arr; dims=2))
    dailymean = NaNMath.mean(meanline)
    poisson_noise_prop = NaNMath.mean(poisson_noise_prop_vec)

    poisson_noise_arr = ensemble_noise_arr .* poisson_noise_prop_vec
    poisson_noise_daily_mean = NaNMath.mean(poisson_noise_arr)
    dynamic_noise_arr = ensemble_noise_arr .- poisson_noise_arr
    dynamic_noise_daily_mean = NaNMath.mean(dynamic_noise_arr)

    fig = Figure()
    ax = Axis(fig[2, 1]; xlabel=xlabel, ylabel=ylabel)

    for noise_sim in eachcol(ensemble_noise_arr)
        lines!(
            ax,
            times,
            noise_sim;
            color=(:gray, 0.2),
        )
    end

    lines!(
        ax,
        times,
        meanline;
        color=:black,
    )

    Label(
        fig[1, :],
        "Noise: $(noisedir), Total Daily Mean: $(round(dailymean, digits = 2)), Poisson Noise Proportion: $(round(poisson_noise_prop, digits = 2))\nDynamic Daily Mean: $(round(dynamic_noise_daily_mean, digits = 2)) Poisson Daily Mean: $(round(poisson_noise_daily_mean, digits = 2))",
    )

    rowsize!(fig.layout, 1, 5)
    colsize!(fig.layout, 1, Relative(1))

    return fig
end
