#%%
using DrWatson
@quickactivate "OutbreakDetection"

using UnPack
using ColorSchemes

using OutbreakDetection

include(srcdir("makie-plotting-setup.jl"))
includet(srcdir("ensemble-parameters.jl"))

#%%
ensemble_single_individual_test_spec = IndividualTestSpecification(0.8, 0.8, 0)

#%%
for noise_specification in ensemble_noise_specification_vec
    @show noise_specification
    scenario_specification = ScenarioSpecification(
        ensemble_specification,
        ensemble_outbreak_specification,
        noise_specification,
        OutbreakDetectionSpecification(5, 7, 0.6, 0.3),
        ensemble_single_individual_test_spec,
    )

    ensemble_solution_dict = get_ensemble_file(
        scenario_specification.ensemble_specification
    )

    incarr = get_ensemble_file(
        ensemble_specification, ensemble_outbreak_specification
    )

    ensemble_single_scenario_detection = get_ensemble_file(
        scenario_specification
    )

    noisearr = create_noise_arr(
        noise_specification,
        incarr["ensemble_inc_arr"];
        ensemble_specification = ensemble_specification,
        seed = 1234,
    )

    fig = Figure()
    ax = Axis(fig[2, 1])

    for noise_sim in eachcol(noisearr)
        lines!(
            ax,
            ensemble_time_specification.trange,
            noise_sim;
            color = (:gray, 0.2),
        )
    end

    meanline = vec(mean(noisearr; dims = 2))
    dailymean = mean(meanline)
    noisedir = getdirpath(noise_specification)
    noisepath = replace(noisedir, "/" => "_")
    plotpath = joinpath(
        plotsdir(),
        "ensemble",
        "single-scenario",
        "noise-visualizations"
    )
    mkpath(plotpath)

    lines!(
        ax,
        ensemble_time_specification.trange,
        meanline;
        color = :black,
    )

    Label(
        fig[1, :],
        "Noise: $(noisedir), Daily Mean: $(round(dailymean, digits = 2))",
    )

    rowsize!(fig.layout, 1, 5)
    colsize!(fig.layout, 1, Relative(1))

    save(
        joinpath(
            plotpath,
            "$(noisepath).png",
        ),
        fig; resolution = (2200, 1600),
    )
end
