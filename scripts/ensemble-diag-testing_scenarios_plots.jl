#%%
using DrWatson
@quickactivate "OutbreakDetection"

using GLMakie

include("ensemble-diag-testing_scenarios.jl")

include(srcdir("makie-plotting-setup.jl"))

#%%
sensitivity_vec = collect(0.8:0.1:1.0)
specificity_vec = collect(0.8:0.1:1.0)
detectthreshold_vec = collect(5:5:20)

#%%
ensemble_chars_vec = Vector(
    undef, length(sensitivity_vec) * length(detectthreshold_vec)
)

ensemble_static_noise_arr = create_static_NoiseSpecification(
    [10.0],
    SimTimeParameters(;
        tmin = 0.0, tmax = 365.0 * 100, tstep = 1.0
    ),
    0.0,
    0.1,
    1_000,
)

#%%
prog = Progress(length(sensitivity_vec) * length(detectthreshold_vec))
for (i, ((sens, spec), detectthrehold)) in enumerate(
    Iterators.product(
        zip(sensitivity_vec, specificity_vec),
        detectthreshold_vec
    ),
)
    ensemble_spec =
        let time_params = SimTimeParameters(;
                tmin = 0.0, tmax = 365.0 * 100, tstep = 1.0
            )
            ScenarioSpecification(
                EnsembleSpecification(
                    ("seasonal-infectivity-import", "tau-leaping"),
                    StateParameters(
                        500_000,
                        Dict(
                            :s_prop => 0.1,
                            :e_prop => 0.01,
                            :i_prop => 0.01,
                            :r_prop => 0.88,
                        ),
                    ),
                    DynamicsParameters(500_000, 5, 0.2),
                    time_params,
                    1_000,
                ),
                OutbreakSpecification(5, 30, 500),
                ensemble_static_noise_arr,
                OutbreakDetectionSpecification(10, 7, 0.3, 0.3, 3),
                IndividualTestSpecification(sens, spec),
            )
        end

    ensemble_chars_file = get_scenario_file("scenario", ensemble_spec)

    ensemble_chars_vec[i] = ensemble_chars_file["OT_chars"]

    next!(prog)
end

#%%
OT_sens_spec_fig = Figure()
for (i, ((x, (sens, spec)), (y, detectthrehold))) in enumerate(
    Iterators.product(
        enumerate(zip(sensitivity_vec, specificity_vec)),
        enumerate(detectthreshold_vec),
    ),
)
    gl = OT_sens_spec_fig[x, y] = GridLayout()
    ax = Axis(
        gl[1, 1];
        xlabel = "Characteristic Value",
        ylabel = "Density",
    )

    hist!(
        ax,
        ensemble_chars_vec[i].sensitivity;
        bins = 0.0:0.01:1.01,
        color = (:blue, 0.5),
        normalization = :pdf,
    )

    hist!(
        ax,
        ensemble_chars_vec[i].specificity;
        bins = 0.0:0.01:1.01,
        color = (:red, 0.5),
        normalization = :pdf,
    )

    Label(
        gl[2, :],
        "Individual Sensitivity/Specificity: $(sens), Detection Threshold: $(detectthrehold)";
        word_wrap = true
    )
    colsize!(gl, 1, Relative(1))
end

OT_sens_spec_fig
