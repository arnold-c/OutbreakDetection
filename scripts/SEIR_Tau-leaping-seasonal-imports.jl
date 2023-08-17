"""
This is a simulation of an SIR model that uses Tau-leaping, with commuter
imports. All jumps are manually defined.
"""
#%%
using DrWatson
@quickactivate "OutbreakDetection"

using DataFrames, DataFramesMeta, LinearAlgebra
using WGLMakie, AlgebraOfGraphics, ColorSchemes, Colors
using BenchmarkTools, JLD2, Random, ProgressMeter, StatsBase, Distributions
using IterTools, FLoops, FreqTables, ThreadsX, ProtoStructs

WGLMakie.activate!()
#= CairoMakie.activate!(type = "pdf") =#
set_aog_theme!()
# Set depending on size of screen
update_theme!(; resolution = (1300, 900))

#%%
# Revise will keep track of file change_arr and reload functions as necessary
includet(srcdir("Julia/DrWatson-helpers.jl"))
includet(funsdir("transmission-functions.jl"))
includet(funsdir("plotting-functions.jl"))
includet(funsdir("cleaning-functions.jl"))
includet(funsdir("SEIR-model.jl"))


#%%
################################################################################
########################## Ensemble Analysis ###################################
################################################################################


create_sir_quantiles_plot(
    ensemble_seir_summary; labels = state_labels, colors = seircolors,
    annual = true, caption = caption, tstep = param_dict[:tstep], xlims = (80, 100),
    ylims = (0, 1_000),
)

#%%
################################################################################
########################## Above-Below Analysis ################################
################################################################################
function calculate_outbreak_thresholds(outbreakrle)
    # Calculate upper and lower indices of consecutive days of infection
    outbreakaccum = accumulate(+, outbreakrle[2])
    outbreakuppers = outbreakaccum[findall(==(1), outbreakrle[1])]
    outbreaklowers = filter(
        x -> x <= maximum(outbreakuppers),
        outbreakaccum[findall(==(0), outbreakrle[1])] .+ 1,
    )

    return (outbreaklowers, outbreakuppers)
end

#%%
outbreakthreshold = 5
minoutbreakdur = 30
minoutbreaksize = 500

function create_inc_infec_arr!(
    incarr, ensemblejumparr, outbreakthreshold, minoutbreakdur, minoutbreaksize
)
    prog = Progress(size(ensemblejumparr, 3))
    @floop for sim in axes(ensemblejumparr, 3)
        # Copy new infections to array
        incarr[:, 1, sim] = @view(ensemblejumparr[1, :, sim])
        # Calculate if new infection is above or below threshold
        incarr[:, 2, sim] =
            @view(incarr[:, 1, sim]) .>= outbreakthreshold

        # Calculate the total number of infections above threshold in a consecutive string of days
        # Calculate the number of consecutive days of infection above or below threshold
        above5rle = rle(@view(incarr[:, 2, sim]))

        ## Calculate upper and lower indices of consecutive days of infection
        above5lowers, above5uppers = calculate_outbreak_thresholds(above5rle)

        for (lower, upper) in zip(above5lowers, above5uppers)
            # Calculate number of infections between lower and upper indices
            period_sum = sum(@view(incarr[lower:upper, 1, sim]))
            incarr[lower:upper, 3, sim] .= period_sum

            # Determine if there is an outbreak between lower and upper indices
            if upper - lower >= minoutbreakdur && period_sum >= minoutbreaksize
                incarr[lower:upper, 4, sim] .= 1
            end
        end

        next!(prog)
    end
end

function create_inc_infec_arr(
    ensemble_jump_arr, outbreakthreshold, minoutbreakdur, minoutbreaksize
)
    incarr = zeros(
        Int64, size(ensemble_jump_arr, 2), 4, size(ensemble_jump_arr, 3)
    )

    create_inc_infec_arr!(
        incarr,
        ensemble_jump_arr,
        outbreakthreshold,
        minoutbreakdur,
        minoutbreaksize,
    )

    return incarr
end

inc_infec_arr = create_inc_infec_arr(
    ensemble_jump_arr, outbreakthreshold, minoutbreakdur, minoutbreaksize
)

#%%
outbreakcols = [ColorSchemes.magma[i] for i in (200, 20)]

above5fig = Figure()
above5ax_prev = Axis(above5fig[1, 1]; ylabel = "Prevalence")
above5ax_inc = Axis(above5fig[2, 1]; ylabel = "Incidence")
above5ax_periodsum = Axis(
    above5fig[3, 1]; xlabel = "Time (years)", ylabel = "Period Sum"
)

linkxaxes!(above5ax_prev, above5ax_inc, above5ax_periodsum)

times = collect(0:param_dict[:tstep]:tmax) ./ 365

lines!(above5ax_prev, times, ensemble_seir_arr[2, :, 1])
lines!(above5ax_inc, times, inc_infec_arr[:, 1, 1])
outbreak_fig = barplot!(
    above5ax_periodsum,
    times,
    inc_infec_arr[:, 3, 1];
    color = inc_infec_arr[:, 4, 1],
    colormap = outbreakcols,
)

map(hidexdecorations!, [above5ax_prev, above5ax_inc])

map(
    ax -> xlims!(ax, (92, 94)),
    [above5ax_prev, above5ax_inc, above5ax_periodsum],
)
ylims!(above5ax_periodsum, (0, 10_000))
ylims!(above5ax_inc, (0, 300))

axislegend(
    above5ax_periodsum,
    [PolyElement(; color = col) for col in outbreakcols],
    ["Not Outbreak", "Outbreak"],
)

above5fig

#%%
################################################################################
########################## Background Noise ####################################
################################################################################
background_ode!(du, u, p, t) = (du .= 0.0)
background_noise!(du, u, p, t) = (du .= 0.1)

sde_condition(u, t, integrator) = true
function sde_affect!(integrator)
    if integrator.u[1] < 0.0
        integrator.u[1] = -integrator.u[1]
    end
end

sde_cb = DiscreteCallback(
    sde_condition, sde_affect!; save_positions = (false, false)
)

#%%
# Noise should be incidence, not prevalence
init_noise = [10.0]
tspan = (tmin, tmax)
noise_prob = SDEProblem(background_ode!, background_noise!, init_noise, tspan, p)
noise_sol = solve(
    noise_prob, SRIW1(); callback = sde_cb, tstep = param_dict[:tstep],
    adaptive = false,
)
noise_df = rename(DataFrame(noise_sol), [:time, :noise])

lines(noise_df[:, :time], noise_df[:, :noise])

#%%
# Set noise arr to 3D array (even though not necessary), so it has the same 
# dimensions as the other arrays
noise_arr = zeros(
    Float64, size(ensemble_jump_arr, 2), 1, size(ensemble_jump_arr, 3)
)
@floop for sim in axes(ensemble_jump_arr, 3)
    noise_prob = SDEProblem(
        background_ode!,
        background_noise!,
        init_noise,
        tspan,
        p
    )

    noise_sol = solve(
        noise_prob, SRIW1(); callback = sde_cb, tstep = param_dict[:tstep],
        adaptive = false,
    )

    noise_arr[:, 1, sim] = @view(noise_sol[1, :])
    # Set first noise incidence to 0 as no new noise in the first time step
    noise_arr[1, 1, sim] = 0.0
end

#%% 
# noise_fig = Figure()
# noise_ax = Axis(
#     noise_fig[1, 1]; xlabel = "Time (years)", ylabel = "Noise Incidence"
# )

# for sim in axes(noise_arr, 3)
#     lines!(noise_ax, times, noise_arr[:, 1, sim]; color = (:red, 0.1))
# end

# noise_fig

#%%
################################################################################
############################### Testing ########################################
################################################################################
testlag = 3
moveavglag = 7
detectthreshold = 10
perc_clinic = 0.3
perc_clinic_test = 0.8
perc_tested = perc_clinic * perc_clinic_test
testsens = 0.9
testspec = 0.9

testing_arr = zeros(Int64, tlength, 8, size(inc_infec_arr, 3));
post_odds_arr = zeros(Float64, tlength, 2, size(inc_infec_arr, 3));

function calculate_tested!(outarr, outarr_ind, inarr, perc_tested, sim)
    @. outarr[:, outarr_ind, sim] = round(@view(inarr[:, 1, sim]) * perc_tested)
end

function calculate_pos!(
    npos_vec,
    tested_vec,
    ntested,
    lag,
    sens,
    spec;
    noise = false
)
    if noise
        for day in eachindex(tested_vec)
            if day + lag <= ntested
                npos_vec[day + lag] = Int64(
                    round(tested_vec[day] * (1.0 - spec))
                )
            end
        end
    else
        for day in eachindex(tested_vec)
            if day + lag <= ntested
                npos_vec[day + lag] = Int64(
                    round(tested_vec[day] * sens)
                )
            end
        end
    end

    return nothing
end

function calculate_pos(
    tested_vec,
    lag,
    sens,
    spec;
    noise = false,
)
    ntested = length(tested_vec)
    npos = zeros(Int64, ntested)

    calculate_pos!(
        npos,
        tested_vec,
        ntested,
        lag,
        sens,
        spec;
        noise = noise,
    )

    return npos
end

#%%
function calculate_movingavg!(invec, outvec, testlag, avglag; Float = true)
    if Float
        avgfunc =
            (invec, day, avglag) -> mean(@view(invec[(day - avglag + 1):day]))
    else
        avgfunc =
            (invec, day, avglag) -> Int64(round(
                mean(@view(invec[(day - avglag + 1):day]))
            ))
    end
    for day in eachindex(invec)
        if day >= testlag + avglag + 1
            outvec[day] = avgfunc(invec, day, avglag)
        end
    end
end

function calculate_movingavg(invec, testlag, avglag)
    outvec = zeros(Float64, size(invec, 1), 1)

    calculate_movingavg!(invec, outvec, testlag, avglag)

    return outvec
end

function detectoutbreak!(outbreakvec, incvec, avgvec, threshold, avglag)
    @. outbreakvec = ifelse(incvec >= threshold || avgvec >= threshold, 1, 0)

    return nothing
end

function detectoutbreak(incvec, avgvec, threshold, avglag)
    outbreak = zeros(Int64, length(incvec))

    detectoutbreak!(outbreak, incvec, avgvec, threshold, avglag)

    return outbreak
end

#%%
function create_testing_arr!(
    testarr, incarr, noisearr, posoddsarr, perc_tested, testlag, testsens,
    testspec,
    detectthreshold, moveavglag,
)
    ntested = size(testarr, 1)

    # prog = Progress(size(incarr, 3))
    @floop for sim in 1:size(incarr, 3)
        # Number of infectious individuals tested
        calculate_tested!(testarr, 1, incarr, perc_tested, sim)

        # Number of noise individuals tested
        calculate_tested!(testarr, 2, noisearr, perc_tested, sim)

        # Number of test positive INFECTED individuals
        calculate_pos!(
            @view(testarr[:, 3, sim]),
            @view(testarr[:, 1, sim]),
            ntested,
            testlag,
            testsens,
            testspec;
            noise = false,
        )

        # Number of test positive NOISE individuals
        calculate_pos!(
            @view(testarr[:, 4, sim]),
            @view(testarr[:, 2, sim]),
            ntested,
            testlag,
            testsens,
            testspec;
            noise = true,
        )

        # Number of test positive TOTAL individuals
        @. testarr[:, 5, sim] =
            @view(testarr[:, 3, sim]) + @view(testarr[:, 4, sim])

        # Calculate moving average of TOTAL test positives
        calculate_movingavg!(
            @view(testarr[:, 5, sim]),
            @view(testarr[:, 6, sim]),
            testlag, moveavglag;
            Float = false,
        )

        # TOTAL Test positive individuals trigger outbreak response 
        detectoutbreak!(
            @view(testarr[:, 7, sim]),
            @view(testarr[:, 5, sim]),
            @view(testarr[:, 6, sim]),
            detectthreshold, moveavglag,
        )

        # # Posterior prob of infectious / total test positive
        @. @view(posoddsarr[:, 1, sim]) =
            @view(testarr[:, 3, sim]) / @view(testarr[:, 5, sim])
        calculate_movingavg!(
            @view(posoddsarr[:, 1, sim]),
            @view(posoddsarr[:, 2, sim]),
            testlag, moveavglag,
        )

        # Triggered outbreak equal to actual outbreak status
        @. testarr[:, 8, sim] =
            @view(testarr[:, 7, sim]) == @view(incarr[:, 4, sim])

        # next!(prog)
    end

    return nothing
end

function create_testing_arr(
    incarr, noisearr, perc_tested, testlag, testsens, testspec, detectthreshold,
    moveavglag,
)
    testarr = zeros(Int64, size(incarr, 1), 6, size(incarr, 3))
    posoddsarr = zeros(Float64, size(incarr, 1), 2, size(incarr, 3))

    create_testing_arr!(
        testarr, incarr, noisearr, posoddsarr, perc_tested, testlag, testsens,
        testspec,
        detectthreshold, moveavglag,
    )

    return testarr
end

#%%
# @benchmark create_testing_arr!(
#     $testing_arr,
#     $inc_infec_arr,
#     $noise_arr,
#     $post_odds_arr,
#     $perc_tested,
#     $testlag,
#     $testsens,
#     $testspec,
#     $detectthreshold,
#     $moveavglag,
# )

create_testing_arr!(
    testing_arr,
    inc_infec_arr,
    noise_arr,
    post_odds_arr,
    perc_tested,
    testlag,
    testsens,
    testspec,
    detectthreshold,
    moveavglag,
)

#%%
inc_test_fig = Figure()
inc_test_ax1 = Axis(inc_test_fig[1, 1]; ylabel = "Incidence")
inc_test_ax2 = Axis(inc_test_fig[2, 1]; ylabel = "Test Positive")
inc_test_ax3 = Axis(
    inc_test_fig[3, 1];
    xlabel = "Time (years)",
    ylabel = "7d Avg Test Positive"
)

lines!(
    inc_test_ax1, times, inc_infec_arr[:, 1, 1];
    color = inc_infec_arr[:, 4, 1],
    colormap = outbreakcols,
)
lines!(
    inc_test_ax2, times, testing_arr[:, 3, 1];
    color = testing_arr[:, 7, 1],
    colormap = outbreakcols,
)
lines!(
    inc_test_ax3, times, testing_arr[:, 6, 1];
    color = testing_arr[:, 7, 1],
    colormap = outbreakcols,
)

linkxaxes!(inc_test_ax1, inc_test_ax2, inc_test_ax3)

map(hidexdecorations!, [inc_test_ax1, inc_test_ax2])

# map(
#     ax -> xlims!(ax, (1750 / 365, 1850 / 365)),
#     [inc_test_ax1, inc_test_ax2, inc_test_ax3],
# )
map(ax -> ylims!(ax, (0, 20)), [inc_test_ax1, inc_test_ax2, inc_test_ax3])

hlines!(
    inc_test_ax1, 5;
    color = :black,
    linestyle = :dash,
    linewidth = 2,
)

map(
    ax -> hlines!(
        ax,
        detectthreshold;
        color = :black,
        linestyle = :dash,
        linewidth = 2,
    ),
    [inc_test_ax2, inc_test_ax3],
)

Legend(
    inc_test_fig[1, 2],
    [PolyElement(; color = col) for col in outbreakcols],
    ["Not Outbreak", "Outbreak"],
    "True\nOutbreak Status",
)

Legend(
    inc_test_fig[2:3, 2],
    [PolyElement(; color = col) for col in outbreakcols],
    ["Not Outbreak", "Outbreak"],
    "Detected\nOutbreak Status",
)

inc_test_fig

#%%
testing_fig = Figure()
plot_test_sims = 4
plot_test_coords = 1:(plot_test_sims ÷ 2)
for (sim, ax) in
    zip(1:plot_test_sims, IterTools.product(plot_test_coords, plot_test_coords))
    row = ax[1]
    col = ax[2]

    fig_ax = Symbol("testing_ax_" * string(sim))

    @eval $(fig_ax) = Axis(
        testing_fig[$row, $col]; xlabel = "Time (years)",
        ylabel = "Tested"
    )

    @eval lines!(
        $(fig_ax), times, testing_arr[:, 1, $sim];
        color = :red,
        label = "Infectious",
    )
    @eval lines!(
        $(fig_ax), times, testing_arr[:, 2, $sim];
        color = :blue,
        label = "Noise",
    )
    @eval lines!(
        $(fig_ax), times, testing_arr[:, 5, $sim];
        color = :black,
        label = "Total Positive",
    )
end

linkxaxes!(testing_ax_1, testing_ax_2)
linkxaxes!(testing_ax_3, testing_ax_4)

linkyaxes!(testing_ax_1, testing_ax_3)
linkyaxes!(testing_ax_2, testing_ax_4)

map(hidexdecorations!, [testing_ax_1, testing_ax_3])
map(hideydecorations!, [testing_ax_3, testing_ax_4])

Legend(
    testing_fig[3, :],
    testing_ax_1,
    "Type of Individual";
    orientation = :horizontal,
)

testing_fig

#%%
@proto struct OutbreakThresholdChars{A,B,C,D}
    crosstab::A
    tp::B
    tn::B
    fp::B
    fn::B
    sens::C
    spec::C
    ppv::C
    npv::C
    noutbreaks::B
    ndetectoutbreaks::B
    outbreakbounds::D
    detectoutbreakbounds::D
end

#%%
function calculate_ot_characterstics(test_arr, infec_arr, ind)
    crosstab = freqtable(testing_arr[:, 5, ind], inc_infec_arr[:, 4, ind])

    tp = crosstab[2, 2]
    tn = crosstab[1, 1]
    fp = crosstab[2, 1]
    fn = crosstab[1, 2]

    sens = tp / (tp + fn)
    spec = tn / (tn + fp)

    ppv = tp / (tp + fp)
    npv = tn / (tn + fn)

    return crosstab, tp, tn, fp, fn, sens, spec, ppv, npv
end

function calculate_noutbreaks(outbreakrle)
    return length(findall(==(1), outbreakrle[1]))
end

#%%
OT_chars = ThreadsX.map(
    axes(inc_infec_arr, 3)
) do sim
    outbreakrle = rle(@view(inc_infec_arr[:, 4, sim]))
    detectrle = rle(@view(testing_arr[:, 7, sim]))

    OutbreakThresholdChars(
        calculate_ot_characterstics(testing_arr, inc_infec_arr, sim)...,
        calculate_noutbreaks(outbreakrle),
        calculate_noutbreaks(detectrle),
        reduce(hcat, collect(calculate_outbreak_thresholds(outbreakrle))),
        reduce(hcat, collect(calculate_outbreak_thresholds(detectrle))),
    )
end

#%%
OT_chars[1].crosstab
OT_chars[1].outbreakbounds
OT_chars[1].detectoutbreakbounds
OT_chars[1].noutbreaks
OT_chars[1].ndetectoutbreaks

# Note that an outbreak isn't detected continously!
testing_arr[80:100, :, 1]

#%%
otchars_vec = zeros(Float64, length(OT_chars), 6);

@floop for sim in eachindex(OT_chars)
    otchars_vec[sim, 1] = OT_chars[sim].sens
    otchars_vec[sim, 2] = OT_chars[sim].spec
    otchars_vec[sim, 3] = OT_chars[sim].ppv
    otchars_vec[sim, 4] = OT_chars[sim].npv
    otchars_vec[sim, 5] = OT_chars[sim].noutbreaks
    otchars_vec[sim, 6] = OT_chars[sim].ndetectoutbreaks
end

#%%
outbreak_dist_fig = Figure()
outbreak_dist_ax = Axis(
    outbreak_dist_fig[1, 1]; xlabel = "Proportion of Time Series with Outbreak"
)

hist!(
    outbreak_dist_ax,
    vec(sum(@view(inc_infec_arr[:, 4, :]); dims = 1)) ./ size(inc_infec_arr, 1);
    bins = 0.0:0.01:0.7,
    color = (:blue, 0.5),
    strokecolor = :black,
    strokewidth = 1,
    normalization = :pdf,
    label = "True Outbreaks",
)

hist!(
    outbreak_dist_ax,
    vec(sum(@view(testing_arr[:, 7, :]); dims = 1)) ./ size(testing_arr, 1);
    bins = 0.0:0.01:0.7,
    color = (:red, 0.5),
    strokecolor = :black,
    strokewidth = 1,
    normalization = :pdf,
    label = "Tested Outbreaks",
)

Legend(outbreak_dist_fig[1, 2], outbreak_dist_ax, "Outbreak Proportion")

outbreak_dist_fig

#%%
noutbreaks_fig = Figure()
noutbreaks_ax = Axis(noutbreaks_fig[1, 1]; xlabel = "Number of Outbreaks")

hist!(
    noutbreaks_ax,
    @view(otchars_vec[:, 5]);
    bins = 0.0:10.0:450.0,
    color = (:blue, 0.5),
    strokecolor = :black,
    strokewidth = 1,
    normalization = :pdf,
    label = "True Outbreaks",
)

hist!(
    noutbreaks_ax,
    @view(otchars_vec[:, 6]);
    bins = 0.0:10.0:450.0,
    color = (:red, 0.5),
    strokecolor = :black,
    strokewidth = 1,
    normalization = :pdf,
    label = "Tested Outbreaks",
)

Legend(noutbreaks_fig[1, 2], noutbreaks_ax, "# Outbreaks")

noutbreaks_fig

#%%
sens_spec_fig = Figure()
sens_spec_ax = Axis(sens_spec_fig[1, 1]; xticks = 0.0:0.1:1.0)

hist!(
    sens_spec_ax,
    @view(otchars_vec[:, 1]);
    bins = 0.0:0.01:1.01,
    color = (:blue, 0.5),
    strokecolor = :black,
    strokewidth = 1,
    label = "Sensitivity",
    normalization = :pdf,
)

hist!(
    sens_spec_ax,
    @view(otchars_vec[:, 2]);
    bins = 0.0:0.01:1.01,
    color = (:red, 0.5),
    strokecolor = :black,
    strokewidth = 1,
    label = "Specificity",
    normalization = :pdf,
)

vlines!(
    sens_spec_ax,
    [mean(@view(otchars_vec[:, i])) for i in 1:2];
    color = :black,
    linestyle = :dash,
    linewidth = 4,
)

Legend(sens_spec_fig[1, 2], sens_spec_ax, "Characteristic")

sens_spec_fig

#%%
ppv_npv_fig = Figure()
ppv_npv_ax = Axis(ppv_npv_fig[1, 1]; xticks = 0.0:0.1:1.0)

hist!(
    ppv_npv_ax,
    @view(otchars_vec[:, 3]);
    bins = 0.0:0.01:1.01,
    color = (:green, 0.5),
    strokecolor = :black,
    strokewidth = 1,
    label = "PPV",
    normalization = :pdf,
)

hist!(
    ppv_npv_ax,
    @view(otchars_vec[:, 4]);
    bins = 0.0:0.01:1.01,
    color = (:purple, 0.5),
    strokecolor = :black,
    strokewidth = 1,
    label = "NPV",
    normalization = :pdf,
)

vlines!(
    ppv_npv_ax,
    [mean(@view(otchars_vec[:, i])) for i in 3:4];
    color = :black,
    linestyle = :dash,
    linewidth = 4,
)

Legend(ppv_npv_fig[1, 2], ppv_npv_ax, "Characteristic")

ppv_npv_fig

#%%
# TODO: Calculate threshold chars over range of test lags and detection thresholds. Use constant R0 and ind test chars for the moment.
@proto struct MeanOutbreakThresholdChars{A,B}
    testlag::A
    incthreshold::A
    tp::A
    tn::A
    fp::A
    fn::A
    sens::B
    spec::B
    ppv::B
    npv::B
    noutbreaks::A
    ndetectoutbreaks::A
end

testlag_vec = 0:1:2
detectthreshold_vec = [1, collect(5:5:10)...]

function create_OTchars_struct(incarr, testarr)
    ThreadsX.map(
        axes(incarr, 3)
    ) do sim
        outbreakrle = rle(@view(incarr[:, 4, sim]))
        return detectrle = rle(@view(testarr[:, 5, sim]))

        OutbreakThresholdChars(
            calculate_ot_characterstics(testarr, incarr, sim)...,
            calculate_noutbreaks(outbreakrle),
            calculate_noutbreaks(detectrle),
            reduce(hcat, collect(calculate_outbreak_thresholds(outbreakrle))),
            reduce(hcat, collect(calculate_outbreak_thresholds(detectrle))),
        )
    end
end

function calculate_mean_ot_chars(
    testlag,
    detectthreshold;
    ensemblejumparr = ensemble_jump_arr,
    outbreakthreshold = outbreakthreshold,
    minoutbreakdur = minoutbreakdur,
    minoutbreaksize = minoutbreaksize,
    noisearr = noise_arr,
    perctested = perc_tested,
    testsens = testsens,
    testspec = testspec,
    moveavglag = moveavglag,
)
    @info "Creating Incidence Array"
    incarr = create_inc_infec_arr(
        ensemble_jump_arr,
        outbreakthreshold,
        minoutbreakdur,
        minoutbreaksize
    )

    @info "Creating Testing Array"
    testing_arr = zeros(
        Int64, size(incarr, 1), 6, size(incarr, 3)
    )

    create_testing_arr!(
        testing_arr,
        incarr,
        noise_arr,
        perc_tested,
        testlag,
        testsens,
        testspec,
        detectthreshold,
        moveavglag,
    )

    @info "Calculating OT characteristics"
    OT_chars = create_OTchars_struct(incarr, testing_arr)

    return OT_chars
end

#%%
test = calculate_mean_ot_chars(7, 10)

#%%
for (testlag, detectthreshold) in
    IterTools.product(testlag_vec, detectthreshold_vec)
end
