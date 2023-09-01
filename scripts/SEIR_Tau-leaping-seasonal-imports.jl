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
############################### Testing ########################################
################################################################################



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
