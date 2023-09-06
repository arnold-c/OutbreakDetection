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
