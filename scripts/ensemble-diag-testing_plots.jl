#%%
using DrWatson
@quickactivate "OutbreakDetection"

includet(srcdir("Julia/DrWatson-helpers.jl"))
includet(scriptsdir("ensemble-diag-testing-sim.jl"))
includet(funsdir("plotting-functions.jl"))

#%%
incidence_testing_plot(
    inc_infec_arr,
    testing_arr,
    time_p,
    detectthreshold;
    sim = 1
)

#%%
testing_plot(testing_arr, time_p)

#%%
ensemble_outbreak_distribution_plot(testing_arr, inc_infec_arr)

#%%
ensemble_OTChars_plot(
    OT_Chars,
    :sensitivity,
    :specificity;
    bins = 0.0:0.01:1.01,
    char1_label = "Sensitivity",
    char2_label = "Specificity",
    xlabel = "Proportion",
    legendlabel = "Characterstic",
)

#%%
ensemble_OTChars_plot(OT_Chars, :noutbreaks, :ndetectoutbreaks)

#%%
ensemble_OTChars_plot(
    OT_Chars,
    :ppv,
    :npv;
    bins = 0.0:0.01:1.01,
    char1_label = "PPV",
    char2_label = "NPV",
    char1_color = :green,
    char2_color = :purple,
    xlabel = "Proportion",
    legendlabel = "Characterstic",
)
