export run_OutbreakThresholdChars_creation, OutbreakThresholdChars_creation,
    get_ensemble_file

"""
    run_define_outbreaks(dict_of_outbreak_spec_params)

Define outbreaks for multiple outbreak specifications.
"""
function run_define_outbreaks(dict_of_outbreak_spec_params)
    for outbreak_spec_params in dict_of_outbreak_spec_params
        DrWatson.@produce_or_load(
            define_outbreaks,
            outbreak_spec_params,
            "$(outbreak_spec_params[:dirpath])";
            filename = "ensemble-incidence-array",
            loadfile = false,
            force = true
        )
    end
    return
end

"""
    define_outbreaks(incidence_param_dict)

Define outbreaks and create scenario specifications.
"""
function define_outbreaks(incidence_param_dict)
    UnPack.@unpack ensemble_spec,
        ensemble_inc_vecs,
        outbreak_spec,
        noise_spec_vec,
        outbreak_detection_spec_vec,
        test_spec_vec, seed =
        incidence_param_dict

    ensemble_inc_arr, ensemble_thresholds_vec = create_inc_infec_arr(
        ensemble_inc_vecs, outbreak_spec
    )

    non_clinical_case_test_spec_vec = filter(
        spec -> !(spec in CLINICAL_TEST_SPECS),
        test_spec_vec,
    )

    non_clinical_case_outbreak_detection_spec_vec = filter(
        spec -> spec.percent_clinic_tested !== 1.0,
        outbreak_detection_spec_vec,
    )

    non_clinical_case_ensemble_scenarios = create_combinations_vec(
        ScenarioSpecification,
        (
            [ensemble_spec],
            [outbreak_spec],
            noise_spec_vec,
            non_clinical_case_outbreak_detection_spec_vec,
            non_clinical_case_test_spec_vec,
        ),
    )

    clinical_case_outbreak_detection_spec_vec = filter(
        spec -> spec.percent_clinic_tested == 1.0,
        outbreak_detection_spec_vec,
    )

    clinical_case_ensemble_scenarios = create_combinations_vec(
        ScenarioSpecification,
        (
            [ensemble_spec],
            [outbreak_spec],
            noise_spec_vec,
            clinical_case_outbreak_detection_spec_vec,
            CLINICAL_TEST_SPECS,
        ),
    )

    ensemble_scenarios = vcat(
        non_clinical_case_ensemble_scenarios, clinical_case_ensemble_scenarios
    )

    scenario_param_dict = DrWatson.dict_list(
        DrWatson.@dict(
            scenario_spec = ensemble_scenarios,
            ensemble_inc_arr,
            thresholds_vec = [ensemble_thresholds_vec],
            seed = seed
        )
    )

    run_OutbreakThresholdChars_creation(scenario_param_dict)

    return DrWatson.@strdict ensemble_inc_arr ensemble_thresholds_vec
end

"""
    run_OutbreakThresholdChars_creation(dict_of_OTchars_params)

Create outbreak threshold characteristics for multiple scenarios.
"""
function run_OutbreakThresholdChars_creation(
        dict_of_OTchars_params
    )
    return FLoops.@floop for OTChars_params in dict_of_OTchars_params
        DrWatson.@produce_or_load(
            OutbreakThresholdChars_creation,
            OTChars_params,
            "$(OTChars_params[:scenario_spec].dirpath)";
            filename = "ensemble-scenario",
            loadfile = false
        )
    end
end

"""
    OutbreakThresholdChars_creation(OT_chars_param_dict)

Create outbreak threshold characteristics for a single scenario.
"""
function OutbreakThresholdChars_creation(OT_chars_param_dict)
    UnPack.@unpack scenario_spec, ensemble_inc_arr, thresholds_vec, seed =
        OT_chars_param_dict
    UnPack.@unpack noise_specification,
        outbreak_specification,
        outbreak_detection_specification,
        individual_test_specification = scenario_spec

    noise_array, noise_means = create_noise_arr(
        noise_specification,
        ensemble_inc_arr;
        ensemble_specification = scenario_spec.ensemble_specification,
        seed = seed,
    )

    testarr, test_movingvg_arr = create_testing_arrs(
        ensemble_inc_arr,
        noise_array,
        outbreak_detection_specification,
        individual_test_specification,
    )

    OT_chars = calculate_OutbreakThresholdChars(
        testarr, ensemble_inc_arr, thresholds_vec, noise_means
    )

    return DrWatson.@strdict OT_chars
end

"""
    get_ensemble_file()

Get ensemble file from disk. Multiple dispatch for different specification types.
"""
function get_ensemble_file() end

function get_ensemble_file(
        ensemble_spec::EnsembleSpecification, outbreak_spec::OutbreakSpecification
    )
    dirpath = joinpath(ensemble_spec.dirpath, outbreak_spec.dirpath)

    return JLD2.load(collect_ensemble_file("incidence-array", dirpath)...)
end

function get_ensemble_file(spec::EnsembleSpecification)
    return JLD2.load(collect_ensemble_file("solution", spec.dirpath)...)
end

function get_ensemble_file(spec::EnsembleSpecification, quantile::Int64)
    return JLD2.load(
        collect_ensemble_file("quantiles_$(quantile)", spec.dirpath)...
    )
end

function get_ensemble_file(spec::ScenarioSpecification)
    return JLD2.load(collect_ensemble_file("scenario", spec.dirpath)...)
end

function collect_ensemble_file(type, dirpath)
    filecontainer = []
    for f in readdir(dirpath)
        match_ensemble_file!(type, dirpath, filecontainer, f)
    end
    if length(filecontainer) != 1
        println("Matched $(length(filecontainer)) files, when should be 1")
    end
    return filecontainer
end

function match_ensemble_file!(criteria, dirpath, container, file)
    return if occursin(criteria, file)
        push!(container, joinpath(dirpath, file))
    end
end
