export get_test_description,
    table_test_type,
    plot_test_description

"""
    get_test_description(test_specification::IndividualTestSpecification)

Get a human-readable description of the test specification.
"""
function get_test_description(test_specification::IndividualTestSpecification)
    sens = test_specification.sensitivity
    spec = test_specification.specificity
    lag = test_specification.test_result_lag

    description = if sens == 1.0 && spec == 0.0 && lag == 0
        "Clinical case definition"
    elseif (sens == 0.98 && spec == 0.98) || (sens == 0.95 && spec == 0.98)
        "ELISA-like ($(sens * 100)% sens, $(spec * 100)% spec, $(lag) day lag)"
    elseif sens == 1.0 && spec == 1.0
        "Perfect test ($(lag) day lag)"
    elseif sens == spec && sens < 1.0 && lag == 0
        "RDT-like ($(sens * 100)% sens/spec)"
    else
        error("Unknown test specification pattern: sensitivity=$(sens), specificity=$(spec), lag=$(lag)")
    end

    return description
end

"""
    table_test_type(sensitivity, specificity, test_lag)

Get a table-formatted test type description.
"""
function table_test_type(
        sensitivity::Float64,
        specificity::Float64,
        test_lag::Int64
    )
    return if sensitivity == 1.0 && specificity == 0.0 && test_lag == 0
        "Clinical Case Definition"
    elseif sensitivity == specificity && sensitivity < 1.0 && test_lag == 0
        "Imperfect Test ($(Int64(round(sensitivity * 100; digits = 0)))%)"
    elseif sensitivity == 1.0 && specificity == 1.0
        "Perfect Test"
    else
        error("Unknown test type pattern: sensitivity=$(sensitivity), specificity=$(specificity), test_lag=$(test_lag)")
    end
end

"""
    plot_test_description(test_specification::IndividualTestSpecification)

Get a plot-formatted test description.
"""
function plot_test_description(test_specification::IndividualTestSpecification)
    sens = test_specification.sensitivity
    spec = test_specification.specificity
    lag = test_specification.test_result_lag

    return if sens == 1.0 && spec == 0.0 && lag == 0
        "Clinical Case Definition"
    elseif sens == spec && sens < 1.0 && lag == 0
        "Imperfect Test ($(Int64(round(sens * 100; digits = 0)))%)"
    elseif sens == 1.0 && spec == 1.0
        "Perfect Test ($lag-day lag)"
    else
        error("Unknown test specification pattern: sensitivity=$(sens), specificity=$(spec), lag=$(lag)")
    end
end
