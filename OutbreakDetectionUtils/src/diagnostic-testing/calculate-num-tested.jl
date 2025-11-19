export calculate_tested!

"""
    calculate_tested!(outvec, invec, perc_tested)

Calculate number of individuals tested based on percentage.
"""
function calculate_tested!(outvec, invec, perc_tested)
    return @. outvec = round(invec * perc_tested)
end
