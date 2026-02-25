export AlertOutbreakMatchingStrategy,
    SingleOutbreakPerAlert,
    MultipleOutbreaksPerAlert

"""
    AbstractAlertOutbreakMatchingStrategy

Abstract type for strategies that define how alerts are matched to outbreaks.
"""
abstract type AbstractAlertOutbreakMatchingStrategy end

"""
    SingleOutbreakPerAlert <: AbstractAlertOutbreakMatchingStrategy

Match each alert to at most one outbreak (the first overlapping outbreak).

This is the historical and default behavior for outbreak-alert matching.
"""
struct SingleOutbreakPerAlert <: AbstractAlertOutbreakMatchingStrategy end

"""
    MultipleOutbreaksPerAlert <: AbstractAlertOutbreakMatchingStrategy

Allow a single alert to match multiple overlapping outbreaks.

This strategy is useful for long alerts that overlap more than one outbreak.
"""
struct MultipleOutbreaksPerAlert <: AbstractAlertOutbreakMatchingStrategy end

"""
    AlertOutbreakMatchingStrategy

Sum type for alert-outbreak matching behavior.

Variants:
- `SingleOutbreakPerAlert`: Each alert can match at most one outbreak.
- `MultipleOutbreaksPerAlert`: A single alert can match multiple outbreaks.
"""
LightSumTypes.@sumtype AlertOutbreakMatchingStrategy(
    SingleOutbreakPerAlert,
    MultipleOutbreaksPerAlert,
) <: AbstractAlertOutbreakMatchingStrategy
