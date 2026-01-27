export AlertFilteringStrategy, AllAlerts, PostOutbreakStartAlerts

"""
    AbstractAlertFilteringStrategy

Abstract type for alert filtering strategies used in outbreak detection matching.

Subtypes define different rules for which alerts should be considered when matching
alerts to outbreaks during threshold optimization.
"""
abstract type AbstractAlertFilteringStrategy end

"""
    AllAlerts <: AbstractAlertFilteringStrategy

Include all alerts when matching to outbreaks, regardless of timing.

This is the default behavior that includes alerts even if they start before
the outbreak begins (resulting in negative detection delays).

# Example
```julia
strategy = AlertFilteringStrategy(AllAlerts())
matched = match_outbreak_detection_bounds(outbreaks, alerts, strategy)
```

# See Also
- [`PostOutbreakStartAlerts`](@ref): Alternative that excludes pre-outbreak alerts
- [`AlertFilteringStrategy`](@ref): Sum type containing filtering strategies
"""
struct AllAlerts <: AbstractAlertFilteringStrategy end

"""
    PostOutbreakStartAlerts <: AbstractAlertFilteringStrategy

Only include alerts that start on or after the outbreak start time.

This filtering strategy excludes alerts that begin before the outbreak starts,
ensuring all detection delays are non-negative. This is useful when you want
to evaluate detection performance only for alerts that could plausibly be
responses to the outbreak rather than coincidental pre-outbreak alerts.

# Matching Rule
An alert is included only if `alert_lower >= outbreak_lower`.

# Example
```julia
strategy = AlertFilteringStrategy(PostOutbreakStartAlerts())
matched = match_outbreak_detection_bounds(outbreaks, alerts, strategy)
# All detection delays in matched results will be >= 0
```

# See Also
- [`AllAlerts`](@ref): Alternative that includes all alerts
- [`AlertFilteringStrategy`](@ref): Sum type containing filtering strategies
"""
struct PostOutbreakStartAlerts <: AbstractAlertFilteringStrategy end

"""
    AlertFilteringStrategy

Sum type for alert filtering strategies in outbreak detection matching.

This is a LightSumTypes-based sum type that can hold one of two variants:
- `AllAlerts`: Include all alerts (default, backward compatible)
- `PostOutbreakStartAlerts`: Only include alerts starting on or after outbreak start

# Usage
```julia
# Default: include all alerts
strategy = AlertFilteringStrategy(AllAlerts())

# Only post-outbreak alerts
strategy = AlertFilteringStrategy(PostOutbreakStartAlerts())

# Use in matching
matched = match_outbreak_detection_bounds(outbreaks, alerts, strategy)
```

# See Also
- [`AllAlerts`](@ref): Include all alerts
- [`PostOutbreakStartAlerts`](@ref): Exclude pre-outbreak alerts
- [`match_outbreak_detection_bounds`](@ref): Function that uses this strategy
"""
LightSumTypes.@sumtype AlertFilteringStrategy(
    AllAlerts, PostOutbreakStartAlerts
) <: AbstractAlertFilteringStrategy
