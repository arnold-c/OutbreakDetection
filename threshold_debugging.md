# Threshold Debugging Log

## Goal
Match the exact threshold values from the old code version at revision lwxvloll, which used the previous DataFrame based approach of setting up the optimizations.
It also handled a number of steps differently, which will be explained one-by-one.
The values in the summary DataFrame below demonstrate the optimal thresholds at 8x dynamical noise.

```
4×13 DataFrame
 Row │ Noise Type            Test Type             Test Lag  10%       20%       30%       40%       50%       60%       70%       80%       90%       100%
     │ Cat…                  String                Int64     Float64?  Float64?  Float64?  Float64?  Float64?  Float64?  Float64?  Float64?  Float64?  Float64?
─────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │ Dynamical noise       Imperfect Test (85%)         0     1.172     2.734     3.516     5.859     6.641     0.391     0.391     0.391    11.328    12.109
   2 │ Dynamical noise       Imperfect Test (90%)         0     1.172     2.734     3.516     4.297     6.641     7.422     8.984    10.547    11.328     0.391
   3 │ All noise structures  Perfect Test                14     1.172     2.734     3.516     5.859     6.641     7.422     8.984    10.547    11.328    12.109
   4 │ All noise structures  Perfect Test                 0     1.172     2.734     3.516     5.859     6.641     7.422     8.984    10.547    11.328    12.109
```

## Things that had to be reverted in this JJ revision

To match the exact optimization values using the new StructVector approach that also included some other refactoring (jj change id `olzstqvm`, git commit id `4648eb3c`), some decisions and changes had to be reversed.
This document summarizes those reversals.
By far the largest impact came from changing the calculation of the alert-outbreak matching logic.

### Moving average calculation

`./OutbreakDetectionCore/src/utilities/calculate-moving-average.jl`

The core calculations are the same, but instead of returning a Float64 for the moving average (of the test positives, principally), round and convert the moving average to an integer.

### Noise dynamics parameters

`./OutbreakDetectionCore/src/noise/noise-dynamics-parameters.jl`

In the previous version, the noise dynamics parameters struct was creating using the target disease's `beta_mean` values i.e,. measles' `beta_mean`, not the value calculated for rubella.

### Parameter sampling

In the previous version, there were not that many parameters that were sampled in each of the simulations of the ensemble.
In the refactored version, most are sampled.

`./OutbreakDetectionCore/src/diagnostic-testing/calculate-num-tested.jl`
`./OutbreakDetectionCore/src/diagnostic-testing/calculate-num-positive.jl`
`./OutbreakDetectionCore/src/types/dynamics-parameters.jl`

Instead of sampling from a Binomial distribution to calculate the numbers tested and positive, just multiple the integer number by the percentage, round, and then convert to an integer.
This was changed to sampling to reduce the impact of small numbers all being rounded to 0 equally.

### Beta value calculation

`./OutbreakDetectionCore/src/simulation/transmission-functions.jl`

The previous version calculated `beta` using the SIR equation, not SEIR equation that includes the movement out of the latent state.

`./OutbreakDetectionCore/src/types/dynamics-parameters.jl`

The `DynamicsParameterSpecification` constructor also used to normalize the value of `beta_mean` by the initial population size `N`, rather than using a strictly frequency-dependent approach where the normalization is included in the simulation of infections.

```julia
function DynamicsParameterSpecification(
        state_specification::StateParameters,
        target_disease_dynamics_params::TargetDiseaseDynamicsParameters,
        common_disease_dynamics_params::CommonDiseaseDynamicsParameters,
    )
...
    beta_mean = calculate_beta(
        target_disease_dynamics_params.R_0,
        sigma,
        gamma,
        mu
    ) / state_specification.init_states.N
    ...
end
```

### Threshold bounds

`./scripts/optimal-thresholds_optims.jl`
`./scripts/optimal-thresholds_plots.jl`

The optimisation of the threshold previously used bounds between 0 and 50.
When running the new version with values between 0 and 20, different Sobol' points are created, so slightly different optimum are found.

### Fixed noise vaccination coverage

`./OutbreakDetectionCore/src/threshold-optimization/evaluate-missing-optimizations.jl`

As part of the refactoring, the new implementation uses MultistartOptimization to optimize the vaccination level in the dynamical noise simulations to achieve X times noise level.
The benefit of the refactor is that if any other scenarios need to be run, or a change to the target or noise disease etc. then the optimization will handle those changes.

### Endemic noise initialization

`./OutbreakDetectionCore/src/noise/noise-recreation.jl`

During the refactoring, calculations of the endemic state for the noise dynamics was added to seed the noise simulations in an endemic state to avoid always starting with a large outbreak.
This wasn't present in the original implementation.


### SEIR model beta value index

`./OutbreakDetectionCore/src/simulation/seir-model.jl`

In the refactored version, the SEIR model loop uses the previous value of the beta to determine how many individuals are infected at the current time step (also using the previous values of the states).
This is correct as the otherwise the individuals are being affected by the state of a dynamical system that doesn't (quite) exist yet.
However, the previous implementation doesn't do this, so roll back for the testing.
Makes a very negligible difference.

```julia
function seir_mod!(
        state_vec::ASV,
        inc_vec::AI,
        Reff_vec::AF1,
        beta_vec::AF2,
        states::StaticArrays.SVector{5, Int64},
        dynamics_params::DynamicsParameters,
        time_params::SimTimeParameters,
    ) where {
        ASV <: AbstractVector{StaticArrays.SVector{5, Int64}},
        AI <: AbstractVector{Int64},
        AF1 <: AbstractVector{Float64},
        AF2 <: AbstractVector{Float64},
    }
    ...

    @inbounds for i in 2:(tlength)
        vaccination_coverage = dynamics_params.vaccination_coverage

        state_vec[i], inc_vec[i] = seir_mod_loop(
            state_vec[i - 1],
            beta_vec[i], # <- change here from beta_vec[i-1]
            mu_timestep,
            epsilon_timestep,
            sigma_timestep,
            gamma_timestep,
            R_0,
            vaccination_coverage,
            timestep,
        )
        ...
    end
    return nothing
end
```

### SEIR model imported infections

`./OutbreakDetectionCore/src/simulation/seir-model.jl`

In the refactored version, we calculated the number of importer individuals as a binomial sample from the number of remaining susceptible individuals (after removing those from direct contact infections).
The previous version approximated this as a function of the total population and $R_0$.
The refactored version could be considered more correct as only susceptible individuals can be infected, but ultimately that's a framing question as to whether they are commuter-style imports or not (the refactored version likely results in slightly fewer infections when prevalence is high, and more infections when prevalence is low).

```julia
function seir_mod_loop(...)::Tuple{StaticArrays.SVector{5, Int64}, Int64}

    @inbounds begin
        ...
        import_inf = Random.rand(
            Distributions.Poisson(epsilon_timestep * N / R_0) # <- previously epsilon_timestep * remaining_S
        ) # Import: S -> E
        ...
    end
    ...
end
```

### SEIR model Poisson transitions

`./OutbreakDetectionCore/src/simulation/seir-model.jl`

In the refactored version, we use the `_smart_transition()` function to use Poisson samples when the size of the compartments is large relative to the rate, and Binomial when small.
This is done for computational efficiency, and makes a negligible difference to the number of individuals moving between compartments.

### Alert-outbreak matching

`./OutbreakDetectionCore/src/detection/match-alert-outbreak-thresholds.jl`

By far the biggest reason for the difference between the refactored and original implementation results comes from the change to  the alert-outbreak matching.
As can be seen by the plots below, updating how the matching of outbreaks and alerts is performed has a drastic impact on imperfect tests at high dynamical noise.
In the previous implementation, a single alert can be matched correctly to multiple outbreaks, assuming the alert is long enough to span the gap between successive outbreaks.
In the newer implementation, an alert is only matched to the first outbreak that it overlaps with.
This refactor also affects perfect tests with a lag, albeit to a much smaller degree.
In one respect, this accounts for the reality that once an alert has triggered and been investigated, it probably is not going to get much follow up if there is not a cessation in alert status.

#### Results including all reversions (inc. matching)

![](./plots/optimal-thresholds_alert-threshold-plot_all-revisions.svg)
![](./plots/optimal-thresholds_accuracy-plot_all-revisions.svg)

#### Results including all reversions (exc. matching)

![](./plots/optimal-thresholds_alert-threshold-plot_new-matching.svg)
![](./plots/optimal-thresholds_accuracy-plot_new-matching.svg)
