using OutbreakDetectionUtils: StatsBase
#%%
using DrWatson
@quickactivate "OutbreakDetection"

using UnPack: @unpack
using StatsBase: StatsBase
using OutbreakDetectionUtils

#%%
ensemble_model_type = ("seasonal-infectivity-import", "tau-leaping")

nyears = 100
ensemble_time_specification = SimTimeParameters(;
    tmin = 0.0, tmax = 365.0 * nyears, tstep = 1.0
)

ensemble_state_specification = StateParameters(
    500_000,
    Dict(:s_prop => 0.05, :e_prop => 0.0, :i_prop => 0.0, :r_prop => 0.95),
)

ensemble_dynamics_specification = DynamicsParameters(
    ensemble_state_specification.init_states.N,
    27,
    0.2,
    SIGMA,
    GAMMA,
    16.0,
    0.8,
)

ensemble_nsims = 100

ensemble_specification = EnsembleSpecification(
    ensemble_model_type,
    ensemble_state_specification,
    ensemble_dynamics_specification,
    ensemble_time_specification,
    ensemble_nsims,
)

#%%
dynamical_noise_R0 = [5.0]
dynamical_noise_latent_period = [7]
dynamical_noise_duration_infection = [14]
dynamical_noise_correlation = [
    "in-phase"
    # "out-of-phase",
    # "none",
]
dynamical_noise_mean_scaling_vec = [0.15]

#%%
function calculate_dynamic_vaccination_coverage(
    target_scaling,
    measles_daily_incidence,
    dynamical_noise_specification_parameters,
    ensemble_specification;
    vaccination_range = [0.0, 1.0],
    maxiters = 10,
    atol = 0.02,
    showprogress = false,
)
    @unpack R0,
    latent_period,
    duration_infection,
    correlation,
    poisson_component = dynamical_noise_specification_parameters

    @assert length(vaccination_range) == 2

    target_noise = target_scaling * measles_daily_incidence

    mean_noise = map(
        coverage -> calculate_mean_dynamical_noise(
            R0,
            latent_period,
            duration_infection,
            correlation,
            poisson_component,
            coverage,
            ensemble_specification,
        ),
        vaccination_range,
    )

    distance = map(noise -> noise - target_noise, mean_noise)
    distance_order_indices = sortperm(distance)
    vaccination_range = vaccination_range[distance_order_indices]
    mean_noise = mean_noise[distance_order_indices]
    sort!(distance)
    abs_distance = abs.(distance)

    # if doesn't have a negative and positive distance end and warning
    # for each step, calculate bisect coverage, mean noise distance, and if
    # absolute distance less than that of the negative/positve bounds (depending on
    # the sign of the bisect distance)

    i = 1
    while i <= maxiters
        @assert distance[1] < 0 && distance[2] > 0

        if showprogress
            println(
                "Run $i: vaccination ranges = $(vaccination_range), noise range = $mean_noise distance range = $(distance)"
            )
        end

        bisect_vaccination = round(
            (vaccination_range[1] + vaccination_range[2]) / 2; digits = 4
        )

        new_noise = calculate_mean_dynamical_noise(
            R0,
            latent_period,
            duration_infection,
            correlation,
            poisson_component,
            bisect_vaccination,
            ensemble_specification,
        )
        new_distance = new_noise - target_noise

        if abs(new_distance) <= atol
            println("Converged in $i iterations")
            return bisect_vaccination, new_noise
        end

        if new_distance < 0
            # @assert abs(new_distance) < abs_distance[1]

            distance[1] = new_distance
            abs_distance[1] = abs(new_distance)
            vaccination_range[1] = bisect_vaccination
            mean_noise[1] = new_noise
        end

        if new_distance > 0
            # @assert abs(new_distance) < abs_distance[2]

            distance[2] = new_distance
            abs_distance[2] = abs(new_distance)
            vaccination_range[2] = bisect_vaccination
            mean_noise[2] = new_noise
        end

        if showprogress
            println(
                "bisect vaccination = $(bisect_vaccination), bisect noise = $new_noise, bisect distance = $new_distance\n"
            )
        end
        i += 1
    end

    min_abs_distance_index = partialsortperm(abs_distance, 1)
    min_distance = distance[min_abs_distance_index]
    min_vaccination = vaccination_range[min_abs_distance_index]
    min_noise = mean_noise[min_abs_distance_index]

    error(
        "Warning: Minimum distance not met\nVaccination that produced the closest value: $min_vaccination\n Associated distance and daily noise: $min_distance $min_noise"
    )
    return nothing
end

function calculate_mean_dynamical_noise(
    R0,
    latent_period,
    duration_infection,
    correlation,
    poisson_component,
    vaccination_coverage,
    ensemble_specification,
)
    return create_noise_arr(
        DynamicalNoiseSpecification(
            "dynamical",
            R0,
            latent_period,
            duration_infection,
            correlation,
            poisson_component,
            vaccination_coverage,
        ),
        nothing;
        ensemble_specification = ensemble_specification,
    )[2].mean_noise
end

dynamical_noise_spec = (;
    R0 = 5.0,
    latent_period = 7,
    duration_infection = 14,
    correlation = "in-phase",
    poisson_component = 0.15,
)

#%%
ensemble_inc_arr = get_ensemble_file(
    ensemble_specification, OutbreakSpecification(5, 30, 500)
)["ensemble_inc_arr"]

mean_measles = StatsBase.mean(ensemble_inc_arr[:, 1, :])

#%%
dynamical_noise_coverages = map(
    ((scale, coverage_range),) -> calculate_dynamic_vaccination_coverage(
        scale,
        mean_measles,
        dynamical_noise_spec,
        ensemble_specification;
        maxiters = 20,
        vaccination_range = coverage_range,
        atol = 0.01,
        showprogress = true,
    ),
    zip(
        [1, 2, 4, 6, 8],
        [[0.8, 0.9], [0.7, 0.8], [0.4, 0.6], [0.2, 0.4], [0.0, 0.1]],
    ),
)

#%%
map(
    ((t, c),) -> t[2] - c * mean_measles,
    zip(dynamical_noise_coverages, [1, 2, 4, 6, 8]),
)
