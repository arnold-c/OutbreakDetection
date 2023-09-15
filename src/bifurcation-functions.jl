using ProgressMeter

function birth_rate_bifurcation_simulation!(
    birth_rate_seir_arr,
    birth_rate_change_arr,
    birth_rate_jump_arr,
    beta_arr,
    initial_states,
    rates,
    birth_rate_vec,
    dynamics_parameters,
    time_parameters,
)
    for (k, birth_rate_run) in pairs(birth_rate_vec)
        mu_run = calculate_mu(birth_rate_run)

        epsilon_run = calculate_import_rate(
            mu_run, dynamics_parameters.R_0, initial_states.N
        )

        dynamics_parameters = DynamicsParameters(
            dynamics_parameters.beta_mean,
            dynamics_parameters.beta_force,
            dynamics_parameters.sigma,
            dynamics_parameters.gamma,
            mu_run,
            birth_rate_run,
            epsilon_run,
            dynamics_parameters.R_0,
        )

        seir_mod!(
            @view(birth_rate_seir_arr[:, :, k]),
            @view(birth_rate_change_arr[:, :, k]),
            @view(birth_rate_jump_arr[:, :, k]),
            beta_arr,
            initial_states,
            rates,
            dynamics_parameters,
            time_parameters;
            type = "det",
        )
    end
    return nothing
end

function birth_rate_birucation_summary(
    birth_rate_state_arr,
    birth_rate_vec,
    years;
    state_labels = seir_state_labels,
)
    annual_summary = zeros(
        Float64, length(years), length(state_labels), length(birth_rate_vec)
    )

    for birth_rate in eachindex(birth_rate_vec),
        state in eachindex(state_labels),
        (year, day) in pairs(years)

        annual_summary[year, state, birth_rate] = maximum(
            birth_rate_state_arr[
                day:(day + 364), state, birth_rate
            ]
        )
    end

    return annual_summary
end
