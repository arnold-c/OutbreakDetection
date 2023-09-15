using ProgressMeter
using FLoops

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

        run_dynamics_parameters = DynamicsParameters(
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
            run_dynamics_parameters,
            time_parameters;
            type = "det",
        )
    end
    return nothing
end

function bifurcation_summary(
    state_arr,
    variable_vec,
    years;
    state_labels = seir_state_labels,
)
    annual_summary = zeros(
        Float64, length(years), length(state_labels), length(variable_vec)
    )

    @floop for birth_rate in eachindex(variable_vec),
        state in eachindex(state_labels),
        (year, day) in pairs(years)

        annual_summary[year, state, birth_rate] = maximum(
            state_arr[
                day:(day + 364), state, birth_rate
            ]
        )
    end

    return annual_summary
end

function beta_force_bifurcation_simulation!(
    beta_force_seir_arr,
    beta_force_change_arr,
    beta_force_jump_arr,
    beta_arr,
    initial_states,
    rates,
    beta_force_vec,
    dynamics_parameters,
    time_parameters;
    birth_rate = 50,)
    mu = calculate_mu(birth_rate)
    epsilon = calculate_import_rate(
        mu, dynamics_parameters.R_0, initial_states.N
    )

    for (k, beta_force_run) in pairs(beta_force_vec)
        run_dynamics_parameters = DynamicsParameters(
            dynamics_parameters.beta_mean,
            beta_force_run,
            dynamics_parameters.sigma,
            dynamics_parameters.gamma,
            mu,
            birth_rate,
            epsilon,
            dynamics_parameters.R_0,
        )

        seir_mod!(
            @view(beta_force_seir_arr[:, :, k]),
            @view(beta_force_change_arr[:, :, k]),
            @view(beta_force_jump_arr[:, :, k]),
            beta_arr,
            initial_states,
            rates,
            run_dynamics_parameters,
            time_parameters;
            type = "det",
        )
    end
    return nothing
end
