export create_combinations_vec, create_ensemble_spec_combinations

"""
    create_combinations_vec(custom_function, combinations)

Create a vector of all combinations by applying custom_function to each combination.
"""
function create_combinations_vec(custom_function, combinations)
    combs = Iterators.product(combinations...)

    return vec(map(combination -> custom_function(combination...), combs))
end

"""
    create_ensemble_spec_combinations(beta_force_vec, seasonality_vec, sigma_vec, 
                                     gamma_vec, annual_births_per_k_vec, R_0_vec, 
                                     vaccination_coverage_vec, N_vec, 
                                     init_states_prop_dict, model_types_vec, 
                                     time_p_vec, nsims_vec)

Create all combinations of ensemble specifications from parameter vectors.
"""
function create_ensemble_spec_combinations(
        beta_force_vec,
        seasonality_vec,
        sigma_vec,
        gamma_vec,
        annual_births_per_k_vec,
        R_0_vec,
        vaccination_coverage_vec,
        N_vec,
        init_states_prop_dict,
        model_types_vec,
        time_p_vec,
        nsims_vec,
    )
    ensemble_spec_combinations = Iterators.product(
        beta_force_vec,
        seasonality_vec,
        sigma_vec,
        gamma_vec,
        annual_births_per_k_vec,
        R_0_vec,
        vaccination_coverage_vec,
        N_vec,
        init_states_prop_dict,
        model_types_vec,
        time_p_vec,
        nsims_vec,
    )

    ensemble_spec_vec = Vector(undef, length(ensemble_spec_combinations))

    for (
            i,
            (
                beta_force,
                seasonality,
                sigma,
                gamma,
                annual_births_per_k,
                R_0,
                vaccination_coverage,
                N,
                init_states_prop,
                model_type,
                time_p,
                nsims,
            ),
        ) in enumerate(ensemble_spec_combinations)
        mu = calculate_mu(annual_births_per_k)
        beta_mean = calculate_beta(R_0, gamma, mu, 1, N)
        epsilon = calculate_import_rate(mu, R_0, N)

        ensemble_spec_vec[i] = EnsembleSpecification(
            model_type,
            StateParameters(
                N, init_states_prop
            ),
            DynamicsParameters(
                beta_mean,
                beta_force,
                seasonality,
                sigma,
                gamma,
                mu,
                annual_births_per_k,
                epsilon,
                R_0,
                vaccination_coverage,
            ),
            time_p,
            nsims,
        )
    end

    return ensemble_spec_vec
end
