function create_sir_df(sol, states = [:S, :I, :R])
    @chain DataFrame(sol) begin
        rename!([:time, states...])
        transform!(_, states => (+) => :N)
        stack(_, [states..., :N]; variable_name = :State, value_name = :Number)
    end
end

function create_sir_df(sol::RODESolution, states = [:S, :I, :R])
    @chain DataFrame(sol) begin
        rename!(s -> replace(s, "(t)" => ""), _)
        rename!(:timestamp => :time)
        transform!(_, states => (+) => :N)
        stack(_, [states..., :N]; variable_name = :State, value_name = :Number)
    end
end

function create_sir_df(sol::ODESolution, states = [:S, :I, :R])
    @chain DataFrame(sol) begin
        rename!(s -> replace(s, "(t)" => ""), _)
        rename!(:timestamp => :time)
        transform!(_, states => (+) => :N)
        stack(_, [states..., :N]; variable_name = :State, value_name = :Number)
    end
end

function create_sir_beta_dfs(sol, states = [:S, :I, :R])
    state_df = create_sir_df(sol, states)
    
    beta_df = select(state_df, [:time, :β])
    unique!(beta_df)

    select!(state_df, Not(:β))

    return state_df, beta_df
end


function create_sir_sim_array!(; jump_sol)
    sir_array[1:3, :] = Array(jump_sol)
    sir_array[4, :] = sum(sir_array[1:3, :]; dims = 1)

    return nothing
end

function create_sir_all_sims_array!(; nsims, prob, alg, δt)
    for i in 1:nsims
        jump_sol = solve(prob, alg; saveat = δt)
        create_sir_sim_array!(; jump_sol = jump_sol)

        all_sims_array[:, :, i] = sir_array
    end
end

function create_sir_all_sims_array!(ensemble_sol::EnsembleSolution, nsims::Int)
    for i in 1:nsims
        if size(ensemble_sol.u[i], 2) != size(all_sims_array, 2)
            skip
        else
            all_sims_array[1:3, :, i] = Array(ensemble_sol.u[i])
        end
    end

    all_sims_array[4, :, :] = sum(all_sims_array[1:3, :, :]; dims = 1)

    return nothing
end

function create_sir_all_sims_array!(ensemble_sol::EnsembleSolution, nsims::Int; β = true)
    for i in 1:nsims
        if size(ensemble_sol.u[i], 2) != size(all_sims_array, 2)
            skip
        else
            all_sims_array[1:3, :, i] = Array(ensemble_sol.u[i])[1:3, :]
            all_sims_array[5, :, i] = Array(ensemble_sol.u[i])[4, :]
        end
    end

    all_sims_array[4, :, :] = sum(all_sims_array[1:3, :, :]; dims = 1)

    return nothing
end

function create_sir_all_sim_quantiles!(; quantiles)
    for time in 1:size(all_sims_array, 2)
        for state in 1:size(all_sims_array, 1)
            sim_quantiles[:, time, state] = quantile(
                skipmissing(replace(all_sims_array[state, time, :], NaN => missing)),
                quantiles,
            )
        end
    end
end