function create_sir_df(sol)
    @chain Tables.table(sol') begin
        DataFrame(_)
        rename!(:Column1 => :S, :Column2 => :I, :Column3 => :R)
        @rtransform! :N = :S + :I + :R
        hcat(_, DataFrame(; time=sol.t))
        stack(_, [:S, :I, :R, :N]; variable_name=:State, value_name=:Number)
    end
end

function create_sir_sim_array!(; jump_sol)
    sir_array[1:3, :] = Array(jump_sol)
    sir_array[4, :] = sum(sir_array[1:3, :]; dims=1)
    sir_array[5, :] .= jump_sol.t

    return nothing
end

function create_sir_all_sims_array!(; nsims, prob, alg, δt)
    for i in 1:nsims
        jump_sol = solve(prob, alg; saveat=δt)
        create_sir_sim_array!(; jump_sol=jump_sol)

        all_sims_array[:, :, i] = sir_array
    end
end

function create_sir_all_sims_array!(ensemble_sol::EnsembleSolution, nsims::Int)
    for i in 1:nsims
        all_sims_array[1:3, :, i] = Array(ensemble_sol.u[i])
    end

    all_sims_array[4, :, :] = sum(all_sims_array[1:3, :, :], dims = 1)
end

function create_sir_all_sim_quantiles!(; quantiles)
    for time in 1:size(all_sims_array, 2)
        # sim_means[:, time] = mean(all_sims_array[1:4, time, :], dims = 2)

        for state in 1:4
            sim_quantiles[:, time, state] = quantile(
                skipmissing(all_sims_array[state, time, :]), quantiles
            )
        end
    end
end

