using DrWatson
@quickactivate "OutbreakDetection"

using DataFrames, DataFramesMeta
using ModelingToolkit, DifferentialEquations
using FLoops

function create_sir_df(sir_array::Matrix, trange, states = [:S, :I, :R, :N])
    if size(sir_array, 1) == length(states)
        sir_array = sir_array'
    end
    return create_sir_df_inner(
        hcat(trange, DataFrame(Tables.table(sir_array))), states
    )
end

function create_sir_df_inner(sir_df::DataFrame, states)
    @chain sir_df begin
        rename!([:time, states...])
        if :N in states
            stack(_, [states...]; variable_name = :State, value_name = :Number)
        else
            transform!(_, states => (+) => :N)
            stack(
                _, [states..., :N]; variable_name = :State, value_name = :Number
            )
        end
    end
end

function create_sir_df(sol, states)
    return create_sir_df_inner(DataFrame(sol), states)
end

function create_sir_df(sol::RODESolution, states = [:S, :I, :R])
    @chain DataFrame(sol) begin
        if "value1" in names(_)
            rename!(_, [:time, states...])
        else
            rename!(s -> replace(s, "(t)" => ""), _)
            rename!(:timestamp => :time)
        end
        if :N in states
            stack(_, [states...]; variable_name = :State, value_name = :Number)
        else
            transform!(_, states => (+) => :N)
            stack(
                _, [states..., :N]; variable_name = :State, value_name = :Number
            )
        end
    end
end

function create_sir_df(sol::ODESolution, states = [:S, :I, :R])
    @chain DataFrame(sol) begin
        if "value1" in names(DataFrame(_))
            rename!(_, [:time, states...])
        else
            rename!(s -> replace(s, "(t)" => ""), _)
            rename!(:timestamp => :time)
        end
        if :N in states
            stack(_, [states...]; variable_name = :State, value_name = :Number)
        else
            transform!(_, states => (+) => :N)
            stack(
                _, [states..., :N]; variable_name = :State, value_name = :Number
            )
        end
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

function create_sir_all_sims_array_multithread!(prob, nsims, alg, saveat)
    @floop for i in 1:nsims
        all_sims_array[1:3, :, i] = Array(solve(prob, alg; saveat = saveat))
        all_sims_array[4, :, i] = sum(all_sims_array[1:3, :, i]; dims = 1)
    end
end

function create_sir_all_sims_array!(
    ensemble_sol::EnsembleSolution, nsims::Int, array = all_sims_array
)
    @floop for i in 1:nsims
        if size(ensemble_sol.u[i], 2) != size(array, 2)
            skip
        else
            array[1:4, :, i] = Array(ensemble_sol.u[i])
        end
    end

    return nothing
end

function create_sir_all_sims_array(ensemble_sol::EnsembleSolution, nsims::Int)
    all_sims_array = zeros(
        size(ensemble_sol.u[1], 1), size(ensemble_sol.u[1], 2), nsims
    )

    create_sir_all_sims_array!(ensemble_sol, nsims, all_sims_array)

    return all_sims_array
end

function create_sir_all_sims_array!(
    ensemble_sol::EnsembleSolution, nsims::Int, β::Bool
)
    if β == false
        return create_sir_all_sims_array!(ensemble_sol, nsims)
    end

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

function create_sir_all_sim_quantiles!(all_sims_array, sim_quantiles; quantiles)
    @floop for state in 1:size(all_sims_array, 1)
        for time in 1:size(all_sims_array, 2)
            sim_quantiles[:, time, state] = quantile(
                skipmissing(
                    replace(all_sims_array[state, time, :], NaN => missing)
                ),
                quantiles,
            )
        end
    end
end

function create_sir_all_sim_quantiles(all_sims_array; quantiles)
    quantile_array = zeros(
        length(quantiles), size(all_sims_array, 2), size(all_sims_array, 1)
    )

    create_sir_all_sim_quantiles!(
        all_sims_array, quantile_array; quantiles = quantiles
    )

    return quantile_array
end
