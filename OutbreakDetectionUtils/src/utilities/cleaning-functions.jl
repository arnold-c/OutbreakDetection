using DataFrames: DataFrames
using DataFramesMeta: DataFramesMeta
using Chain: Chain
using Statistics: Statistics
using Tables: Tables

export create_sir_df, create_sir_beta_dfs,
    create_sir_all_sim_quantiles, create_sir_all_sim_quantiles!

function create_sir_df(sir_array::Matrix, trange, states = [:S, :I, :R, :N])
    if size(sir_array, 1) == length(states)
        sir_array = sir_array'
    end
    return create_sir_df_inner(
        hcat(trange, DataFrames.DataFrame(Tables.table(sir_array))), states
    )
end

function create_sir_df_inner(sir_df::DataFrames.DataFrame, states)
    return Chain.@chain sir_df begin
        DataFrames.rename!([:time, states...])
        if :N in states
            DataFrames.stack(
                _, [states...]; variable_name = :State, value_name = :Number
            )
        else
            DataFrames.transform!(_, states => (+) => :N)
            DataFrames.stack(
                _, [states..., :N]; variable_name = :State, value_name = :Number
            )
        end
    end
end

function create_sir_df(sol, states)
    return create_sir_df_inner(DataFrames.DataFrame(sol), states)
end

function create_sir_beta_dfs(sol, states = [:S, :I, :R])
    state_df = create_sir_df(sol, states)

    beta_df = DataFrames.select(state_df, [:time, :beta])
    unique!(beta_df)

    DataFrames.select!(state_df, DataFrames.Not(:beta))

    return state_df, beta_df
end

function create_sir_all_sim_quantiles!(all_sims_array, sim_quantiles; quantiles)
    for state in axes(all_sims_array, 2), time in axes(all_sims_array, 1)
        sim_quantiles[:, time, state] = Statistics.quantile(
            skipmissing(
                replace(@view(all_sims_array[time, state, :]), NaN => missing)
            ),
            quantiles,
        )
    end
    return
end

function create_sir_all_sim_quantiles(all_sims_array; quantiles)
    quantile_array = zeros(
        length(quantiles), size(all_sims_array, 1), size(all_sims_array, 2)
    )

    create_sir_all_sim_quantiles!(
        all_sims_array, quantile_array; quantiles = quantiles
    )

    return quantile_array
end
