function create_sir_df(sol)
    @chain Tables.table(sol') begin
        DataFrame(_)
        rename!(:Column1 => :S, :Column2 => :I, :Column3 => :R)
        @rtransform! :N = :S + :I + :R
        hcat(_, DataFrame(time = sol.t))
        stack(_, [:S, :I, :R, :N], variable_name = :State, value_name = :Number)
    end
end