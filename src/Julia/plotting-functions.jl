function create_sir_plot(
        sol_df;
        colors = ["dodgerblue4", "firebrick3", "chocolate2", "purple"]
    )
    sir_plot = data(sol_df) *
        mapping(
            :time => "Time (days)", :Number,
            color = :State => sorter("S", "I", "R", "N")
        ) *
        visual(Lines, linewidth = 4)

    draw(
        sir_plot;
        palettes = (; color = colors)
        )
end