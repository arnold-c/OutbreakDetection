using GLMakie
using AlgebraOfGraphics

set_aog_theme!()
update_theme!(; size = (1300, 800)) # Laptop monitor
# update_theme!(; size = (2100, 1200)) # Home monitor
# update_theme!(; size = (950, 550)) # Work monitor
# CairoMakie.activate!(; type = "png", px_per_unit = 1.0)
GLMakie.activate!(; float = true)
