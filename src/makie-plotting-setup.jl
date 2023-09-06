using GLMakie
using AlgebraOfGraphics
using ColorSchemes

set_aog_theme!()
# Set depending on size of screen
update_theme!(; resolution = (2200, 1300))
#= update_theme!(; resolution = (850, 600)) =#
GLMakie.activate!(; float = true)
