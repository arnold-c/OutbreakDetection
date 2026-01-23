using OutbreakDetection
using Test

@testset "OutbreakDetection.jl" begin
    include("./test_sort_test_specifications.jl")
    include("./test_plotting_helpers.jl")
end
