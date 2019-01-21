using SimpleANOVA
using Test
using DataFrames

@testset "SimpleANOVA Tests" begin
include("test_anova.jl")
include("test_show.jl")
end
