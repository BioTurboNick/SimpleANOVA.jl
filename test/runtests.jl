using CategoricalArrays
using DataFrames
using Documenter
using SimpleANOVA
using Test

@testset "SimpleANOVA Tests" begin
    DocMeta.setdocmeta!(
        SimpleANOVA,
        :DocTestSetup,
        :(using SimpleANOVA);
        recursive=true
    )

    # Ensures that jldoctests also pass.
    doctest(SimpleANOVA)

    include("test_anova.jl")
    include("test_pretests.jl")
    include("test_contrasts.jl")
    include("test_show.jl")
end
