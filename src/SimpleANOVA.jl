module SimpleANOVA

#=

Workflow:
results = anova() # get table for all effects and detect interactions
plot(results)     # assess whether interactions are trivial or large

# if interactions are trivial, can immediately do posthoc test of choice on the independent main effects
tukey(results)

# if interactions are nontrivial, do subanovas with one less factor, at each level of the removed factor, with α/k for k levels
subanova(results)

# finally, do pairwise (e.g. tukey) tests among levels of one factor
tukey(results, [1,2])
=#


using Distributions, Requires, InvertedIndices

include("data/AnovaEffect.jl")
include("data/AnovaValue.jl")
include("data/AnovaFactor.jl")
include("data/AnovaResult.jl")
include("data/AnovaData.jl")
include("data/FactorType.jl")
include("anova.jl")
include("pretests.jl")

function __init__()
    @require DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0" include("anova_dataframes.jl")
    @require DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0" include("pretests_dataframes.jl") # can these be combined?
    @require Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80" include("anova_plots.jl")
end


export anova, ftest, plot, levene

end
