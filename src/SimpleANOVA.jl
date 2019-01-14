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

#=
function anovasubjectskernel()
    N = ncells * nreplicates
    nfactors = nnestedfactors + ncrossedfactors

    # collapse replicate dimension
    cellsums = eltype(observations) <: Number && nreplicates == 1 ? observations : sumfirstdim(observations)
    C = sum(cellsums) ^ 2 / N
    total = totalcalc(observations, N, C)

    # for moment, assuming first dimension is subject
    subjectsss = sum(sum(cellsums, dims = 2) .^ 2) / 3 -  C # 3 is number of elements for the subject, dims = 2 is all but the subject dimension

    #factors and interactions calculated as normal

    withinsubjectsss = totalss - subjectsss

    withinsubjectsdf = nsubjects * nfactoroutside * factorinsidedf
    subjectswithinfactorsss = subjectss - factorsss  #guess
    subjectswithinfactorsdf = subjectdf - factorsdf

    withinsubjectinteractionsss = withinsubjectsss - withinsubjectsfactorss - interactionss
    withinsubjectinteractionsdf = withisubjectsdf - withinsubjectsfactordf - interactiondf


    # test for B; B / withinsubjectinteractions
    #test for A; A / subjects within factor A
    # test for interaction; interaction / withinsubjectinteraction
end
=#

end
