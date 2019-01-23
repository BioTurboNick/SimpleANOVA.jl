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

function rmanovakernel(observations, ncells, nreplicates)
    N = ncells * nreplicates

    cellsums = eltype(observations) <: Number && nreplicates == 1 ? observations : sumfirstdim(observations)
    C = sum(cellsums) ^ 2 / N
    total = totalcalc(observations, N, C)
    amongallnested, nestedsums, ncrossedfactorlevels, nnestedfactorlevels = amongnestedfactorscalc(cellsums, nfactorlevels, nnestedfactors, nreplicates, C)
    cells = cellscalc(cellsums, nreplicates, ncells, C)

    subjectsss =
    subjectsdf = prod(nsubjectsfactorlevels) * nreplicates - 1
    subjects = AnovaValue(subjectsss, subjectsdf)

    # assumes all among factors and within factors are fixed, subjects is random

    # assumes that amongfactors are topmost
    if namongfactors == 1
        amongfactor = factorscalc(nestedsums, namongfactors, namongfactorlevels, N, C, amongfactornames) |> first
        subjectswithin_amongfactor_ss = subjects.ss - amongfactor.ss
        subjectswithin_amongfactor_df = subjects.df - amongfactor.df
        subjectswithin_amongfactor = AnovaFactor(subjectswithin_amongfactor_ss, subjectswithin_amongfactor_df)
    elseif namongfactors == 2
        amongfactors = factorscalc(nestedsums, namongfactors, namongfactorlevels, N, C, amongfactornames)
        amonginteractions, amonginteractionsmap = interactionscalc(cells, nestedsums, amongfactors, namongfactors, namongfactorlevels, nnestedfactorlevels, nreplicates, C, amongfactornames) # 9 kb allocated here!
        subjectswithin_amongfactor_ss = subjects.ss - sum(f -> f.ss, amongfactors) - amonginteractionsmap[(1,2)].ss
        subjectswithin_amongfactor_df = subjects.df - sum(f -> f.df, amongfactors) - amonginteractionsmap[(1,2)].df
        subjectswithin_amongfactor = AnovaFactor(subjectswithin_amongfactor_ss, subjectswithin_amongfactor_df)
    end

    withinsubjectsss = total.ss - subjects.ss
    withinsubjectsdf = total.df - subjects.df
    withinsubjects = AnovaValue(withinsubjectsss, withinsubjectsdf)

    # needs to adjust for not being topmost factors
    if nwithinfactors == 1
        withinfactor = factorscalc(nestedsums, nwithinfactors, nwithinfactorlevels, N, C, withinfactornames)
        factor_subjectswithin_interaction_df = subjectswithin_amongfactor_df * withinfactor.df
        factor_subjectswithin_interaction_ss = ?????????
    elseif nwithinfactors == 2
        withinfactors = factorscalc(nestedsums, nwithinfactors, nwithinfactorlevels, N, C, withinfactornames)
        withininteractions, withininteractionsmap = interactionscalc(cells, nestedsums, withinfactors, nwithinfactors, nwithinfactorlevels, nnestedfactorlevels, nreplicates, C, withinfactornames) # 9 kb allocated here!
        factor_subjectswithin_interaction_df = subjectswithin_amongfactor_df .* [f.df for f in withinfactors]
        factor_subjectswithin_interaction_ss = ?????????
    end

    if namongfactors + nwithinfactors == 3
        interaction_subjectswithin_df = subjectswithin_amongfactor.df .* prod(f -> f.df, withinfactors)
        interaction_subjectswithin_ss = ?????????
    end

    amongsubjectsdenominator = subjectswithin_amongfactor
    withinsubjectsdenominators .= pairinteractionsubjectswithinms, threeinteractionsubjectswithinms
end


end
