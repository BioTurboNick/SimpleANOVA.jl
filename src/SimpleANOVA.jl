module SimpleANOVA

using Distributions
include("InvertedIndices.jl")
include("AnovaEffect.jl")
include("AnovaValue.jl")
include("AnovaFactor.jl")
include("AnovaResult.jl")
include("AnovaData.jl")
include("FactorType.jl")

const totalname = "Total"
const cellsname = "Cells"
const errorname = "Error"
const remaindername = "Remainder"

"""
    anova(observations, [factortypes, factornames])

`observations` - Matrix containing the values. Each dimension is a factor level, such that observations[2,5,3] indicates
the 2nd level of the first factor, the 5th level of the second factor, and the 3rd level of the third factor. May
contain values or vectors of values, where the vector contains replicates.

`factortypes` - Vector indicating the `FactorType` for each factor. If present, `replicates` must appear first, any
`nested` after, and then `random` or `fixed` in any order. Specify `replicates` if the first dimension of the
`observations` matrix contains replicate values (vs. contained in vectors). If too few values are provided, remaining
are assumed to be `fixed`.

`factornames` - Vector of names for each factor, excluding the replicate factor. If empty, will be automatically
populated.

Notes: Requires balanced data. The last index will be your "topmost" factor.

Output: `AnovaData` structure containing the test results for each factor.

Examples:
N-way fixed-effects ANOVA with replicates in vectors: anova(observations)

N-way fixed-effects ANOVA with replicates in first dimension: anova(observations)

N-way fixed-effects ANOVA without replicates in first dimension: anova(observations, replicates = false)

2-way ANOVA with Factor A random and Factor B fixed with replicates in vectors: anova(observations, [random])

2-way ANOVA with Factor A fixed and Factor B random with replicates in vectors: anova(observations, [fixed, random])

2-way fixed-effects ANOVA with 2 random nested factors with replicates in first dimension:
anova(observations, [nested, nested])
"""
function anova(observations::AbstractArray{T}, factortypes::Vector{FactorType} = FactorType[]; factornames::Vector{<:AbstractString} = String[], hasreplicates = true, withinsubjects = false) where {T <: Union{Number, AbstractVector{<:Number}}}
    length(observations) > 0 || return

    firstlevelreplicates = eltype(observations) <: Number ? hasreplicates : false
    nfactors = ndims(observations) - (firstlevelreplicates ? 1 : 0)

    # defaults to assuming all unspecified factors are fixed
    if length(factortypes) < nfactors
        nremaining = nfactors - length(factortypes)
        append!(factortypes, repeat([fixed], nremaining))
    end

    validate(factortypes, factornames, nfactors)

    # automatically assigns alphabetical names if not provided.
    if isempty(factornames)
        factornames = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z"][1:nfactors]
        reverse!(factornames)
    end

    nnestedfactors = count(f -> f == nested, factortypes)
    crossedfactortypes = filter(f -> f ∈ [fixed, random], factortypes)
    reverse!(crossedfactortypes)
    ncrossedfactors = length(crossedfactortypes)
    ncrossedfactors < 4 || error("ANOVA with 4 or more crossed factors is not supported.")
    nreplicates = firstlevelreplicates ? size(observations, 1) : length(observations[1])
    firstlevelreplicates || all(c -> length(c) == nreplicates, observations) || throw(ErrorException("All cells must have the same number of replicates."))
    ncells = Int.(length(observations) / (firstlevelreplicates ? nreplicates : 1))
    nfactorlevels = firstlevelreplicates ? [size(observations)...][Not(1)] : [size(observations)...]

    crossedfactornames = factornames[factortypes .≠ nested]
    nestedfactornames = factornames[factortypes .== nested]

    #10kb allocated before this point
    anovakernel(observations, nreplicates, ncells, nnestedfactors, ncrossedfactors, nfactorlevels, crossedfactortypes, crossedfactornames, nestedfactornames)
end

function validate(factortypes::Vector{FactorType}, factornames::Vector{<:AbstractString}, nfactors)
    if !isempty(factortypes)
        length(factortypes) == nfactors || error("factortypes must have an entry for each factor.")
        nested ∉ factortypes || factortypes[1:count(t -> t == nested, factortypes)] |> unique |> length == 1 || error("nested entries must come before crossed factors")
    end

    if !isempty(factornames)
        nfactors ≤ 26 || error("Can only automatically name up to 26 factors. Provide names explicitly.")
        length(factornames) == nfactors || error("factornames must have an entry for each factor.")
    end
end

function anovakernel(observations, nreplicates, ncells, nnestedfactors, ncrossedfactors, nfactorlevels, crossedfactortypes, crossedfactornames, nestedfactornames)
    N = ncells * nreplicates
    nfactors = nnestedfactors + ncrossedfactors

    # collapse replicate dimension
    cellsums = eltype(observations) <: Number && nreplicates == 1 ? observations : sumfirstdim(observations)
    C = sum(cellsums) ^ 2 / N
    total = totalcalc(observations, N, C)
    amongallnested, nestedsums, ncrossedfactorlevels, nnestedfactorlevels = amongnestedfactorscalc(cellsums, nfactorlevels, nnestedfactors, nreplicates, C)
    cells = cellscalc(cellsums, nreplicates, ncells, C)

    crossedfactors = factorscalc(nestedsums, ncrossedfactors, ncrossedfactorlevels, N, C, crossedfactornames) # 3kb allocated here, possibly can't be avoided
    interactions, interactionsmap = interactionscalc(cells, nestedsums, crossedfactors, ncrossedfactors, ncrossedfactorlevels, nnestedfactorlevels, nreplicates, C, crossedfactornames) # 9 kb allocated here!
    nestedfactors = nestedfactorscalc(amongallnested, nnestedfactors, crossedfactors, interactions, nestedfactornames)
    error = errorcalc(total, amongallnested, cells, [crossedfactors; interactions[1:end-1]], nnestedfactors, nreplicates)

    reverse!(crossedfactors)

    numerators = getnumerators(crossedfactors, ncrossedfactors, nnestedfactors, nestedfactors, interactions)
    crossedbasedenominator = nnestedfactors > 0 ? nestedfactors[end] : error;
    denominators = getdenominators(nnestedfactors, nestedfactors, nreplicates, crossedbasedenominator, error, total, crossedfactors, ncrossedfactors, crossedfactortypes, interactionsmap)# 2-3kb allocated, 50 allocations

    # drop least significant test if nreplicates == 1; either the lowest interaction level, or lowest nesting level if present
    if nreplicates == 1
        droppedfactor = pop!(numerators)
        pop!(denominators)
    end

    # perform test
    results = ftest.(numerators, denominators)

    data = AnovaData([total; results])
    nnestedfactors > 0 && nreplicates == 1 && push!(data.effects, droppedfactor)
    error.df > 0 && push!(data.effects, error)

    return data
end

function sumfirstdim(observations::T) where {T <: AbstractArray{<:AbstractVector}}
    map(sumfirstdim, observations)
end

function sumfirstdim(observations::T) where {T <: AbstractArray{<:Number}}
    dropdims(sum(observations, dims = 1), dims = 1)
end

function sumfirstdim(observations::T) where {T <: AbstractVector{<:Number}}
    sum(observations, dims = 1) |> first
end

function totalcalc(observations, N, C)
    ss = sum(c -> sum(c.^2), observations) - C
    df = N - 1
    AnovaValue(totalname, ss, df)
end

function factorscalc(cellsums, nfactors, nfactorlevels, N, C, factornames)
    factorindices = 1:nfactors
    ss = map(i -> sum(sum(cellsums, dims = factorindices[Not(i)]) .^ 2) / (N / nfactorlevels[i]), factorindices) .- C
    df = nfactorlevels .- 1
    AnovaFactor.(factornames, ss, df)
end

function cellscalc(cellsums, nreplicates, ncells, C)
    ss = sum(cellsums .^ 2) / nreplicates - C
    df = ncells - 1
    AnovaValue(cellsname, ss, df)
end

function errorcalc(total, amongallnested, cells, otherfactors, nnestedfactors, nreplicates)
    if nnestedfactors > 0
        otherfactor = amongallnested[1]
        name = errorname
    elseif nreplicates > 1
        otherfactor = cells
        name = errorname
    else
        otherfactor = AnovaValue("", sum(f -> f.ss, otherfactors), sum(f -> f.df, otherfactors))
        name = remaindername
    end

    ss = total.ss - otherfactor.ss
    df = total.df - otherfactor.df
    AnovaFactor(name, ss, df)
end

function amongnestedfactorscalc(cellsums, nfactorlevels, nnestedfactors, nreplicates, C)
    nestedsums = cellsums
    nlowerfactorlevels = 1
    nupperfactorlevels = nfactorlevels
    amongallnested = Vector{AnovaValue}(undef, nnestedfactors)
    for i ∈ 1:nnestedfactors
        # collapse each nested factor
        ss = sum(nestedsums .^ 2 ./ (nreplicates * prod(nlowerfactorlevels))) - C
        df = prod(nupperfactorlevels) - 1
        amongallnested[i] = AnovaValue(ss, df)
        nlowerfactorlevels = nfactorlevels[1:i]
        nupperfactorlevels = nfactorlevels[(i+1):end]
        nestedsums = sumfirstdim(nestedsums)
    end
    amongallnested, nestedsums, nupperfactorlevels, nlowerfactorlevels
end

function interactionscalc(cells, nestedsums, crossedfactors, ncrossedfactors, ncrossedfactorlevels, nnestedfactorlevels, nreplicates, C, crossedfactorlabels)
    interactionsmap = Dict{Any,AnovaFactor}()
    interactions = Vector{AnovaFactor}()
    if ncrossedfactors > 1
        pairwise = pairwisecalc(nestedsums, crossedfactors, ncrossedfactors, ncrossedfactorlevels, nnestedfactorlevels, nreplicates, C, crossedfactorlabels)
        interactionsmap[(1,2)] = pairwise[1,2]
        interactionsmap[(2,1)] = pairwise[1,2]
        push!(interactions, interactionsmap[(1,2)])
        if ncrossedfactors > 2
            threewise = threewisecalc(cells, crossedfactors, pairwise, crossedfactorlabels)
            interactionsmap[(1,3)] = pairwise[1,3]
            interactionsmap[(3,1)] = pairwise[1,3]
            interactionsmap[(2,3)] = pairwise[2,3]
            interactionsmap[(3,2)] = pairwise[2,3]
            interactionsmap[(1,2,3)] = threewise
            push!(interactions, interactionsmap[(1,3)], interactionsmap[(2,3)])
            reverse!(interactions)
            push!(interactions, interactionsmap[(1,2,3)])
        end
    end
    return interactions, interactionsmap
end

function pairwisecalc(nestedsums, crossedfactors, ncrossedfactors, ncrossedfactorlevels, nnestedfactorlevels, nreplicates, C, crossedfactorlabels)
    pairwise = Array{AnovaFactor,2}(undef, ncrossedfactors, ncrossedfactors)
    factorindexes = 1:ncrossedfactors
    for i ∈ factorindexes
        for j ∈ (i+1):ncrossedfactors
            otherfactorindexes = intersect(factorindexes[Not(i)], factorindexes[Not(j)])
            ss = sum(sum(nestedsums, dims = otherfactorindexes) .^ 2 ./ (prod(ncrossedfactorlevels[otherfactorindexes]) * prod(nnestedfactorlevels) * nreplicates)) - C - crossedfactors[i].ss - crossedfactors[j].ss
            df = crossedfactors[i].df * crossedfactors[j].df
            pairwise[i,j] = pairwise[j,i] = AnovaFactor("$(crossedfactorlabels[j]) × $(crossedfactorlabels[i])", ss, df)
        end
    end
    return pairwise
end

function threewisecalc(cells, factors, pairwise, crossedfactorlabels)
    ss = cells.ss - sum(f -> f.ss, factors) - pairwise[1,2].ss - pairwise[1,3].ss - pairwise[2,3].ss
    df = prod(f -> f.df, factors)
    AnovaFactor("$(crossedfactorlabels[3]) × $(crossedfactorlabels[2]) × $(crossedfactorlabels[1])", ss, df)
end

function nestedfactorscalc(amongallnested, nnestedfactors, crossedfactors, interactions, nestedfactorlabels)
    nestedfactors = []
    if nnestedfactors > 0
        otherfactors = [crossedfactors; interactions]
        otherss = sum(f -> f.ss, otherfactors)
        otherdf = sum(f -> f.df, otherfactors)
        for i ∈ nnestedfactors:-1:1
            ss = amongallnested[i].ss - otherss
            df = amongallnested[i].df - otherdf
            otherss += ss
            otherdf += df
            push!(nestedfactors, AnovaFactor(nestedfactorlabels[i], ss, df))
        end
    end
    return nestedfactors
end

function getnumerators(crossedfactors, ncrossedfactors, nnestedfactors, nestedfactors, interactions)
    numerators = copy(crossedfactors)

    if ncrossedfactors > 1
        push!(numerators, interactions[1])
        if ncrossedfactors > 2
            push!(numerators, interactions[2], interactions[3], interactions[4])
        end
    end

    if nnestedfactors > 0
        append!(numerators, nestedfactors)
    end

    return numerators
end

function getdenominators(nnestedfactors, nestedfactors, nreplicates, crossedbasedenominator, error, total, crossedfactors, ncrossedfactors, crossedfactortypes, interactionsmap)
    if ncrossedfactors == 1
        denominators = [crossedbasedenominator]
    elseif ncrossedfactors == 2
        if all(f -> f == fixed, crossedfactortypes)
            crosseddenominators = repeat([crossedbasedenominator], ncrossedfactors)
        elseif all(f -> f == random, crossedfactortypes)
            crosseddenominators = repeat([interactionsmap[(1,2)]], ncrossedfactors)
        else
            crosseddenominators = map(f -> f == fixed ? interactionsmap[(1,2)] : crossedbasedenominator, crossedfactortypes)
        end

        denominators = [crosseddenominators; crossedbasedenominator]
    elseif ncrossedfactors == 3
        if all(f -> f == fixed, crossedfactortypes)
            crosseddenominators = repeat([crossedbasedenominator], ncrossedfactors)
            pairwiseinteractiondenominators = repeat([crossedbasedenominator], ncrossedfactors)
        elseif all(f -> f == random, crossedfactortypes)
            crosseddenominators = Vector{AnovaFactor}(undef, ncrossedfactors)
            for i ∈ 1:ncrossedfactors
                otherfactors = (1:ncrossedfactors)[Not(i)]
                j = otherfactors[1]
                k = otherfactors[2]
                crosseddenominators[i] = threewayinteraction(interactionsmap[(i,j)], interactionsmap[(i,k)], interactionsmap[(1,2,3)])
            end
            pairwiseinteractiondenominators = repeat([interactionsmap[(1,2,3)]], ncrossedfactors)
        elseif count(f -> f == random, crossedfactortypes) == 1
            i = findfirst(f -> f == random, crossedfactortypes)

            crosseddenominators = Vector{AnovaFactor}(undef, ncrossedfactors)
            crosseddenominators[i] = crossedbasedenominator

            fixedindexes = (1:ncrossedfactors)[Not(i)]
            for j ∈ fixedindexes
                crosseddenominators[j] = interactionsmap[(i,j)]
            end

            fixedinteractionindex = sum(fixedindexes) - 2
            pairwiseinteractiondenominators = Vector{AnovaFactor}(undef, ncrossedfactors)
            pairwiseinteractiondenominators[fixedinteractionindex] = interactionsmap[(1,2,3)]
            pairwiseinteractiondenominators[Not(fixedinteractionindex)] .= crossedbasedenominator
        elseif count(f -> f == random, crossedfactortypes) == 2
            i = findfirst(f -> f == fixed, crossedfactortypes)
            otherfactors = (1:ncrossedfactors)[Not(i)]
            j = otherfactors[1]
            k = otherfactors[2]

            crosseddenominators = Vector{AnovaFactor}(undef, ncrossedfactors)
            crosseddenominators[i] = threewayinteraction(interactionsmap[(i,j)], interactionsmap[(i,k)], interactionsmap[(1,2,3)])
            crosseddenominators[otherfactors] .= interactionsmap[(j,k)]

            randominteractionindex = sum(otherfactors) - 2
            pairwiseinteractiondenominators = Vector{AnovaFactor}(undef, ncrossedfactors)
            pairwiseinteractiondenominators[randominteractionindex] = crossedbasedenominator
            pairwiseinteractiondenominators[Not(randominteractionindex)] .= interactionsmap[(1,2,3)]
        end

        denominators = [crosseddenominators; pairwiseinteractiondenominators; crossedbasedenominator]
    end

    # determine correct denominators for nested factors
    if nnestedfactors > 0
        nesteddenominators = Vector{AnovaFactor}(undef, nnestedfactors)
        nesteddenominators[1:end-1] = nestedfactors[2:end]
        nesteddenominators[end] = error
        append!(denominators, nesteddenominators)
    end

    denominators
end

function threewayinteraction(interaction_ab, interaction_bc, interaction_abc)
    reducedmeansquare(factor::AnovaFactor) = factor.ms ^ 2 / factor.df

    ms = interaction_ab.ms + interaction_bc.ms - interaction_abc.ms
    df = ms ^ 2 / (reducedmeansquare(interaction_ab) + reducedmeansquare(interaction_bc) + reducedmeansquare(interaction_abc))
    AnovaFactor("", ms * df, df, ms)
end

function ftest(x, y)
    f = x.ms / y.ms
    fdist = FDist(x.df, y.df)
    p = ccdf(fdist, f)
    AnovaResult(x, f, p)
end

export anova, ftest

end
