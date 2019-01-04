module SimpleAnova

using Statistics
using Distributions
include("InvertedIndices.jl")

import Main.InvertedIndices.Not

@enum FactorType fixed random nested replicate

abstract type AnovaEffect
end

struct AnovaValue <: AnovaEffect
    ss::Float64
    df::Float64
end

struct AnovaFactor <: AnovaEffect
    ss::Float64
    df::Float64
    ms::Float64
end

Broadcast.broadcastable(a::AnovaFactor) = (a,) # workaround for current behavior

AnovaFactor(ss, df) = AnovaFactor(ss, df, ss / df)


#=
Index indicates the value of that factor level. E.g. [3,4] specifies Factor A level 3 and Factor B level 4

Examples

1-way ANOVA:
observations = Array{Vector{Float64}, 1}(undef, 4)
observations[1] = [60.8, 57.0, 65.0, 58.6, 61.7]
observations[2] = [68.7, 67.7, 74.9, 66.3, 69.8]
observations[3] = [102.6, 102.1, 100.2, 96.5, 100.4]
observations[4] = [87.9, 84.2, 83.1, 85.7, 90.3]

-or-

observations = [60.8 68.7 102.6 87.9;
                57.0 67.7 102.1 84.2;
                65.0 74.9 100.2 83.1;
                58.6 66.3 96.5 85.7;
                61.7 69.8 100.4 90.3]

Factor
    1      2      3      4
 60.8   68.7  102.6   87.9
 57.0   67.7  102.1   84.2
 65.0   74.9  100.2   83.1
 58.6   66.3   96.5   85.7
 61.7   69.8  100.4   90.3

            SS  DF      MS      F       p
Total   4823.0  19
Factor  4685.0   3  1561.7  181.8   1e-12
Error    137.5  16     8.6

2-way ANOVA without replication:
observations = Array{Vector{Float64}, 2}(undef, 3, 4)
observations[1,1] = [123]
observations[1,2] = [138]
observations[1,3] = [110]
observations[1,4] = [151]
observations[2,1] = [145]
observations[2,2] = [165]
observations[2,3] = [140]
observations[2,4] = [167]
observations[3,1] = [156]
observations[3,2] = [176]
observations[3,3] = [185]
observations[3,4] = [175]

-or-

observations = [123 138 110 151; 145 165 140 167; 156 176 185 175]

                Factor B
                  1   2   3   4
Factor A    1   123 138 110 151
            2   145 165 140 167
            3   156 176 185 175

                SS  DF      MS     F        p
Total       5594.9  11
Factor A    3629.2   2  1814.6  12.8    0.007
Factor B    1116.9   3   372.3   2.6    0.145
Remainder    848.8   6   141.5

2-way ANOVA with replication

Specified in Array with cells nested as vectors
observations = Array{Vector{Float64}, 2}(undef, 2, 2)
observations[1,1] = [16.5, 18.4, 12.7, 14.0, 12.8]
observations[1,2] = [14.5, 11.0, 10.8, 14.3, 10.0]
observations[2,1] = [39.1, 26.2, 21.3, 35.8, 40.2]
observations[2,2] = [32.0, 23.8, 28.8, 25.0, 29.3]

-or-

Specified in multidimensional array with 1st dimension as replicate
observations = cat(hcat([16.5, 18.4, 12.7, 14.0, 12.8], [39.1, 26.2, 21.3, 35.8, 40.2]),
                   hcat([14.5, 11.0, 10.8, 14.3, 10.0], [32.0, 23.8, 28.8, 25.0, 29.3]), dims = 3)
    Note: must specify that first dimension is a replicate

                Factor B
                   1       2
Factor A    1   16.5    14.5
                18.4    11.0
                12.7    10.8
                14.0    14.3
                12.8    10.0

            2   39.1    32.0
                26.2    23.8
                21.3    28.8
                35.8    25.0
                40.2    29.3

                SS  DF      MS     F        p
Total       1827.7  11
Cells       1461.3   3
Factor A    1386.1   2  1386.1  60.5     8e-7
Factor B      70.3   3    70.3   3.1      0.1
Factor AxB      4.9  1     4.9   0.2      0.6
Error        366.4   6   366.4


1-way ANOVA with 1 nested factor, replicates

Specified in Array with cells nested as vectors
observations = Array{Vector{Float64}, 2}(undef, 2, 3)
observations[1,1] = [102, 104]
observations[1,2] = [108, 110]
observations[1,3] = [104, 106]
observations[2,1] = [103, 104]
observations[2,2] = [109, 108]
observations[2,3] = [105, 107]

-or-

Specified in multidimensional array with 1st dimension as replicate
observations = cat(hcat([102, 104], [103, 104]),
                   hcat([108, 110], [109, 108]),
                   hcat([104, 106], [105, 107]), dims = 3)
    Note: must specify that first dimension is a replicate

                        Factor B
                          1   2   3
Nested Factor A     1   102 108 104
                        104 110 106

                    2   103 109 105
                        104 108 107

                      SS    DF    MS     F
Total               71.7    11
Across Factor A     62.7     5
Factor B            61.2     2  30.6
Factor A             1.5     3   0.5
Error                9.0     6   1.5


2-way ANOVA with 1 nested factor, no replicates
Specified in multidimensional array with 1st dimension as replicate
observations = cat(hcat([16.5, 18.4, 12.7, 14.0, 12.8], [39.1, 26.2, 21.3, 35.8, 40.2]),
                   hcat([14.5, 11.0, 10.8, 14.3, 10.0], [32.0, 23.8, 28.8, 25.0, 29.3]), dims = 3)
    Note: must specify that first dimension is a replicate

                Factor B
                   1       2
Factor A    1   16.5    14.5
                18.4    11.0
                12.7    10.8
                14.0    14.3
                12.8    10.0

            2   39.1    32.0
                26.2    23.8
                21.3    28.8
                35.8    25.0
                40.2    29.3


3-way ANOVA
Specified in Array with cells nested as vectors
observations = Array{Vector{Float64}, 3}(undef, 2, 3, 3)
observations[1,1,1] = [1.9, 1.8, 1.6, 1.4]
observations[1,1,2] = [2.1, 2.0, 1.8, 2.2]
observations[1,1,3] = [1.1, 1.2, 1.0, 1.4]
observations[1,2,1] = [2.3, 2.1, 2.0, 2.6]
observations[1,2,2] = [2.4, 2.6, 2.7, 2.3]
observations[1,2,3] = [2.0, 2.1, 1.9, 2.2]
observations[1,3,1] = [2.9, 2.8, 3.4, 3.2]
observations[1,3,2] = [3.6, 3.1, 3.4, 3.2]
observations[1,3,3] = [2.9, 2.8, 3.0, 3.1]
observations[2,1,1] = [1.8, 1.7, 1.4, 1.5]
observations[2,1,2] = [2.3, 2.0, 1.9, 1.7]
observations[2,1,3] = [1.4, 1.0, 1.3, 1.2]
observations[2,2,1] = [2.4, 2.7, 2.4, 2.6]
observations[2,2,2] = [2.0, 2.3, 2.1, 2.4]
observations[2,2,3] = [2.4, 2.6, 2.3, 2.2]
observations[2,3,1] = [3.0, 3.1, 3.0, 2.7]
observations[2,3,2] = [3.1, 3.0, 2.8, 3.2]
observations[2,3,3] = [3.2, 2.9, 2.8, 2.9]

-or-

Specified in multidimensional array with 1st dimension as replicate
observations = cat(cat(hcat([1.9, 1.8, 1.6, 1.4], [1.8, 1.7, 1.4, 1.5]),
                       hcat([2.3, 2.1, 2.0, 2.6], [2.4, 2.7, 2.4, 2.6]),
                       hcat([2.9, 2.8, 3.4, 3.2], [3.0, 3.1, 3.0, 2.7]), dims = 3),
                   cat(hcat([2.1, 2.0, 1.8, 2.2], [2.3, 2.0, 1.9, 1.7]),
                       hcat([2.4, 2.6, 2.7, 2.3], [2.0, 2.3, 2.1, 2.4]),
                       hcat([3.6, 3.1, 3.4, 3.2], [3.1, 3.0, 2.8, 3.2]), dims = 3),
                   cat(hcat([1.1, 1.2, 1.0, 1.4], [1.4, 1.0, 1.3, 1.2]),
                       hcat([2.0, 2.1, 1.9, 2.2], [2.4, 2.6, 2.3, 2.2]),
                       hcat([2.9, 2.8, 3.0, 3.1], [3.2, 2.9, 2.8, 2.9]), dims = 3), dims = 4)

Factor A      1                             2                             3
Factor B      1         2         3         1         2         3         1         2         3
Factor C      1    2    1    2    1    2    1    2    1    2    1    2    1    2    1    2    1    2
            1.9  1.8  2.3  2.4  2.9  3.0  2.1  2.3  2.4  2.0  3.6  3.1  1.1  1.4  2.0  2.4  2.9  3.2
            1.8  1.7  2.1  2.7  2.8  3.1  2.0  2.0  2.6  2.3  3.1  3.0  1.2  1.0  2.1  2.6  2.8  2.9
            1.6  1.4  2.0  2.4  3.4  3.0  1.8  1.9  2.7  2.1  3.4  2.8  1.0  1.3  1.9  2.3  3.0  2.8
            1.4  1.5  2.6  2.6  3.2  2.7  2.2  1.7  2.3  2.4  3.2  3.2  1.4  1.2  2.2  2.2  3.1  2.9

Factor Types                                 FFF         FFR       FRF       RFF       FRR       RFR       RRF     RRR
                            SS    DF    MS     F     p     F         F         F         F         F         F       F
Total                     30.4    71
Cells                     28.4    17
    Factors
        Factor A           1.8     2   0.9  24.5  3e-8  4.9       3.3      24.5       3.3       4.9       2.2      2.2
        Factor B          24.7     2  12.3   332 5e-31  141       332      44.8      44.8      40.0       141     40.0
        Factor C          9e-3     1  9e-3   0.2   0.6  0.2       0.1      5e-2      4e-2      5e-2       0.1     5e-2
    All Interactions
        Pair Interactions
            Factor AxB     1.1     4   0.3   7.4  8e-5  5.0       7.4       7.4       7.4       5.0       5.0      5.0
            Factor AxC     0.4     2   0.2   5.0  1e-2  5.0       3.4       5.0       3.4       5.0       3.4      3.4
            Factor BxC     0.2     2  9e-2   2.4   0.1  2.4       2.4       1.6       1.6       1.6       2.4      1.6
        Factor AxBxC       0.2     4  6e-2   1.5   0.2  1.5       1.5       1.5       1.5       1.5       1.5      1.5
Error                      2.0    54  4e-2

Currently only works for 1-way, 2-way, and 3-way ANOVAs
Next: expand to fully nested 2-way ANOVAs
=#

function anova(observations::AbstractArray{T}, factortypes::Vector{FactorType} = [fixed]) where {T <: Union{Number, AbstractVector{<:Number}}}
    length(observations) > 0 || return
    validate(factortypes, ndims(observations))
    firstlevelreplicates = first(factortypes) == replicate

    nnestedfactors = count(f -> f == nested, factortypes)
    crossedfactortypes = filter(f -> f ∈ [fixed, random], factortypes)
    ncrossedfactors = length(crossedfactortypes)
    ncrossedfactors < 4 || throw(ErrorException("ANOVA with 4 or more crossed factors is not supported."))
    nreplicates = firstlevelreplicates ? size(observations, 1) : length(observations[1])
    firstlevelreplicates || all(c -> length(c) == nreplicates, observations) || throw(ErrorException("All cells must have the same number of replicates."))
    ncells = Int.(length(observations) / (firstlevelreplicates ? nreplicates : 1))
    nfactorlevels = firstlevelreplicates ? [size(observations)...][Not(1)] : [size(observations)...]

    anovakernel(observations, nreplicates, ncells, nnestedfactors, ncrossedfactors, nfactorlevels, crossedfactortypes)
end

function validate(factortypes::Vector{FactorType}, ndims; noreplicates = false)
    length(factortypes) == ndims || throw(ErrorException("factortypes must have an entry for each factor."))
    if noreplicates
        replicate ∉ factortypes || throw(ErrorException("replicates are not valid for this structure."))
    else
        replicate ∉ factortypes || first(factortypes) == replicate || throw(ErrorException("replicate must be the first entry if present"))
    end
    factortypes = filter(t -> t ≠ replicate, factortypes)
    nested ∉ factortypes || length(unique(factortypes[1:count(t -> t == nested, factortypes)])) == 1 || throw(ErrorException("nested entries must come before crossed factors"))
end

function anovakernel(observations, nreplicates, ncells, nnestedfactors, ncrossedfactors, nfactorlevels, crossedfactortypes)
    N = ncells * nreplicates
    nfactors = nnestedfactors + ncrossedfactors

    # collapse replicate dimension
    cellsums = calccellsums(observations, nfactors, nfactorlevels)
    C = sum(cellsums) ^ 2 / N
    total = totalcalc(observations, N, C)
    amongallnested, nestedsums, ncrossedfactorlevels, nnestedfactorlevels = amongnestedfactorscalc(cellsums, nfactorlevels, nnestedfactors, C)

    cells = cellscalc(cellsums, nreplicates, ncells, C)
    nonerror = nnestedfactors > 0 ? amongallnested[1] : cells
    error = errorcalc(total, nonerror, nfactorlevels, nreplicates)

    crossedfactors = factorscalc(nestedsums, ncrossedfactors, ncrossedfactorlevels, N, C)
    interactions, interactionsmap = interactionscalc(cells, nestedsums, crossedfactors, ncrossedfactors, ncrossedfactorlevels, nnestedfactorlevels, nreplicates, C)
    nestedfactors = nestedfactorscalc(amongallnested, nnestedfactors, crossedfactors, interactions)

    numerators = getnumerators(crossedfactors, ncrossedfactors, nnestedfactors, interactions)
    denominators = getdenominators(nnestedfactors, nestedfactors, nreplicates, error, total, crossedfactors, ncrossedfactors, crossedfactortypes, interactionsmap)

    # drop least significant term if nreplicates == 1
    if nreplicates == 1
        pop!(numerators)
        pop!(denominators)
    end

    # perform test
    ftest.(numerators, denominators)
end

function calccellsums(observations::T, nfactors, nfactorlevels) where {T <: AbstractArray{<:AbstractVector{<:Number}}}
    map(c -> sum(c), observations)
end

function calccellsums(observations::T, nfactors, nfactorlevels) where {T <: AbstractArray{<:Number}}
    ndims(observations) > nfactors || return observations
    reshape(sum(observations, dims = 1), (nfactorlevels...)) # check dropdims
end

function totalcalc(observations, N, C)
    ss = sum(c -> sum(c.^2), observations) - C
    df = N - 1
    AnovaValue(ss, df)
end

function factorscalc(cellsums, nfactors, nfactorlevels, N, C)
    factorindices = 1:nfactors
    ss = map(i -> sum(sum(cellsums, dims = factorindices[Not(i)]) .^ 2) / (N / nfactorlevels[i]), factorindices) .- C
    df = nfactorlevels .- 1
    AnovaFactor.(ss, df)
end

function cellscalc(cellsums, nreplicates, ncells, C)
    ss = sum(cellsums .^ 2) / nreplicates - C
    df = ncells - 1
    AnovaValue(ss, df)
end

function errorcalc(total, cells, nfactorlevels, nreplicates)
    ss = total.ss - cells.ss
    df = prod(nfactorlevels) * (nreplicates - 1)
    AnovaFactor(ss, df)
end

function remaindercalc(total, factors)
    ss = total.ss - sum(f -> f.ss, factors)
    df = prod(f -> f.df, factors)
    AnovaFactor(ss,df)
end

function amongnestedfactorscalc(cellsums, nfactorlevels, nnestedfactors, C)
    nestedsums = cellsums
    nlowerfactorlevels = 1
    nupperfactorlevels = nfactorlevels
    amongallnested = Vector{AnovaValue}(undef, nnestedfactors)
    for i ∈ 1:nnestedfactors
        # collapse each nested factor

        amongallnestedSS = sum(nestedsums .^ 2 ./ (nreplicates * prod(nlowerfactorlevels))) - C
        amongallnestedDF = prod(nupperfactorlevels) - 1

        amongallnested[i] = AnovaValue(amongallnestedSS, amongallnestedDF)

        nlowerfactorlevels = nfactorlevels[1:i]
        nupperfactorlevels = nfactorlevels[(i+1):end]
        nestedsums = reshape(sum(nestedsums, dims = i), (nupperfactorlevels...))
    end

    amongallnested, nestedsums, nupperfactorlevels, nlowerfactorlevels
end

function interactionscalc(cells, nestedsums, crossedfactors, ncrossedfactors, ncrossedfactorlevels, nnestedfactorlevels, nreplicates, C)
    interactionsmap = Dict{Any,AnovaFactor}()
    interactions = Vector{AnovaFactor}()
    if ncrossedfactors > 1
        pairwise = pairwisecalc(nestedsums, crossedfactors, ncrossedfactors, ncrossedfactorlevels, nnestedfactorlevels, nreplicates, C)
        interactionsmap[(1,2)] = pairwise[1,2]
        interactionsmap[(2,1)] = pairwise[1,2]
        push!(interactions, interactionsmap[(1,2)])

        if ncrossedfactors > 2
            threewise = threewisecalc(cells, crossedfactors, pairwise)
            interactionsmap[(1,3)] = pairwise[1,3]
            interactionsmap[(3,1)] = pairwise[1,3]
            interactionsmap[(2,3)] = pairwise[2,3]
            interactionsmap[(3,2)] = pairwise[2,3]
            interactionsmap[(1,2,3)] = threewise
            push!(interactions, interactionsmap[(1,3)], interactionsmap[(2,3)], interactionsmap[(1,2,3)])
        end
    end

    return interactions, interactionsmap
end

function pairwisecalc(nestedsums, crossedfactors, ncrossedfactors, ncrossedfactorlevels, nnestedfactorlevels, nreplicates, C)
    pairwise = Array{AnovaFactor,2}(undef, ncrossedfactors, ncrossedfactors)
    factorindexes = 1:ncrossedfactors
    for i ∈ factorindexes
        for j ∈ (i+1):ncrossedfactors
            otherfactorindexes = intersect(factorindexes[Not(i)], factorindexes[Not(j)])

            ss = sum(sum(nestedsums, dims = otherfactorindexes) .^ 2 ./ (prod(ncrossedfactorlevels[otherfactorindexes]) * prod(nnestedfactorlevels) * nreplicates)) - C - crossedfactors[i].ss - crossedfactors[j].ss
            df = crossedfactors[i].df * crossedfactors[j].df

            pairwise[i,j] = pairwise[j,i] = AnovaFactor(ss, df)
        end
    end
    pairwise
end

function threewisecalc(cells, factors, pairwise)
    ss = cells.ss - sum(f -> f.ss, factors) - pairwise[1,2].ss - pairwise[1,3].ss - pairwise[2,3].ss
    df = prod(f -> f.df, factors)
    AnovaFactor(ss, df)
end

function nestedfactorscalc(amongallnested, nnestedfactors, crossedfactors, interactions)
    nestedfactors = []
    interactions
    if nnestedfactors > 0
        nestedss = amongallnested[1].ss - sum(f -> f.ss, [crossedfactors; interactions])
        nesteddf = amongallnested[1].df - sum(f -> f.df, [crossedfactors; interactions])
        nestedfactors = [AnovaFactor(nestedss, nesteddf)]
    end
    nestedfactors
end

function getnumerators(crossedfactors, ncrossedfactors, nnestedfactors, interactions)
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

function getdenominators(nnestedfactors, nestedfactors, nreplicates, error, total, crossedfactors, ncrossedfactors, crossedfactortypes, interactionsmap)
    basedenominator = nnestedfactors > 0 ? nestedfactors[end] :
                         nreplicates > 1 ? error : remaindercalc(total, crossedfactors)

    if ncrossedfactors == 1
        denominators = [basedenominator]
    elseif ncrossedfactors == 2
        if all(f -> f == fixed, crossedfactortypes)
            crosseddenominators = repeat([basedenominator], ncrossedfactors)
        elseif all(f -> f == random, crossedfactortypes)
            crosseddenominators = repeat([interactionsmap[(1,2)]], ncrossedfactors)
        else
            crosseddenominators = map(f -> f == fixed ? basedenominator : interactionsmap[(1,2)], crossedfactortypes)
        end

        denominators = [crosseddenominators; basedenominator]
    elseif ncrossedfactors == 3
        if all(f -> f == fixed, crossedfactortypes)
            crosseddenominators = repeat([basedenominator], ncrossedfactors)
            pairwiseinteractiondenominators = repeat([basedenominator], ncrossedfactors)
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
            crosseddenominators[i] = basedenominator

            fixedindexes = (1:ncrossedfactors)[Not(i)]
            for j ∈ fixedindexes
                crosseddenominators[j] = interactions[(i,j)]
            end

            fixedinteractionindex = sum(fixedindexes) - 2
            pairwiseinteractiondenominators = Vector{AnovaFactor}(undef, ncrossedfactors)
            pairwiseinteractiondenominators[fixedinteractionindex] = interactionsmap[(1,2,3)]
            pairwiseinteractiondenominators[Not(fixedinteractionindex)] .= basedenominator
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
            pairwiseinteractiondenominators[randominteractionindex] = basedenominator
            pairwiseinteractiondenominators[Not(randominteractionindex)] .= interactionsmap[(1,2,3)]
        end

        denominators = [crosseddenominators; pairwiseinteractiondenominators; basedenominator]
    end

    # determine correct denominators for nested factors
    if nnestedfactors > 0
        nesteddenominators = Vector{AnovaFactor}(undef, nnestedfactors)
        nesteddenominators[1] = error
        for i ∈ 2:nnestedfactors
            nesteddenominators[i] = nestedfactors[i - 1]
        end
        append!(denominators, nesteddenominators)
    end

    denominators
end

function threewayinteraction(interaction_ab, interaction_bc, interaction_abc)
    reducedmeansquare(factor::AnovaFactor) = factor.ms ^ 2 / factor.df

    ms = interaction_ab.ms + interaction_bc.ms - interaction_abc.ms
    df = ms ^ 2 / (reducedmeansquare(interaction_ab) + reducedmeansquare(interaction_bc) + reducedmeansquare(interaction_abc))
    AnovaFactor(ms * df, df, ms)
end

function ftest(x, y)
    f = x.ms / y.ms
    fdist = FDist(x.df, y.df)
    p = ccdf(fdist, f)
    (f,p)
end

export anova, FactorType

end
