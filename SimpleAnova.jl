module SimpleAnova

using Statistics
using Distributions
include("InvertedIndices.jl")
include("AnovaValue.jl")
include("AnovaFactor.jl")
include("AnovaResult.jl")
include("AnovaData.jl")
include("FactorType.jl")

import Main.InvertedIndices.Not

const totalname = "Total"
const cellsname = "Cells"
const errorname = "Error"
const remaindername = "Remainder"


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

1-way ANOVA with 2 nested factors, replicates
Specified in multidimensional array with 1st dimension as replicate
observations = cat(cat(hcat([3.17, 4.41, 1.81, 1.74], [2.81, 4.98, 2.62, 2.53]),
                       hcat([3.0, 3.02, 4.73, 1.77], [4.13, 0.71, 3.18, 3.34]), dims = 3),
                   cat(hcat([2.42, 1.28, 1.4, 2.56], [1.26, 1.08, 1.42, 0.85]),
                       hcat([0.36, 0.35, 2.64, 3.75], [3.86, 2.53, 3.97, 3.03]), dims = 3),
                   cat(hcat([5.24, 2.24, 0.18, 3.06], [4.04, 4.14, 0.33, 4.61]),
                       hcat([6.06, 1.61, 2.25, 2.44], [0.02, 3.95, 0.87, 2.0]), dims = 3), dims = 4)
    Note: must specify that first dimension is a replicate






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


1-way ANOVA within subjects
observations = Array{Vector{Float64}, 2}(undef, 7, 3)
observations[1,1] = [164]
observations[1,2] = [152]
observations[1,3] = [178]
observations[2,1] = [202]
observations[2,2] = [181]
observations[2,3] = [222]
observations[3,1] = [143]
observations[3,2] = [136]
observations[3,3] = [132]
observations[4,1] = [210]
observations[4,2] = [194]
observations[4,3] = [216]
observations[5,1] = [228]
observations[5,2] = [219]
observations[5,3] = [245]
observations[6,1] = [173]
observations[6,2] = [159]
observations[6,3] = [182]
observations[7,1] = [161]
observations[7,2] = [157]
observations[7,3] = [165]

-or-

observations = [164 152 178; 202 181 222; 143 136 132; 210 194 216; 228 219 245; 173 159 182; 161 157 165]

                Factor
                  1   2   3
Subjects    1   164 152 178
            2   202 181 222
            3   143 136 132
            4   210 194 216
            5   228 219 245
            6   173 159 182
            7   161 157 165
=#

"""
    anova(observations, [factortypes, factorlabels])

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

N-way fixed-effects ANOVA with replicates in first dimension: anova(observations, [replicates])

2-way ANOVA with Factor 1 random and Factor B fixed with replicates in vectors: anova(observations, [random])

2-way ANOVA with Factor 1 fixed and Factor B random with replicates in vectors: anova(observations, [fixed, random])

2-way fixed-effects ANOVA with 2 random nested factors with replicates in first dimension:
anova(observations, [replicates, nested, nested])
"""
function anova(observations::AbstractArray{T}, factortypes::Vector{FactorType} = FactorType[], factornames::Vector{<:AbstractString} = String[], withinsubjects = false) where {T <: Union{Number, AbstractVector{<:Number}}}
    length(observations) > 0 || return

    # if empty, defaults to assuming ndims = nfactors and all factors are fixed
    if isempty(factortypes) || length(factortypes) < ndims(observations)
        nremaining = ndims(observations) - length(factortypes)
        factortypes = [factortypes; repeat([fixed], nremaining)]
    end

    firstlevelreplicates = first(factortypes) == replicates

    if isempty(factornames)
        factornames = ["A", "B", "C", "D", "E", "F"][1:(ndims(observations) - (firstlevelreplicates ? 1 : 0))]
        reverse!(factorlabels)
    end

    validate(factortypes, ndims(observations))

    # for randomized blocks, only real differences is that "blocks" factor isn't tested. Option?
    # for repeated measures, only difference is that "subjects" factor isn't tested. Option?
    #     if 2 factors, treated as 3-way ANOVA with subjects as random factor. Where does remainder come in with a 3-way?
    # for latin squares, test as a 3-factor ANOVA with 1 fixed and 2 random factors
    # multiway blocked/repeated measures with within-and-among factors will need another look, but may just be
    #     the same as nesting? No, looks like not because the factor could be fixed

    # If missing data in without-replicate design, need to calculate and provide the bias. Can provide that with below function?
    # Although, idea is for user to have to make few decisions...

    # Get replicate counts per cell.
    # If equal, proceed as normal
    # If proportional, proceed with matrix of nreplicates (question: how to make generalized?)
    # The function will not handle replacement of mising if disproportional, but Shearer's function will be available to allow user
    # to create such estimates for themselves.
    # Also provide a random deletion function?
    # Note: GLM can deal with high unequal replication

    # can provide power equations too.

    #  Different type Sums of Squares only apply for unbalanced designs
    #  Type I is bad and should never be done.
    #  Type II is only valid if there is no interaction. This is the same advice of my stats book.
    #  Type III is valid if there are interactions, but usually it isn't useful to interpret if there are interactions.
    #  I'll stick with Type II and document it.

    nnestedfactors = count(f -> f == nested, factortypes)
    crossedfactortypes = filter(f -> f ∈ [fixed, random], factortypes)
    ncrossedfactors = length(crossedfactortypes)
    ncrossedfactors < 4 || throw(ErrorException("ANOVA with 4 or more crossed factors is not supported."))
    nreplicates = firstlevelreplicates ? size(observations, 1) : length(observations[1])
    firstlevelreplicates || all(c -> length(c) == nreplicates, observations) || throw(ErrorException("All cells must have the same number of replicates."))
    ncells = Int.(length(observations) / (firstlevelreplicates ? nreplicates : 1))
    nfactorlevels = firstlevelreplicates ? [size(observations)...][Not(1)] : [size(observations)...]

    nonreplicatefactortypes = filter(t -> t ≠ replicates, factortypes)
    crossedfactornames = factornames[nonreplicatefactortypes .≠ nested]
    nestedfactornames = factornames[nonreplicatefactortypes .== nested]

    if withinsubjects
        anovasubjectskernel(observations, nreplicates)
    else
        anovakernel(observations, nreplicates, ncells, nnestedfactors, ncrossedfactors, nfactorlevels, crossedfactortypes, crossedfactornames, nestedfactornames)
    end
end

function anova(observations::AbstractVector{<:Number}, factorassignments::AbstractVector{AbstractVector{<:Int}}, factortypes::Vector{FactorType} = FactorType[], factornames::Vector{<:AbstractString} = String[])
    # take a vector of observations and a vector containing a vector for each factor assigning the observations to a factor level of that factor.

    # ensure observations are balanced
    # reorganize observations into matrix and complete as normal
    anova(observationsmatrix, factortypes, factornames)
end

function anova(observations::AbstractVector{<:Number}, factorassignments::AbstractVector{AbstractVector{<:Int}}, factortypes::Vector{FactorType} = FactorType[], factornames::Vector{<:AbstractString} = String[])
    # extract data from DataFame, place into matrix, and then proceed
    anova(observations, factorassignments, factortypes, factornames)
end

function validate(factortypes::Vector{FactorType}, ndims; noreplicates = false)
    length(factortypes) == ndims || throw(ErrorException("factortypes must have an entry for each factor."))
    if noreplicates
        replicates ∉ factortypes || throw(ErrorException("replicates are not valid for this structure."))
    else
        replicates ∉ factortypes || first(factortypes) == replicates || throw(ErrorException("replicates must be the first entry if present"))
    end
    nonreplicatefactortypes = filter(t -> t ≠ replicates, factortypes)
    nested ∉ factortypes || nonreplicatefactortypes[1:count(t -> t == nested, factortypes)] |> unique |> length == 1 || throw(ErrorException("nested entries must come before crossed factors"))
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
    nonerror = nnestedfactors > 0 ? amongallnested[1] : cells
    error = errorcalc(nreplicates > 1 ? errorname : remaindername, total, nonerror)

    crossedfactors = factorscalc(nestedsums, ncrossedfactors, ncrossedfactorlevels, N, C, crossedfactornames)
    interactions, interactionsmap = interactionscalc(cells, nestedsums, crossedfactors, ncrossedfactors, ncrossedfactorlevels, nnestedfactorlevels, nreplicates, C, crossedfactornames)
    nestedfactors = nestedfactorscalc(amongallnested, nnestedfactors, crossedfactors, interactions, nestedfactornames)

    numerators = getnumerators(crossedfactors, ncrossedfactors, nnestedfactors, nestedfactors, interactions)

    crossedbasedenominator = nnestedfactors > 0 ? nestedfactors[end] : error;

    denominators = getdenominators(nnestedfactors, nestedfactors, nreplicates, crossedbasedenominator, error, total, crossedfactors, ncrossedfactors, crossedfactortypes, interactionsmap)

    # drop least significant test if nreplicates == 1; either the lowest interaction level, or lowest nesting level if present
    if nreplicates == 1
        pop!(numerators)
        pop!(denominators)
    end

    # perform test
    results = ftest.(numerators, denominators) # need to order them from highest to lowest

    AnovaData([total; results; error])
end

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


    #test for A; A / subjects within factor A
    # test for B; B / withinsubjectinteractions
    # test for interaction; interaction / withinsubjectinteraction
end

function sumfirstdim(observations::T) where {T <: AbstractArray{<:AbstractVector}}
    map(sumfirstdim, observations)
end

function sumfirstdim(observations::T) where {T <: AbstractArray{<:Number}}
    if (ndims(observations) > 1)
        dropdims(sum(observations, dims = 1), dims = 1)
    else
        sum(observations)
    end
end

function sumfirstdim(observations::T) where {T <: AbstractVector{<:Number}}
    first(sum(observations, dims = 1))
end

function totalcalc(observations, N, C)
    ss = sum(c -> sum(c.^2), observations) - C
    df = N - 1
    AnovaValue(totalname, ss, df)
end

function factorscalc(cellsums, nfactors, nfactorlevels, N, C, factorlabels)
    factorindices = 1:nfactors
    ss = map(i -> sum(sum(cellsums, dims = factorindices[Not(i)]) .^ 2) / (N / nfactorlevels[i]), factorindices) .- C
    df = nfactorlevels .- 1
    AnovaFactor.(factorlabels, ss, df)
end

function cellscalc(cellsums, nreplicates, ncells, C)
    ss = sum(cellsums .^ 2) / nreplicates - C
    df = ncells - 1
    AnovaValue(cellsname, ss, df)
end

function errorcalc(name, total, nonerror)
    ss = total.ss - nonerror.ss
    df = total.df - nonerror.df
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
            push!(interactions, interactionsmap[(1,3)], interactionsmap[(2,3)], interactionsmap[(1,2,3)])
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

            pairwise[i,j] = pairwise[j,i] = AnovaFactor("$(crossedfactorlabels[i]) × $(crossedfactorlabels[j])", ss, df)
        end
    end
    return pairwise
end

function threewisecalc(cells, factors, pairwise, crossedfactorlabels)
    ss = cells.ss - sum(f -> f.ss, factors) - pairwise[1,2].ss - pairwise[1,3].ss - pairwise[2,3].ss
    df = prod(f -> f.df, factors)
    AnovaFactor("$(crossedfactorlabels[1]) × $(crossedfactorlabels[2]) × $(crossedfactorlabels[3])", ss, df)
end

function nestedfactorscalc(amongallnested, nnestedfactors, crossedfactors, interactions, nestedfactorlabels)
    nestedfactors = []
    if nnestedfactors > 0
        otherfactors = [crossedfactors; interactions]
        otherss = sum(f -> f.ss, otherfactors)
        otherdf = sum(f -> f.df, otherfactors)

        ss = map(f -> f.ss - otherss, amongallnested)
        df = map(f -> f.df - otherdf, amongallnested)
        append!(nestedfactors, AnovaFactor.(nestedfactorlabels, ss, df))
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
            crosseddenominators = map(f -> f == fixed ? crossedbasedenominator : interactionsmap[(1,2)], crossedfactortypes)
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
    AnovaResult(x, f, p)
end

export anova, AnovaEffect, AnovaValue, AnovaFactor, AnovaResult, FactorType

end
