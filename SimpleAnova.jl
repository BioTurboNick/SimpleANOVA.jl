module SimpleAnova

using Statistics
using Distributions
using InvertedIndices

"""
    anova

    measurements = N-element vector of each measurement
    factors = M-element vector of CategoricalArrays assigning each measurement to a level of each of M factors
    factortype = M-element vector of FactorType

    Requires equal replication, assumes no missing data
"""

@enum FactorType fixed random

#=
function anova(measurements::Vector{T}, factors::Vector{CategoricalArray{Int}}, factortype::Vector{FactorType})
    N = length(measurements)
    factors = length(factors)

    all(length.(factors) .== N) || throw(ErrorException("Factor arrays must each have $N entries."))
    factortype == nfactors || throw(ErrorException("Factortype must have an entry for each factor."))

    nfactorlevels = length.(levels.(factors))

        squared = measurements .^ 2

end
=#

abstract type AnovaEffect
end

struct AnovaValue <: AnovaEffect
    ss::Float64
    df::Int
end

struct AnovaFactor <: AnovaEffect
    ss::Float64
    df::Int
    ms::Float64
end


AnovaFactor(ss, df) = AnovaFactor(ss, df, ss / df)

"""
    anova(observations)

observations - multidimensional array containing observations. Each dimension of
the array is a crossed factor. The elements of the array may be:
    - numbers or `missing` (only 1 observation per combination)
    - vectors (of vectors) of numbers or `missing` (nested random factors)

Attempts to fill missing values
"""

function ftest(x, y)
    f = x.ms / y.ms
    fdist = FDist(x.df, y.df)
    p = ccdf(fdist, f)
    (f,p)
end


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
Factor B      70.3   3    70.3   3.1    0.099
Error        366.4   6   366.4


Currently only works for 1-way and 2-way ANOVAs
Next: expand to fully nested 2-way ANOVAs
=#

function anova(observations::T, factortypes::Vector{FactorType} = [fixed]; firstlevelreplicates::Bool = false) where {T <: AbstractArray{<:Number}}
    length(observations) > 0 || return
    nfactors = ndims(observations) - (firstlevelreplicates ? 1 : 0)
    if nfactors > 1
        length(factortypes) == nfactors || throw(ErrorException("factortypes must have an entry for each factor."))
    end

    nreplicates = firstlevelreplicates ? size(observations, 1) : length(observations[1])
    nreplicates > 0 || return

    ncells = Int.(length(observations) / (firstlevelreplicates ? nreplicates : 1))


    nfactorlevels = firstlevelreplicates ? [size(observations)...][Not(1)] : [size(observations)...]

    anovakernel(observations, nreplicates, ncells, nfactors, nfactorlevels, factortypes)
end


function anova(observations::T, factortypes::Vector{FactorType} = [fixed]) where {T <: AbstractArray{<:AbstractVector{<:Number}}}
    length(observations) > 0 || return
    nfactors = ndims(observations)
    if nfactors > 1
        length(factortypes) == nfactors || throw(ErrorException("factortypes must have an entry for each factor."))
    end

    nreplicates = length(observations[1])
    nreplicates > 0 || return

    ncells = length(observations)

    all(c -> length(c) == nreplicates, observations) || throw(ErrorException("All cells must have the same number of replicates."))


    nfactorlevels = [size(observations)...]

    anovakernel(observations, nreplicates, ncells, nfactors, nfactorlevels, factortypes)
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

function totalcalc(observations, N, C)
    ss = sum(c -> sum(c.^2), observations) - C
    df = N - 1
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

function pairwiseinteractionscalc(cells, factors)
    factors_ss = map(f -> f.ss, factors)
    factors_df = map(f -> f.df, factors)

    ss = cells.ss .- (factors_ss .+ factors_ss')  # symmetric matrix of interaction terms, diagonal is meaningless
    df = factors_df .+ factors_df'
    AnovaFactor.(ss, df)
end

function calccellsums(observations::T, nfactors, nfactorlevels) where {T <: AbstractArray{<:AbstractVector{<:Number}}}
    map(c -> sum(c), observations)
end

function calccellsums(observations::T, nfactors, nfactorlevels) where {T <: AbstractArray{<:Number}}
    ndims(observations) > nfactors || return observations
    reshape(sum(observations, dims = 1), (nfactorlevels...))
end

function anovakernel(observations, nreplicates, ncells, nfactors, nfactorlevels, factortypes)
    N = ncells * nreplicates

    # collapse replicate dimension
    cellsums = calccellsums(observations, nfactors, nfactorlevels)

    C = sum(cellsums) ^ 2 / N

    total = totalcalc(observations, N, C)
    factors = factorscalc(cellsums, nfactors, nfactorlevels, N, C)

    if nreplicates > 1
        if nfactors > 1
            cells = nfactors == 1 ? factors[1] : cellscalc(cellsums, nreplicates, ncells, C)
            error = errorcalc(total, cells, nfactorlevels, nreplicates)
            pairwiseinteractions = pairwiseinteractionscalc(cells, factors)

            denominators = map(f -> f == fixed ? error : pairwiseinteractions[1,2], factortypes)

            f, p = ftest.(factors, denominators)
        else
            factor = factors[1]
            error = errorcalc(total, factor, nfactorlevels, nreplicates)
            f, p = ftest(factor, error)
        end
    else
        remainder = remaindercalc(total, factors)

        f, p = ftest.(factors, Ref(remainder))
    end
end

export anova, FactorType

end
