"""
    levene(observations::Array{Union{Number, Vector{Number}}})
    levene(observations::Vector{Number}, factorassignments::Vector{Vector{Any}})
    levene(df::DataFrame, observationscolumn::Symbol, factorcolumns::Vector{Symbol})

Performs Levene's homogeniety of variance test.

Operates on at least 3 crossed factors (fixed or random) and arbitrarily many random nested factors on balanced data. Levene's test is identical to performing an ANOVA on the absolute deviations of each cell.

If Levene's test is significant, the data is heteroscedastic (variances differ). If it is not, variances could differ or the sample size could be too small to tell.

# Arguments
- `observations`: Array containing the values to test. For the array, each dimension is a factor level, such that observations[2,5,3] indicates the 2nd level of the first factor, the 5th level of the second factor, and the 3rd level of the third factor. May contain values or vectors of values, where the vector contains replicates. Factors should be ordered with least significant first. For the vector, must provide `factorassignments` to specify factor levels.
- `factorassignments`: Vector of vectors of integers specifying how each observation is assigned to a factor level. Provide this when `observations` is given as a vector. Factor levels do not have to be consecutive or ordered. Nested factors must reuse factor levels currently.

Notes: Throws an exception if it can tell there are no replicates. Do not use with on any data structure you would pass into `anova()` using `hasrepliates=false`.

Output: `AnovaData` structure containing the test results for each factor.

# Examples
```julia
levene(observations)
```
"""
function levene(observations::AbstractArray{T}) where {T <: Union{Number, AbstractVector{<:Number}}}
    length(observations) > 0 || return

    firstlevelreplicates = eltype(observations) <: Number ? true : false
    nfactors = ndims(observations) - (firstlevelreplicates ? 1 : 0)

    nreplicates = firstlevelreplicates ? size(observations, 1) : length(observations[1])
    nreplicates > 1 || throw(ErrorException("Levene's test is invalid if there are no replicates."))
    firstlevelreplicates || all(c -> length(c) == nreplicates, observations) || throw(ErrorException("All cells must have the same number of replicates."))

    ncells = Int.(length(observations) / (firstlevelreplicates ? nreplicates : 1))
    nfactorlevels = firstlevelreplicates ? [size(observations)...][Not(1)] : [size(observations)...]

    observations = firstlevelreplicates ? reshape(observations, nreplicates, ncells) : reshape(observations, ncells)
    levenekernel(observations, nreplicates, ncells)
end

function levene(observations::AbstractVector{T}, factorassignments::AbstractVector{<:AbstractVector}) where {T <: Number}
    length(observations) > 0 || return
    nfactors = length(factorassignments)
    N = length(observations)
    N % nfactors == 0 || error("Design is unbalanced.")
    all(length.(factorassignments) .== N) || error("Each observation must have an assignment for each factor.")

    factorlevels = factorassignments .|> unique .|> sort
    nfactorlevels = length.(factorlevels)
    all(N .& nfactorlevels .== 0) || error("Design is unbalanced.")
    factorlevelcounts = [[count(l -> l == factorlevels[i][j], factorassignments[i]) for j ∈ 1:nfactorlevels[i]] for i ∈ 1:nfactors]
    nperfactorlevel = factorlevelcounts .|> unique
    all(nperfactorlevel .|> length .== 1) || error("Design is unbalanced.")
    nperfactorlevel = nperfactorlevel .|> first

    if !(isa(factorassignments, Number)) || any(maximum.(factorlevels) .> nfactorlevels)
        compressedfactorlevels = [1:i for i ∈ nfactorlevels]
        factorlevelremapping = [factorlevels[i] .=> compressedfactorlevels[i] for i ∈ 1:nfactors]
        factorassignments = [replace(factorassignments[i], factorlevelremapping[i]...) for i ∈ 1:nfactors]
    end

    #deduplicate above

    ncells = prod(nfactorlevels)
    nreplicates = Int(N / ncells)

    nlevels = [nreplicates; nfactorlevels]
    sortorder = sortperm(repeat(1:nreplicates, Int(N / nreplicates)) .+
                         sum([factorassignments[i] .* prod(nlevels[1:i]) for i ∈ 1:nfactors]))
    observationsmatrix = reshape(observations[sortorder], nreplicates, ncells)

    levenekernel(observationsmatrix, nreplicates, ncells)
end

function levenekernel(observations, nreplicates, ncells)
    cellsums = sumfirstdim(observations)
    cellmeans = cellsums ./ nreplicates
    if eltype(observations) <: Number
        absdevs = abs.(observations .- cellmeans')
    else
        absdevs = [abs.(observations[i] .- cellmeans[i]) for i ∈ eachindex(observations)]
    end
    anovakernel(absdevs, nreplicates, ncells, 0, 1, [ncells], [fixed], ["Groups"], String[])
end

export levene
