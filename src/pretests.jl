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

    nlevels = [nreplicates; ncells]
    sortorder = sortperm(repeat(1:nreplicates, Int(N / nreplicates)) .+
                         sum([factorassignments[i] .* prod(nlevels[1:i]) for i ∈ 1:nfactors]))
    observationsmatrix = reshape(observations[sortorder], nlevels...)



    levenekernel(observationsmatrix, nreplicates, ncells)
end

# need median version

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
