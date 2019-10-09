const totalname = "Total"
const cellsname = "Cells"
const errorname = "Error"
const remaindername = "Remainder"



"""
    anova(observations::Array{Union{Number, Vector{Number}}}, factortypes = FactorType[]; factornames = String[], hasreplicates = true)
    anova(observations::Vector{Number}, factorassignments::Vector{Vector{Any}}, factortypes = FactorType[]; factornames = String[], hasreplicates = true)
    anova(df::DataFrame, observationscolumn::Symbol, factorcolumns::Vector{Symbol}, factortypes = FactorType[]; factornames = String[])

Performs an Analysis of Variance (ANOVA) calculation.

Operates on up to 3 crossed factors (fixed or random) and arbitrarily many random nested factors, with or without
replicates, on balanced data.

# Arguments
- `observations`: Array containing the values to test. For the array, each dimension is a factor level, such that observations[2,5,3] indicates the 2nd level of the first factor, the 5th level of the second factor, and the 3rd level of the third factor. May contain values or vectors of values, where the vector contains replicates. Factors should be ordered with least significant first. For the vector, must provide `factorassignments` to specify factor levels.
- `factorassignments`: Vector of vectors of integers specifying how each observation is assigned to a factor level. Provide this when `observations` is given as a vector. Factor levels do not have to be consecutive or ordered. Nested factors must reuse factor levels currently.
- `factortypes`: Vector indicating the `FactorType` for each factor. If present, `replicates` must appear first, any `nested` after, and then `random` or `fixed` in any order. Specify `replicates` if the first dimension of the `observations` matrix contains replicate values (vs. contained in vectors). If too few values are provided, remaining are assumed to be `fixed`.
- `factornames`: Vector of names for each factor, excluding the replicate factor. If empty, will be automatically populated alphabetically.

Notes: The last index will be the top factor in the table.

Output: `AnovaData` structure containing the test results for each factor.

# Examples
```julia
anova(observations)                        # N-way fixed-effects ANOVA with replicates (vectors or first dimension)
anova(observations, hasreplicates = false) # N-way fixed-effects ANOVA without replicates (first dimension)
anova(observations, [random])              # N-way ANOVA with lower random factor and 1 or 2 upper fixed factors
anova(observations, [random])              # N-way ANOVA with lower random factor and 1 or 2 upper fixed factors
anova(observations, [fixed, random])       # N-way ANOVA with 1 lower fixed factor, 1 random factor, and 0 or 1 upper fixed factor
anova(observations, [nested, random])      # N-way fixed-effects ANOVA with 1 random nested factor, 1 random factor, and 1-2 fixed factors
```

# Glossary
- observation: The dependent variable.
- factor: An independent variable.
- factor level: A value of a factor.
- balanced: All combinations of factor levels have the same number of observations.
- crossed factor: A factor with levels that combine with the levels of all other crossed factors.
- fixed factor: A factor with fixed effects (e.g. treatment, concentration, exposure time).
- random factor: A factor with random effects (e.g. location, individual).
- nested factor: A random factor where the levels are unique to a combination of crossed factor levels (e.g. replicate).
- sum of squares (SS): A measure of variance that is dependent on sample size. Also called "sum of squared deviations."
- degrees of freedom (DF, ν): The number of bins in which the values could have been moved, if random.
- mean square (MS): SS / DF. Corrects for the larger variance expected if random values can be assigned to more bins. Also called "mean squared error" or "mean squared deviation."
- F-statistic: The division of MS values produce a result belonging to the "F distribution", the shape of which depends on the DF of the numerator and denominator. The location of this value on the distribution provides the p-value.
- p-value: The probability that, if all measurements had been drawn from the same population, you would obtain data at least as extreme as contained in your observations.
- effect size: The standardized difference in the measurement caused by the factor.
"""
function anova(observations::AbstractArray{T}, factortypes::Vector{FactorType} = FactorType[]; factornames::Vector{<:AbstractString} = String[], hasreplicates = true) where {T <: Union{Number, AbstractVector{<:Number}}}
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

#=
function anova(data::AnovaData, crossedfactors::Vector{Int}, )
    # performs a subtest of the specified crossed factors within level of the remaing crossed factors, using the original denominators
end

# possible bug: function anova(data::AnovaData, crossedfactors::Vector{Int}, ) with normal functions ==> hang when called?
=#

function anova(observations::AbstractVector{T}, factorassignments::AbstractVector{<:AbstractVector}, factortypes::Vector{FactorType} = FactorType[]; factornames::Vector{<:AbstractString} = String[]) where {T <: Number}
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

    nreplicates = Int(N / prod(nfactorlevels))

    nlevels = [nreplicates; nfactorlevels]
    sortorder = sortperm(repeat(1:nreplicates, Int(N / nreplicates)) .+
                         sum([factorassignments[i] .* prod(nlevels[1:i]) for i ∈ 1:nfactors]))
    observationsmatrix = reshape(observations[sortorder], nlevels...)
    anova(observationsmatrix, factortypes, factornames = factornames, hasreplicates = nreplicates > 1)
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
    constant = sum(cellsums) ^ 2 / N
    total = totalcalc(observations, N, constant)
    amongallnested, crossedcellsums, ncrossedfactorlevels, nnestedfactorlevels = amongnestedfactorscalc(cellsums, nfactorlevels, nnestedfactors, nreplicates, constant)
    cells = cellscalc(cellsums, nreplicates, ncells, constant)

    npercrossedcell = nreplicates * prod(nnestedfactorlevels)

    crossedfactors = factorscalc(crossedcellsums, ncrossedfactors, ncrossedfactorlevels, N, constant, crossedfactornames) # 3kb allocated here, possibly can't be avoided
    interactions = interactionscalc(cells, crossedcellsums, crossedfactors, ncrossedfactors, ncrossedfactorlevels, nnestedfactorlevels, nreplicates, constant, crossedfactornames, npercrossedcell) # 9 kb allocated here!
    nestedfactors = nestedfactorscalc(amongallnested, nnestedfactors, crossedfactors, interactions, nestedfactornames)

    reverse!(crossedfactors)
    reverse!(ncrossedfactorlevels)
    if ncrossedfactors == 3 && nreplicates > 1
        interactions[1:3] = reverse(interactions[1:3])
    end

    error = errorcalc(total, amongallnested, cells, [crossedfactors; interactions[1:end-1]], nnestedfactors, nreplicates) # could probably condense the otherfactor here
    numerators = getnumerators(crossedfactors, ncrossedfactors, nnestedfactors, nestedfactors, interactions)
    denominators = getdenominators(interactions, nestedfactors, error, crossedfactortypes, ncrossedfactors, nnestedfactors)

    # drop least significant test if nreplicates == 1; lowest nesting level. Lowest interaction previously dropped if no nested factors. Should look to do same for nesting.
    if nnestedfactors > 0 && nreplicates == 1
        droppedfactor = pop!(numerators)
        pop!(denominators)
    end

    # perform test
    results = ftest.(numerators, denominators)

    crossedcellmeans = crossedcellsums ./ npercrossedcell
    results = effectsizescalc(results, denominators, total, ncrossedfactors, npercrossedcell, ncrossedfactorlevels, crossedfactortypes, nnestedfactors, nnestedfactorlevels, nreplicates) # note: effect size doesn't account for nesting
    data = AnovaData([total; results], total, ncrossedfactors, ncrossedfactorlevels, npercrossedcell, crossedfactors, denominators[1:ncrossedfactors], crossedcellmeans)
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

function interactionscalc(cells, nestedsums, crossedfactors, ncrossedfactors, ncrossedfactorlevels, nnestedfactorlevels, nreplicates, C, crossedfactorlabels, npercrossedcell)
    interactions = Vector{AnovaFactor}()
    if (ncrossedfactors == 2 && npercrossedcell > 1) || ncrossedfactors == 3
        pairwise = pairwisecalc(nestedsums, crossedfactors, ncrossedfactors, ncrossedfactorlevels, nnestedfactorlevels, nreplicates, C, crossedfactorlabels)
        push!(interactions, pairwise[1,2])
        if ncrossedfactors == 3
            push!(interactions, pairwise[1,3])
            push!(interactions, pairwise[2,3])
            if npercrossedcell > 1
                threewise = threewisecalc(cells, crossedfactors, pairwise, crossedfactorlabels)
                push!(interactions, threewise)
            end
        end
    end
    interactions
end

function pairwisecalc(nestedsums, crossedfactors, ncrossedfactors, ncrossedfactorlevels, nnestedfactorlevels, nreplicates, C, crossedfactorlabels)
    pairwise = Array{AnovaFactor,2}(undef, ncrossedfactors, ncrossedfactors)
    factorindexes = 1:ncrossedfactors
    for i ∈ factorindexes
        for j ∈ (i+1):ncrossedfactors
            otherfactorindexes = factorindexes[Not([i, j])]
            ss = sum(sum(nestedsums, dims = otherfactorindexes) .^ 2 ./ (prod(ncrossedfactorlevels[otherfactorindexes]) * prod(nnestedfactorlevels) * nreplicates)) - C - crossedfactors[i].ss - crossedfactors[j].ss
            df = crossedfactors[i].df * crossedfactors[j].df
            pairwise[i,j] = pairwise[j,i] = AnovaFactor("$(crossedfactorlabels[j]) × $(crossedfactorlabels[i])", ss, df)
        end
    end
    pairwise
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
    nestedfactors
end

function getnumerators(crossedfactors, ncrossedfactors, nnestedfactors, nestedfactors, interactions)
    numerators = copy(crossedfactors)
    if ncrossedfactors > 1 && length(interactions) > 0
        push!(numerators, interactions[1])
        ncrossedfactors > 2 && push!(numerators, interactions[2], interactions[3], interactions[4])
    end
    nnestedfactors > 0 && append!(numerators, nestedfactors)
    numerators
end

function getdenominators(interactions, nestedfactors, error, crossedfactortypes, ncrossedfactors, nnestedfactors)
    nesteddenominators = Vector{AnovaFactor}(undef, nnestedfactors)

    if nnestedfactors > 0
        crossedbasedenominator = nestedfactors[1]
        nesteddenominators[1:end-1] = nestedfactors[2:end]
        nesteddenominators[end] = error
    else
        crossedbasedenominator = error
    end

    if ncrossedfactors == 1
        denominators = [crossedbasedenominator]
    elseif length(interactions) == 0
        denominators = repeat([crossedbasedenominator], ncrossedfactors)
    else # 2 or 3 crossed factors
        ninteractions = ncrossedfactors == 2 ? 1 : 4
        crosseddenominators = Vector{AnovaFactor}(undef, ncrossedfactors)
        pairinteractiondenominators = Vector{AnovaFactor}(undef, ninteractions - 1)

        if all(f -> f == fixed, crossedfactortypes)
            crosseddenominators = repeat([crossedbasedenominator], ncrossedfactors)
            pairinteractiondenominators = repeat([crossedbasedenominator], ninteractions - 1)

        elseif all(f -> f == random, crossedfactortypes)
            if ncrossedfactors == 2
                crosseddenominators = repeat(interactions, ncrossedfactors)

            else # 3 crossed factors
                crosseddenominators = Vector{AnovaFactor}(undef, 3)
                for i ∈ 1:3
                    j1, j2 = (1:3)[Not(i)]
                    crosseddenominators[i] = threewayinteraction(interactions[i + j1 - 2], interactions[i + j2 - 2], interactions[4])
                end
                pairinteractiondenominators = repeat([interactions[4]], ninteractions - 1)

            end

        else # has both random and fixed factors
            if ncrossedfactors == 2
                crosseddenominators = map(f -> f == fixed ? interactions[1] : crossedbasedenominator, crossedfactortypes)

            else # 3 crossed factors
                randomcount = count(f -> f == random, crossedfactortypes)
                if randomcount == 1
                    randomi = findfirst(f -> f == random, crossedfactortypes)
                    fixedi = (1:3)[Not(randomi)]
                    crosseddenominators = Vector{AnovaFactor}(undef, 3)
                    crosseddenominators[fixedi] = [interactions[i + randomi - 2] for i ∈ fixedi]
                    crosseddenominators[randomi] = crossedbasedenominator

                    fixedinteractioni = sum(fixedi) - 2
                    pairinteractiondenominators = Vector{AnovaFactor}(undef, 3)
                    pairinteractiondenominators[fixedinteractioni] = interactions[4]
                    pairinteractiondenominators[(1:3)[Not(fixedinteractioni)]] .= crossedbasedenominator

                else # 2 random factors
                    fixedi = findfirst(f -> f == fixed, crossedfactortypes)
                    randomi = (1:3)[Not(fixedi)]
                    crosseddenominators = Vector{AnovaFactor}(undef, 3)
                    crosseddenominators[fixedi] = threewayinteraction(interactions[fixedi + randomi[1] - 2], interactions[fixedi + randomi[2] - 2], interactions[4])
                    crosseddenominators[randomi] .= interactions[randomi[1] + randomi[2] - 2]

                    randominteractioni = sum(randomi) - 2
                    pairinteractiondenominators = Vector{AnovaFactor}(undef, 3)
                    pairinteractiondenominators[(1:3)[Not(randominteractioni)]] .= interactions[4]
                    pairinteractiondenominators[randominteractioni] = crossedbasedenominator

                end
            end
        end
        denominators = [crosseddenominators; pairinteractiondenominators; crossedbasedenominator]
    end
    append!(denominators, nesteddenominators)
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

function effectsizescalc(results, denominators, total, ncrossedfactors, npercrossedcell, ncrossedfactorlevels, crossedfactortypes, nnestedfactors, nnestedfactorlevels, nreplicates)
    differences = [results[i].ms - denominators[i].ms for i ∈ eachindex(results)]
    crossedfactordfs = [r.df for r ∈ results[1:ncrossedfactors]]

    if nreplicates == 1 && nnestedfactors > 0
        nnestedfactors -= 1
        nnestedfactorlevels = nnestedfactorlevels[1:(end-1)]
    end

    if ncrossedfactors == 1
        if nnestedfactors == 0
            ω² = [(results[1].ss - results[1].df * denominators[1].ms) / (total.ss + denominators[1].ms)]
        else
            effectdenominators = repeat([nreplicates], nnestedfactors + 1)
            nfactorlevels = [ncrossedfactorlevels; nnestedfactorlevels]
            effectdenominators[1] *= prod(nfactorlevels)
            factors = ones(Int, nnestedfactors + 1)
            factors[1] = crossedfactordfs[1]
            for i ∈ 2:nnestedfactors
                effectdenominators[2:(end - i + 1)] .*= nfactorlevels[end - i + 2]
            end
            σ² = factors .* differences ./ effectdenominators
            σ²total = sum(σ²) + denominators[end].ms
            ω² = σ² ./ σ²total
        end
    else
        if ncrossedfactors == 2
            if npercrossedcell > 1
                interactions = [[1,2]]
                imax = 3
            else
                interactions = []
                imax = 2
            end
        else
            if npercrossedcell > 1
                interactions = [[1,2], [1,3], [2,3], [1,2,3]]
                imax = 7
            else
                interactions = [[1,2], [1,3], [2,3]]
                imax = 6
            end
        end

        icrossed = 1:ncrossedfactors
        iother = ncrossedfactors < imax ? ((ncrossedfactors + 1):imax) : []
        factors = Vector{Int}(undef, imax)
        factors[icrossed] = [crossedfactortypes[i] == fixed ? crossedfactordfs[i] : 1 for i ∈ icrossed]
        factors[iother] = [prod(factors[x]) for x ∈ interactions]
        effectsdenominators = repeat([npercrossedcell], imax)

        israndom = crossedfactortypes .== random
        isfixed = crossedfactortypes .== fixed
        crossedeffectsdenominators = effectsdenominators[icrossed]
        crossedeffectsdenominators[isfixed] .*= prod(ncrossedfactorlevels)
        crossedeffectsdenominators[israndom] .*= [prod(ncrossedfactorlevels[Not(i)]) for i ∈ icrossed[israndom]]
        effectsdenominators[icrossed] = crossedeffectsdenominators
        effectsdenominators[iother] .*= [prod(ncrossedfactorlevels[Not(icrossed[israndom] ∩ x)]) for x ∈ interactions]

        σ² = factors .* differences[1:imax] ./ effectsdenominators

        if nnestedfactors > 0
            nestedrange = (length(results) .- nnestedfactors .+ 1):length(results)
            nestedeffectdenominators = repeat([nreplicates], nnestedfactors)
            for i ∈ 1:(nnestedfactors - 1)
                nestedeffectdenominators[1:(end - i + 1)] .*= nnestedfactorlevels[end - i + 2]
            end
            σ²nested = differences[nestedrange] ./ nestedeffectdenominators
            σ² = [σ²; σ²nested]
        end

        σ²total = sum(σ²) + denominators[end].ms
        ω² = σ² ./ σ²total
    end
    AnovaResult.(results, ω²)
end
