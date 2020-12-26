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

    isrepeatedmeasures = within ∈ factortypes
    isrepeatedmeasures && hasreplicates && throw(ErrorException("Repeated measures design cannot have replicates; set hasreplicates = false."))

    firstlevelreplicates = eltype(observations) <: Number ? hasreplicates : false
    nfactors = ndims(observations) - (firstlevelreplicates ? 1 : 0)

    # defaults to assuming all unspecified factors are fixed
    if length(factortypes) < nfactors
        nremaining = nfactors - length(factortypes)
        append!(factortypes, repeat([fixed], nremaining))
    end

    replace!(factortypes, block => subject)

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

    crossedfactornames = factornames[factortypes .∉ Ref((nested, within))]
    nestedfactornames = factornames[factortypes .== nested]

    if isrepeatedmeasures
        anovakernel(observations, nreplicates, ncells, nnestedfactors, ncrossedfactors, nfactorlevels, crossedfactortypes, crossedfactornames, nestedfactornames, true)
    else
        anovakernel(observations, nreplicates, ncells, nnestedfactors, ncrossedfactors, nfactorlevels, crossedfactortypes, crossedfactornames, nestedfactornames, false)
    end
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
        within ∉ factortypes || (count(t -> t == within, factortypes) == 1 && factortypes[1] == within) || error("maximum of one within entry and must be first")
        subject ∉ factortypes || (count(t -> t == subject, factortypes) == 1 && factortypes[2] == subject) || error("maximum of one subject entry and must be second")
    end
    if !isempty(factornames)
        nfactors ≤ 26 || error("Can only automatically name up to 26 factors. Provide names explicitly.")
        length(factornames) == nfactors || error("factornames must have an entry for each factor.")
    end
end

function anovakernel(observations, nreplicates, ncells, nnestedfactors, ncrossedfactors, nfactorlevels, crossedfactortypes, crossedfactornames, nestedfactornames, isrepeatedmeasures)
    N = ncells * nreplicates
    nfactors = nnestedfactors + ncrossedfactors

    # collapse replicate dimension
    cellmeans = eltype(observations) <: Number && nreplicates == 1 ? Float64.(observations) : meanfirstdim(observations)
    total = totalcalc(observations)
    cells = cellscalc(cellmeans, nreplicates)

    if isrepeatedmeasures
        withinsubjects, nestedmeans = withinsubjectscalc(cellmeans, nfactorlevels, "Subjects")
        withinfactor = withincalc(cellmeans, nfactorlevels, factornames[2])
        nreplicates = nfactorlevels[1]
        nfactorlevels = nfactorlevels[2:end]
    else
        withinsubjects = ()
        withinfactor = ()
        nestedmeans = cellmeans
    end

    amongallnested, crossedcellmeans, ncrossedfactorlevels, nnestedfactorlevels = amongnestedfactorscalc(nestedmeans, nfactorlevels, nnestedfactors, nreplicates)

    npercrossedcell = nreplicates * prod(nnestedfactorlevels)
    crossedfactors = factorscalc(crossedcellmeans, ncrossedfactors, ncrossedfactorlevels, N, crossedfactornames)
    interactions = interactionscalc(cells, crossedcellmeans, crossedfactors, crossedfactornames, npercrossedcell)
    nestedfactors = nestedfactorscalc(amongallnested, crossedfactors, interactions, nestedfactornames)

    reverse!(crossedfactors)
    reverse!(ncrossedfactorlevels)
    if ncrossedfactors == 3 && nreplicates > 1
        interactions[1:3] = reverse(interactions[1:3])
    end

    

    error = if isrepeatedmeasures
        AnovaFactor(remaindername, withinsubjects - withinfactor)
    elseif length(amongallnested) > 0
        AnovaFactor(errorname, total - amongallnested[1])
    elseif nreplicates > 1
        AnovaFactor(errorname, total - cells)
    else
        AnovaFactor(remaindername, total - sum([crossedfactors; interactions[1:end-1]]))
    end

    numerators = getnumerators(withinfactor, crossedfactors, nestedfactors, interactions)
    denominators = getdenominators(isrepeatedmeasures, interactions, nestedfactors, error, crossedfactortypes)

    # drop least significant test if nreplicates == 1; lowest nesting level. Lowest interaction previously dropped if no nested factors. Should look to do same for nesting.
    if nnestedfactors > 0 && nreplicates == 1
        droppedfactor = pop!(numerators)
        pop!(denominators)
    end

    # perform test
    results = ftest.(numerators, denominators)

    results = effectsizescalc(results, denominators, total, ncrossedfactors, npercrossedcell, ncrossedfactorlevels, crossedfactortypes, nnestedfactors, nnestedfactorlevels, nreplicates)
    data = AnovaData([total; withinsubjects; results], total, ncrossedfactors, ncrossedfactorlevels, npercrossedcell, crossedfactors, denominators[1:ncrossedfactors], crossedcellmeans)
    nnestedfactors > 0 && nreplicates == 1 && push!(data.effects, droppedfactor)
    error.df > 0 && push!(data.effects, error)

    return data
end
#=
function rmanovakernel(observations, nreplicates, ncells, nnestedfactors, ncrossedfactors, nfactorlevels, crossedfactortypes, crossedfactornames, nestedfactornames, withinfactorname)
    N = ncells * nreplicates

    cellsums = eltype(observations) <: Number && nreplicates == 1 ? observations : sumfirstdim(observations)
    C = sum(cellsums) ^ 2 / N
    total = totalcalc(observations, N, C)

    withinnested = squaresum()


    amongallnested, nestedsums, ncrossedfactorlevels, nnestedfactorlevels = amongnestedfactorscalc(cellsums, nfactorlevels, nnestedfactors, nreplicates, C)
    cells = cellscalc(cellsums, nreplicates, ncells, C)

    
end
=#

#=
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
=#

meanfirstdim(observations::AbstractArray{<:AbstractVector}) = map(meanfirstdim, observations)
meanfirstdim(observations::AbstractArray{<:Number}) = dropdims(mean(observations, dims = 1), dims = 1)
meanfirstdim(observations::AbstractVector{<:Number}) = mean(observations, dims = 1) |> first

upcat(x::AbstractArray) = x
upcat(x::AbstractArray{<:AbstractVector}) = reshape(vcat(x...), (size(x[1])..., size(x)...))

function totalcalc(observations)
    observations = upcat(observations)
    df = length(observations) - 1
    ss = var(observations) * df
    AnovaValue(totalname, ss, df)
end

function cellscalc(cellmeans, nreplicates)
    df = length(cellmeans) - 1
    ss = var(cellmeans) * nreplicates * df
    AnovaValue(cellsname, ss, df)
end

function factorscalc(cellmeans, nfactors, nfactorlevels, N, factornames)
    df = nfactorlevels .- 1
    ss = [var(mean(cellmeans, dims=(1:nfactors)[Not(i)]), corrected = false) * N for i ∈ 1:nfactors]
    AnovaFactor.(factornames, ss, df)
end

function errorcalc(total, amongallnested, cells, otherfactors, nreplicates)
    if length(amongallnested) > 0
        otherfactor = amongallnested[1]
        name = errorname
    elseif nreplicates > 1
        otherfactor = cells
        name = errorname
    else
        otherfactor = sum(otherfactors)
        name = remaindername
    end

    AnovaFactor(name, total - otherfactor)
end

function amongnestedfactorscalc(nestedmeans, nfactorlevels, nnestedfactors, nreplicates)
    nlowerfactorlevels = 1
    nupperfactorlevels = nfactorlevels
    amongallnested = Vector{AnovaValue}(undef, nnestedfactors)
    @views for i ∈ 1:nnestedfactors
        # collapse each nested factor
        df = prod(nupperfactorlevels) - 1
        ss = var(nestedmeans) * nreplicates * df * prod(nlowerfactorlevels)
        amongallnested[i] = AnovaValue(ss, df)
        nlowerfactorlevels = nfactorlevels[1:i]
        nupperfactorlevels = nfactorlevels[(i+1):end]
        nestedmeans = meanfirstdim(nestedmeans)
    end
    amongallnested, nestedmeans, collect(nupperfactorlevels), collect(nlowerfactorlevels)
end

function withinsubjectscalc(cellmeans, nfactorlevels, withinfactorname)
    df = nfactorlevels[1] - 1
    ss = sum(var(cellmeans, dims=1)) * df
    AnovaValue(withinfactorname, ss, df * nfactorlevels[2]), meanfirstdim(cellmeans)
end

function withincalc(cellmeans, nfactorlevels, withinfactorname)
    df = nfactorlevels[1] - 1
    ss = var(mean(cellmeans, dims=2)) * (prod(nfactorlevels[2:end]) * df)
    AnovaFactor(withinfactorname, ss, df)
end

function interactionscalc(cells, crossedcellmeans, crossedfactors, crossedfactorlabels, npercrossedcell)
    ncrossedfactors = length(crossedfactors)
    if (ncrossedfactors == 2 && npercrossedcell > 1) || ncrossedfactors == 3
        npairs = ncrossedfactors == 2 ? 1 : 3
        ntriples = ncrossedfactors == 3 && npercrossedcell > 1 ? 1 : 0
        interactions = Vector{AnovaFactor}(undef, npairs + ntriples)
        interactions[1] = pairwisecalc(crossedcellmeans, crossedfactors, npercrossedcell, crossedfactorlabels, 1, 2)
        if ncrossedfactors == 3
            interactions[2:3] = [pairwisecalc(crossedcellmeans, crossedfactors, npercrossedcell, crossedfactorlabels, i, 3) for i ∈ 1:2]
            if npercrossedcell > 1
                interactions[4] = threewisecalc(cells, crossedfactors, interactions, crossedfactorlabels)
            end
        end
        return interactions
    else
        return AnovaFactor[]
    end
end

function pairwisecalc(crossedcellmeans, crossedfactors, npercrossedcell, crossedfactorlabels, i, j)
    otherfactorindexes = (1:length(crossedfactors))[Not([i, j])]
    ss = (var(mean(crossedcellmeans, dims = otherfactorindexes), corrected=false) * npercrossedcell * length(crossedcellmeans) - crossedfactors[i].ss - crossedfactors[j].ss)
    df = crossedfactors[i].df * crossedfactors[j].df
    AnovaFactor("$(crossedfactorlabels[j]) × $(crossedfactorlabels[i])", ss, df)
end

function threewisecalc(cells, factors, pairwise, crossedfactorlabels)
    ss = cells.ss - sum(f.ss for f ∈ factors) - sum(p.ss for p ∈ pairwise[1:3])
    df = prod(f -> f.df, factors)
    AnovaFactor("$(crossedfactorlabels[3]) × $(crossedfactorlabels[2]) × $(crossedfactorlabels[1])", ss, df)
end

function nestedfactorscalc(amongallnested, crossedfactors, interactions, nestedfactorlabels)
    nestedfactors = []
    nnestedfactors = length(amongallnested)
    if nnestedfactors > 0
        otherfactors = sum([crossedfactors; interactions])
        for i ∈ nnestedfactors:-1:1
            nested = amongallnested[i] - otherfactors
            otherfactors += nested
            push!(nestedfactors, AnovaFactor(nestedfactorlabels[i], nested))
        end
    end
    nestedfactors
end

function getnumerators(withinfactor, crossedfactors, nestedfactors, interactions)
    numerators = AnovaEffect[]
    if !isnothing(withinfactor)
        push!(numerators, withinfactor)
    end
    append!(numerators, crossedfactors)
    ncrossedfactors = length(crossedfactors)
    if ncrossedfactors > 1 && length(interactions) > 0
        push!(numerators, interactions[1])
        ncrossedfactors > 2 && push!(numerators, interactions[2], interactions[3], interactions[4])
    end
    length(nestedfactors) > 0 && append!(numerators, nestedfactors)
    numerators
end

function getdenominators(isrepeatedmeasures, interactions, nestedfactors, error, crossedfactortypes)
    nnestedfactors = length(nestedfactors)
    ncrossedfactors = length(crossedfactortypes)
    nesteddenominators = Vector{AnovaFactor}(undef, nnestedfactors)

    if isrepeatedmeasures
        subjectsdenominators = [error]
    else
        subjectsdenominators = []
    end

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
    append!(denominators, subjectsdenominators)
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
    differences = [results[i].ms - denominators[i].ms for i ∈ eachindex(results)] # 1 kb between this line and next
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
        if ncrossedfactors == 2 # this whole block not quite 1 kb
            if npercrossedcell > 1
                interactionindexes = ([1,2],)
                imax = 3
            else
                interactionindexes = ()
                imax = 2
            end
        else
            if npercrossedcell > 1
                interactionindexes = ([1,2], [1,3], [2,3], [1,2,3])
                imax = 7
            else
                interactionindexes = ([1,2], [1,3], [2,3])
                imax = 6
            end
        end

        icrossed = 1:ncrossedfactors # this whole block 1 kb
        iother = ncrossedfactors < imax ? ((ncrossedfactors + 1):imax) : []
        factors = Vector{Int}(undef, imax)
        factors[icrossed] = [crossedfactortypes[i] == fixed ? crossedfactordfs[i] : 1 for i ∈ icrossed]
        factors[iother] = [prod(factors[x]) for x ∈ interactionindexes]

        effectsdenominators = repeat([npercrossedcell], imax)
        israndom = [x == random for x ∈ crossedfactortypes] # Originally used broadcasted equality (.==) but causes high allocations as of 1.3.0-rc3
        isfixed = [x == fixed for x ∈ crossedfactortypes]
        crossedeffectsdenominators = effectsdenominators[icrossed]
        crossedeffectsdenominators[isfixed] .*= prod(ncrossedfactorlevels)
        crossedeffectsdenominators[israndom] .*= [prod(ncrossedfactorlevels[Not(i)]) for i ∈ icrossed[israndom]]
        effectsdenominators[icrossed] = crossedeffectsdenominators
        effectsdenominators[iother] .*= [prod(ncrossedfactorlevels[Not(icrossed[israndom] ∩ x)]) for x ∈ interactionindexes] # 7kb - set intersection is 1kb, has to be done for each interaction

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
