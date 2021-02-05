const totalname = "Total"
const cellsname = "Cells"
const errorname = "Error"
const remaindername = "Remainder"

"""
    anova(observations::Array{Union{Number, Vector{Number}}}, factortypes = FactorType[]; factornames = String[], hasreplicates = true)
    anova(observations::Vector{Number}, factorassignments::Vector{Vector{Any}}, factortypes = FactorType[]; factornames = String[], hasreplicates = true)
    anova(df::DataFrame, observationscolumn::Symbol, factorcolumns::Vector{Symbol}, factortypes = FactorType[]; factornames = String[])

Performs an Analysis of Variance (ANOVA) calculation.

Operates on fixed or random factors, subject/block factors, and nested factors, with or without replicates, on balanced data.

Limitations:
- If any factors are `random` type, limited to 3-way
- Repeated measures ANOVA (with subject/block factor) limited to 3 `fixed` factors which may be partitioned as within or among subjects

Operates on up to 3 crossed factors (fixed or random) and arbitrarily many random nested factors, with or without
replicates, on balanced data.

# Arguments
- `observations`: Array containing the values to test. For the array, each dimension is a factor level, such that observations[2,5,3] indicates the 2nd level of the first factor, the 5th level of the second factor, and the 3rd level of the third factor. May contain values or vectors of values, where the vector contains replicates. Factors should be ordered with least significant first. For the vector, must provide `factorassignments` to specify factor levels.
- `factorassignments`: Vector of vectors of integers specifying how each observation is assigned to a factor level. Provide this when `observations` is given as a vector. Factor levels do not have to be consecutive or ordered. Nested factors must reuse factor levels currently.
- `factortypes`: Vector indicating the `FactorType` for each factor. Factors must be ordered with `nested` first, then `random` or `fixed` in any order. For repeated measures, `subject` or `block` must appear after within-subject `fixed` and before among-subject `fixed`. If too few values are provided, remaining are assumed to be `fixed`.
- `factornames`: Vector of names for each factor, excluding the replicate factor. If empty, will be automatically populated alphabetically.
- `hasreplicates`: Boolean to specify if the first level should be considered 

Notes: The last index will be the top factor in the table.

Output: `AnovaData` structure containing the test results for each factor.

# Examples
```julia
anova(observations)                           # N-way fixed-effects ANOVA with replicates (vectors or first dimension)
anova(observations, hasreplicates = false)    # N-way fixed-effects ANOVA without replicates (first dimension)
anova(observations, [random])                 # N-way ANOVA with lower random factor and 1 or 2 upper fixed factors
anova(observations, [random])                 # N-way ANOVA with lower random factor and 1 or 2 upper fixed factors
anova(observations, [fixed, random])          # N-way ANOVA with 1 lower fixed factor, 1 random factor, and 0 or 1 upper fixed factor
anova(observations, [nested, random])         # N-way fixed-effects ANOVA with 1 random nested factor, 1 random factor, and 1-2 fixed factors
anova(observations, [fixed, subject])         # N-way repeated measures ANOVA with 1 within-subjects fixed factor
anova(observations, [fixed, block])           # N-way repeated measures ANOVA with 1 within-block fixed factor
anova(observations, [nested, fixed, subject]) # N-way repeated measures ANOVA with 1 within-subjects fixed factor and 1 nested factor
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
- subject/block factor: A nested factor that is subjected to multiple levels of another factor.
- sum of squares (SS): A measure of variance that is dependent on sample size. Also called "sum of squared deviations."
- degrees of freedom (DF, ν): The number of bins in which the values could have been moved, if random.
- mean square (MS): SS / DF. Corrects for the larger variance expected if random values can be assigned to more bins. Also called "mean squared error" or "mean squared deviation."
- F-statistic: The division of MS values produce a result belonging to the "F distribution", the shape of which depends on the DF of the numerator and error. The location of this value on the distribution provides the p-value.
- p-value: The probability that, if all measurements had been drawn from the same population, you would obtain data at least as extreme as contained in your observations.
- effect size: The standardized difference in the measurement caused by the factor.
"""
function anova(observations::AbstractVector{T}, factorassignments::AbstractVector{<:AbstractVector}, factortypes::Vector{FactorType} = FactorType[]; factornames::Vector{<:AbstractString} = String[]) where {T <: Number}
    # convert vector arguments into a multidimensional matrix

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

function anova(observations::AbstractArray{T}, factortypes::Vector{FactorType} = FactorType[]; factornames::Vector{<:AbstractString} = String[], hasreplicates = true) where {T <: Union{Number, AbstractVector{<:Number}}}
    length(observations) > 0 || return

    isrepeatedmeasures = subject ∈ factortypes

    eltype(observations) <: Number ||
        all(c -> length(c) == length(observations[1]), observations) ||
        error("All cells must have the same number of replicates.")

    if eltype(observations) <: AbstractVector && length(observations[1]) == 1
        observations = first.(observations)
        hasreplicates = false
    else
        observations = observations |> upcat
    end

    observations = observations .|> Float64

    nfactors = ndims(observations) - (hasreplicates ? 1 : 0)
    
    # defaults to assuming all unspecified factors are fixed
    if length(factortypes) < nfactors
        nremaining = nfactors - length(factortypes)
        append!(factortypes, repeat([fixed], nremaining))
    end
    replace!(factortypes, block => subject)

    # automatically assigns alphabetical names if not provided.
    if isempty(factornames)
        factornames = string.('A':'Z')[1:min(nfactors, 26)]
        reverse!(factornames)
    end

    validate(observations, factortypes, factornames, nfactors, hasreplicates)

    anovakernel(observations, factornames, factortypes, hasreplicates)
end

function validate(observations, factortypes::Vector{FactorType}, factornames::Vector{<:AbstractString}, nfactors, hasreplicates)
    length(factortypes) == nfactors ||
            error("`factortypes` must have an entry for each factor.")

    nested ∉ factortypes ||
        factortypes[1:count(isnested, factortypes)] |> unique |> length == 1 ||
        error("Nested factors must come before any other factors.")

    if subject ∈ factortypes
        count(issubject, factortypes) == 1 ||
            error("Maximum of one subject/block factor.")

        notnestedfactortypes = filter(f -> !isnested(f), factortypes)
        subject ∈ notnestedfactortypes[1:min(4, end)] ||
            error("Subject/block factor must be in the first four factors.")
        
        length(notnestedfactortypes) < 5 ||
            error("Maximum of 3 within-subjects or among-subjects factors.")

        random ∉ notnestedfactortypes ||
            error("Random factors are not supported with subject/block factor.")
    end

    crossedfactortypes = filter(f -> f ∈ [fixed, random], factortypes)
    length(crossedfactortypes) < 4 ||
        all(isfixed, crossedfactortypes) ||
        error("ANOVA with 4 or more crossed factors is not supported if any are random.")

    length(factornames) == nfactors || error("factornames must have an entry for each factor.")
end



#=
function anova(data::AnovaData, crossedfactors::Vector{Int}, )
    # performs a subtest of the specified crossed factors within level of the remaing crossed factors, using the original errors
end

# possible bug: function anova(data::AnovaData, crossedfactors::Vector{Int}, ) with normal functions ==> hang when called?
=#

function anovakernel(observations::AbstractArray{<:Number}, factornames, factortypes, hasreplicates)
    isrepeatedmeasures = subject ∈ factortypes

    totaldf = length(observations) - 1
    totalvar = anovavalue(totalname, var(observations), totaldf)

    if hasreplicates
        nreplicates = size(observations, 1)
        cellmeans = meanfirstdim(observations)
    else
        nreplicates = 1
        cellmeans = observations
    end
    
    cellsdf = length(cellmeans) - 1
    cellsvar = anovavalue(cellsname, var(cellmeans) * nreplicates, cellsdf)

    errorvar = AnovaFactor(errorname, totalvar - cellsvar)

    nnested = count(isnested, factortypes)
    nestedfactornames = @view factornames[1:nnested]
    nnestedfactorlevels = size(cellmeans)[1:nnested]
    notnestedfactornames = @view factornames[(nnested + 1):end]
    notnestedfactortypes = @view factortypes[(nnested + 1):end]

    nestedvars, nestederrorvars, cellmeans, nnestedreplicates = nestedfactors(cellmeans, nreplicates, nestedfactornames, errorvar)

    if isrepeatedmeasures
        si = findfirst(x -> x == subject, notnestedfactortypes)
        factorvars, factorvarsdict = anovafactors(dropdims(mean(cellmeans, dims = si), dims = si), nnestedreplicates * size(cellmeans, si), notnestedfactornames[Not(si)])
        subjectinteractionvarsdict = subjectinteractions(cellmeans, si, nnestedreplicates, notnestedfactornames)
        merge!(factorvarsdict, subjectinteractionvarsdict)
        factorerrorvars = anovasubjecterrors(factorvarsdict, notnestedfactortypes)
    else
        factorvars, factorvarsdict = anovafactors(cellmeans, nnestedreplicates, notnestedfactornames)
        factorerrorvars = nnested > 0 ? anovaerrors(factorvarsdict, notnestedfactortypes, nestedvars[1]) :
                          hasreplicates ? anovaerrors(factorvarsdict, notnestedfactortypes, errorvar) :
                          anovaerrors(factorvarsdict, notnestedfactortypes, factorvars[end])
    end

    append!(factorvars, nestedvars)
    append!(factorerrorvars, nestederrorvars)

    if !hasreplicates
        if isrepeatedmeasures
            remaindervar = factorerrorvars[end]
        else
            remaindervar = pop!(factorvars)
            pop!(factorerrorvars)
        end
        factorresults = ftest.(factorvars, factorerrorvars)
        errorvar = AnovaFactor(remaindername, remaindervar)
    else
        factorresults = ftest.(factorvars, factorerrorvars)
    end

    # not yet refactored, or accomodating for subject factors
    npercrossedcell = length(observations) ÷ length(cellmeans)
    nfactorlevels = Int[size(cellmeans)...]
    effectsizes = isrepeatedmeasures ? repeat([NaN], length(factorresults)) :
                                       effectsizescalc(factorresults, factorerrorvars, totalvar, npercrossedcell, nfactorlevels, notnestedfactortypes, nnested, nnestedfactorlevels, nreplicates)

    factorresults = AnovaResult.(factorresults, effectsizes)

    allresults = [totalvar; factorresults; errorvar]

    return AnovaData(allresults, factorerrorvars, nfactorlevels, nreplicates, cellmeans)
end

function anovafactors(cellmeans, nreplicates, factornames)
    # calculate all simple factors and all interactions

    N = length(cellmeans) * nreplicates

    nfactors = ndims(cellmeans)
    factors = 1:nfactors

    allfactors = []
    allfactorvars = AnovaFactor[]
    allfactorsdict = Dict{Any, AnovaEffect}()
    for i ∈ factors
        ifactors = collect(combinations(reverse(factors), i))
        iotherfactors = [factors[Not(i...)] for i ∈ ifactors]
        iupperfactorvars = [allfactorvars[findall(x -> x ⊆ i, allfactors)] for i ∈ ifactors]

        ifactornames = [makefactorname(factornames[i]) for i ∈ ifactors]
        ifactorss = [var(mean(cellmeans, dims = iotherfactors[i]), corrected = false) * N - sum(iupperfactorvars[i]).ss for i ∈ eachindex(iotherfactors)]
        ifactordf = isempty(iupperfactorvars[1]) ? [size(cellmeans, invertfactorindex(i, nfactors)) - 1 for i ∈ factors] :
                                                   [prod(f.df for f ∈ iupperfactorvars[j][1:i]) for j ∈ eachindex(iotherfactors)]
        ifactorvars = AnovaFactor.(ifactornames, ifactorss, ifactordf)
        
        append!(allfactors, ifactors)
        append!(allfactorvars, ifactorvars)
        merge!(allfactorsdict, Dict(invertfactorindex.((tuple.(i...) for i ∈ ifactors), nfactors) .=> ifactorvars))
    end

    return allfactorvars, allfactorsdict
end

function anovaerrors(factorvarsdict, factortypes, errorvar)
    # assign proper error terms for each factor

    length(factortypes) == 1 && return [errorvar]
    all(isfixed, factortypes) && return repeat([errorvar], length(factorvarsdict))

    factortypes = reverse(factortypes)

    if length(factortypes) == 2
        if factortypes[1] == factortypes[2]
            factor1error = factor2error = factorvarsdict[(1,2)]
        else
            if factortypes[1] == fixed
                factor1error = factorvarsdict[(1,2)]
                factor2error = errorvar
            else
                factor1error = errorvar
                factor2error = factorvarsdict[(1,2)]
            end
        end
        factorerrorvars = [factor1error; factor2error; errorvar]

    elseif length(factortypes) == 3
        if factortypes[1] == factortypes[2] == factortypes[3]
            factor1error = threeway_random_error(factorvarsdict[(1,2)], factorvarsdict[(1,3)], factorvarsdict[(1,2,3)])
            factor2error = threeway_random_error(factorvarsdict[(1,2)], factorvarsdict[(2,3)], factorvarsdict[(1,2,3)])
            factor3error = threeway_random_error(factorvarsdict[(1,3)], factorvarsdict[(2,3)], factorvarsdict[(1,2,3)])
            interaction12error = interaction13error = interaction23error = factorvarsdict[(1,2,3)]

        elseif factortypes[1] == factortypes[2]
            if factortypes[1] == fixed
                factor1error = factorvarsdict[(1,3)]
                factor2error = factorvarsdict[(2,3)]
                factor3error = errorvar
                interaction12error = factorvarsdict[(1,2,3)]
                interaction13error = interaction23error = errorvar

            else
                factor1error = factor2error = factorvarsdict[(1,2)]
                factor3error = threeway_random_error(factorvarsdict[(1,3)], factorvarsdict[(2,3)], factorvarsdict[(1,2,3)])
                interaction12error = errorvar
                interaction13error = interaction23error = factorvarsdict[(1,2,3)]

            end
        elseif factortypes[1] == factortypes[3]
            if factortypes[1] == fixed
                factor1error = factorvarsdict[(1,2)]
                factor2error = errorvar
                factor3error = factorvarsdict[(2,3)]
                interaction12error = interaction23error = errorvar
                interaction13error = factorvarsdict[(1,2,3)]

            else
                factor1error = factor3error = factorvarsdict[(1,3)]
                factor2error = threeway_random_error(factorvarsdict[(1,2)], factorvarsdict[(2,3)], factorvarsdict[(1,2,3)])
                interaction12error = interaction23error = factorvarsdict[(1,2,3)]
                interaction13error = errorvar

            end
        else
            if factortypes[2] == fixed
                factor1error = errorvar
                factor2error = factorvarsdict[(1,2)]
                factor3error = factorvarsdict[(1,3)]
                interaction12error = interaction13error = errorvar
                interaction23error = factorvarsdict[(1,2,3)]

            else
                factor1error = threeway_random_error(factorvarsdict[(1,2)], factorvarsdict[(1,3)], factorvarsdict[(1,2,3)])
                factor2error = factor3error = factorvarsdict[(2,3)]
                interaction12error = interaction13error = factorvarsdict[(1,2,3)]
                interaction23error = errorvar

            end
        end
        factorerrorvars = [factor1error; factor2error; factor3error; interaction12error; interaction13error; interaction23error; errorvar]

    else
        error("More than 3 factors with any random are not supported.")
    end

    return factorerrorvars
end

function anovasubjecterrors(factorvarsdict, factortypes)
    factortypes = reverse(factortypes)
    si = findfirst(x -> x == subject, factortypes) # could just hardcode, but in case I generalize later

    nfactors = length(factortypes)

    if nfactors == 2
        if si == 2
            # one among-subject factor
            factorerrorvars = [factorvarsdict[(1, si)]]
        else
            # one within-subject factor    
            factorerrorvars = [factorvarsdict[(si, 2)]]
        end

    elseif nfactors == 3
        if si == 3
            # two among-subject factors
            factorerrorvars = [factorvarsdict[(1, si)], factorerrorvars[(2, si)],
                               factorerrorvars[(1, 2, si)]]
        elseif si == 2
            # one among-subject factor and one within-subject factor
            factorerrorvars = [factorvarsdict[(1, si)]; factorvarsdict[(1, si, 3)]; 
                               factorvarsdict[(1, si, 3)]]
        else
            # two within-subject factors
            factorerrorvars = [factorvarsdict[(si, 2)]; factorvarsdict[(si, 3)];
                               factorvarsdict[(si, 2, 3)]]
        end

    elseif nfactors == 4
        if si == 4
            # three among factors
            factorerrorvars = [factorvarsdict[(1, 2, 3, si)]; factorvarsdict[(1, 2, 3, si)]; factorvarsdict[(1, 2, 3, si)];
                               factorvarsdict[(1, 2, 3, si)]; factorvarsdict[(1, 2, 3, si)]; factorvarsdict[(1, 2, 3, si)];
                               factorvarsdict[(1, 2, 3, si)]]
        elseif si == 3
            # two among factors and one within-subject factor
            factorerrorvars = [factorvarsdict[(1, 2, si)]; factorvarsdict[(1, 2, si)];
                               factorvarsdict[(1, 2, si)];
                               factorvarsdict[(1, 2, si, 4)]; factorvarsdict[(1, 2, si, 4)]; factorvarsdict[(1, 2, si, 4)];
                               factorvarsdict[(1, 2, si, 4)]]
        elseif si == 2
            # one among factor and two within-subject factors
            factorerrorvars = [factorvarsdict[(1, si)];
                               factorvarsdict[(1, si, 3)]; factorvarsdict[(1, si, 4)];
                               factorvarsdict[(1, si, 3)]; factorvarsdict[(1, si, 4)]; factorvarsdict[(1, si, 3, 4)];
                               factorvarsdict[(1, si, 3, 4)]]
        else
            # three within-subject factors
            factorerrorvars = [factorvarsdict[(si, 2)]; factorvarsdict[(si, 3)]; factorvarsdict[(si, 4)];
                               factorvarsdict[(si, 2, 3)]; factorvarsdict[(si, 2, 4)]; factorvarsdict[(si, 3, 4)];
                               factorvarsdict[(si, 2, 3, 4)]]
        end

    else
        error("More than 3 non-subject factors are not supported.")
    end

    return AnovaFactor.(factorerrorvars)
end

function nestedfactors(cellmeans, nreplicates, factornames, errorvar)
    # compute nested factors and collapse nested levels

    nnested = length(factornames)

    nestedvars = AnovaFactor[]
    for i ∈ 1:nnested
        nestedmeans = mean(cellmeans, dims = 1)
        ss = sum((cellmeans .- nestedmeans) .^ 2) * nreplicates
        df = (size(cellmeans, 1) - 1) * prod(size(nestedmeans))
        push!(nestedvars, AnovaFactor(factornames[i], ss, df))
        nreplicates *= size(cellmeans, 1)
        cellmeans = dropdims(nestedmeans, dims = 1)
    end
    nestederrorvars = nnested > 0 ? [errorvar; nestedvars[1:(end - 1)]] :
                                    []
    return reverse!(nestedvars), reverse!(nestederrorvars), cellmeans, nreplicates
end

function subjectfactor(cellmeans, si, nreplicates, factornames)
    withinfactors = 1:(si - 1)
    nreplicates *= prod(size(cellmeans)[withinfactors])
    cellmeans = dropdims(mean(cellmeans, dims = withinfactors), dims = tuple(withinfactors...))
    subjectmeans = mean(cellmeans, dims = 1)
    ss = sum((cellmeans .- subjectmeans) .^ 2) * nreplicates
    df = (size(cellmeans, 1) - 1) * prod(size(subjectmeans))
    AnovaFactor(factornames[si], ss, df)
end

function subjectinteractions(cellmeans, si, nreplicates, factornames)
    cellmeans = reshape(cellmeans, (size(cellmeans)[1:si]..., prod(size(cellmeans)[(si + 1):end])))
    namongfactors = length(factornames) - si
    
    factorvarsdict = Dict{Any, AnovaEffect}()
    for i ∈ axes(cellmeans, si + 1)
        _, ifactorvarsdict = @views anovafactors(selectdim(cellmeans, si + 1, i), nreplicates, factornames[1:si])

        if isempty(factorvarsdict)
            for k ∈ keys(ifactorvarsdict)
                1 ∈ k || continue

                newk = ((1:namongfactors)..., (k .+ namongfactors)...)
                ifactor = ifactorvarsdict[k]
                factorvarsdict[newk] = AnovaValue(makefactorname(factornames[invertfactorindex.(collect(k[2:end]), length(factornames))], factornames[si], factornames[(si + 1):end]), ifactor)
            end
        else
            for k ∈ keys(ifactorvarsdict)
                1 ∈ k || continue

                newk = ((1:namongfactors)..., (k .+ namongfactors)...)
                factorvarsdict[newk] += ifactorvarsdict[k]
            end
        end
    end
    return factorvarsdict
end

anovavalue(name, variance, df) = AnovaValue(name, variance * df, df)

meanfirstdim(observations::AbstractArray{<:Number}) = dropdims(mean(observations, dims = 1), dims = 1)
meanfirstdim(observations::AbstractArray{<:Array}) = dropdims(mean(observations, dims = 1), dims = 1)

upcat(x::AbstractArray) = x
upcat(x::AbstractArray{<:AbstractVector}) = reshape(vcat(x...), (size(x[1])..., size(x)...))

invertfactorindex(index, nfactors) = nfactors .- index .+ 1

makefactorname(factorname::AbstractString) = factorname

function makefactorname(factornames::AbstractVector{<:AbstractString})
    length(factornames) > 0 || return ""
    
    factorname = factornames[1]
    
    for i ∈ 2:length(factornames)
        factorname *= " × " * factornames[i]
    end

    return factorname
end

function makefactorname(withinfactornames::AbstractVector{<:AbstractString}, subjectfactorname, amongfactornames::AbstractVector{<:AbstractString})
    subjectname = subjectfactorname 

    if !isempty(amongfactornames)
        subjectname *= "[" * makefactorname(amongfactornames) * "]"
    end

    makefactorname([withinfactornames; subjectname])
end

function threeway_random_error(interaction_ab, interaction_bc, interaction_abc)
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


effectsizescalc(factorresults, factorerrorvars, totalvar, npercrossedcell, nfactorlevels, notnestedfactortypes, nnested, nnestedfactorlevels, nreplicates)
function effectsizescalc(factorresults, factorerrorvars, totalvar, npercrossedcell, ncrossedfactorlevels, crossedfactortypes, nnestedfactors, nnestedfactorlevels, nreplicates)
    ncrossedfactors = length(crossedfactortypes)
    differences = [factorresults[i].ms - factorerrorvars[i].ms for i ∈ eachindex(factorresults)] # 1 kb between this line and next
    ncrossedfactorlevels = reverse(ncrossedfactorlevels)
    nnestedfactorlevels = reverse(nnestedfactorlevels)
    crossedfactortypes = reverse(crossedfactortypes)

    if nreplicates == 1 && nnestedfactors > 0
        nnestedfactors -= 1
        nnestedfactorlevels = nnestedfactorlevels[1:(end-1)]
    end

    if ncrossedfactors == 1
        if nnestedfactors == 0
            ω² = [(factorresults[1].ss - factorresults[1].df * factorerrorvars[1].ms) / (totalvar.ss + factorerrorvars[1].ms)]
        else
            effectdenominators = repeat([nreplicates], nnestedfactors + 1)
            nfactorlevels = [ncrossedfactorlevels...; nnestedfactorlevels...]
            effectdenominators[1] *= prod(nfactorlevels)
            factors = ones(Int, nnestedfactors + 1)
            factors[1] = factorresults[1].df
            for i ∈ 2:nnestedfactors
                effectdenominators[2:(end - i + 1)] .*= nfactorlevels[end - i + 2]
            end
            σ² = factors .* differences ./ effectdenominators
            σ²total = sum(σ²) + factorerrorvars[end].ms
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
        
        factors = calculate_effect_numerator_factors(imax, factorresults, ncrossedfactors, crossedfactortypes, interactionindexes)

        function f(imax, ncrossedfactors, ncrossedfactorlevels, npercrossedcell, crossedfactortypes, interactionindexes)
            icrossed = 1:ncrossedfactors
            iother = (ncrossedfactors + 1):imax

            effectsdenominators = fill(npercrossedcell, imax)
            israndomtype = israndom.(crossedfactortypes) # Originally used broadcasted equality (.==) but causes high allocations as of 1.3.0-rc3
            crossedeffectsdenominators = effectsdenominators[icrossed] # responsible for 2 allocations as a view, 1 allocation without a view?
            crossedeffectsdenominators .*= prod(ncrossedfactorlevels)
            crossedeffectsdenominators[israndomtype] .÷= @view ncrossedfactorlevels[israndomtype] # responsible for 3 allocations
            effectsdenominators[iother] .*= [prod(@view ncrossedfactorlevels[Not(icrossed[israndomtype] ∩ x)]) for x ∈ interactionindexes] # 7kb - set intersection is 1kb, has to be done for each interaction
            effectsdenominators
        end
        effectsdenominators = f(imax, ncrossedfactors, ncrossedfactorlevels, npercrossedcell, crossedfactortypes, interactionindexes)

        σ² = factors .* differences[1:imax] ./ effectsdenominators

        if nnestedfactors > 0
            nestedrange = (length(factorresults) .- nnestedfactors .+ 1):length(factorresults)
            nestedeffectdenominators = repeat([nreplicates], nnestedfactors)
            for i ∈ 1:(nnestedfactors - 1)
                nestedeffectdenominators[1:(end - i + 1)] .*= nnestedfactorlevels[end - i + 2]
            end
            σ²nested = differences[nestedrange] ./ nestedeffectdenominators
            σ² = [σ²; σ²nested]
        end

        σ²total = sum(σ²) + factorerrorvars[end].ms
        ω² = σ² ./ σ²total
    end
    return ω²
end

function calculate_effect_numerator_factors(imax, factorresults, ncrossedfactors, crossedfactortypes, interactionindexes)
    icrossed = 1:ncrossedfactors
    iother = (ncrossedfactors + 1):imax
    factors = ones(Int, imax)
    for i ∈ icrossed
        isfixed(crossedfactortypes[i]) && (factors[i] = factorresults[i].df)
    end
    for i ∈ eachindex(interactionindexes)
        factors[iother[i]] = prod(@view factors[interactionindexes[i]])
    end
    factors
end
