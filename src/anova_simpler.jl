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
anova(observations, [fixed, subject])      # N-way repeated measures ANOVA with 1 within-subjects fixed factor
anova(observations, [fixed, block])        # N-way repeated measures ANOVA with 1 within-block fixed factor
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
function anova1(observations::AbstractArray{T}, factortypes::Vector{FactorType} = FactorType[]; factornames::Vector{<:AbstractString} = String[], hasreplicates = true) where {T <: Union{Number, AbstractVector{<:Number}}}
    length(observations) > 0 || return

    isrepeatedmeasures = subject ∈ factortypes

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
        factornames = string.('A':'Z')[1:nfactors]
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

    crossedfactornames = factornames[factortypes .!= nested]
    nestedfactornames = factornames[factortypes .== nested]

    if isrepeatedmeasures
        anovakernel(observations, nreplicates, ncells, nnestedfactors, ncrossedfactors, nfactorlevels, crossedfactortypes, crossedfactornames, nestedfactornames, true)
    else
        anovakernel(observations, nreplicates, ncells, nnestedfactors, ncrossedfactors, nfactorlevels, crossedfactortypes, crossedfactornames, nestedfactornames, false)
    end
end

#=
function anova(data::AnovaData, crossedfactors::Vector{Int}, )
    # performs a subtest of the specified crossed factors within level of the remaing crossed factors, using the original errors
end

# possible bug: function anova(data::AnovaData, crossedfactors::Vector{Int}, ) with normal functions ==> hang when called?
=#

function anova1(observations::AbstractVector{T}, factorassignments::AbstractVector{<:AbstractVector}, factortypes::Vector{FactorType} = FactorType[]; factornames::Vector{<:AbstractString} = String[]) where {T <: Number}
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
        subject ∉ factortypes || (count(t -> t == subject, factortypes) == 1 && (factortypes[2] == subject || factortypes[3] == subject)) || error("maximum of one subject entry and must be second or third")
    end
    if !isempty(factornames)
        nfactors ≤ 26 || error("Can only automatically name up to 26 factors. Provide names explicitly.")
        length(factornames) == nfactors || error("factornames must have an entry for each factor.")
    end
end

function anovakernel(observations, nreplicates, ncells, nnestedfactors, ncrossedfactors, nfactorlevels, crossedfactortypes, crossedfactornames, nestedfactornames, isrepeatedmeasures)
    N = ncells * nreplicates
    nfactors = nnestedfactors + ncrossedfactors

    
end

anovavalue(name, variance, df) = AnovaValue(name, variance * df, df)
anovafactor(name, variance, df) = AnovaFactor(anovavalue(name, variance, df))
meanfirstdim(observations::AbstractArray{<:Number}) = dropdims(mean(observations, dims = 1), dims = 1)
upcat(x::AbstractArray) = x
upcat(x::AbstractArray{<:AbstractVector}) = reshape(vcat(x...), (size(x[1])..., size(x)...))

function ftest(x, y)
    f = x.ms / y.ms
    fdist = FDist(x.df, y.df)
    p = ccdf(fdist, f)
    AnovaResult(x, f, p)
end

mutable struct AnovaData2
    effects::Vector{AnovaEffect}
end

function anova1(observations, factornames, factortypes, isrepeatedmeasures)
    observations = upcat(observations)

    nreplicates = size(observations, 1)

    totaldf = length(observations) - 1
    totalvar = anovavalue(totalname, var(observations), totaldf)

    cellmeans = meanfirstdim(observations)
    cellsdf = length(cellmeans) - 1
    cellsvar = anovavalue(cellsname, var(cellmeans) * nreplicates, cellsdf)

    errorvar = AnovaFactor(errorname, totalvar - cellsvar)

    if isrepeatedmeasures
        subjectindex = findfirst(x -> x == subject, factortypes)
        nfactors = ndims(observations)
        otherfactors = (1:nfactors)[Not(subjectindex)]
        notherfactorlevels = collect(nfactorlevels)[otherfactors]
        subjectsvar = AnovaValue("Subjects", sum((mean(observations, dims=otherfactors) .- mean(observations)) .^ 2) * prod(notherfactorlevels), nfactorlevels[subjectindex] - 1)
        withinsubjectsvar = totalvar - subjectsvar
        cellmeans = dropdims(mean(observations, dims = subjectindex), dims = subjectindex)
        nreplicates *= nfactorlevels[subjectindex]

        factorvars, factorerrorvars = anovakernel(cellmeans, nreplicates, withinsubjectsvar, errorvar, factornames, factortypes)

        # subject interactions UNGENERALIZED
        sa = AnovaFactor("S/A", sum((mean(observations, dims=2) .- mean(observations)) .^ 2) * nfactorlevels[3] - factorvars[1].ss - subjectsvar.ss, subjectsvar.df * factorvars[1].df)
        sb = AnovaFactor("S/B", sum((mean(observations, dims=3) .- mean(observations)) .^ 2) * nfactorlevels[2] - factorvars[2].ss - subjectsvar.ss, subjectsvar.df * factorvars[2].df)
        sab = withinsubjectsvar - sum(factorvars) - sa - sb

        # can I have anovakernel do the subject interactions?

        factorerrorvar = subjectsvar - sum(factorvars)
        replace!(factorerrorvars, errorvar => factorerrorvar)
        
        

        remaindervar = withinsubjectsvar - sum(factorvars)
        replace!(withinfactorerrorvars, errorvar => remaindervar)
    else
        amongallnestedvars, cellmeans, nnestedfactorlevels, factornames, factortypes = amongnestedfactorscalc!(cellmeans, factornames, factortypes)
        nreplicates *= prod(nnestedfactorlevels)

        factorvars, factorerrorvars = anovakernel(cellmeans, nreplicates, cellsvar, errorvar, factornames, factortypes)

        nestedvars = nestedfactorscalc(amongallnestedvars, sum(factorvars))
        replace!(factorerrorvars, errorvar => nestedvars[1])
    end

    factorresults = ftest.(factorvars, factorerrorvars)

    if isrepeatedmeasures
        withinfactorresults = ftest.(withinfactorvars, withinfactorerrorvars)
        append!(factorresults, withinfactorresults)
    else
        nestedresults = ftest.(nestedvars, [nestedvars[2:end]; errorvar])
        append!(factorresults, nestedresults)
    end

    return AnovaData2([totalvar; factorresults; errorvar])
end

#=
For one-way top level: factorA error = A within Subjects
For two-way top level: factorA/B error = AxB within Subjects

For within factors and interaction with top level: interaction between the factor and A within Subjects
For interaction among within factors: interaction between the lowest factor and A within Subjects (doesn't seem right?)
For three-way interaction: interaction bteween B, C, and A within Subjcects


Subjects within A = Subjects - Factor A
=#




#=


Generally speaking, a 2-way interaction is calculated as:

ss = sum((mean(observations, dims=otherfactors) .- mean(observations)) .^ 2)) * prod(nfactorlevels[otherfactors]) * nreplicates - factorAss - factorBss

and has df = factorAdf * factorBdf

a 3-way interaction is calculated as:

ss = sum((mean(observations, dims=otherfactors) .- mean(observations)) .^ 2) * prod(nfactorlevels[otherfactors]) * nreplicates - sum(factorsABC) - sum(pairwiseinteractionsABC)

and has df = factorAdf * factorBdf * factorCdf

=#



makefactorname(factorname::AbstractString) = factorname

function makefactorname(factornames::AbstractVector{<:AbstractString})
    length(factornames) > 0 || return ""
    
    factorname = factornames[1]
    
    for i ∈ 2:length(factornames)
        factorname *= " × " * factornames[i]
    end

    return factorname
end


# nested levels
function amongnestedfactorscalc(cellmeans, factornames, factortypes)
    nnestedfactors = count(x -> x == nested, factortypes)
    nlowerfactorlevels = (1,)
    nupperfactorlevels = size(cellmeans)
    amongallnestedvars = Vector{AnovaFactor}(undef, nnestedfactors)
    @views for i ∈ 1:nnestedfactors
        # collapse each nested factor
        df = prod(nupperfactorlevels) - 1
        ss = var(cellmeans) * nreplicates * df * prod(nlowerfactorlevels)
        amongallnestedvars[i] = AnovaFactor(factornames[i], ss, df)
        nlowerfactorlevels = (nlowerfactorlevels..., nupperfactorlevels[1])
        nupperfactorlevels = nupperfactorlevels[2:end]
        cellmeans = meanfirstdim(cellmeans)
    end
    return amongallnestedvars, cellmeans, collect(nlowerfactorlevels), factornames[nnestedfactors + 1:end], factortypes[nnestedfactors + 1:end]
end

function nestedfactorscalc(amongallnestedvars, otherfactorvars)
    nestedfactors = []
    nnestedfactors = length(amongallnestedvars)
    if nnestedfactors > 0
        for i ∈ nnestedfactors:-1:1
            nestedvar = amongallnestedvars[i] - otherfactorvars
            otherfactorvars += nestedvar
            push!(nestedfactors, AnovaFactor(amongallnestedvars[i].name, nestedvar))
        end
    end
    return nestedfactors
end

function anovafactors(cellmeans, nreplicates, factornames)
    # calculate all simple factors and all interactions

    N = length(cellmeans) * nreplicates

    factors = 1:ndims(cellmeans)

    allfactors = []
    allfactorvars = AnovaFactor[]
    for i ∈ factors
        ifactors = collect(combinations(factors, i))
        iotherfactors = [factors[Not(i...)] for i ∈ ifactors]
        iupperfactorvars = [allfactorvars[findall(x -> x ⊆ i, allfactors)] for i ∈ ifactors]

        ifactornames = [makefactorname(reverse(factornames[i])) for i ∈ ifactors]
        ifactorss = [var(mean(cellmeans, dims = iotherfactors[i]), corrected = false) * N - sum(iupperfactorvars[i]).ss for i ∈ eachindex(iotherfactors)]
        ifactordf = isempty(iupperfactorvars[1]) ? [size(cellmeans, i) for i ∈ factors] :
                                                   [prod(f.df for f ∈ iupperfactorvars[i]) for i ∈ eachindex(iotherfactors)]
        ifactorvars = AnovaFactor.(ifactornames, ifactorss, ifactordf)
        
        append!(allfactors, reverse!(ifactors))
        append!(allfactorvars, reverse!(ifactorvars))
    end

    return allfactorvars
end

# one-way
function anovakernel(cellmeans::AbstractVector, nreplicates, _1, errorvar, factornames, _2)
    factordf = length(cellmeans) - 1
    factorvar = anovafactor(factornames[1], var(cellmeans) * nreplicates, factordf)
    return [factorvar], [errorvar]
end

# two-way
function anovakernel(cellmeans::AbstractMatrix, nreplicates, _, errorvar, factornames, factortypes)
    N = length(cellmeans) * nreplicates

    factor1var, factor2var, 

    if factortypes[1] == factortypes[2]
        if factortypes[1] == fixed
            factor1error = errorvar
            factor2error = errorvar
        else
            factor1error = interaction12var
            factor2error = interaction12var
        end
    else
        if factortypes[1] == fixed
            factor1error = interaction12var
            factor2error = errorvar
        else
            factor1error = errorvar
            factor2error = interaction12var
        end
    end

    return [factor2var; factor1var; interaction12var],
           [factor2error; factor1error; interaction12error]
end

# three-way
function anovakernel(cellmeans::AbstractArray{T, 3}, nreplicates, cellsvar, errorvar, factornames, factortypes) where T
    N = length(cellmeans) * nreplicates

    factor1df = size(cellmeans, 1) - 1
    factor1var = AnovaFactor(factornames[1], var(mean(cellmeans, dims=(2,3)), corrected = false) * N, factor1df)
    
    factor2df = size(cellmeans, 2) - 1
    factor2var = AnovaFactor(factornames[2], var(mean(cellmeans, dims=(1,3)), corrected = false) * N, factor2df)

    factor3df = size(cellmeans, 3) - 1
    factor3var = AnovaFactor(factornames[3], var(mean(cellmeans, dims=(1,2)), corrected = false) * N, factor3df)

    interaction12df = factor1df * factor2df
    interaction12var = AnovaFactor("$(factornames[2]) × $(factornames[1])", var(mean(cellmeans, dims = 3), corrected = false) * N - factor1var.ss - factor2var.ss, interaction12df)

    interaction13df = factor1df * factor3df
    interaction13var = AnovaFactor("$(factornames[3]) × $(factornames[1])", var(mean(cellmeans, dims = 2), corrected = false) * N - factor1var.ss - factor3var.ss, interaction13df)

    interaction23df = factor2df * factor3df
    interaction23var = AnovaFactor("$(factornames[3]) × $(factornames[2])", var(mean(cellmeans, dims = 1), corrected = false) * N - factor2var.ss - factor3var.ss, interaction23df)

    interaction123df = factor1df * factor2df * factor3df
    interaction123var = AnovaFactor("$(factornames[3]) × $(factornames[2]) × $(factornames[1])", cellsvar.ss - factor1var.ss - factor2var.ss - factor3var.ss - interaction12var.ss - interaction13var.ss - interaction23var.ss, interaction123df)

    if factortypes[1] == factortypes[2] == factortypes[3]
        if factortypes[1] == fixed
            factor1error = factor2error = factor3error = errorvar
            interaction12error = interaction13error = interaction23error = errorvar
        else
            factor1error = threeway_random_error(interaction12var, interaction13var, interaction123var)
            factor2error = threeway_random_error(interaction12var, interaction23var, interaction123var)
            factor3error = threeway_random_error(interaction13var, interaction23var, interaction123var)
            interaction12error = interaction13error = interaction23error = interaction123var
        end
    elseif factortypes[1] == factortypes[2]
        if factortypes[1] == fixed
            factor1error = interaction13var
            factor2error = interaction23var
            factor3error = errorvar
            interaction12error = interaction123var
            interaction13error = interaction23error = errorvar
        else
            factor1error = factor2error = interaction12var
            factor3error = threeway_random_error(interaction13var, interaction23var, interaction123var)
            interaction12error = errorvar
            interaction13error = interaction23error = interaction123var
        end
    elseif factortypes[1] == factortypes[3]
        if factortypes[1] == fixed
            factor1error = interaction12var
            factor2error = errorvar
            factor3error = interaction23var
            interaction12error = interaction23error = errorvar
            interaction13error = interaction123var
        else
            factor1error = factor3error = interaction13var
            factor2error = threeway_random_error(interaction12var, interaction23var, interaction123var)
            interaction12error = interaction23error = interaction123var
            interaction13error = errorvar
        end
    else
        if factortypes[2] == fixed
            factor1error = errorvar
            factor2error = interaction12var
            factor3error = interaction13var
            interaction12error = interaction13error = errorvar
            interaction23error = interaction123var
        else
            factor1error = threeway_random_error(interaction12var, interaction13var, interaction123var)
            factor2error = factor3error = interaction23var
            interaction12error = interaction13error = interaction123var
            interaction23error = errorvar
        end
    end

    interaction123error = errorvar

    return [factor3var; factor2var; factor1var; interaction23var; interaction13var; interaction12var; interaction123var],
           [factor3error; factor2error; factor1error; interaction23error; interaction13error; interaction12error; interaction123error]
end

function threeway_random_error(interaction_ab, interaction_bc, interaction_abc)
    reducedmeansquare(factor::AnovaFactor) = factor.ms ^ 2 / factor.df
    ms = interaction_ab.ms + interaction_bc.ms - interaction_abc.ms
    df = ms ^ 2 / (reducedmeansquare(interaction_ab) + reducedmeansquare(interaction_bc) + reducedmeansquare(interaction_abc))
    AnovaFactor("", ms * df, df, ms)
end


