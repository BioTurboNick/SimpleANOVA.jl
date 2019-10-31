
pvalue(dist, x) = min(2 * min(cdf(dist, x), ccdf(dist, x)), 1.0)

"""
    contrast(anovaresult::AnovaResult, groupassignment::Vector{Int})

Calculate a single-factor contrast against the groups specified in `groupassignment`.
Valid groups are "0", "1", and "2", where "0" excludes the group from the comparison.

Note: If you do nonorthogonal contrasts, use the Bonferroni or Šidák correction to adjust the
α level (the level at which you consider a result likely to be true):
  Bonferroni: α′ = α/c for c contrasts
  Šidák:      α′ = 1 - (1 - α)^(1 / c) (slightly better)

Note: Effect size is calcluated using the overall error term. Other choices are possible,
including average of each group error; or the error associated with a control.
"""
function contrast(anovaresult::AnovaData, groupassignment::Vector{Int}, factorindex::Int = 1)
    0 < factorindex ≤ anovaresult.ncrossedfactors                            || error("factor index must be a valid factor")
    all(0 .≤ groupassignment .≤ 2)                                           || error("valid groups are 0, 1, and 2")
    length(groupassignment) == anovaresult.ncrossedfactorlevels[factorindex] || error("each level must be assigned to a group")

    lowerfactorlevels = prod(anovaresult.ncrossedfactorlevels[1:(factorindex-1)])
    upperfactorlevels = prod(anovaresult.ncrossedfactorlevels[(factorindex+1):anovaresult.ncrossedfactors])

    # factorindex is in the order the factors appear in the result
    # but the crossedcellmeans dimensions are in reverse order (as input)
    group1levels = repeat(groupassignment .== 1, outer = upperfactorlevels, inner = lowerfactorlevels)
    group2levels = repeat(groupassignment .== 2, outer = upperfactorlevels, inner = lowerfactorlevels)
    group1count = count(group1levels)
    group2count = count(group2levels)

    contrastcoefficients = zeros(anovaresult.ncrossedfactorlevels[factorindex] * lowerfactorlevels * upperfactorlevels)
    contrastcoefficients[group1levels] .= 1 / group1count
    contrastcoefficients[group2levels] .= -1 / group2count

    errorfactor = anovaresult.effects[end]
    ψ = sum(contrastcoefficients .* anovaresult.crossedcellmeans)
    contrast = anovaresult.npercrossedcell * ψ ^ 2
    error = sum(contrastcoefficients .^ 2) * errorfactor.ms

    f = contrast / error
    fdist = FDist(1, errorfactor.df)
    p = ccdf(fdist, f)

    effectsize = sqrt(f / (f  + errorfactor.df)) # t^2 == f
    #cohen's d: effectsize = ψ / sqrt(errorfactor.ms) ; preferred when group sizes very different

    AnovaContrastResult(contrast, errorfactor.df, f, p, effectsize)
end



helmertcontrasts(anovaresult, factorindex = 1) = differencecontrast(anovaresult, factorindex, true)

"""
    differencecontrast(anovaresult::AnovaData, factorindex = 1, reverseorder = false)

Compute orthogonal contrasts on the factor levels in the original order. Forward direction also known as
a "Helmert" contrast; revere direction may also be called "Difference" (as in SPSS). See `contrast`
function for more.
"""
function differencecontrasts(anovaresult::AnovaData, factorindex = 1, reverseorder = false)
    0 < factorindex ≤ anovaresult.ncrossedfactors || error("factor index must be a valid factor")

    levels = anovaresult.ncrossedfactorlevels[1]
    groupassignments = [[repeat([0], i - 1); 1; repeat([2], levels - i)] for i ∈ 1:(levels - 1)]
    reverseorder && reverse!.(groupassignments)

    AnovaContrastResults(contrast.(Ref(anovaresult), groupassignments, factorindex))
end

"""
    repeatedcontrast(anovaresult::AnovaData, factorindex = 1)

Compute contrasts between neighboring levels. Non-orthogonal. See `contrast` function for more.
"""
function repeatedcontrasts(anovaresult::AnovaData, factorindex = 1)
    0 < factorindex ≤ anovaresult.ncrossedfactors || error("factor index must be a valid factor")

    levels = anovaresult.ncrossedfactorlevels[factorindex]
    groupassignments = [[repeat([0], i - 1); 1; 2; repeat([0], levels - i - 1)] for i ∈ 1:(levels - 1)]

    AnovaContrastResults(contrast.(Ref(anovaresult), groupassignments, factorindex))
end

"""
    simplecontrast(anovaresult::AnovaData, factorindex = 1, controlindex = 1)

Compute contrasts of each level to a single level (control). Non-orthogonal. See `contrast` function for more.
"""
function simplecontrasts(anovaresult::AnovaData, factorindex = 1, controlindex = 1)
    0 < factorindex ≤ anovaresult.ncrossedfactors || error("factor index must be a valid factor")
    0 < controlindex ≤ anovaresult.ncrossedfactorlevels[factorindex] || error("index must be for a valid factor level")

    levels = anovaresult.ncrossedfactorlevels[factorindex]
    groupassignments = [zeros(Int, levels) for i ∈ 1:(levels - 1)]
    otherlevels = (1:levels)[Not(controlindex)]
    for i ∈ 1:(levels - 1)
        groupassignments[i][controlindex] = 1
        groupassignments[i][otherlevels[i]] = 2
    end

    AnovaContrastResults(contrast.(Ref(anovaresult), groupassignments), factorindex)
end

#= To do this one, need to be able to code a factor level as in both groups
"""
    deviationcontrasts(anovaresult::AnovaData, factorindex = 1, controlindex = 1)

Compute contrasts of each level except one (control) to all levels. Non-orthogonal.
See `contrast` function for more.
"""
function deviationcontrast(anovaresult::AnovaData, factorindex = 1, controlindex = 1)
    0 < factorindex ≤ anovaresult.ncrossedfactors || error("factor index must be a valid factor")
    0 < controlindex ≤ anovaresult.ncrossedfactorlevels[1] || error("index must be for a valid factor level")

    levels = anovaresult.ncrossedfactorlevels[1]
    groupassignments = [zeros(Int, levels) for i ∈ 1:(levels - 1)]
    otherlevels = (1:levels)[Not(controlindex)]
    for i ∈ 1:(levels - 1)
        groupassignments[i][controlindex] = 1
        groupassignments[i][otherlevels[i]] = 2
    end

    [groupassignments[i][i] = 1 for i ∈ (1:levels)[Not(controlindex)]]
    [groupassignments[i][i] = 1 for i ∈ (1:levels)[Not(controlindex)]]

    contrast.(Ref(anovaresult), groupassignments)
end
=#
