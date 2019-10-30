
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
function contrast(anovaresult::AnovaData, groupassignment::Vector{Int})
    length(anovaresult.effects) == 3                               || error("1-factor only")
    all(0 .≤ groupassignment .≤ 2)                                 || error("valid groups are 0, 1, and 2")
    length(groupassignment) == anovaresult.ncrossedfactorlevels[1] || error("each level must be assigned to a group")

    group1levels = groupassignment .== 1
    group2levels = groupassignment .== 2
    group1count = count(group1levels)
    group2count = count(group2levels)

    contrastcoefficients = zeros(anovaresult.ncrossedfactorlevels[1])
    contrastcoefficients[group1levels] .= 1 / group1count
    contrastcoefficients[group2levels] .= -1 / group2count

    errorfactor = anovaresult.effects[end]
    ψ = sum(contrastcoefficients .* anovaresult.crossedcellmeans)
    contrast = anovaresult.npercrossedcell * ψ ^ 2
    error = sum(contrastcoefficients .^ 2) * errorfactor.ms

    f = contrast / error
    fdist = FDist(1, errorfactor.df)
    p = ccdf(fdist, f)

    effectsize = abs(ψ) / sqrt(errorfactor.ms)

    (f, p, effectsize)
end



helmertcontrast(args...) = differencecontrast(args...)

# SPSS uses Difference to just mean reverse order
function differencecontrast(anovaresult::AnovaResult, levelorder = [])


end
