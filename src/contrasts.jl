
pvalue(dist, x) = min(2 * min(cdf(dist, x), ccdf(dist, x)), 1.0)

function contrast(anovaresult::AnovaResult, groupassignment)
    # for 1 factor only
    group1larger = length(group1levels) ≥ length(group2levels)

    unitcontrast =

    weightedmeans = sum(contrastcoeffs_0vs10_20_40 .* anovaresult.cellmeans)
    error = sqrt(anovaresult.effects[end].ms * sum(contrastcoeffs_0vs10_20_40 .^ 2 ./ 10))
    weightedmeans / error
end



helmertcontrast(args...) = differencecontrast(args...)

# SPSS uses Difference to just mean reverse order
function differencecontrast(anovaresult::AnovaResult, levelorder = [])


end
