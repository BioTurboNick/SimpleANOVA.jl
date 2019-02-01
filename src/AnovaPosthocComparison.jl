"""
    AnovaPosthocComparison

A set of values for an Anova effect which has a mean square.

`levels` - the levels being compared

`difference` - the difference between the level means

`df` - degrees of freedom

`se` - standard error

`q` - Studentized Range test statistic
"""
struct AnovaPosthocComparison
    levels::Tuple{Int, Int}
    difference::Float64
    df::Float64
    se::Float64
    q::Float64
    p::Float64
end

import Base.isapprox
isapprox(x::AnovaPosthocComparison, y::AnovaPosthocComparison) =
    x.levels == y.levels &&
    x.df == y.df &&
    x.se ≈ y.se &&
    (isnan(x.q) && isnan(y.q) || x.q ≈ y.q) &&
    (isnan(x.p) && isnan(y.p) || x.p ≈ y.p)

export AnovaPosthocComparison
