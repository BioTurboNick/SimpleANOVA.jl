"""
    AnovaPosthocComparison

A set of values for a pairwise comparison.

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
    x.difference ≈ y.difference &&
    x.df == y.df &&
    x.se ≈ y.se &&
    x.q ≈ y.q &&
    x.p ≈ y.p

export AnovaPosthocComparison
