"""
    AnovaResult <: AnovaEffect

A set of values for an Anova factor which has been tested

`name` - the name of this factor

`ss` - sum of squares

`df` - degrees of freedom

`ms` - mean square

`f` - the F statistic

`p` - the probability of a Type I error (incorrect rejection of null hypothesis)

`ω²` - the effect size
"""
struct AnovaResult <: AnovaEffect
    name::String
    ss::Float64
    df::Float64
    ms::Float64
    f::Float64
    p::Float64
    ω²::Float64
end

AnovaResult(factor, f, p) = AnovaResult(factor.name, factor.ss, factor.df, factor.ms, f, p, 0)

AnovaResult(result, ω²) = AnovaResult(result.name, result.ss, result.df, result.ms, result.f, result.p, ω²)

import Base.isapprox
isapprox(x::AnovaResult, y::AnovaResult) =
    x.name == y.name &&
    x.ss ≈ y.ss &&
    x.df == y.df &&
    x.ms ≈ y.ms &&
    x.f ≈ y.f &&
    x.p ≈ y.p &&
    x.ω² ≈ y.ω²
