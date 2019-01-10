"""
    AnovaResult <: AnovaEffect

A set of values for an Anova factor which has been tested

`name` - the name of this factor

`ss` - sum of squares

`df` - degrees of freedom

`ms` - mean square

`f` - the F statistic

`p` - the probability of a Type I error (imcorrect rejection of null hypothesis)
"""
struct AnovaResult <: AnovaEffect
    name::AbstractString
    ss::Float64
    df::Float64
    ms::Float64
    f::Float64
    p::Float64
end

AnovaResult(factor, f, p) = AnovaResult(factor.name, factor.ss, factor.df, factor.ms, f, p)

import Base.show
function show(io::IO, x::AnovaResult)
    println(io, "$(x.name)    $(round(x.ss, sigdigits = 5))    $(round(x.df, sigdigits = 2))    $(round(x.ms, sigdigits = 5))    $(round(x.f, sigdigits = 3))    $(round(x.p, sigdigits = 3))")
end
