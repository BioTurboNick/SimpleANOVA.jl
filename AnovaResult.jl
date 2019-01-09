include("AnovaEffect.jl")

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
    println(io, "$(x.name)    $(x.ss)    $(x.df)    $(x.ms)    $(x.f)    $(x.p)")
end
