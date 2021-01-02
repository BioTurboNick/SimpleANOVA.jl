"""
    AnovaFactor <: AnovaEffect

A set of values for an Anova effect which has a mean square.

`name` - the name of this factor

`ss` - sum of squares

`df` - degrees of freedom

`ms` - mean square
"""
struct AnovaFactor <: AnovaEffect
    name::String
    ss::Float64
    df::Float64
    ms::Float64
end

AnovaFactor(name, ss, df) = AnovaFactor(name, ss, df, ss / df)
AnovaFactor(name, av::AnovaEffect) = AnovaFactor(name, av.ss, av.df)
AnovaFactor(av::AnovaEffect) = AnovaFactor(av.name, av.ss, av.df)

import Base.isapprox
isapprox(x::AnovaFactor, y::AnovaFactor) =
    x.name == y.name &&
    x.ss ≈ y.ss &&
    x.df == y.df &&
    (isnan(x.ms) && isnan(y.ms) || x.ms ≈ y.ms)

import Base.zero
zero(::Type{AnovaFactor}) = AnovaFactor("", 0, 1, 0)
