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

import Base.isapprox
isapprox(x::AnovaFactor, y::AnovaFactor) =
    x.name == y.name &&
    x.ss ≈ y.ss &&
    x.df == y.df &&
    (isnan(x.ms) && isnan(y.ms) || x.ms ≈ y.ms)

export AnovaFactor
