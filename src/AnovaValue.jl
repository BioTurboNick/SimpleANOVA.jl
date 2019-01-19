"""
    AnovaValue <: AnovaEffect

A set of values for an Anova item for which a mean square is not required.

`name` - the name of this value

`ss` - sum of squares

`df` - degrees of freedom
"""
struct AnovaValue <: AnovaEffect
    name::String
    ss::Float64
    df::Float64
end

AnovaValue(ss, df) = AnovaValue("", ss, df)

import Base.isapprox
isapprox(x::AnovaValue, y::AnovaValue) =
    x.name == y.name &&
    x.ss ≈ y.ss &&
    x.df == y.df

export AnovaValue
