"""
    AnovaValue <: AnovaEffect

A set of values for an Anova item for which a mean square is not required.

`name` - the name of this value

`ss` - sum of squares

`df` - degrees of freedom
"""
struct AnovaValue <: AnovaEffect
    name::AbstractString
    ss::Float64
    df::Float64
end

AnovaValue(ss, df) = AnovaValue("", ss, df)

import Base.show
function show(io::IO, x::AnovaValue)
    println(io, "$(x.name)    $(round(x.ss, sigdigits = 5))    $(round(x.df, sigdigits = 2))")
end
