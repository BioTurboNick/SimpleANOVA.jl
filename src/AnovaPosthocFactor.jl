"""
    AnovaPosthocFactor

A set of values for an Anova effect which has a mean square.

`name` - the name of this factor

`ss` - sum of squares

`df` - degrees of freedom

`ms` - mean square
"""
struct AnovaPosthocFactor
    name::String
    comparisons::Vector{AnovaPosthocComparison}
end

Broadcast.broadcastable(a::T) where {T <: AnovaPosthocFactor} = (a,) # workaround for current behavior

export AnovaPosthocFactor
