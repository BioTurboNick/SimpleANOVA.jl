"""
    AnovaFactor <: AnovaEffect

A set of values for an Anova effect which has a mean square.

`name` - the name of this factor

`ss` - sum of squares

`df` - degrees of freedom

`ms` - mean square
"""
struct AnovaFactor <: AnovaEffect
    name::AbstractString
    ss::Float64
    df::Float64
    ms::Float64
end

AnovaFactor(name, ss, df) = AnovaFactor(name, ss, df, ss / df)
