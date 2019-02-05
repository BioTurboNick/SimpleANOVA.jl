"""
    AnovaPosthocFactor

Contains all comparisons between levels of one factor.
"""
struct AnovaPosthocFactor
    name::String
    comparisons::Vector{AnovaPosthocComparison}
end

Broadcast.broadcastable(a::T) where {T <: AnovaPosthocFactor} = (a,) # workaround for current behavior

import Base.isapprox
isapprox(x::AnovaPosthocFactor, y::AnovaPosthocFactor) =
    x.name == y.name &&
    all(x.comparisons .≈ y.comparisons)

export AnovaPosthocFactor
