"""
    AnovaEffect

Supertype for all ANOVA data items: `AnovaValue`, `AnovaFactor`, `AnovaResult`
"""
abstract type AnovaEffect end

Broadcast.broadcastable(a::T) where {T <: AnovaEffect} = (a,) # workaround for current behavior

export AnovaEffect
