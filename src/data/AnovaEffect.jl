"""
    AnovaEffect

Supertype for all ANOVA data items: `AnovaValue`, `AnovaFactor`, `AnovaResult`
"""
abstract type AnovaEffect end

import Base.-, Base.+
-(x::AnovaEffect, y::AnovaEffect) = AnovaValue(x.name, x.ss - y.ss, x.df - y.df)
+(x::AnovaEffect, y::AnovaEffect) = AnovaValue(x.name, x.ss + y.ss, x.df + y.df)

Broadcast.broadcastable(a::T) where {T <: AnovaEffect} = (a,) # workaround for current behavior
