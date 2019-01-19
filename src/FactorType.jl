"""
    FactorType

`fixed` - A crossed fixed-effects factor

`random` - A crossed random-effects factor

`nested` - A random factor fully nested within another factor
"""
@enum FactorType fixed random nested

Broadcast.broadcastable(a::FactorType) = (a,) # workaround for current behavior

export FactorType, fixed, random, nested
