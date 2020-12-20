"""
    FactorType

`fixed` - A crossed fixed-effects factor

`random` - A crossed random-effects factor

`nested` `subject` `block` - A random factor fully nested within another factor

`within` - A fixed-effects factor applied across a nested factor
"""
@enum FactorType fixed random nested subject block within

Broadcast.broadcastable(a::FactorType) = (a,) # workaround for current behavior
