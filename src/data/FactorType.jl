"""
    FactorType

`fixed` - A crossed fixed-effects factor

`random` - A crossed random-effects factor

`nested` - A random factor fully nested within another factor

`subject` `block` - A random factor subjected to multiple levels of another factor
"""
@enum FactorType fixed random nested subject block

Broadcast.broadcastable(a::FactorType) = (a,) # workaround for current behavior
