"""
    FactorType

`fixed` `manipulated` - A crossed fixed-effects factor

`random` `measured` - A crossed random-effects factor

`nested` - A random factor fully nested within another factor

`subject` `block` - A random factor subjected to multiple levels of another factor
"""
@enum FactorType fixed manipulated random measured nested subject block

Broadcast.broadcastable(a::FactorType) = (a,) # workaround for current behavior

isnested(x) = x == nested
isfixed(x) = x == fixed || x == manipulated
israndom(x) = x == random || x == measured
issubject(x) = x == subject || x == block
