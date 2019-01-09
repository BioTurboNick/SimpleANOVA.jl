@enum FactorType fixed random nested replicates

Broadcast.broadcastable(a::FactorType) = (a,) # workaround for current behavior
