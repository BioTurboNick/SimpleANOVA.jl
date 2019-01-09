include("AnovaEffect.jl")

"""
    AnovaData

Container for the complete results of an ANOVA test.
"""
mutable struct AnovaData
    effects::Vector{AnovaEffect}
end

import Base.show
function show(io::IO, x::AnovaData)
    println(io, "Analysis of Variance results")
    println(io, "Name    SS    DF    MS")
    for effect ∈ x.effects
        show(io, effect)
    end
end
