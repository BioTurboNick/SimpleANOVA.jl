"""
    AnovaData

Container for the complete results of an ANOVA test.
"""
mutable struct AnovaData
    effects::Vector{AnovaEffect}
end

import Base.show
function show(io::IO, ad::AnovaData)
    println(io, "Analysis of Variance results")

    colnames = ["SS", "DF", "MS", "F", "p"]
    rownames = [e.name for e ∈ ad.effects]

    rownamewidth = max(7, maximum([length(name) for name ∈ rownames])) + 1
    columnwidth = max(12, maximum([length(name) for name ∈ colnames]))

    paddedrownames = [rpad(name, rownamewidth) for name ∈ rownames]

    compactshow(x) = lpad(sprint(show, x, context=:compact=>true), columnwidth)

    headerrow = rpad("Effect", rownamewidth) * join(lpad.(colnames, columnwidth))

    rows = [compactshow(e.ss) * compactshow(e.df) *
            (typeof(e) ∈ [AnovaFactor, AnovaResult] ?
                compactshow(e.ms) : lpad("", columnwidth)) *
            (e isa AnovaResult ?
                compactshow(e.f) * compactshow(e.p) : lpad("", columnwidth) * lpad("", columnwidth)) for e ∈ ad.effects]

   println(io, headerrow)
   foreach(i -> println(io, paddedrownames[i] * rows[i]), 1:length(rows))
end
