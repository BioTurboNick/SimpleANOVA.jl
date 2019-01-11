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

    colnames = ["Effect", "SS", "DF", "MS", "F", "p"]
    rownames = [e.name for e ∈ ad.effects]
    nrows = length(ad.effects)

    compactshow(x) = sprint(show, x, context=:compact=>true)

    ss = [e.ss |> compactshow for e ∈ ad.effects]
    df = [e.df |> Int |> compactshow for e ∈ ad.effects]
    ms = [typeof(e) ∈ [AnovaFactor, AnovaResult] ? e.ms |> compactshow  : "" for e ∈ ad.effects]
    f = [e isa AnovaResult ? e.f |> compactshow : "" for e ∈ ad.effects]
    p = [e isa AnovaResult ? e.p |> compactshow : "" for e ∈ ad.effects]

    columnwidths = [length.(values) |> maximum for values ∈ [rownames, ss, df, ms, f, p]] .+ 2

    headerrow = join(lpad.(colnames, columnwidths))
    rows = [lpad(rownames[i], columnwidths[1]) *
            lpad(ss[i], columnwidths[2]) *
            lpad(df[i], columnwidths[3]) *
            lpad(ms[i], columnwidths[4]) *
            lpad(f[i], columnwidths[5]) *
            lpad(p[i], columnwidths[6]) for i ∈ 1:nrows]

    println(io, headerrow)
    foreach(i -> println(io, rows[i]), 1:length(rows))
end
