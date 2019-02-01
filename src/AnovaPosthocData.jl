"""
    AnovaPosthocData

Container for the complete results of a posthoc analysis following ANOVA.
"""
mutable struct AnovaPosthocData
    anova::AnovaData
    comparisons::Vector{AnovaPosthocFactor}
end

Broadcast.broadcastable(a::T) where {T <: AnovaPosthocData} = (a,) # workaround for current behavior

import Base.show
function show(io::IO, apd::AnovaPosthocData)
    colnames = ["Effect", "Δ", "DF", "SE", "q", "p"]
    for factor ∈ apd.comparisons
        println(factor.name)
        for comparison ∈ factor.comparisons
            println(comparison)
        end
    end
    #=
    rownames = [e.name for e ∈ ad.effects]
    nrows = length(ad.effects)

    compactshow(x) = sprint(show, x, context=:compact=>true)

    ss = [e.ss |> compactshow for e ∈ ad.effects]
    df = [e.df |> Int |> compactshow for e ∈ ad.effects]
    ms = [typeof(e) ∈ [AnovaFactor, AnovaResult] ? e.ms |> compactshow  : "" for e ∈ ad.effects]
    f = [e isa AnovaResult ? e.f |> compactshow : "" for e ∈ ad.effects]
    p = [e isa AnovaResult ? e.p |> compactshow : "" for e ∈ ad.effects]

    ndec = [length.(last.(split.(values, "."))) for values ∈ [ss, ms, f, p]]
    maxndec = maximum.(ndec)
    rpadding = [maxndec[i] .- ndec[i] for i ∈ 1:4]

    ss = [rpad(ss[i], rpadding[1][i] + length(ss[i])) for i ∈ 1:nrows]
    ms = [rpad(ms[i], rpadding[2][i] + length(ms[i])) for i ∈ 1:nrows]
    f = [rpad(f[i], rpadding[3][i] + length(f[i])) for i ∈ 1:nrows]
    p = [rpad(p[i], rpadding[4][i] + length(p[i])) for i ∈ 1:nrows]

    colwidths = [length.(values) |> maximum for values ∈ [[colnames[1]; rownames],
                                                          [colnames[2]; ss],
                                                          [colnames[3]; df],
                                                          [colnames[4]; ms],
                                                          [colnames[5]; f],
                                                          [colnames[6]; p]]]
    colwidths[2:end] .+= 2

    headerrow = join(lpad.(colnames, colwidths))
    separator = repeat("-", sum(colwidths))
    rows = [lpad(rownames[i], colwidths[1]) *
            lpad(ss[i], colwidths[2]) *
            lpad(df[i], colwidths[3]) *
            lpad(ms[i], colwidths[4]) *
            lpad(f[i], colwidths[5]) *
            lpad(p[i], colwidths[6]) for i ∈ 1:nrows]

    println(io)
    println(io, "Analysis of Variance Results")
    println(io)
    println(io, headerrow)
    println(io, separator)
    foreach(i -> println(io, rows[i]), 1:length(rows))
    =#
end

export AnovaPosthocData
