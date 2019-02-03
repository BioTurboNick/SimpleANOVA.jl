"""
    AnovaPosthocData

Container for the complete results of a posthoc analysis following ANOVA.
"""
mutable struct AnovaPosthocData
    anova::AnovaData
    factorcomparisons::Vector{AnovaPosthocFactor}
end

Broadcast.broadcastable(a::T) where {T <: AnovaPosthocData} = (a,) # workaround for current behavior

import Base.show
function show(io::IO, apd::AnovaPosthocData)
    colnames = ["Levels", "Δ", "DF", "SE", "q", "p"]

    println(io)
    println(io, "Posthoc Comparison Results")
    println(io)

    for factorcomparisons ∈ apd.factorcomparisons
        println(factorcomparisons.name)
        rownames = [join(c.levels, "×") for c ∈ factorcomparisons.comparisons]
        nrows = length(factorcomparisons.comparisons)

        compactshow(x) = sprint(show, x, context=:compact=>true)

        differences = [c.difference |> compactshow for c ∈ factorcomparisons.comparisons]
        df = [c.df |> Int |> compactshow for c ∈ factorcomparisons.comparisons]
        se = [c.se |> compactshow for c ∈ factorcomparisons.comparisons]
        q = [c.q |> compactshow for c ∈ factorcomparisons.comparisons]
        p = [c.p |> compactshow for c ∈ factorcomparisons.comparisons]

        ndec = [length.(last.(split.(values, "."))) for values ∈ [differences, se, q, p]]
        maxndec = maximum.(ndec)
        rpadding = [maxndec[i] .- ndec[i] for i ∈ 1:4]

        differences = [rpad(differences[i], rpadding[1][i] + length(differences[i])) for i ∈ 1:nrows]
        se = [rpad(se[i], rpadding[2][i] + length(se[i])) for i ∈ 1:nrows]
        q = [rpad(q[i], rpadding[3][i] + length(q[i])) for i ∈ 1:nrows]
        p = [rpad(p[i], rpadding[4][i] + length(p[i])) for i ∈ 1:nrows]

        colwidths = [length.(values) |> maximum for values ∈ [[colnames[1]; rownames],
                                                              [colnames[2]; differences],
                                                              [colnames[3]; df],
                                                              [colnames[4]; se],
                                                              [colnames[5]; q],
                                                              [colnames[6]; p]]]
        colwidths[2:end] .+= 2

        headerrow = join(lpad.(colnames, colwidths))
        separator = repeat("-", sum(colwidths))
        rows = [lpad(rownames[i], colwidths[1]) *
                lpad(differences[i], colwidths[2]) *
                lpad(df[i], colwidths[3]) *
                lpad(se[i], colwidths[4]) *
                lpad(q[i], colwidths[5]) *
                lpad(p[i], colwidths[6]) for i ∈ 1:nrows]

        println(io, headerrow)
        println(io, separator)
        foreach(i -> println(io, rows[i]), 1:length(rows))
    end
end

export AnovaPosthocData
