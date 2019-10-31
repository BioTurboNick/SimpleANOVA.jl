"""
    AnovaValue <: AnovaEffect

A set of values for an Anova item for which a mean square is not required.

`contrast` - the contrast value

`df` - degrees of freedom

`f` - the F statistic

`p` - the probability of a Type I error (incorrect rejection of null hypothesis)

`r` - the effect size
"""
struct AnovaContrastResult <: AnovaEffect
    contrast::Float64
    df::Float64
    f::Float64
    p::Float64
    r::Float64
end

import Base.isapprox
isapprox(x::AnovaContrastResult, y::AnovaContrastResult) =
    x.contrast ≈ y.contrast &&
    x.df ≈ y.df &&
    x.f ≈ y.f &&
    x.p ≈ y.p &&
    x.r ≈ y.r

import Base.show
show(io::IO, acr::AnovaContrastResult) = show(io, [acr])
function show(io::IO, acr::Vector{AnovaContrastResult})
    colnames = ["Contrast", "DF", "F", "p", "r"]
    nrows = length(acr)

    compactshow(x) = sprint(show, x, context=:compact=>true)

    contrast = [e.contrast |> compactshow for e ∈ acr]
    df = [e.df |> Int |> compactshow for e ∈ acr]
    f = [e.f |> compactshow for e ∈ acr]
    p = [e.p |> compactshow for e ∈ acr]
    r = [e.r |> compactshow for e ∈ acr]

    ndec = [length.(last.(split.(values, "."))) for values ∈ [contrast, f, p, r]]
    maxndec = maximum.(ndec)
    rpadding = [maxndec[i] .- ndec[i] for i ∈ 1:4]

    ss = [rpad(contrast[i], rpadding[1][i] + length(contrast[i])) for i ∈ 1:nrows]
    f = [rpad(f[i], rpadding[2][i] + length(f[i])) for i ∈ 1:nrows]
    p = [rpad(p[i], rpadding[3][i] + length(p[i])) for i ∈ 1:nrows]
    r = [rpad(r[i], rpadding[4][i] + length(r[i])) for i ∈ 1:nrows]

    colwidths = [length.(values) |> maximum for values ∈ [[colnames[1]; contrast],
                                                          [colnames[2]; df],
                                                          [colnames[3]; f],
                                                          [colnames[4]; p],
                                                          [colnames[5]; r]]]
    colwidths[2:end] .+= 2

    headerrow = join(lpad.(colnames, colwidths))
    separator = repeat("-", sum(colwidths))
    rows = [lpad(contrast[i], colwidths[1]) *
            lpad(df[i], colwidths[2]) *
            lpad(f[i], colwidths[3]) *
            lpad(p[i], colwidths[4]) *
            lpad(r[i], colwidths[5]) for i ∈ 1:nrows]

    println(io)
    println(io, "Analysis of Variance Contrast Results")
    println(io)
    println(io, headerrow)
    println(io, separator)
    foreach(i -> println(io, rows[i]), 1:length(rows))
end
