include("AnovaPosthocComparison.jl")
include("AnovaPosthocFactor.jl")
include("AnovaPosthocData.jl")

"""
    tukey(anova::AnovaData, α)
    hsd(anova::AnovaData, α)
    honestlysignificantdifference(anova::AnovaData, α)
    multiplecomparison(anova::AnovaData, α)

Performs the Tukey multiple comparisons posthoc test.

NOTE: While this method calculates the p-values for each comparison, the Tukey procecure requires ignoring p-values from further tests between intermediate pairs once a non-significant one is found.
For example:
If you have 5 levels ranked by mean, and you test 5 vs 1, 5 vs 2, 5 vs 3, and find that 5 vs 3 is non-significant, this is sufficient to conclude that 5 = 4 = 3, and thus do not test 5 vs 4 or 4 vs 3.

If the test concludes that a level is equal to two adjacent levels but those levels are not equal to each other, it is considered to be ambiguous as to which group it belongs to.

As this test is less powerful than ANOVA, it is possible for ANOVA to find a significant effect but this test finds that all levels are equal.
"""
tukey(args...) = multiplecomparison(args...)
hsd(args...) = multiplecomparison(args...)
honestlysignificantdifference(args...) = multiplecomparison(args...)
multiplecomparison(anova::AnovaData) = multiplecomparisonkernel(anova, tukeygroups, "Tukey HSD")

"""
    snk(anova::AnovaData)
    studentnewmankeuls(anova::AnovaData)
    newmankeuls(anova::AnovaData)
    multiplerange(anova::AnovaData)

Performs the (Student-)Newman-Keuls multiple range posthoc test.

Tends to be more powerful than Tukey, but some suggest it is more likely to lead to Type I error than intended.
"""
snk(args...) = newmankeulsmultiplerange(args...)
studentnewmankeuls(args...) = newmankeulsmultiplerange(args...)
newmankeuls(args...) = newmankeulsmultiplerange(args...)
newmankeulsmultiplerange(anova::AnovaData) = multiplecomparisonkernel(anova, snkgroups, "Newman-Keuls")

"""
    wsd(anova::AnovaData)
    whollysignificantdifference(anova::AnovaData)
    multiplecomparisonandrange(anova::AnovaData)

Performs a compromise test between Tukey and (Student-)Newman-Keuls.

Computes the critical values for both tests and uses the mean.

Note: Tukey's test has sometimes been referred to as "wholly significant difference test" as well.
"""
wsd(args...) = multiplecomparisonandrange(args...)
whollysignificantdifference(args...) = multiplecomparisonandrange(args...)
multiplecomparisonandrange(anova::AnovaData) = multiplecomparisonkernel(anova, wsdgroups, "WSD")

function multiplecomparisonkernel(anova::AnovaData, ngroupsfunc, testtype::String)
    nfactors = anova.ncrossedfactors
    nfactorlevels = anova.ncrossedfactorlevels |> reverse
    ngroups = ngroupsfunc(nfactorlevels)

    factormeans = mean.(anova, 1:nfactors) |> reverse!
    df = [f.df for f ∈ anova.crossedfactorsdenominators]
    ms = [f.ms for f ∈ anova.crossedfactorsdenominators]
    se = sqrt.(ms ./ anova.npercrossedcell)
    diffs = [abs.(f .- f') for f ∈ factormeans]
    q = diffs ./ se

    p = [srdistccdf.(df[i], ngroups[i], q[i]) for i ∈ 1:nfactors]

    comparisons = [[AnovaPosthocComparison((j,i), diffs[k][i,j], df[k], se[k], q[k][i,j], p[k][i,j]) for j ∈ 1:nfactorlevels[k]
                                                                                                     for i ∈ (j + 1):nfactorlevels[k]]
                                                                                                     for k ∈ 1:nfactors]
    factorcomparisons = [AnovaPosthocFactor(anova.crossedfactors[i].name, comparisons[i]) for i ∈ 1:nfactors]

    return AnovaPosthocData(anova, factorcomparisons, testtype)
end

tukeygroups(nfactorlevels::Vector{Int}) = nfactorlevels
function snkgroups(nfactorlevels::Vector{Int})
    factorlevels = range.(1, nfactorlevels)
    return [abs.(l .- l') .+ 1 for l ∈ factorlevels]
end
function wsdgroups(nfactorlevels::Vector{Int})
    ngroups1 = snkgroups(nfactorlevels)
    ngroups2 = tukeygroups(nfactorlevels)
    return [(ngroups1[i] .+ ngroups2[i]) ./ 2 for i ∈ 1:length(nfactorlevels)]
end

import Rmath: libRmath
srdistccdf(ν, k, x) = ccall((:ptukey, libRmath), Float64, (Float64, Float64, Float64, Float64, Int, Int), x, 1, k, ν, 0, 0)
srdistinvccdf(ν, k, x) = ccall((:qtukey, libRmath), Float64, (Float64, Float64, Float64, Float64, Int, Int), x, 1, k, ν, 0, 0)







"""
    pooled(anova::AnovaData)

Pools the mean squares for the interaction term(s) and error term and recalculates the test.

This should only be used if no interaction is detected. Increases the power of the test, but the chance of committing a Type  I error is less certain.

Generally, the conservative unpooled approach has been shown to be superior.
"""
function pooled(anova::AnovaData)

end

"""
    dmr(anova::AnovaData)
    duncan(anova::AnovaData)
    duncanmultiplerange(anova::AnovaData)


"""
dmr(args...) = duncanmultiplerange(args...)
duncan(args...) = duncanmultiplerange(args...)
function duncanmultiplerange(anova::AnovaData)
    error("Not implemented")
end

"""
    dunnett(anova::AnovaData)
    controltoall(anova::AnovaData)

Performs the Dunnett method for comparing each factor level to a control group, and not all pairwise comparisons.
"""
dunnett(args...) = controltoall(args...)
function controlcomparison(anova::AnovaData)
    # distribution is complicated and not standard
    error("Not implemented")
end

"""
    bonferroni(anova::AnovaData)
    dunn(anova::AnovaData)
    bonferronicorrection(anova::AnovaData)

Tests multiple comparisons using the Bonferroni p-value correction.
"""
bonferroni(args...) = bonferronicorrection(args...)
dunn(args...) = bonferronicorrection(args...)
function bonferronicorrection(anova::AnovaData)
    # use t distribution but alpha / m where m = number of comparisons as the critical factor, I think.
    error("Not implemented")
end

"""
    holmbonferroni(anova::AnovaData)
    holm(anova::AnovaData)
    holmbonferronicorrection(anova::AnovaData)

Tests multiple comparisons using the Bonferroni p-value correction.
"""
holmbonferroni(args...) = holmbonferronicorrection(args...)
holm(args...) = holmbonferronicorrection(args...)
function holmbonferronicorrection(anova::AnovaData)
    error("Not implemented")
end

"""
    lsd(anova::AnovaData)
    leastsignificantdifference(anova::AnovaData)



"""
lsd(args...) = leastsignificantdifference(args...)
function leastsignificantdifference(anova::AnovaData)
end

"""
    s(anova::AnovaData)
    scheffé(anova::AnovaData)
    scheffe(anova::AnovaData)
    multiplecontrasts(anova::AnovaData)

Performs Scheffé multiple contrasts.

Less powerful than Tukey in the typical case, but allows factor levels to be grouped prior to comparison.
"""
s(args...) = s(args...)
scheffé(args...) = multiplecontrasts(args...)
scheffe(args...) = multiplecontrasts(args...)
function multiplecontrasts(anova::AnovaData)
end

"""
    benjaminihochberg(anova::AnovaData)
    fdr(anova::AnovaData)
    falsediscoveryratecorrection(anova::AnovaData)
"""
benjaminihochberg(args...) = falsediscoveryratecorrection(args...)
fdr(args...) = falsediscoveryratecorrection(args...)
function falsediscoveryratecorrection(anova::AnovaData)
end


export pooled
export tukey, hsd, honestlysignificantdifference, multiplecomparison
export snk, studentnewmankeuls, newmankeuls, newmankeulsmultiplerange
export dmr, duncan, duncanmultiplerange
export wsd, whollysignificantdifference, multiplecomparisonandrange
export dunnett, controltoall
export bonferroni, dunn, bonferronicorrection
export holmbonferroni, holm, holmbonferronicorrection
export lsd, leastsignificantdifference
export s, scheffé, scheffe, multiplecontrasts
export benjaminihochberg, fdr, falsediscoveryratecorrection
