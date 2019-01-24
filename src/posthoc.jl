"""
    pooled(anova::AnovaData)

Pools the mean squares for the interaction term(s) and error term and recalculates the test.

This should only be used if no interaction is detected. Increases the power of the test, but the chance of committing a Type  I error is less certain.

Generally, the conservative unpooled approach has been shown to be superior.
"""
function pooled(anova::AnovaData)

end

"""
    tukey(anova::AnovaData)
    hsd(anova::AnovaData)
    honestlysignificantdifference(anova::AnovaData)
    multiplecomparison(anova::AnovaData)

Performs the Tukey multiple comparisons posthoc test.
"""
tukey(args...) = multiplecomparison(args...)
hsd(args...) = multiplecomparison(args...)
honestlysignificantdifference(args...) = multiplecomparison(args...)
function multiplecomparison(anova::AnovaData)
end

"""
    snk(anova::AnovaData)
    studentnewmankeuls(anova::AnovaData)
    newmankeuls(anova::AnovaData)
    multiplerange(anova::AnovaData)

Performs the (Student-)Newman-Keuls multiple range posthoc test.

Tends to be more powerful than Tukey, but some suggest it is more likely to lead to Type I error than intended.
"""
snk(args...) = multiplerange(args...)
studentnewmankeuls(args...) = multiplerange(args...)
newmankeuls(args...) = multiplerange(args...)
function multiplerange(anova::AnovaData)
end

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
function multiplecomparisonandrange(anova::AnovaData)
end

"""
    dunnett(anova::AnovaData)
    controltoall(anova::AnovaData)

Performs the Dunnett method for comparing each factor level to a control group, and not all pairwise comparisons.
"""
dunnett(args...) = controltoall(args...)
function controltoall(anova::AnovaData)
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

export pooled
export tukey, hsd, honestlysignificantdifference, multiplecomparison
export snk, studentnewmankeuls, newmankeuls, multiplerange
export wsd, whollysignificantdifference, multiplecomparisonandrange
export dunnett, controltoall
export s, scheffé, scheffe, multiplecontrasts
