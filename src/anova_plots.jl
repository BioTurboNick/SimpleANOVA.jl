using .Plots

import .Plots.plot
"""
    plot(anova::AnovaData)

Creates factor plots for ANOVA data.

Enables inspection of ANOVA for the effect size of interactions and relationships between factor levels.

Import Plots package to use this.
"""
function plot(anova::AnovaData)
    nfactors = length(anova.crossedfactors)
    if nfactors == 1
        error("Not useful to plot just one factor.")
    end
    nfactorlevels = size(anova.cellmeans)
    clibrary(:cmocean)
    plots = Array{Plots.Plot, 2}(undef, nfactors, nfactors)
    for i in 1:nfactors
        otherfactorindexes = (1:nfactors)[Not(i)]
        factormidpoint = (1 + nfactorlevels[i]) / 2
        factormeans = [mean(anova.cellmeans, dims = (1:nfactors)[Not(i)]) |> vec for i = 1:nfactors]
        for jindex in 1:(nfactors - 1)
            j = otherfactorindexes[jindex]
            pairfactormeans = mean(anova.cellmeans, dims = otherfactorindexes[Not(jindex)])
            plots[i,j] = plot([selectdim(pairfactormeans, j, k) |> vec for k = 1:nfactorlevels[j]], legend = false)
            scatter!(factormeans[i], markershape = :+, markercolor = :gray)
            scatter!(repeat([factormidpoint], nfactorlevels[j]), factormeans[j], markershape = :x, markercolor = :gray)
        end
        plots[i,i] = plot()

    end
    plot(plots..., layout = (nfactors, nfactors))
end
