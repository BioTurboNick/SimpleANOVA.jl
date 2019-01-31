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
    plots = Array{Plots.Plot, 2}(undef, nfactors, nfactors)
    for i in 1:nfactors
        otherfactorindexes = (1:nfactors)[Not(i)]
        for jindex in 1:(nfactors - 1)
            j = otherfactorindexes[jindex]
            factormeans = mean(anova.cellmeans, dims = otherfactorindexes[Not(jindex)])
            plots[i,j] = plot([selectdim(factormeans, j, k) |> vec for k = 1:nfactorlevels[j]], legend = false)
        end
        plots[i,i] = plot()
        factormidpoint = (1 + nfactorlevels[i]) / 2
    end
    plot(plots..., layout = (nfactors, nfactors))
#=
    factormeans = [mean(anova.cellmeans, dims = (1:nfactors)[Not(i)]) |> vec for i = 1:nfactors]


    plot!(factormeans, layout = (nfactors, nfactors), seriestype=:scatter)
    plot!(repeat([factormidpoint], nfactorlevels), factormeans[2]', seriestype=:scatter)

    plot!(factormeans[1], seriestype=:scatter)
    plot!(repeat([factorAmidpoint], nfactorlevels[2]), factormeans[2], seriestype=:scatter)
    =#
end
