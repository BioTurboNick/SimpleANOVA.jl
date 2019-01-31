using .Plots

import .Plots.plot
function plot(anova::AnovaData)
    nfactors = length(anova.crossedfactors)
    if nfactors > 2
        error("Not yet able to handle more than 2 factors.")
    end
    nfactorlevels = size(anova.cellmeans)
    plot(legend = false)
    foreach(i -> plot!(results.cellmeans[:,i], seriestype=:line), 1:nfactorlevels[2])
    factorAmidpoint = (1 + nfactorlevels[1]) / 2
    factormeans = [mean(anova.cellmeans, dims = (1:nfactors)[Not(i)]) for i = 1:nfactors]
    plot!(factormeans[1], seriestype=:scatter)
    plot!(repeat([factorAmidpoint], nfactorlevels[2]), factormeans[2]', seriestype=:scatter)
end
