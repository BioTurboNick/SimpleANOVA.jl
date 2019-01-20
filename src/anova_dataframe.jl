function anova(df::DataFrame, observationscolumn::Symbol, factorcolumns::Vector{Symbol}, factortypes::Vector{FactorType} = FactorType[]; factornames::Vector{<:AbstractString} = String[])
    observations = df[observationscolumn]
    length(observations) > 0 || return
    eltype(observations) <: Number || error("Obervations must be numeric")
    isempty(factornames) && (factornames = [String(col) for col ∈ factorcolumns])
    anova(df[observationscolumn], df[factorcolumns], factortypes, factornames = factornames)
end
