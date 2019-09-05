using .DataFrames

function levene(df::DataFrame, observationscolumn::Symbol, factorcolumns::Vector{Symbol})
    observations = df[observationscolumn]
    length(observations) > 0 || return
    eltype(observations) <: Number || error("Obervations must be numeric")
    isempty(factornames) && (factornames = [String(col) for col ∈ factorcolumns])
    levene(df[observationscolumn], [df[x] for x ∈ factorcolumns])
end
