var documenterSearchIndex = {"docs":
[{"location":"api/#API","page":"API","title":"API","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Modules = [SimpleANOVA]\nOrder = [:function, :type]","category":"page"},{"location":"api/#SimpleANOVA.anova-Union{Tuple{T}, Tuple{AbstractVector{T}, AbstractVector{<:AbstractVector}}, Tuple{AbstractVector{T}, AbstractVector{<:AbstractVector}, Vector{FactorType}}} where T<:Number","page":"API","title":"SimpleANOVA.anova","text":"anova(observations::Array{Union{Number, Vector{Number}}}, factortypes = FactorType[]; factornames = String[], hasreplicates = true)\nanova(observations::Vector{Number}, factorassignments::Vector{Vector{Any}}, factortypes = FactorType[]; factornames = String[], hasreplicates = true)\nanova(df::DataFrame, observationscolumn::Symbol, factorcolumns::Vector{Symbol}, factortypes = FactorType[]; factornames = String[])\n\nPerforms an Analysis of Variance (ANOVA) calculation.\n\nOperates on fixed or random factors, subject/block factors, and nested factors, with or without replicates, on balanced data.\n\nLimitations:\n\nIf any factors are random type, limited to 3-way\nRepeated measures ANOVA (with subject/block factor) limited to 3 fixed factors which may be partitioned as within or among subjects\n\nOperates on up to 3 crossed factors (fixed or random) and arbitrarily many random nested factors, with or without replicates, on balanced data.\n\nArguments\n\nobservations: Array containing the values to test. For the array, each dimension is a factor level, such that observations[2,5,3] indicates the 2nd level of the first factor, the 5th level of the second factor, and the 3rd level of the third factor. May contain values or vectors of values, where the vector contains replicates. Factors should be ordered with least significant first. For the vector, must provide factorassignments to specify factor levels.\nfactorassignments: Vector of vectors of integers specifying how each observation is assigned to a factor level. Provide this when observations is given as a vector. Factor levels do not have to be consecutive or ordered. Nested factors must reuse factor levels currently.\nfactortypes: Vector indicating the FactorType for each factor. Factors must be ordered with nested first, then random/measured or fixed/manipulated in any order. For repeated measures, subject or block must appear after within-subject fixed and before among-subject fixed. If too few values are provided, remaining are assumed to be fixed.\nfactornames: Vector of names for each factor, excluding the replicate factor. If empty, will be automatically populated alphabetically.\nhasreplicates: Boolean to specify if the first level should be considered \n\nNotes: The last index will be the top factor in the table.\n\nOutput: AnovaData structure containing the test results for each factor.\n\nExamples\n\nanova(observations)                           # N-way fixed-effects ANOVA with replicates (vectors or first dimension)\nanova(observations, hasreplicates = false)    # N-way fixed-effects ANOVA without replicates (first dimension)\nanova(observations, [random])                 # N-way ANOVA with lower random factor and 1 or 2 upper fixed factors\nanova(observations, [random])                 # N-way ANOVA with lower random factor and 1 or 2 upper fixed factors\nanova(observations, [fixed, random])          # N-way ANOVA with 1 lower fixed factor, 1 random factor, and 0 or 1 upper fixed factor\nanova(observations, [nested, random])         # N-way fixed-effects ANOVA with 1 random nested factor, 1 random factor, and 1-2 fixed factors\nanova(observations, [fixed, subject])         # N-way repeated measures ANOVA with 1 within-subjects fixed factor\nanova(observations, [fixed, block])           # N-way repeated measures ANOVA with 1 within-block fixed factor\nanova(observations, [nested, fixed, subject]) # N-way repeated measures ANOVA with 1 within-subjects fixed factor and 1 nested factor\n\nGlossary\n\nobservation: The dependent variable.\nfactor: An independent variable.\nfactor level: A value of a factor.\nbalanced: All combinations of factor levels have the same number of observations.\ncrossed factor: A factor with levels that combine with the levels of all other crossed factors.\nfixed/manipulated factor: A factor with fixed effects (e.g. treatment, concentration, exposure time).\nrandom/measured factor: A factor with random effects (e.g. location, individual).\nnested factor: A random factor where the levels are unique to a combination of crossed factor levels (e.g. replicate).\nsubject/block factor: A nested factor that is subjected to multiple levels of another factor.\nsum of squares (SS): A measure of variance that is dependent on sample size. Also called \"sum of squared deviations.\"\ndegrees of freedom (DF, ν): The number of bins in which the values could have been moved, if random.\nmean square (MS): SS / DF. Corrects for the larger variance expected if random values can be assigned to more bins. Also called \"mean squared error\" or \"mean squared deviation.\"\nF-statistic: The division of MS values produce a result belonging to the \"F distribution\", the shape of which depends on the DF of the numerator and error. The location of this value on the distribution provides the p-value.\np-value: The probability that, if all measurements had been drawn from the same population, you would obtain data at least as extreme as contained in your observations.\neffect size (ω²): The standardized difference in the measurement caused by the factor. This is the generalized version which is largely stable with respect to study design.\n\n\n\n\n\n","category":"method"},{"location":"api/#SimpleANOVA.contrast","page":"API","title":"SimpleANOVA.contrast","text":"contrast(anovaresult::AnovaResult, groupassignment::Vector{Int})\n\nCalculate a single-factor contrast against the groups specified in groupassignment. Valid groups are \"0\", \"1\", and \"2\", where \"0\" excludes the group from the comparison.\n\nNote: If you do nonorthogonal contrasts, use the Bonferroni or Šidák correction to adjust the α level (the level at which you consider a result likely to be true):   Bonferroni: α′ = α/c for c contrasts   Šidák:      α′ = 1 - (1 - α)^(1 / c) (slightly better) where c = number of contrasts\n\nNote: Effect size is calcluated using the error term associated with the factor. Other choices are possible, including average of each group error; or the error associated with a control.\n\n\n\n\n\n","category":"function"},{"location":"api/#SimpleANOVA.differencecontrasts","page":"API","title":"SimpleANOVA.differencecontrasts","text":"differencecontrast(anovaresult::AnovaData, factorindex = 1, reverseorder = false)\n\nCompute orthogonal contrasts on the factor levels in the original order. Forward direction also known as a \"Helmert\" contrast; revere direction may also be called \"Difference\" (as in SPSS). See contrast function for more.\n\n\n\n\n\n","category":"function"},{"location":"api/#SimpleANOVA.levene-Union{Tuple{AbstractArray{T}}, Tuple{T}} where T<:Union{Number, AbstractVector{<:Number}}","page":"API","title":"SimpleANOVA.levene","text":"levene(observations::Array{Union{Number, Vector{Number}}})\nlevene(observations::Vector{Number}, factorassignments::Vector{Vector{Any}})\nlevene(df::DataFrame, observationscolumn::Symbol, factorcolumns::Vector{Symbol})\n\nPerforms Levene's homogeniety of variance test.\n\nOperates on at least 3 crossed factors (fixed or random) and arbitrarily many random nested factors on balanced data. Levene's test is identical to performing an ANOVA on the absolute deviations of each cell.\n\nIf Levene's test is significant, the data is heteroscedastic (variances differ). If it is not, variances could differ or the sample size could be too small to tell.\n\nArguments\n\nobservations: Array containing the values to test. For the array, each dimension is a factor level, such that observations[2,5,3] indicates the 2nd level of the first factor, the 5th level of the second factor, and the 3rd level of the third factor. May contain values or vectors of values, where the vector contains replicates. Factors should be ordered with least significant first. For the vector, must provide factorassignments to specify factor levels.\nfactorassignments: Vector of vectors of integers specifying how each observation is assigned to a factor level. Provide this when observations is given as a vector. Factor levels do not have to be consecutive or ordered. Nested factors must reuse factor levels currently.\n\nNotes: Throws an exception if it can tell there are no replicates. Do not use with on any data structure you would pass into anova() using hasrepliates=false.\n\nOutput: AnovaData structure containing the test results for each factor.\n\nExamples\n\nlevene(observations)\n\n\n\n\n\n","category":"method"},{"location":"api/#SimpleANOVA.repeatedcontrasts","page":"API","title":"SimpleANOVA.repeatedcontrasts","text":"repeatedcontrast(anovaresult::AnovaData, factorindex = 1)\n\nCompute contrasts between neighboring levels. Non-orthogonal. See contrast function for more.\n\n\n\n\n\n","category":"function"},{"location":"api/#SimpleANOVA.simplecontrasts","page":"API","title":"SimpleANOVA.simplecontrasts","text":"simplecontrast(anovaresult::AnovaData, controlindex = 1, factorindex = 1)\n\nCompute contrasts of each level to a single level (control). Non-orthogonal. See contrast function for more.\n\n\n\n\n\n","category":"function"},{"location":"api/#SimpleANOVA.AnovaContrastResult","page":"API","title":"SimpleANOVA.AnovaContrastResult","text":"AnovaValue <: AnovaEffect\n\nA set of values for an Anova item for which a mean square is not required.\n\ncontrast - the contrast value\n\ndf - degrees of freedom\n\nf - the F statistic\n\np - the probability of a Type I error (incorrect rejection of null hypothesis)\n\nr - the effect size\n\n\n\n\n\n","category":"type"},{"location":"api/#SimpleANOVA.AnovaData","page":"API","title":"SimpleANOVA.AnovaData","text":"AnovaData\n\nContainer for the complete results of an ANOVA test.\n\n\n\n\n\n","category":"type"},{"location":"api/#SimpleANOVA.AnovaEffect","page":"API","title":"SimpleANOVA.AnovaEffect","text":"AnovaEffect\n\nSupertype for all ANOVA data items: AnovaValue, AnovaFactor, AnovaResult\n\n\n\n\n\n","category":"type"},{"location":"api/#SimpleANOVA.AnovaFactor","page":"API","title":"SimpleANOVA.AnovaFactor","text":"AnovaFactor <: AnovaEffect\n\nA set of values for an Anova effect which has a mean square.\n\nname - the name of this factor\n\nss - sum of squares\n\ndf - degrees of freedom\n\nms - mean square\n\n\n\n\n\n","category":"type"},{"location":"api/#SimpleANOVA.AnovaResult","page":"API","title":"SimpleANOVA.AnovaResult","text":"AnovaResult <: AnovaEffect\n\nA set of values for an Anova factor which has been tested\n\nname - the name of this factor\n\nss - sum of squares\n\ndf - degrees of freedom\n\nms - mean square\n\nf - the F statistic\n\np - the probability of a Type I error (incorrect rejection of null hypothesis)\n\nω² - the effect size\n\n\n\n\n\n","category":"type"},{"location":"api/#SimpleANOVA.AnovaValue","page":"API","title":"SimpleANOVA.AnovaValue","text":"AnovaValue <: AnovaEffect\n\nA set of values for an Anova item for which a mean square is not required.\n\nname - the name of this value\n\nss - sum of squares\n\ndf - degrees of freedom\n\n\n\n\n\n","category":"type"},{"location":"api/#SimpleANOVA.FactorType","page":"API","title":"SimpleANOVA.FactorType","text":"FactorType\n\nfixed manipulated - A crossed fixed-effects factor\n\nrandom measured - A crossed random-effects factor\n\nnested - A random factor fully nested within another factor\n\nsubject block - A random factor subjected to multiple levels of another factor\n\n\n\n\n\n","category":"type"},{"location":"#SimpleANOVA","page":"SimpleANOVA","title":"SimpleANOVA","text":"","category":"section"},{"location":"","page":"SimpleANOVA","title":"SimpleANOVA","text":"Analysis of Variance for Julia, the old-fashioned way.","category":"page"},{"location":"","page":"SimpleANOVA","title":"SimpleANOVA","text":"This is a basic attempt to get a simple ANOVA implementation for Julia that works with data directly - no linear models.","category":"page"},{"location":"","page":"SimpleANOVA","title":"SimpleANOVA","text":"The goal is to allow one function to do as much for you as possible, automatically choosing the right calculations.","category":"page"},{"location":"","page":"SimpleANOVA","title":"SimpleANOVA","text":"Handles ANOVA with up to 3 crossed factors (fixed or random) and arbitrarily many nested factors. Requires equal replication. (If you have missing values, there are techniques to fill them in, within limits; e.g. use the average of the cell for one value). It uses multidimensional arrays to interpret the structure of the data. Replicates should either be indexed along the first dimension or contained in a vector, with Factor B and Factor A the next available indices.","category":"page"},{"location":"","page":"SimpleANOVA","title":"SimpleANOVA","text":"Can also work with multiple vectors and DataFrames.","category":"page"},{"location":"#Basic-example","page":"SimpleANOVA","title":"Basic example","text":"","category":"section"},{"location":"","page":"SimpleANOVA","title":"SimpleANOVA","text":"using Random # hide\nRandom.seed!(1) # hide\nusing CategoricalArrays\nusing DataFrames\nusing SimpleANOVA\n\ndf = DataFrame(\n  id = categorical(repeat(1:2; inner=6)),\n  a = categorical(repeat(1:3, inner=2, outer=2)),\n  b = categorical(repeat(1:2, 6)),\n  value = round.(rand(Float64, 12), digits=1)\n)\n\nanova(df, :value, [:b, :id, :a], [fixed, subject])","category":"page"},{"location":"#Other-examples","page":"SimpleANOVA","title":"Other examples","text":"","category":"section"},{"location":"","page":"SimpleANOVA","title":"SimpleANOVA","text":"data                           # N-dimensional matrix of observations\nlevene(data)                   # test data for homogeniety of variance\nresult = anova(data)           # conduct the test\nplot(result)                   # create pairwise factor plots\ndifferencecontrasts(result, 2) # perform orthogonal series of difference contrasts for the levels of factor 2","category":"page"},{"location":"","page":"SimpleANOVA","title":"SimpleANOVA","text":"data                           # vector of observations\nfactors                        # vector of factor level assignment vectors\nlevene(data)                   # test data for homogeniety of variance\nresult = anova(data, factors)  # conduct the test\nplot(result)                   # create pairwise factor plots\ncontrast(result, [1, 2, 2, 2]) # calculate the contrast between factor level 1 of factor 1 with remaining factor levels","category":"page"},{"location":"","page":"SimpleANOVA","title":"SimpleANOVA","text":"df                                         # DataFrame\nfactors                                    # vector of symbols for factor assignment columns\nlevene(df, :observations, factors)         # test data for homogeniety of variance\nresult = anova(df, :observations, factors) # conduct the test\nplot(result)                               # create pairwise factor plots\nsimplecontrasts(result)                    # calculate the contrast between the first factor level (control) to each other level","category":"page"},{"location":"#Differences-from-SPSS","page":"SimpleANOVA","title":"Differences from SPSS","text":"","category":"section"},{"location":"","page":"SimpleANOVA","title":"SimpleANOVA","text":"Choice of error terms for the F tests in mixed ANOVA follows Zar 1999, Biostatistical Analysis, and Howell 2013, Statistical Methods for Psychology, which differs from SPSS 25 Univariate GLM as follows:","category":"page"},{"location":"","page":"SimpleANOVA","title":"SimpleANOVA","text":"2-way ANOVA with 1 fixed and 1 random factor","category":"page"},{"location":"","page":"SimpleANOVA","title":"SimpleANOVA","text":" A (Fixed) B (Random)\nSPSS Error Error\nSimpleANOVA.jl A×B Error","category":"page"},{"location":"","page":"SimpleANOVA","title":"SimpleANOVA","text":"3-way ANOVA with 2 fixed and 1 random factors","category":"page"},{"location":"","page":"SimpleANOVA","title":"SimpleANOVA","text":" A (Fixed) B (Fixed) C (Random) A×B A×C B×C A×B×C\nSPSS A×C B×C A×C + B×C - A×B×C A×B×C A×B×C A×B×C Error\nSimpleANOVA.jl A×C B×C Error A×B×C Error Error Error","category":"page"},{"location":"","page":"SimpleANOVA","title":"SimpleANOVA","text":"3-way ANOVA with 1 fixed and 2 random factors","category":"page"},{"location":"","page":"SimpleANOVA","title":"SimpleANOVA","text":" A (Fixed) B (Random) C (Random) A×B A×C B×C A×B×C\nSPSS A×C + A×B - A×B×C A×B + B×C - A×B×C A×C + B×C - A×B×C A×B×C A×B×C A×B×C Error\nSimpleANOVA.jl A×C + A×B - A×B×C B×C B×C A×B×C A×B×C Error Error","category":"page"},{"location":"#Changelog","page":"SimpleANOVA","title":"Changelog","text":"","category":"section"},{"location":"","page":"SimpleANOVA","title":"SimpleANOVA","text":"New in v0.5: ω² effect size (Disclaimer: effect size calculations for nested and 3-way mixed ANOVA is inferred and may not be correct.)","category":"page"},{"location":"","page":"SimpleANOVA","title":"SimpleANOVA","text":"New in v0.6: Linear contrasts added for planned post-hoc tests between factor levels. (Disclaimer: Uncertain if contrasts for multiple factors are correct.)","category":"page"},{"location":"","page":"SimpleANOVA","title":"SimpleANOVA","text":"New in v0.7: Repeated measures added. Effect size not yet supported.","category":"page"},{"location":"","page":"SimpleANOVA","title":"SimpleANOVA","text":"New in v0.8: Implemented generalized ω² effect size which is robust to experimental design choice.","category":"page"}]
}
