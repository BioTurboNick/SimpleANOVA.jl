# SimpleANOVA

Analysis of Variance for Julia, the old-fashioned way.

This is a basic attempt to get a simple ANOVA implementation for Julia that works with data directly - no linear models.

The goal is to allow one function to do as much for you as possible, automatically choosing the right calculations.

Handles ANOVA with up to 3 crossed factors (fixed or random) and arbitrarily many nested factors. Requires equal replication. (If you have missing values, there are techniques to fill them in, within limits; e.g. use the average of the cell for one value).
It uses multidimensional arrays to interpret the structure of the data. Replicates should either be indexed along the first dimension or contained in a vector, with Factor B and Factor A the next available indices.

Can also work with multiple vectors and DataFrames.

Basic example
--------
```@example index
using Random # hide
Random.seed!(1) # hide
using CategoricalArrays
using DataFrames
using SimpleANOVA

df = DataFrame(
  id = categorical(repeat(1:2; inner=6)),
  a = categorical(repeat(1:3, inner=2, outer=2)),
  b = categorical(repeat(1:2, 6)),
  value = round.(rand(Float64, 12), digits=1)
)

anova(df, :value, [:b, :id, :a], [fixed, subject])
```

Other examples
--------
```
data                           # N-dimensional matrix of observations
levene(data)                   # test data for homogeniety of variance
result = anova(data)           # conduct the test
plot(result)                   # create pairwise factor plots
differencecontrasts(result, 2) # perform orthogonal series of difference contrasts for the levels of factor 2
```
```
data                           # vector of observations
factors                        # vector of factor level assignment vectors
levene(data)                   # test data for homogeniety of variance
result = anova(data, factors)  # conduct the test
plot(result)                   # create pairwise factor plots
contrast(result, [1, 2, 2, 2]) # calculate the contrast between factor level 1 of factor 1 with remaining factor levels
```
```
df                                         # DataFrame
factors                                    # vector of symbols for factor assignment columns
levene(df, :observations, factors)         # test data for homogeniety of variance
result = anova(df, :observations, factors) # conduct the test
plot(result)                               # create pairwise factor plots
simplecontrasts(result)                    # calculate the contrast between the first factor level (control) to each other level
```

Differences from SPSS
---------------------
Choice of error terms for the F tests in mixed ANOVA follows Zar 1999, _Biostatistical Analysis_, and Howell 2013, _Statistical Methods for Psychology_, which differs from SPSS 25 Univariate GLM as follows:

2-way ANOVA with 1 fixed and 1 random factor

|                | A (Fixed) | B (Random) |
|----------------|-----------|------------|
| SPSS           | Error     | Error      |
| SimpleANOVA.jl | A×B       | Error      |

3-way ANOVA with 2 fixed and 1 random factors

|                | A (Fixed) | B (Fixed) | C (Random)        | A×B   | A×C   | B×C   | A×B×C |
|----------------|-----------|-----------|-------------------|-------|-------|-------|-------|
| SPSS           | A×C       | B×C       | A×C + B×C - A×B×C | A×B×C | A×B×C | A×B×C | Error |
| SimpleANOVA.jl | A×C       | B×C       | Error             | A×B×C | Error | Error | Error |

3-way ANOVA with 1 fixed and 2 random factors

|                | A (Fixed)         | B (Random)        | C (Random)        | A×B   | A×C   | B×C   | A×B×C |
|----------------|-------------------|-------------------|-------------------|-------|-------|-------|-------|
| SPSS           | A×C + A×B - A×B×C | A×B + B×C - A×B×C | A×C + B×C - A×B×C | A×B×C | A×B×C | A×B×C | Error |
| SimpleANOVA.jl | A×C + A×B - A×B×C | B×C               | B×C               | A×B×C | A×B×C | Error | Error |

Changelog
--------
**New in v0.5**: ω² effect size (Disclaimer: effect size calculations for nested and 3-way mixed ANOVA is inferred and may not be correct.)

**New in v0.6**: Linear contrasts added for planned post-hoc tests between factor levels. (Disclaimer: Uncertain if contrasts for multiple factors are correct.)

**New in v0.7**: Repeated measures added. Effect size not yet supported.

**New in v0.8**: Implemented generalized ω² effect size which is robust to experimental design choice.
