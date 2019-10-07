# SimpleANOVA.jl

[![Build Status](https://travis-ci.org/BioTurboNick/SimpleANOVA.jl.svg?branch=master)](https://travis-ci.org/BioTurboNick/SimpleANOVA.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/e2uhgjm8bwbn9sja/branch/master?svg=true)](https://ci.appveyor.com/project/BioTurboNick/simpleanova-jl/branch/master)
[![codecov.io](https://codecov.io/github/BioTurboNick/SimpleANOVA.jl/coverage.svg?branch=master)](https://codecov.io/github/BioTurboNick/SimpleANOVA.jl?branch=master)

Analysis of Variance for Julia, the old-fashioned way.

This is a basic attempt to get a simple ANOVA implementation for Julia that works with data directly - no linear models.

The goal is to allow one function to do as much for you as possible, automatically choosing the right calculations.

Handles ANOVA with up to 3 crossed factors (fixed or random) and arbitrarily many nested factors. Requires equal replication. (If you have missing values, there are techniques to fill them in, within limits; e.g. use the average of the cell for one value).
It uses multidimensional arrays to interpret the structure of the data. Replicates should either be indexed along the first dimension or contained in a vector, with Factor B and Factor A the next available indices.

Can also work with multiple vectors and DataFrames.

**New in v0.5**: ω² effect size (Disclaimer: effect size calculations for nested and 3-way mixed ANOVA is inferred and may not be correct.)

See docstring for usage.

**Experimental, use at own risk!**

Examples
--------
```
data                 # N-dimensional matrix of observations
levene(data)         # test data for homogeniety of variance
result = anova(data) # conduct the test
plot(result)         # create pairwise factor plots
```
```
data                          # vector of observations
factors                       # vector of factor level assignment vectors
levene(data)                  # test data for homogeniety of variance
result = anova(data, factors) # conduct the test
plot(result)                  # create pairwise factor plots
```
```
df                                         # DataFrame
factors                                    # vector of symbols for factor assignment columns
levene(df, :observations, factors)         # test data for homogeniety of variance
result = anova(df, :observations, factors) # conduct the test
plot(result)                               # create pairwise factor plots
```

Differences from SPSS
---------------------
Choice of error terms for the F tests in mixed ANOVA follows Zar 1999, _Biostatistical Analysis_, and Howell 2013, _Statistical Methods for Psychology_, which differs from SPSS 25 Univariate GLM as follows:

2-way ANOVA with 1 fixed and 1 random factor

|                | Fixed | Random |
|----------------|-------|--------|
| SPSS           | Error | Error  |
| SimpleANOVA.jl | A×B   | Error  |

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



