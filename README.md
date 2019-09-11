# SimpleANOVA.jl

[![Build Status](https://travis-ci.org/BioTurboNick/SimpleANOVA.jl.svg?branch=master)](https://travis-ci.org/BioTurboNick/SimpleANOVA.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/e2uhgjm8bwbn9sja/branch/master?svg=true)](https://ci.appveyor.com/project/BioTurboNick/simpleanova-jl/branch/master)
[![codecov.io](https://codecov.io/github/BioTurboNick/SimpleANOVA.jl/coverage.svg?branch=master)](https://codecov.io/github/BioTurboNick/SimpleANOVA.jl?branch=master)

Analysis of Variance for Julia, the old-fashioned way.

This is a basic attempt to get a simple ANOVA implementation for Julia that works with data directly - no linear models.

The goal is to allow one function to do as much for you as possible, automatically choosing the right calculations.

Handles ANOVA with up to 3 crossed factors (fixed or random) and arbitrarily many nested factors. Requires equal replication.

It uses multidimensional arrays to interpret the structure of the data. Replicates should either be indexed along the first dimension or contained in a vector, with Factor B and Factor A the next available indices.

Can also work with multiple vectors and DataFrames.

See docstring for usage.

Examples:
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

Experimental, use at own risk!
