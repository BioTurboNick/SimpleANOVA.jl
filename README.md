# SimpleANOVA.jl

[![Build Status](https://travis-ci.org/BioTurboNick/SimpleANOVA.jl.svg?branch=master)](https://travis-ci.org/BioTurboNick/SimpleANOVA.jl)
[![appveyor][https://img.shields.io/appveyor/ci/BioTurboNick/simpleanova-jl/master.svg?label=Windows
]](https://ci.appveyor.com/project/BioTurboNick/simpleanova-jl)
[![codecov.io](https://codecov.io/github/BioTurboNick/SimpleANOVA.jl/coverage.svg?branch=master)](https://codecov.io/github/BioTurboNick/SimpleANOVA.jl?branch=master)

Analysis of Variance for Julia, the old-fashioned way.

This is a basic attempt to get a simple ANOVA implementation for Julia that works with data directly - no linear models.

The goal is to allow one function to do as much for you as possible, automatically choosing the right calculations.

Handles ANOVA with up to 3 crossed factors (fixed or random) and arbitrarily many nested factors. Requires equal replication.

It uses multidimensional arrays to interpret the structure of the data. Replicates should either be indexed along the first dimension or contained in a vector, with Factor B and Factor A the next available indices.

See docstring for usage.

Note: Uses parts from [InvertedIndices.jl](https://github.com/mbauman/InvertedIndices.jl)


VERY EXPERIMENTAL - NOT READY FOR REAL USE
