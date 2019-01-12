# SimpleANOVA.jl

[![Build Status](https://travis-ci.org/BioTurboNick/SimpleANOVA.jl.svg?branch=master)](https://travis-ci.org/BioTurboNick/SimpleANOVA.jl)

Analysis of Variance for Julia, the old-fashioned way.

This is a basic attempt to get a simple ANOVA implementation for Julia that works with data directly - no linear models.

The goal is to allow one function to do as much for you as possible, automatically choosing the right calculations.

Handles ANOVA with up to 3 crossed factors (fixed or random) and arbitrarily many nested factors. Requires equal replication.

It uses multidimensional arrays to interpret the structure of the data. Replicates should either be indexed along the first dimension or contained in a vector, with Factor A and Factor B the next available indices.

See docstring for usage.


VERY EXPERIMENTAL - NOT READY FOR REAL USE
