# SimpleANOVA.jl
Analysis of Variance for Julia, the old-fashioned way.

This is a basic attempt to get a simple ANOVA implementation for Julia that works with data directly - no linear models.

The goal is to allow one function to do as much for you as possible, automatically choosing the right calculations.

It uses multidimensional arrays to interpret the structure of the data. Replicates should either be indexed along the first dimension or contained in a vector, with Factor A and Factor B the next available indices.

A 1-way ANOVA may be called via `anova(observations)`.

A 2-way ANOVA may be called via e.g. `anova(observations, [fixed, random])` if replicates are in vectors or there are no replicates, or
                                     `anova(observations, [replicate, fixed, random])` if replicates are in the first index.
                                     
A 3-way ANOVA follows the same pattern as 2-way ANOVA.

Nested ANOVA is planned.

VERY EXPERIMENTAL - NOT READY FOR REAL USE
