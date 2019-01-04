# SimpleANOVA.jl
Analysis of Variance for Julia, the old-fashioned way.

This is a basic attempt to get a simple ANOVA implementation for Julia that works with data directly - no linear models.

The goal is to allow one function to do as much for you as possible, automatically choosing the right calculations.

Handles ANOVA with up to 3 crossed factors (fixed or random) and arbitrarily many nested factors. Requires equal replication.

It uses multidimensional arrays to interpret the structure of the data. Replicates should either be indexed along the first dimension or contained in a vector, with Factor A and Factor B the next available indices.

`factortype` argument must be provided to specify the type of each dimension of the matrix. `fixed`, `random`, `nested`, or `replicate`. E.g. a 3-way ANOVA (1 random, 2 fixed) with 2 nested factors and replicates contained in vectors would be specified by `[nested, nested, random, fixed, fixed]`. If the replicates are instead the first dimension of the matrix: `[replicates, nested, nested, random, fixed, fixed]`. `replicates` must come first, followed by `nested`, followed by `random` or `fixed` in any order.


VERY EXPERIMENTAL - NOT READY FOR REAL USE
