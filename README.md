# SimpleANOVA.jl

Analysis of Variance for Julia, the old-fashioned way.

| Documentation | CI Status
|:-----------------:|:------------------:|
| [![][docs-stable-img]][docs-stable-url] [![][docs-latest-img]][docs-latest-url] | [![][ci-img]][ci-url] |

[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: https://BioTurboNick.github.io/SimpleANOVA.jl/dev

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://BioTurboNick.github.io/SimpleANOVA.jl/stable

[ci-img]: https://github.com/BioTurboNick/SimpleANOVA.jl/workflows/CI-stable/badge.svg
[ci-url]: https://github.com/BioTurboNick/SimpleANOVA.jl/actions?query=workflow%3ACI-stable+branch%3Amaster

This is a basic attempt to get a simple ANOVA implementation for Julia that works with data directly - no linear models.

The goal is to allow one function to do as much for you as possible, automatically choosing the right calculations.

Handles ANOVA with up to 3 crossed factors (fixed or random) and arbitrarily many nested factors. Requires equal replication. (If you have missing values, there are techniques to fill them in, within limits; e.g. use the average of the cell for one value).
It uses multidimensional arrays to interpret the structure of the data. Replicates should either be indexed along the first dimension or contained in a vector, with Factor B and Factor A the next available indices.

Can also work with multiple vectors and DataFrames.

**New in v0.5**: ω² effect size (Disclaimer: effect size calculations for nested and 3-way mixed ANOVA is inferred and may not be correct.)

**New in v0.6**: Linear contrasts added for planned post-hoc tests between factor levels. (Disclaimer: Uncertain if contrasts for multiple factors are correct.)

**New in v0.7**: Repeated measures added. Effect size not yet supported.

**New in v0.8**: Implemented generalized ω² effect size which is robust to experimental design choice.

See docstrings for usage.

**Experimental, use at own risk!**
