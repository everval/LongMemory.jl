```@meta
CurrentModule = LongMemory
```

# Long Memory

Documentation for [LongMemory.jl](https://github.com/everval/LongMemory.jl).

## Description

**LongMemory.jl** is a package for time series long memory modelling in [*Julia*](https://julialang.org/).

The package provides functions for *generating long memory*, *estimating the parameters of the models*, and *forecasting*.

Generating methods include *fractional differencing*, *stochastic error duration*, and *cross-sectional aggregation*.

Estimators include *classic* ones used to estimate the Hurst effect, those inspired by the *log-periodogram regression*, and *parametric* ones.

Forecasting is provided for all parametric estimators.

Moreover, the package adds *plotting capabilities* to illustrate long memory dynamics and forecasting.

Finally, the package includes the *Nile River minima* and *Northern Hemisphere Temperature Anomalies* data sets to illustrate the use of the functions.

## Installation

The package is registered in the Julia registry and can be installed at the REPL with `] add LongMemory`.

## Generation

[Long Memory Generation](@ref) contains the documentation for the functions to generate time series long memory models.

## Classic Estimation

[Classic Estimator for the Hurst Effect](@ref) contains the documentation for the functions to estimate the rescaled range (R/S) statistic and the Hurst coefficient.

## Semiparametric Estimators

[Semiparametric Estimators for Long Memory](@ref) contains the documentation for the functions to estimate the long memory parameter based on the log-periodogram regression. Estimators include the Geweke and Porter-Hudak estimators and the Whittle estimator, as well as bias-reduced versions of them.

## Parametric Estimation

[Parametric Estimators for Long Memory](@ref) contains the documentation for the functions to estimate time series long memory models by parametric methods. Of particular interest is the ARFIMA model. Moreover, a method to estimate the HAR model, a specification usually used as a proxy for long memory dynamics, is also available.

## Forecasting

[Forecasting Long Memory Processes](@ref) contains the documentation for the functions to forecast long memory processes using the fractional differencing operator and the cross-sectional aggregation method. Forecasting using the HAR model is also available.

## Structural Changes

[Structural Change Tests for Long Memory Processes](@ref) contains the documentation for the functions to estimate structural changes in the long memory parameter of a time series.

## Data

[Data Available](@ref) contains the list of data sets available in the package.

## Examples

You can follow [this vignette](https://everval.github.io/files/LM_notebook_illustration.html) to learn how to use the package.

## List of Functions

[Index](@ref) contains the list of all the functions in the package.

## Benchmarks

The following [notebook](https://everval.github.io/files/LM_notebook_benchmark.html) contains benchmarks for some of the functions in the package against popular **R** packages: *fracdiff* and *longMemoryTS*.

## Author

[J. Eduardo Vera-Valdés](https://everval.github.io/)

## Citation

If you use this package in your research, please cite it as:

Vera-Valdés, J.E. (2024). "LongMemory.jl: Generating, Estimating, and Forecasting Long Memory Models in Julia". arXiv 2401.14077. [https://arxiv.org/abs/2401.14077](https://arxiv.org/abs/2401.14077)

```bibtex
@article{VERAVALDES2024a,
title = {LongMemory.jl: Generating, Estimating, and Forecasting Long Memory Models in Julia},
year = {2024},
author = {Vera-Valdés, J.E.},
journal = {arXiv preprint arXiv:2401.14077},
url = {https://arxiv.org/abs/2401.14077}
}
```

## Report a bug

Please, report any bugs [here](mailto: eduardo@math.aau.dk).
