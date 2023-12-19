```@meta
CurrentModule = LongMemory
```

# Long Memory

Documentation for [LongMemory.jl](https://github.com/everval/LongMemory.jl).

## Description

This package implements functions to generate and estimate time series long memory models.

## Installation

The package is registered in the Julia registry and can be installed at the REPL with `] add LongMemory`.

## Generation

[Long Memory Generation](@ref) contains the documentation for the functions to generate time series long memory models.

## Semiparametric Estimators

[Semiparametric Estimators for Long Memory](@ref) contains the documentation for the functions to estimate the long memory parameter based on the log-periodogram regression. Estimators include the Geweke and Porter-Hudak estimators and the Whittle estimator, as well as bias-reduced versions of them.

## Parametric Estimation

[Parametric Estimators for Long Memory](@ref) contains the documentation for the functions to estimate time series long memory models by parametric methods. Of particular interest is the ARFIMA model. Moreover, a method to estimate the HAR model, a specification usually used as a proxy for long memory dynamics, is also available.

## Classic Estimation

[Classic Estimator for the Hurst Effect](@ref) contains the documentation for the functions to estimate the rescaled range (R/S) statistic and the Hurst coefficient.

## Forecasting

[Forecasting Long Memory Processes](@ref) contains the documentation for the functions to forecast long memory processes using the fractional differencing operator and the cross-sectional aggregation method. Forecasting using the HAR model is also available.

## List of Functions

[Index](@ref) contains the list of all the functions in the package.

## Data

[Data Available](@ref) contains the list of data sets available in the package.

## Examples

You can follow [this vignette](https://everval.github.io/files\LM_notebook.html) to learn how to use the package.

## Author

[J. Eduardo Vera-Valdés](https://everval.github.io/)

## Citation

If you use this package in your research, please cite it as:

```bibtex
@article{VeraValdes2024longmemory,
  title={LongMemory.jl: Generating, Estimating, and Forecasting Long Memory Models in Julia},
  author={Vera-Vald{\'e}s, J Eduardo},
  year={2024},
}
```

## Report a bug

Please, report any bugs [here](mailto: eduardo@math.aau.dk).
