# LongMemory

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://everval.github.io/LongMemory.jl/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://everval.github.io/LongMemory.jl/dev/)
[![Build Status](https://github.com/everval/LongMemory.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/everval/LongMemory.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/everval/LongMemory.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/everval/LongMemory.jl)
[![DOI](https://zenodo.org/badge/697765094.svg)](https://doi.org/10.5281/zenodo.15096772)


## About

**LongMemory.jl** is a package for time series long memory modelling in [***Julia***](https://julialang.org/).

The package provides functions for *generating long memory*, *estimating the parameters of the models*, and *forecasting*.

Generating methods include *fractional differencing*, *stochastic error duration*, and *cross-sectional aggregation*.

Estimators include *classic* ones used to estimate the Hurst effect, those inspired by the *log-periodogram regression*, and *parametric* ones.

Forecasting is provided for all parametric estimators.

Moreover, the package adds *plotting capabilities* to illustrate long memory dynamics and forecasting.

Finally, the package includes the *Nile River minima* and *Northern Hemisphere Temperature Anomalies* data sets to illustrate the use of the functions.

## Installation

The package is registered in the Julia General registry and can be installed with the Julia package manager.

From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```julia
pkg> add LongMemory
```

Or, equivalently, via the `Pkg` API:

```julia
julia> import Pkg; Pkg.add("LongMemory")
```

## Usage

Once installed, the package can be imported with the command:

```julia
julia> using LongMemory
```

## Documentation

The package documentation is available [here](https://everval.github.io/LongMemory.jl/) or the link below.

## Examples

An illustrative example of the package usage can be found [here.](https://everval.github.io/files/LM_notebook_illustration.html)

## Benchmarks

The following [vignette](https://everval.github.io/files/LM_notebook_benchmark.html) contains benchmarks for some of the functions in the package against popular ***R*** packages: ***fracdiff*** and ***longMemoryTS***.

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

## Contributing 

All types of contributions are encouraged and appreciated.

If you find a bug or have a feature request, please open a new [issue](https://github.com/everval/LongMemory.jl/issues). If you would like to contribute code, please open a [pull request](https://github.com/everval/LongMemory.jl/pulls). I welcome all contributions, including bug fixes, documentation improvements, and new features.

Thank you for considering contributing!
