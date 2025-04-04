# LongMemory

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://everval.github.io/LongMemory.jl/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://everval.github.io/LongMemory.jl/dev/)
[![Build Status](https://github.com/everval/LongMemory.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/everval/LongMemory.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/everval/LongMemory.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/everval/LongMemory.jl)
[![DOI](https://zenodo.org/badge/697765094.svg)](https://doi.org/10.5281/zenodo.15096772)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.07708/status.svg)](https://doi.org/10.21105/joss.07708)


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

The following [notebook](https://everval.github.io/files/LM_notebook_benchmark.html) contains benchmarks for some of the functions in the package against popular ***R*** packages: ***fracdiff*** and ***longMemoryTS***.

## Citation

If you use this package in your research, please cite it as:

Vera-Valdés, J. E., (2025). LongMemory.jl: Generating, Estimating, and Forecasting Long Memory Models in Julia. Journal of Open Source Software, 10(108), 7708, [https://doi.org/10.21105/joss.07708](https://doi.org/10.21105/joss.07708)

```bibtex
@article{Vera-Valdés2025, 
author = {J. Eduardo Vera-Valdés}, 
title = {LongMemory.jl: Generating, Estimating, and Forecasting Long Memory Models in Julia}, 
journal = {Journal of Open Source Software},
doi = {10.21105/joss.07708}, 
url = {https://doi.org/10.21105/joss.07708}, 
year = {2025}, 
publisher = {The Open Journal}, 
volume = {10}, 
number = {108}, 
pages = {7708}
 }
```

## Contributing 

All types of contributions are encouraged and appreciated.

If you find a bug or have a feature request, please open a new [issue](https://github.com/everval/LongMemory.jl/issues). If you would like to contribute code, please open a [pull request](https://github.com/everval/LongMemory.jl/pulls). I welcome all contributions, including bug fixes, documentation improvements, and new features.

Thank you for considering contributing!
