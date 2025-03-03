---
title: 'LongMemory.jl: Generating, Estimating, and Forecasting Long Memory Models in Julia'
tags:
  - Julia
  - long memory
  - long-range dependence
  - fractional difference
  - ARFIMA
  - strong persistence
authors:
  - name: J. Eduardo Vera-Valdés
    orcid: 0000-0002-0337-8055
    affiliation: 1
affiliations:
 - name: Aalborg University, Department of Mathematical Sciences, Aalborg, Denmark
   index: 1
date: 14 January 2025
bibliography: paper.bib
---

# Summary

Long memory in time series analysis deals with the notion that data may have a strong dependence on past values. @Hurst1956 is one of the pioneering works on long memory. The author analyzed the flow of the Nile River and noted that water reservoirs that do not account for its long-term dynamics are at risk of overflowing. Long memory models are used in climate, finance, biology, economics, and many other fields. See @Beran2013 for a textbook on the subject.

We say that a time series $x_t$ has long memory if:
$$\gamma_x(k) \approx C_x k^{2d-1}\quad \textnormal{as}\quad k\to\infty, \label{def:cov}$$
where $\gamma_x(k)$ is the autocovariance function and $C_x$ a constant, or if:
$$f_x(\lambda)\approx C_f\lambda^{-2d}\quad \textnormal{as}\quad \lambda\to 0, \label{def:spectral}$$
where $f_x(\lambda)$ is the spectral density function and $C_f$ is a constant. Above, $g(x)\approx h(x)$ as $x\to x_0$ means that $g(x)/h(x)$ converges to $1$ as $x$ tends to $x_0$.

Both properties above can be analyzed graphically by plotting the autocorrelation and periodogram (an estimator of the spectral density), respectively.

As an example, the figure below shows the autocorrelation and periodogram (in logs) for the Nile River minima data. The data are available in `LongMemory.jl` through `NileData()`.

![Nile River minima (top), its autocorrelation function (bottom left), and log-periodogram (bottom right)](NileRiverMin.png)

As the figure shows, the autocorrelation function decays slowly and the periodogram diverges towards infinity near the origin. These are the features of long memory processes described in the definitions above. The `LongMemory.jl` package is concerned with methods for modeling data with this type of behavior.

The following code generates the figure. The `autocorrelation_plot` and `periodogram_plot` functions are part of the `LongMemory.jl` package.

```julia
using LongMemory, Plots
NileMin = NileData()
theme(:ggplot2)
p1 = plot( NileMin.Year , NileMin.NileMin, xlabel="Year", 
      ylabel = "Nile River minima" , legend = false )
p2 = autocorrelation_plot(NileMin.NileMin, 50)
p3 = periodogram_plot(NileMin.NileMin, slope = true)
l = @layout [a; b c]
plot(p1, p2, p3, layout = l, size = (700, 500) )
```

# Statement of need

`LongMemory.jl` is a package for time series long memory modeling in Julia [@juliaprimer]. The package provides functions to generate long memory, estimate model parameters, and forecast. Generating methods include fractional differencing, stochastic error duration, and cross-sectional aggregation. Estimators include the classic ones used to estimate the Hurst effect, those inspired by log-periodogram regression, and parametric ones. Forecasting is provided for all parametric estimators. Moreover, the package adds plotting capabilities to illustrate long memory dynamics and forecasting. For some of the theoretical developments, `LongMemory.jl` provides the first publicly available implementation in any programming language. A notable feature of this package is that all functions are implemented in the same programming language, taking advantage of the ease of use and speed provided by Julia. Therefore, all code is accessible to the user. Multiple dispatch, a novel feature of the language, is used to speed computations and provide consistent calls to related methods. 


# Comparison to existing packages

This section presents benchmarks contrasting different implementations to show the computational efficiency of `LongMemory.jl`. The package is related to the R [@R:lang] packages `LongMemoryTS` [@R:LongMemoryTS] and `fracdiff` [@R:fracdiff]. The former is the most similar to `LongMemory.jl` in terms of functionality, while the latter is the most popular package for long memory time series in R, measured by the number of CRAN downloads. 

The following table presents the summary of the benchmarks. The code for the benchmarks is available on the [package's website](https://everval.github.io/files/LM_notebook_benchmark.html). The benchmarks were made using `BenchmarkTools.jl` [@BenchmarkToolsJL]. 

| Function          | Mean       | Median     | Language:Package     |
|-------------------|------------|------------|----------------------|
| `fi_gen`          | 7.048E+05  | 5.812E+05  | Julia:LongMemory.jl  |
| `fracdiff.sim`    | 1.017E+08  | 1.010E+08  | R:fracdiff           |
| `fracdiff`        | 7.101E+05  | 5.824E+05  | Julia:LongMemory.jl  |
| `diffseries`      | 2.124E+06  | 1.892E+06  | R:fracdiff           |
| `fdiff`           | 4.150E+06  | 3.950E+06  | R:LongMemoryTS       |
| `gph_est`         | 4.165E+04  | 3.330E+04  | Julia:LongMemory.jl  |
| `gph`             | 1.662E+05  | 1.471E+05  | R:LongMemoryTS       |
| `fdGPH`           | 7.058E+06  | 6.022E+06  | R:fracdiff           |
: Comparison of function performance. All sample sizes are $10^4$. 

For long memory generation, the table shows the benchmarks for the function `fi_gen` in `LongMemory.jl` and the function `fracdiff.sim` in the package `fracdiff`. `LongMemoryTS` does not provide a function to directly generate processes with long memory. The results show that `LongMemory.jl` is more than $10^2$ times faster than `fracdiff` at this sample size.

Regarding fractional differencing, the table shows the benchmarks for the function `fracdiff` in `LongMemory.jl` and the functions `diffseries` and `fdiff` in packages `fracdiff` and `LongMemoryTS`, respectively. Note that `LongMemory.jl` is the fastest implementation by a large margin.

Finally, for long memory estimation, the table shows the benchmarks for the function `gph_est` in `LongMemory.jl` and the functions `fdGPH` and `gph` in packages `fracdiff` and `LongMemoryTS`, respectively. The benchmarks show that `LongMemory.jl` is significantly faster, taking advantage of the speed of Julia.

Interestingly, there does not seem to be a package in Python that provides the same functionality as `LongMemory.jl`. 

# Examples of Research Conducted with `LongMemory.jl`

The package has been used in the following research:

- Vera-Valdés, J. Eduardo, and Olivia Kvist. 2024. “Breaching 1.5°C: Give Me the Odds.” arXiv, December. https://doi.org/10.48550/arXiv.2412.13855. [@vera-valdés2024]

- Vera-Valdés, J.E., and Kvist, O. (2024). “Effects of the Paris Agreement and the COVID-19 Pandemic on Volatility Persistence of Stocks Associated with the Climate Crisis: A Multiverse Analysis.” [@VERAVALDES2024c]

Furthermore, see the accompanying examples of the package in use at the [package's website](https://everval.github.io/files/LM_notebook_illustration.html).

# References