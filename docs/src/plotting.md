```@meta
CurrentModule = LongMemory
```

# Plots for Long Memory

Functions to plot the autocorrelation function, log-variance, log-periodogram, rescaled range, and detrended fluctuation analysis. The plots help to visualize the presence of long memory in the data. The `LMPlot` function plots the time series, autocorrelation function, log-periodogram, and log-variance in a single figure.

Since version 2.0, the plotting functions are provided through a package extension: they are available once [Plots.jl](https://docs.juliaplots.org/stable/) is loaded. Install Plots once with `import Pkg; Pkg.add("Plots")` and then use

```julia
using LongMemory, Plots
LMPlot(fi_gen(1000, 0.3))
```

The forecasting plots, such as `fi_forecast_plot`, work in the same way.

```@docs
autocorrelation_plot
log_variance_plot
periodogram_plot
rescaled_range_plot
dfa_plot
LMPlot
```

Documentation for [LongMemory.jl](https://github.com/everval/LongMemory.jl).
