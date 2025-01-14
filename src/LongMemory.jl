""" 
# LongMemory.jl

This package provides functions for generating, estimating, and forecasting long memory time series models. Moreover, the package provides plotting capabilities for the different estimators and forecasters, including classical diagnostic plots. 

## Author

[J. Eduardo Vera-Valdés](https://everval.github.io/)

## Documentation

The documentation of the package can be found [here](https://everval.github.io/LongMemory.jl/).

## Quick start

The package includes the classical Nile River and Northern Hemisphere Temperature datasets. They can be accessed by typing:

```julia
NileData()
```

or

```julia
NHTempData()
```

respectively. The package also includes the functions *NileDataPlot()* and *NHTempDataPlot()* to plot the data along with diagnostic plots. The same can be done for an arbitrary data using the *LMPlot()* function.

## Report bugs and suggestions

Please report any bugs or suggestions [here](mailto: eduardo@math.aau.dk)

"""
module LongMemory

include("GeneratingFunctions.jl")
using .GeneratingFunctions
export fracdiff, csadiff, csa_gen, edm_gen, sds_gen, fi, fi_gen, arfi_gen, arfima_gen, fi_survival_probs

include("DataExamples.jl")
using .DataExamples
export NileData, NHTempData, NileDataPlot, NHTempDataPlot, LMPlot

include("LogPeriodEstimators.jl")
using .LogPeriodEstimators
export gph_est, gph_est_variance, whittle_est, exact_whittle_est, whittle_est_variance, exact_whittle_est_variance, periodogram, periodogram_plot, periodogram

include("ParametricEstimators.jl")
using .ParametricEstimators
export fi_mle_est, csa_mle_est, har_est, fi_var_vals, csa_var_vals, fi_cor_vals, csa_cor_vals, fi_var_matrix, csa_var_matrix, fi_llk, csa_llk, my_toeplitz

include("ClassicEstimators.jl")
using .ClassicEstimators
export smean, sstd, autocovariance, autocorrelation, autocorrelation_plot, sstdk, rescaled_range_est, rescaled_range,  rescaled_range_plot, log_variance_plot, log_variance_est 

include("Forecasters.jl")
using .Forecasters
export fi_ar_coefs, fi_forecast, fi_forecast_plot, csa_forecast, csa_forecast_plot, har_forecast, har_forecast_plot

include("StructuralChanges.jl")
using .StructuralChanges
export lm_change_test



end
