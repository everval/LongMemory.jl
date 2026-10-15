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
The plotting functions require Plots.jl to be loaded:

```julia
using LongMemory, Plots
NileDataPlot()
```

## Report bugs and suggestions

Please report any bugs or suggestions [here](mailto:eduardo@math.aau.dk)

"""
module LongMemory

# The plotting functions are defined in the package extension LongMemoryPlotsExt, loaded with `using Plots`.
# Without Plots, they fall back to this function, which explains how to enable them.
function plots_required(f::Function, args::Tuple)
    if Base.get_extension(LongMemory, :LongMemoryPlotsExt) === nothing
        error("`$(nameof(f))` requires Plots.jl. Install it with `import Pkg; Pkg.add(\"Plots\")` and load it with `using Plots`.")
    end
    throw(MethodError(f, args))
end

# Submodules are included once, in dependency order; they import each other via `using ..Module`.
include("GeneratingFunctions.jl")
using .GeneratingFunctions
export fracdiff, csadiff, csa_gen, edm_gen, sds_gen, fi_gen, arfi_gen, arfima_gen, fi_survival_probs, hwfilter, hwp_gen, davies_harte

include("ClassicEstimators.jl")
using .ClassicEstimators
export smean, sstd, autocovariance, autocorrelation, autocorrelation_plot, sstdk, rescaled_range_est, rescaled_range, rescaled_range_plot, log_variance_plot, log_variance_est, lo_rs_test, vs_test, dfa_est, dfa_plot, dfa_fluctuation

include("LogPeriodEstimators.jl")
using .LogPeriodEstimators
export gph_est, gph_est_variance, whittle_est, exact_whittle_est, whittle_est_variance, exact_whittle_est_variance, periodogram, periodogram_plot, feasible_exact_whittle_est

include("DataExamples.jl")
using .DataExamples
export NileData, NHTempData, NileDataPlot, NHTempDataPlot, LMPlot

include("ParametricEstimators.jl")
using .ParametricEstimators
export fi_mle_est, csa_mle_est, har_est, fi_var_vals, csa_var_vals, fi_cor_vals, csa_cor_vals, fi_var_matrix, csa_var_matrix, fi_llk, csa_llk, durbin_levinson, fi_exact_gen, csa_exact_gen, arfima_whittle_est, arfima_var_vals, arfima_cor_vals, arfima_var_matrix, arfima_llk, arfima_mle_est, arfima_exact_gen

include("Forecasters.jl")
using .Forecasters
export fi_ar_coefs, fi_ma_coefs, fi_forecast, fi_forecast_plot, csa_ma_coefs, csa_forecast, csa_forecast_plot, har_forecast, har_forecast_plot, gregory_coeffs, hitransform, arfima_forecast, arfima_forecast_plot, arfima_ar_coefs, arfima_ma_coefs

include("StructuralChanges.jl")
using .StructuralChanges
export lm_change_test, qu_test


end
