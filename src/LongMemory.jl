""" 
# LongMemory.jl

This package provides functions for generating and estimating the long memory parameter d of a time series.

## Author
[J. Eduardo Vera-Valdés](https://everval.github.io/)

"""
module LongMemory

include("GeneratingFunctions.jl")
using .GeneratingFunctions
export fracdiff, csadiff, csa_gen, edm_gen, fi, fi_gen, arfi_gen, arfima_gen, fi_survival_probs


include("LogPeriodEstimators.jl")
using .LogPeriodEstimators
export gph_est, gph_est_variance, whittle_est, exact_whittle_est, whittle_est_variance, periodogram

include("ParametricEstimators.jl")
using .ParametricEstimators
export fi_mle_est, csa_mle_est, har_est, fi_var_vals, csa_var_vals, fi_cor_vals, csa_cor_vals

include("ClassicEstimators.jl")
using .ClassicEstimators
export rescaled_range_est, rescaled_range, my_std

include("Forecasters.jl")
using .Forecasters
export fi_ar_coefs, fi_forecast, csa_forecast, har_forecast

end
