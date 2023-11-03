""" 
# LongMemory.jl

This package provides functions for generating and estimating the long memory parameter d of a time series.

## Author
[J. Eduardo Vera-Valdés](https://everval.github.io/)

"""
module LongMemory

include("GeneratingFunctions.jl")
using .GeneratingFunctions
export fracdiff, csadiff, csagen, edmgen, fi, figen, arfigen, arfimagen

include("LogPeriodEstimators.jl")
using .LogPeriodEstimators
export gph_est, gph_est_variance, whittle_est, exact_whittle_est, whittle_est_variance, periodogram

include("ParametricEstimators.jl")
using .ParametricEstimators
export fimle_est, csamle_est

end
