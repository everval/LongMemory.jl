"""
# StructuralChanges

This module contains functions to estimate structural changes in models contaning long memory errors.

## Author

[J. Eduardo Vera-Valdés](https://everval.github.io/)

"""
module StructuralChanges

include("LogPeriodEstimators.jl")
import .LogPeriodEstimators: exact_whittle_est

include("GeneratingFunctions.jl")
import .GeneratingFunctions: fracdiff

export lm_change_test


function lm_change_test(y::Array)
    T = length(y)

    # Estimate the long memory parameter
    d0 = exact_whittle_est(y)

    # Integrating the series
    x = fracdiff(y, d0)

    return x

end



end # module