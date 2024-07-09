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

using FFTW

export lm_change_test

"""
    lm_change_test(y::Array; uplim::Real = 0.15, lowlim::Real = 0.15)

Estimates the location of a structural change in a long memory model.

# Arguments
- `y::Array`: The series to be tested.

# Optional arguments
- `uplim::Real = 0.15`: The upper limit of the fraction of the series to be tested.
- `lowlim::Real = 0.15`: The lower limit of the fraction of the series to be tested.

# Output
- `τ::Int`: The estimated location of the structural change.

# Examples
```julia-repl
julia> lm_change_test(randn(100))
```

# Notes
The function estimates the location of a structural change in a long memory model. The function uses the Whittle estimator to estimate the long memory parameter. Then, it integrates the series and computes the t-statistics of the OLS regression of the integrated series on the starred series. The function computes the forward and backward t-statistics and returns the location of the maximum squared t-statistic. Hence, it is robust to the direction of the change.

# References
Martins and Rodrigues (2014), "Testing for persistence change in fractionally integrated models: An application to world inflation rates", Computational Statistics and Data Analysis, 76.
"""
function lm_change_test(y::Array; uplim::Real = 0.15, lowlim::Real = 0.15)
    T = length(y)
    Λₗ = round(Int, T*lowlim)
    Λᵤ = round(Int, T*(1-uplim))
    
    # Estimate the long memory parameter
    d0 = exact_whittle_est(y)

    # Integrating the series
    x = fracdiff(y, d0)
    xstar = starred_var(x)

    # Forward direction
    tvalsup = zeros(Λᵤ - Λₗ + 1, 1)
    for ii in Λₗ:Λᵤ
        tvalsup[ii-Λₗ+1] = simple_ols(x[2:ii], xstar[1:(ii-1)])
    end

    # Backward direction
    w = reverse(x)
    wstar = starred_var(w)
    tvalsdown = zeros(Λᵤ - Λₗ + 1, 1)
    for ii in Λₗ:Λᵤ
        tvalsdown[ii-Λₗ+1] = simple_ols(w[2:ii], wstar[1:(ii-1)])
    end

    tvals = [tvalsup.^2 tvalsdown.^2]
 
    return findmax(tvals)[2][1] + Λₗ - 1
end

## Internal functions
"""
    starred_var(x::Array)

Computes the padded series `x` with detrending coefficients.

# Arguments
- `x::Array`: The series to be detrended.

# Output
- `dx::Array`: The detrended series.
"""
function starred_var(x::Array)
    T = length(x)

    np2 = nextpow(2, 2 * T - 1)

    coefs = collect(1 ./ (1:(T-1)))

    padcoefs = [coefs; zeros(np2 - T + 1, 1)]
    padx = [x[1:(T-1)]; zeros(np2 - T + 1, 1)]

    dx = irfft(rfft(padx) .* rfft(padcoefs), np2)
    return dx[1:(T-1)]
end

"""
    simple_ols(y::Array, x::Array)

Computes the t-statistics of the OLS regression of `y` on `x`.

# Arguments
- `y::Array`: The dependent variable.
- `x::Array`: The independent variable.

# Output
- `tstat::Array`: The t-statistic of the regression.

# Notes
The function should only be used internally. It is not exported. It is used to compute the t-statistics of the OLS regression in the `lm_change_test` function. Hence, it does not adds intercepts to the regression nor it returns the coefficient.
"""
function simple_ols(y::Array, x::Array)
    T = length(y)
    β = x\y

    err = y-x*β
    σ² = (err'*err)/(T-1)
    tstat = β/sqrt(σ²/(x'*x))

    return tstat
end

end # module