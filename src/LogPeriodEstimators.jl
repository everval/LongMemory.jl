"""
# LogPeriodEstimators

This module contains functions to estimate the long memory parameter of a time series using the log-periodogram estimator and the Whittle log-likelihood function.

## Author

[J. Eduardo Vera-Valdés](https://everval.github.io/)

"""
module LogPeriodEstimators

include("GeneratingFunctions.jl")
import .GeneratingFunctions: fracdiff

using FFTW, Optim, Plots

export gph_est, gph_est_variance, whittle_est, exact_whittle_est, whittle_est_variance, periodogram, periodogram_plot


"""
    periodogram_plot(x::Array; slope::Bool=true)

Plots the log-periodogram of a time series `x`.

# Arguments
- `x::Vector`: time series

# Output
- `p1::Plots.Plot`: log-periodogram

# Optional arguments
- `slope::Bool`: If true, the slope of the linear regression is displayed.

# Examples
```julia-repl
julia> periodogram_plot(randn(100,1))
```
"""
function periodogram_plot(x::Array; slope::Bool=false)
    T = length(x)

    I_wc, wc = periodogram(x)

    I_w = I_wc[2:end]
    w = wc[2:end]

    if slope == false
        p1 = plot(w, I_w, xlabel="log-frequency", ylabel="log-periodogram", label="", xaxis=:log, yaxis=:log, line=:scatter)
    end

    if slope == true
        p1 = plot(log.(w), log.(I_w), xlabel="log-frequency", ylabel="log-periodogram", label="", line=:scatter)
        oldylims = ylims(p1)
        X = [ones(length(w)) log.(w)]
        beta = X \ log.(I_w)
        plot!(log.(w), X*beta, line=:dash, label=string("Slope = ", round(beta[2],digits=4)), linewidth = 3)
        ylims!(p1, oldylims)
    end

    return p1
end


"""
    periodogram(x::Array)

Compute the periodogram of a time series `x` using the fast Fourier transform.

# Arguments
- `x::Vector`: time series

# Output
- `I_w::Vector`: periodogram
- `w::Vector`: Fourier frequencies

# Examples
```julia-repl
julia> periodogram(randn(100,1))
```
"""
function periodogram(x::Array)
    T = length(x)

    I_w = abs.(rfft(x)) .^ 2 ./ T
    w = 2 * π * (0:T-1) ./ T

    ind = iseven(T) ? round(Int, T / 2 + 1) : ceil(Int, T / 2)
    I_w, w = I_w[1:ind], w[1:ind]
    return I_w, w
end


"""
    gph_est(x::Array; m=0.5, l=0, br=0::Int)

Estimate the long memory parameter of a time series `x` using the log-periodogram estimator. See [Geweke and Porter-Hudak (1983)](https://onlinelibrary.wiley.com/doi/10.1111/j.1467-9892.1983.tb00371.x) and [Andrews and Guggenberger (2003)](https://www.jstor.org/stable/3082070) for details.

# Arguments
- `x::Vector`: time series
- `m∈(0,1)::Float64`: taper final
- `l∈(0,1)::Float64`: taper initial
- `br::Int64`: number of bias reduction terms

# Output
- `d::Float64`: long memory parameter

# Notes
The function considers the periodogram of the time series `x` for frequencies in the interval `[T^l,T^m]`. The zero frequency is always excluded.
The default values of `m` and `l` are 0.5 and 0, respectively.
The condition `m < l` must hold.

The default value of `br` is 0 which returns the original GPH log-periodogram estimator.

# Examples
```julia-repl
julia> gph_est(randn(100,1))
```
"""
function gph_est(x::Array; m=0.5, l=0, br=0::Int)
    T = length(x)

    if m < l
        error("Taper initial is greater than final")
    end

    last = round(Int, T^m)
    first = max(round(Int, T^l), 2)

    I_wc, wc = periodogram(x)

    I_w = I_wc[first:last]
    w = wc[first:last]

    Y = log.(I_w)
    X = [-2 * log.(w) w .^ collect(0:2:(2*br))']
    β = X \ Y

    return β[1]
end


"""
    gph_est_variance(x::Array; m=0.5, l=0, br=0::Int)

Estimate the variance of the long memory parameter of a time series `x` using the log-periodogram estimator. See [Geweke and Porter-Hudak (1983)](https://onlinelibrary.wiley.com/doi/10.1111/j.1467-9892.1983.tb00371.x) and [Andrews and Guggenberger (2003)](https://www.jstor.org/stable/3082070) for details.

# Arguments
- `x::Vector`: time series

# Optional arguments
- `m∈(0,1)::Float64`: taper final
- `br::Int64`: number of bias reduction terms

# Output
- `varb::Float64`: variance of the long memory parameter

# Notes
Multiple dispatch is used for computation. If the first input is an integer, the function interprets it as the sample size; otherwise, it computes the sample size from the length of the time series.

# Examples
```julia-repl
julia> gph_est_variance(fi(100,0.4))
```
"""
function gph_est_variance(x::Array; m=0.5, br=0::Int)
    T = length(x)
    last = round(Int, T^m)

    if br == 0
        cr = 1
    elseif br == 1
        cr = 9 / 4
    elseif br == 2
        cr = 3.52
    elseif br == 3
        cr = 4.79
    elseif br == 4
        cr = 6.06
    elseif br > 4
        cr = 6.06
        @warn "Variance inflation factor for bias reduction terms greater than 4 are not available. Using 4. It will greatly underestimate the variance."
    else
        cr = 1
        @warn "Invalid number of bias reduction terms. Using 0 for the variance inflation factor."
    end

    varb = cr * (π^2 / 24) / last
    return varb
end

"""
    gph_est_variance(T::Int; m=0.5, l=0, br=0::Int)

Estimate the variance of the long memory parameter of a time series of length `T` using the log-periodogram estimator. See [Geweke and Porter-Hudak (1983)](https://onlinelibrary.wiley.com/doi/10.1111/j.1467-9892.1983.tb00371.x) and [Andrews and Guggenberger (2003)](https://www.jstor.org/stable/3082070) for details.

# Arguments
- `T::Int`: length of the time series

# Optional arguments
- `m∈(0,1)::Float64`: taper final
- `br::Int64`: number of bias reduction terms

# Output
- `varb::Float64`: variance of the long memory parameter

# Notes
Multiple dispatch is used for computation. If the first input is an integer, the function interprets it as the sample size; otherwise, it computes the sample size from the length of the time series.

# Examples
```julia-repl
julia> gph_est_variance(100,0.4)
```
"""
function gph_est_variance(T::Int; m=0.5, br=0::Int)
    last = round(Int, T^m)

    if br == 0
        cr = 1
    elseif br == 1
        cr = 9 / 4
    elseif br == 2
        cr = 3.52
    elseif br == 3
        cr = 4.79
    elseif br == 4
        cr = 6.06
    elseif br > 4
        cr = 6.06
        @warn "Variance inflation factor for bias reduction terms greater than 4 are not available. Using 4. It will greatly underestimate the variance."
    else
        cr = 1
        @warn "Invalid number of bias reduction terms. Using 0 for the variance inflation factor."
    end

    varb = cr * (π^2 / 24) / last
    return varb
end

"""
    whittle_llk(d, x::Array; m=0.5, l=0)

Compute the Whittle log-likelihood function of a time series `x` for a given long memory parameter `d`. See Künsch (1987) for details.

# Arguments
- `d::Float64`: long memory parameter
- `x::Vector`: time series
- `m∈(0,1)::Float64`: taper final
- `l∈(0,1)::Float64`: taper initial

# Output
- `Q::Float64`: Whittle log-likelihood function

# Notes
The function considers the periodogram of the time series `x` for frequencies in the interval `[T^l,T^m]`. The zero frequency is always excluded.
The condition `m < l` must hold. 
The default values of `m` and `l` are 0.5 and 0, respectively.

# Examples
```julia-repl
julia> whittle_llk(0.4,randn(100,1))
```
"""
function whittle_llk(d, x::Array; m=0.5, l=0)

    T = length(x)

    if m < l
        error("Taper initial is greater than final")
    end

    first = max(round(Int, T^l), 2)
    last = round(Int, T^m)

    I_w, w = periodogram(x)

    I_w = I_w[first:last]
    w = w[first:last]

    G = sum(I_w .* (w .^ (2 * d))) / length(w)
    Q = log(G) - 2 * d * sum(log.(w)) / length(w)

    return Q
end


"""
    whittle_est(x::Array; m=0.5, l=0)

Estimate the long memory parameter of a time series `x` using the Whittle log-likelihood function. See Künsch (1987) for details.

# Arguments
- `x::Vector`: time series
- `m∈(0,1)::Float64`: taper final
- `l∈(0,1)::Float64`: taper initial

# Output
- `d::Float64`: long memory parameter

# Notes
The function considers the periodogram of the time series `x` for frequencies in the interval `[T^l,T^m]`. The zero frequency is always excluded.
The condition `m < l` must hold.
The default values of `m` and `l` are 0.5 and 0, respectively.

# Examples
```julia-repl
julia> whittle_est(randn(100,1))
```
"""
function whittle_est(x::Array; m=0.5, l=0)
    d0 = gph_est(x; m=m, l=l)
    whittle = optimize(d -> whittle_llk(first(d), x; m=m, l=l), [d0])

    return whittle.minimizer[1]
end

"""
    whittle_est_variance(x::Array; m=0.5)

Estimate the variance of the estimator for the long memory parameter of a time series `x` using the Whittle log-likelihood function. See Künsch (1987) for details.

# Arguments
- `x::Vector`: time series

# Optional arguments
- `m∈(0,1)::Float64`: taper final

# Output
- `varb::Float64`: variance of the estimator

# Notes
Multiple dispatch is used for computation. If the first input is an integer, the function interprets it as the sample size; otherwise, it computes the sample size from the length of the time series.
    The variance is the same as the one from using the exact Whittle log-likelihood function.

# Examples
```julia-repl
julia> whittle_est_variance(fi(100,0.4))
```
"""
function whittle_est_variance(x::Array; m=0.5)
    T = length(x)
    last = round(Int, T^m)

    cr = 1

    varb = cr * (1 / 4) / last
    return varb
end

""" 
    whittle_est_variance(T::Int;m=0.5

Estimate the variance of the estimator for the long memory parameter of a time series of length `T` using the Whittle log-likelihood function. See Künsch (1987) for details.

# Arguments
- `T::Int`: length of the time series

# Optional arguments
- `m∈(0,1)::Float64`: taper final

# Output
- `varb::Float64`: variance of the estimator

# Notes
Multiple dispatch is used for computation. If the first input is an integer, the function interprets it as the sample size; otherwise, it computes the sample size from the length of the time series.
The variance is the same as the one from using the exact Whittle log-likelihood function.

# Examples
```julia-repl
julia> whittle_est_variance(100,0.4)
```
"""
function whittle_est_variance(T::Int; m=0.5)
    last = round(Int, T^m)

    cr = 1

    varb = cr * (1 / 4) / last
    return varb
end


"""
    exact_whittle_llk(d, x::Array; m=0.5, l=0)

Compute the exact Whittle log-likelihood function of a time series `x` for a given long memory parameter `d`. See [Shimotsu and Phillips (2005)](https://doi.org/10.1214/009053605000000309) for details.

# Arguments
- `d::Float64`: long memory parameter
- `x::Vector`: time series
- `m∈(0,1)::Float64`: taper final
- `l∈(0,1)::Float64`: taper initial

# Output
- `Q::Float64`: Whittle log-likelihood function

# Notes
The function considers the periodogram of the time series `x` for frequencies in the interval `[T^l,T^m]`. The zero frequency is always excluded.
The condition `m < l` must hold.
The default values of `m` and `l` are 0.5 and 0, respectively.

# Examples
```julia-repl
julia> exact_whittle_llk(0.4,randn(100,1))
```
"""
function exact_whittle_llk(d, x::Array; m=0.5, l=0)
    T = length(x)

    if m < l
        error("Taper initial is greater than final")
    end

    last = round(Int, T^m)
    first = max(round(Int, T^l), 2)

    dx = fracdiff(x, d)

    I_w, w = periodogram(dx)

    I_w = I_w[first:last]
    w = w[first:last]

    G = sum(I_w) / length(w)
    Q = log(G) - 2 * d * sum(log.(w)) / length(w)

    return Q
end


"""
    exact_whittle_est(x::Array; m=0.5, l=0)

Estimate the long memory parameter of a time series `x` using the exact Whittle log-likelihood function. See [Shimotsu and Phillips (2005)](https://doi.org/10.1214/009053605000000309) for details.

# Arguments
- `x::Vector`: time series
- `m∈(0,1)::Float64`: taper final
- `l∈(0,1)::Float64`: taper initial

# Output
- `d::Float64`: long memory parameter

# Notes
The function considers the periodogram of the time series `x` for frequencies in the interval `[T^l,T^m]`. The zero frequency is always excluded.
The condition `m < l` must hold.
The default values of `m` and `l` are 0.5 and 0, respectively.

# Examples
```julia-repl
julia> exact_whittle_est(randn(100,1))
```
"""
function exact_whittle_est(x::Array; m=0.5, l=0)
    d0 = gph_est(x; m=m, l=l)
    whittle = optimize(d -> exact_whittle_llk(first(d), x; m=m, l=l), [d0])

    return whittle.minimizer[1]
end

end # module