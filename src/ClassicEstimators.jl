"""
# ClassicEstimators

This module contains functions to estimate the Hurst coefficient, closely related to the long memory parameter, of a time series using the rescaled range (R/S) statistic.

## Author

[J. Eduardo Vera-Valdés](https://everval.github.io/)

"""
module ClassicEstimators

import ..plots_required


export smean, sstd, autocovariance, autocorrelation, autocorrelation_plot, sstdk, rescaled_range_est, rescaled_range, rescaled_range_plot, log_variance_plot, log_variance_est, lo_rs_test, vs_test, dfa_est, dfa_plot, dfa_fluctuation

# Plotting functions. Their methods are defined in the package extension LongMemoryPlotsExt, which is loaded with `using Plots`.
function autocorrelation_plot end
function log_variance_plot end
function rescaled_range_plot end
function dfa_plot end

autocorrelation_plot(args...; kwargs...) = plots_required(autocorrelation_plot, args)
log_variance_plot(args...; kwargs...) = plots_required(log_variance_plot, args)
rescaled_range_plot(args...; kwargs...) = plots_required(rescaled_range_plot, args)
dfa_plot(args...; kwargs...) = plots_required(dfa_plot, args)


"""
    autocovariance(x::AbstractArray, k::Int)
    
Computes the autocovariance of a time series at lags 0, 1, ..., `k-1`.

# Arguments
- `x::Array`: The time series.
- `k::Int`: The number of lags (including lag 0).

# Output
- `acv::Array`: The autocovariances at lags 0, ..., `k-1`.

# Examples    
```julia
julia> autocovariance(randn(100), 2)
```
"""
function autocovariance(x::AbstractArray, k::Int)
    T = length(x)
    μ = smean(x)

    if k >= T
        error("k must be less than T")
    end

    acv = zeros(k)

    for ii = 0:k-1
        acv[ii+1] = sum((x[1:T-ii] .- μ) .* (x[ii+1:T] .- μ)) / (T)
    end

    return acv

end


"""
    autocorrelation(x::AbstractArray, k::Int=30)

Computes the autocorrelation function of a time series at lags 0, 1, ..., `k-1`.

# Arguments
- `x::Array`: The time series.
- `k::Int`: The number of lags (including lag 0). Default is 30.

# Output
- `acf::Array`: The autocorrelation function.

# Notes
To plot the autocorrelation function use `autocorrelation_plot`.

# Examples
```julia
julia> autocorrelation(randn(100), 10)
```
"""
function autocorrelation(x::AbstractArray, k::Int=30)
    acv = autocovariance(x, k)
    if acv[1] ≈ 0.0000
        error("The variance of the time series is zero. The autocorrelation function is not defined.")
    end
    acf = acv ./ acv[1]
    return acf
end

"""
    autocorrelation_plot(x::AbstractArray, k::Int)

Computes the autocorrelation function of a time series.

Requires Plots.jl to be loaded, e.g. with `using LongMemory, Plots`.

# Arguments
- `x::Array`: The time series.
- `k::Int`: The lag of the autocorrelation.

# Output
- `p1::Plots.Plot`: The autocorrelation function.

# Examples    
```julia
julia> autocorrelation_plot(randn(100), 10)
```
"""
autocorrelation_plot(x::AbstractArray, k::Int)


"""
    log_variance_plot(x::AbstractArray; m::Int=100, slope::Bool=false, slope2::Bool=false)
    
Produces the log-variance plot of a time series.

Requires Plots.jl to be loaded, e.g. with `using LongMemory, Plots`.

# Arguments
- `x::Array`: The time series.

# Output
- `p1::Plots.Plot`: The log-variance plot.

# Optional arguments
- `m::Int`: The number of partitions of the time series. Default is 100. A value of m ≈ length(x)/2 is recommended.
- `slope::Bool`: If true, the slope of the linear regression is displayed.
- `slope2::Bool`: If true, the theoretical slope for a short memory process is displayed.

# Notes
This function uses the linear regression method on the log of the variance plot to estimate the long memory parameter.
# Examples    
```julia
julia> log_variance_plot(randn(100,1); m = 40)
```
"""
log_variance_plot(x::AbstractArray)


"""
    log_variance_est(x::AbstractArray; m::Int=100)
    
Computes the long memory estimate using the log-variance.

# Arguments
- `x::Array`: The time series.

# Output
- `d_var::Real`: The estimated long memory parameter computed as (beta+1)/2, where beta is the slope of the linear regression of the log-variance plot.

# Optional arguments
- `m::Int`: The number of partitions of the time series. Default is 100. A value of m ≈ length(x)/2 is recommended.

# Notes
This function uses the linear regression method on the log-variance to estimate the long memory parameter.
# Examples    
```julia
julia> log_variance_est(randn(100,1); m = 40)
```
"""
function log_variance_est(x::AbstractArray; m::Int=100)
    T = length(x)
    if m >= T
        error("m must be less than T")
    end
    k = floor(Int, T / 2) - 1
    step = round(Int, T / m)

    nt = collect(1:step:k)
    cr = length(nt)

    if cr - 1 < 3
        error("At least 3 block sizes are needed for the regression; increase the length of the series or m")
    end

    Y = zeros(cr - 1, 1)
    for ii = 2:cr
        Y[ii-1, 1] = sstdk(x, nt[ii] )
    end

    X = [ones(cr-1,1) log.(nt[2:end])]
    Z = log.(Y)
    beta = X \ Z

    d_var = (beta[2] + 1) / 2

    return d_var

end

"""
    sstdk(x::AbstractArray, m::Int)

Computes the sample variance of the means of the non-overlapping blocks of size `m` of a time series.

# Arguments
- `x::Array`: The time series.
- `m::Int`: The block size.

# Output
- `s2k::Real`: The sample variance of the block means.

# Notes
Only complete blocks are used, i.e. the last `length(x) % m` observations are discarded.

# Examples    
```julia
julia> sstdk(randn(100), 20)
```
"""
function sstdk(x::AbstractArray, m::Int)
    T = length(x)

    nb = T ÷ m # number of complete blocks

    Y = zeros(nb, 1)
    for ii = 1:nb
        Y[ii, 1] = smean( x[ (ii-1)*m+1 : ii*m ] )
    end

    μ = smean(Y)

    s2k = sum( (Y .- μ) .^ 2 ) / nb

    return s2k

end


"""
    rescaled_range_est(x::AbstractArray; k::Int = 20)

Estimates the long memory parameter of a time series using the rescaled range (R/S) statistic.

# Arguments
- `x::Array`: The time series.

# Output
- `d::Real`: The estimated long memory parameter.

# Optional arguments
- `k::Int`: The number of partitions of the time series. Default is 20.

# Notes
This function uses the linear regression method on the log of the rescaled range to estimate the long memory parameter.
The Hurst coefficient is related to the long memory parameter d by the formula H = d + 1/2.

# Examples    
```julia
julia> rescaled_range_est(randn(100))
```
"""
function rescaled_range_est(x::AbstractArray; k::Int=20)
    T = length(x)

    if k >= T
        error("k must be less than T")
    end

    RS = rescaled_range(x; k=k)

    step = round(Int, T / k)

    position = collect(1:step:T)
    popfirst!(position)

    Y = log.(RS)

    X = [ones(length(Y), 1) log.(position)]

    beta = X \ Y

    return beta[2] - 1 / 2
end


"""
    rescaled_range_plot(x::AbstractArray; k::Int = 100, slope::Bool = false, slope2::Bool = false)

Produces the rescaled range plot of a time series.

Requires Plots.jl to be loaded, e.g. with `using LongMemory, Plots`.

# Arguments
- `x::Array`: The time series.

# Output
- `p1::Plots.Plot`: The rescaled range plot.

# Optional arguments
- `k::Int`: The number of partitions of the time series. Default is 100.
- `slope::Bool`: If true, the slope of the linear regression is displayed.
- `slope2::Bool`: If true, the theoretical slope for a short memory process is displayed.

# Notes
This function uses the linear regression method on the log of the rescaled range to estimate the long memory parameter.
The Hurst coefficient is related to the long memory parameter d by the formula H = d + 1/2.

# Examples    
```julia
julia> rescaled_range_plot(randn(300,1))
```
"""
rescaled_range_plot(x::AbstractArray)


"""
    rescaled_range(x::AbstractArray; k::Int = 20)

Computes the rescaled range (R/S) statistic of a time series on the subsamples `x[1:n]` for an increasing grid of `n`.

# Arguments
- `x::Array`: The time series.

# Optional arguments
- `k::Int`: The number of partitions of the time series. Default is 20.

# Output
- `RS::Array`: The rescaled range statistic.

# Examples    
```julia
julia> rescaled_range(randn(100))
```
"""
function rescaled_range(x::AbstractArray; k::Int=20)
    T = length(x)

    if k >= T
        error("k must be less than T")
    end

    step = round(Int, T / k)

    position = collect(1:step:T)
    popfirst!(position)

    RS = zeros(length(position))

    for ii = 1:length(position)
        xn = x[1:position[ii]]
        Z = cumsum(vec(xn) .- smean(xn)) # deviations from the subsample mean
        RS[ii] = (maximum(Z) - minimum(Z)) / sstd(xn)
    end

    return RS
end

"""
    sstd(x::AbstractArray; k::Int=0)

Computes the sample standard deviation of a time series.

# Arguments
- `x::Array`: The time series.

# Output
- `std::Real`: The sample standard deviation.

# Optional arguments
- `k::Int`: The bias-correction term.

# Notes
This function divides by T instead of T-1 if no bias-correction term is provided.

# Examples    
```julia
julia> sstd(randn(100))
```
"""
function sstd(x::AbstractArray; k::Int=0)
    T = length(x)
    μ = smean(x)
    return sqrt(sum((x .- μ) .^ 2) / (T - k))
end

"""
    smean(x::AbstractArray)

Computes the sample mean of a time series.

# Arguments
- `x::Array`: The time series.

# Output
- `mean::Real`: The sample mean.

# Examples    
```julia
julia> smean(randn(100))
```
"""
function smean(x::AbstractArray)
    T = length(x)
    return sum(x) / T
end


"""
    lo_rs_test(x::AbstractArray; q::Int = andrews_bandwidth(x))

Computes Lo's (1991) modified rescaled range test of the null hypothesis of short memory.

# Arguments
- `x::Array`: The time series.

# Optional arguments
- `q::Int`: The number of lags in the Bartlett estimator of the long-run variance. By default, it is chosen with the data-dependent rule of Andrews (1991) as suggested by Lo (1991), see `andrews_bandwidth`.

# Output
- `V::Real`: The test statistic, `R / (σ̂(q) √T)`, where `R` is the range of the partial sums of the demeaned series and `σ̂²(q)` the Bartlett estimator of the long-run variance.
- `pval::Real`: The two-sided p-value from the asymptotic distribution under the null, the range of a Brownian bridge.

# Notes
Under the null of short memory, `V` converges to the range of a Brownian bridge, with 2.5% and 97.5% quantiles 0.809 and 1.862. Large values of `V` indicate long memory and small values antipersistence.
In contrast to the classic R/S statistic, the test is robust to short-range dependence through the long-run variance estimator.

# References
Lo, A. W. (1991). Long-term memory in stock market prices. Econometrica, 59(5), 1279-1313.

# Examples
```julia
julia> lo_rs_test(randn(500))
```
"""
function lo_rs_test(x::AbstractArray; q::Int = andrews_bandwidth(x))
    T = length(x)

    S = cumsum(vec(x) .- smean(x))
    R = maximum(S) - minimum(S) # S[T] = 0, so the range includes zero

    V = R / (sqrt(bartlett_lrv(x, q)) * sqrt(T))

    F = brownian_bridge_range_cdf(V)
    pval = 2 * min(F, 1 - F)

    return V, pval
end


"""
    vs_test(x::AbstractArray; q::Int = andrews_bandwidth(x))

Computes the V/S test of Giraitis, Kokoszka, Leipus and Teyssière (2003) of the null hypothesis of short memory.

# Arguments
- `x::Array`: The time series.

# Optional arguments
- `q::Int`: The number of lags in the Bartlett estimator of the long-run variance. By default, it is chosen with the data-dependent rule of Andrews (1991), see `andrews_bandwidth`.

# Output
- `M::Real`: The test statistic, the sample variance of the partial sums of the demeaned series divided by `T σ̂²(q)`.
- `pval::Real`: The p-value from the asymptotic distribution under the null.

# Notes
Under the null of short memory, `M` converges to the variance of a Brownian bridge, with distribution function `1 + 2 Σ_{k≥1} (-1)^k exp(-2k²π²x)`. The 5% critical value is 0.1869; large values indicate long memory.
The test is a modification of the KPSS test that uses the variance of the partial sums instead of their second moment, which makes it more sensitive to shifts in the variance of the partial sums.

# References
Giraitis, L., Kokoszka, P., Leipus, R., and Teyssière, G. (2003). Rescaled variance and related tests for long memory in volatility and levels. Journal of Econometrics, 112(2), 265-294.

# Examples
```julia
julia> vs_test(randn(500))
```
"""
function vs_test(x::AbstractArray; q::Int = andrews_bandwidth(x))
    T = length(x)

    S = cumsum(vec(x) .- smean(x))

    M = (sum(S .^ 2) - sum(S)^2 / T) / (T^2 * bartlett_lrv(x, q))

    pval = 1 - brownian_bridge_variance_cdf(M)

    return M, pval
end


"""
    dfa_est(x::AbstractArray; k::Int = 20)

Estimates the long memory parameter of a time series using detrended fluctuation analysis (DFA).

# Arguments
- `x::Array`: The time series.

# Optional arguments
- `k::Int`: The number of box sizes, log-spaced between 10 and `T/4`. Default is 20.

# Output
- `d::Real`: The estimated long memory parameter.

# Notes
The profile (the cumulative sum of the demeaned series) is split into non-overlapping boxes of size `n`, a linear trend is fitted by OLS in each box, and the fluctuation `F(n)` is the root mean squared residual.
The slope `α` of the regression of `log F(n)` on `log n` estimates the Hurst coefficient `H = d + 1/2`, so `d = α - 1/2`. See `dfa_fluctuation` for the fluctuation function.

# References
Peng, C.-K., Buldyrev, S. V., Havlin, S., Simons, M., Stanley, H. E., and Goldberger, A. L. (1994). Mosaic organization of DNA nucleotides. Physical Review E, 49(2), 1685-1689.

# Examples
```julia
julia> dfa_est(randn(500))
```
"""
function dfa_est(x::AbstractArray; k::Int = 20)
    nt, F = dfa_fluctuation(x; k=k)

    X = [ones(length(nt), 1) log.(nt)]
    Y = log.(F)

    beta = X \ Y

    return beta[2] - 1 / 2
end


"""
    dfa_plot(x::AbstractArray; k::Int = 20, slope::Bool = false, slope2::Bool = false)

Produces the detrended fluctuation analysis (DFA) plot of a time series.

Requires Plots.jl to be loaded, e.g. with `using LongMemory, Plots`.

# Arguments
- `x::Array`: The time series.

# Output
- `p1::Plots.Plot`: The DFA plot.

# Optional arguments
- `k::Int`: The number of box sizes. Default is 20.
- `slope::Bool`: If true, the slope of the linear regression is displayed.
- `slope2::Bool`: If true, the theoretical slope for a short memory process is displayed.

# Notes
The slope of the regression estimates the Hurst coefficient `H = d + 1/2`. See `dfa_est`.

# Examples
```julia
julia> dfa_plot(randn(500))
```
"""
dfa_plot(x::AbstractArray)


"""
    dfa_fluctuation(x::AbstractArray; k::Int = 20)

Computes the fluctuation function of detrended fluctuation analysis (DFA) with linear detrending.

# Arguments
- `x::Array`: The time series.

# Optional arguments
- `k::Int`: The number of box sizes, log-spaced between 10 and `T/4`. Default is 20.

# Output
- `nt::Array`: The box sizes.
- `F::Array`: The fluctuation function at the box sizes.

# Examples
```julia
julia> dfa_fluctuation(randn(500))
```
"""
function dfa_fluctuation(x::AbstractArray; k::Int = 20)
    T = length(x)

    nt = unique(round.(Int, exp.(range(log(10), log(T / 4), length=k))))

    if T < 40 || length(nt) < 3
        error("The series is too short for DFA; at least 3 box sizes between 10 and T/4 are needed")
    end

    Y = cumsum(vec(x) .- smean(x)) # profile

    F = zeros(length(nt))
    for ii in eachindex(nt)
        n = nt[ii]
        nb = T ÷ n
        tt = [ones(n) collect(1:n)]

        res2 = 0.0
        for jj = 1:nb
            yb = Y[(jj-1)*n+1:jj*n]
            β = tt \ yb
            res2 += sum((yb - tt * β) .^ 2)
        end
        F[ii] = sqrt(res2 / (nb * n))
    end

    return nt, F
end


## Internal functions for the short memory tests
"""
    andrews_bandwidth(x::AbstractArray)

Computes the data-dependent number of lags of Andrews (1991) for the Bartlett long-run variance estimator, based on an AR(1) approximation: `floor((3T/2)^{1/3} |2ρ/(1-ρ²)|^{2/3})`, where `ρ` is the first-order autocorrelation. Internal function used by `lo_rs_test` and `vs_test`.

# References
Andrews, D. W. K. (1991). Heteroskedasticity and autocorrelation consistent covariance matrix estimation. Econometrica, 59(3), 817-858.
"""
function andrews_bandwidth(x::AbstractArray)
    T = length(x)
    ρ = autocorrelation(x, 2)[2]

    q = floor(Int, (3 * T / 2)^(1 / 3) * abs(2 * ρ / (1 - ρ^2))^(2 / 3))

    return min(q, T ÷ 2)
end

"""
    bartlett_lrv(x::AbstractArray, q::Int)

Computes the Bartlett (Newey-West) estimator of the long-run variance of `x` with `q` lags. Internal function used by `lo_rs_test` and `vs_test`.
"""
function bartlett_lrv(x::AbstractArray, q::Int)
    acv = autocovariance(x, q + 1)

    s2 = acv[1]
    for jj = 1:q
        s2 += 2 * (1 - jj / (q + 1)) * acv[jj+1]
    end

    return s2
end

"""
    brownian_bridge_range_cdf(v::Real)

Distribution function of the range of a Brownian bridge (Kuiper's distribution), `1 + 2 Σ_{k≥1} (1 - 4k²v²) exp(-2k²v²)`. Internal function used by `lo_rs_test`.
"""
function brownian_bridge_range_cdf(v::Real)
    if v <= 0
        return 0.0
    end

    F = 1.0
    for k = 1:1000
        F += 2 * (1 - 4 * k^2 * v^2) * exp(-2 * k^2 * v^2)
    end

    return min(max(F, 0.0), 1.0)
end

"""
    brownian_bridge_variance_cdf(x::Real)

Distribution function of the variance of a Brownian bridge, `1 + 2 Σ_{k≥1} (-1)^k exp(-2k²π²x)`. Internal function used by `vs_test`.
"""
function brownian_bridge_variance_cdf(x::Real)
    if x <= 0
        return 0.0
    end

    F = 1.0
    for k = 1:1000
        F += 2 * (-1)^k * exp(-2 * k^2 * π^2 * x)
    end

    return min(max(F, 0.0), 1.0)
end


end # module