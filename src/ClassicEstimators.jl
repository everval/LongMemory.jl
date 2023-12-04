"""
# ClassicEstimators

This module contains functions to estimate them Hurst coefficient, closely related to the long memory parameter, of a time series using the rescaled range (R/S) statistic.

## Author

[J. Eduardo Vera-Valdés](https://everval.github.io/)

"""
module ClassicEstimators

using Plots

export smean, sstd, autocovariance, autocorrelation, autocorrelation_plot, sstdk, rescaled_range_est, rescaled_range, variance_plot, variance_plot_est


"""
    autocovariance(x::Array, k::Int)
    
Computes the autocovariance of a time series.

# Arguments
- `x::Array`: The time series.
- `k::Int`: The lag of the autocovariance.

# Output
- `acv::Array`: The autocovariance.

# Examples    
```julia
julia> autocovariance(randn(100), 2)
```
"""
function autocovariance(x::Array, k::Int)
    T = length(x)
    μ = smean(x)

    if k >= T
        error("k must be less than T")
    end

    acv = zeros(k, 1)

    for ii = 0:k-1
        acv[ii+1, 1] = sum((x[1:T-ii] .- μ) .* (x[ii+1:T] .- μ)) / (T)
    end

    return acv

end


"""
    autocorrelation(x::Array, k::Int; flag::Bool=false)

Computes the autocorrelation function of a time series.

# Arguments
- `x::Array`: The time series.
- `k::Int`: The lag of the autocorrelation.

# Output
- `acf::Array`: The autocorrelation function.

# Optional arguments
- `flag::Bool`: If true, the autocorrelation function is displayed.

# Examples    
```julia
julia> autocorrelation(randn(100), 10)
```
"""
function autocorrelation(x::Array, k::Int=30; flag::Bool=false)
    acv = autocovariance(x, k)
    acf = acv ./ acv[1]
    if flag == true
        display(plot(acf, line=:stem, xlabel="Lags", ylabel="Autocorrelation function", legend=false, marker=:circle))
    end
    return acf
end

"""
    autocorrelation_plot(x::Array, k::Int)

Computes the autocorrelation function of a time series.

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
function autocorrelation_plot(x::Array, k::Int=30)
    acv = autocovariance(x, k)
    acf = acv ./ acv[1]

    p1 = plot(acf, line=:stem, xlabel="Lags", ylabel="Autocorrelation function", legend=false, marker=:circle)

    return p1
end


"""
    variance_plot(x::Array; m::Int=20, slope::Bool=false, slope2::Bool=false)
    
Produces the variance plot of a time series.

# Arguments
- `x::Array`: The time series.

# Output
- `d_var::Real`: The estimated long memory parameter computed as (beta+1)/2, where beta is the slope of the linear regression of the log of the variance plot.

# Optional arguments
- `m::Int`: The number of partitions of the time series. Default is 20. A value of m ≈ length(x)/2 is recommended.
- `slope::Bool`: If true, the slope of the linear regression is displayed.
- `slope2::Bool`: If true, the theoretical slope for a short memory process is displayed.

# Notes
This function uses the linear regression method on the log of the variance plot to estimate the long memory parameter.
# Examples    
```julia
julia> variance_plot(randn(100,1); m = 40)
```
"""
function variance_plot(x::Array; m::Int=20, slope::Bool=false, slope2::Bool=false)
    T = length(x)
    if m >= T
        error("m must be less than T")
    end
    k = floor(Int, T / 2) - 1
    step = round(Int, T / m)

    nt = collect(1:step:k)
    cr = length(nt)

    Y = zeros(cr - 1, 1)
    for ii = 2:cr
        Y[ii-1, 1] = sstdk(x, nt[ii] )
    end

    X = [ones(cr-1,1) log.(nt[2:end])]
    Z = log.(Y)
    beta = X \ Z

    p1 = plot(X[:,2], Z, line=:scatter, label="", xlabel="log-sampling", ylabel="log-variance")
    if slope == true
        oldylims = ylims(p1)
        plot!(X[:,2], X * beta, line=:dash, label=string("Slope = ", round(beta[2], digits=4)), linewidth=3)
        ylims!(p1, oldylims)
    end
    if slope2 == true
        oldylims = ylims(p1)
        plot!(X[:,2], X *  [beta[1]; -1], line=:dashdot, label="Slope = -1", linewidth=3)
        ylims!(p1, oldylims)
    end
    return p1

end


"""
    variance_plot_est(x::Array; m::Int=20)
    
Computes the long memory estimate using the variance plot.

# Arguments
- `x::Array`: The time series.

# Output
- `d_var::Real`: The estimated long memory parameter computed as (beta+1)/2, where beta is the slope of the linear regression of the log of the variance plot.

# Optional arguments
- `m::Int`: The number of partitions of the time series. Default is 20. A value of m ≈ length(x)/2 is recommended.

# Notes
This function uses the linear regression method on the log of the variance plot to estimate the long memory parameter.
# Examples    
```julia
julia> variance_plot_est(randn(100,1); m = 40)
```
"""
function variance_plot_est(x::Array; m::Int=20)
    T = length(x)
    if m >= T
        error("m must be less than T")
    end
    k = floor(Int, T / 2) - 1
    step = round(Int, T / m)

    nt = collect(1:step:k)
    cr = length(nt)

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
    sstdk(x::Array, m::Int)

Computes the sample standard deviation of a time series using the m-th order sample mean.

# Arguments
- `x::Array`: The time series.
- `m::Int`: The order of the sample mean.

# Output
- `std::Real`: The sample standard deviation.

# Examples    
```julia
julia> sstdk(randn(100), m = 20)
```
"""
function sstdk(x::Array, m::Int)
    T = length(x)

    nt = collect(1:m:T)
    cr = length(nt)

    Y = zeros(cr - 1, 1)
    for ii = 1:(cr-1)
        Y[ii, 1] = smean( x[ nt[ii] : (nt[ii+1]-1) ] )
    end

    μ = smean(Y)

    s2k = sum( (Y .- μ) .^ 2 ) / (cr - 1)

    return s2k

end


"""
    rescaled_range_est(x::Array; k::Int = 20)

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
function rescaled_range_est(x::Array; k::Int=20)
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
    rescaled_range(x::Array; k::Int = 20)

Computes the rescaled range (R/S) statistic of a time series.

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
function rescaled_range(x::Array; k::Int=20)
    T = length(x)

    if k >= T
        error("k must be less than T")
    end

    μ = smean(x)

    Y = x .- μ
    Z = cumsum(Y, dims=1)

    step = round(Int, T / k)

    position = collect(1:step:T)
    popfirst!(position)

    RS = zeros(length(position), 1)

    for ii = 1:length(position)
        RS[ii, 1] = (maximum(Z[1:position[ii]]) - minimum(Z[1:position[ii]]))
    end

    RS = RS ./ sstd(x; k=1)

    return RS
end

"""
    sstd(x::Array; k::Int=0)

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
function sstd(x::Array; k::Int=0)
    T = length(x)
    μ = smean(x)
    return sqrt(sum((x .- μ) .^ 2) / (T - k))
end

"""
    smean(x::Array)

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
function smean(x::Array)
    T = length(x)
    return sum(x) / T
end


end # module