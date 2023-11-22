"""
# ClassicEstimators

This module contains functions to estimate them Hurst coefficient, closely related to the long memory parameter, of a time series using the rescaled range (R/S) statistic.

## Author
[J. Eduardo Vera-Valdés](https://everval.github.io/)

"""
module ClassicEstimators

export rescaled_range_est, rescaled_range, sstd, smean, variance_plot


"""
    variance_plot(x::Array; flag::Bool=true)

Computes the variance plot of a time series and estimates the Hurst coefficient.

# Arguments
- `x::Array`: The time series.

# Output
- `beta::Real`: The estimated Hurst coefficient.

# Optional arguments
- `flag::Bool`: If true, the variance plot is displayed.

# Notes
This function uses the linear regression method on the log of the variance plot to estimate the Hurst coefficient.

# Examples    
```julia
julia> variance_plot(randn(100))
```
"""
function variance_plot(x::Array; flag::Bool=true)
    T = length(x)
    k = floor(Int, T / 2)

    Y = zeros(k-1, 1)

    for ii = 2:k
        Y[ii-1, 1] = sstdk(x, ii)
    end

    X = log.(collect(2:k))
    beta = X \ log.(Y)

    if flag == true
        p1 = plot(X, log.(Y), seriestype=:scatter, label="", title="Variance Plot", xlabel="log(k)", ylabel="log(s2(k))")
        plot!(X, beta .* X, line=:dash, label=string("Slope = ", beta) )

        display(p1)
    end

    return beta

end


"""
    sstdk(x::Array, k::Int)

Computes the sample standard deviation of a time series using the k-th order sample mean.

# Arguments
- `x::Array`: The time series.
- `k::Int`: The order of the sample mean.

# Output
- `std::Real`: The sample standard deviation.

# Examples    
```julia
julia> sstdk(randn(100), 2)
```
"""
function sstdk(x::Array, k::Int=1)
    T = length(x)
    μ = smean(x)

    if k >= T
        error("k must be less than T")
    end

    nt = collect(1:k:T)

    m = length(nt)

    Y = zeros(m - 1, 1)
    for ii = 1:m-1
        Y[ii, 1] = smean(x[nt[ii]:nt[ii+1]])
    end

    s2k = sum((Y .- μ) .^ 2) / (m - 1)

    return s2k

end


"""
    rescaled_range_est(x::Array)

Estimates the Hurst coefficient of a time series using the rescaled range (R/S) statistic.

# Arguments
- `x::Array`: The time series.

# Output
- `H::Real`: The estimated Hurst coefficient.

# Notes
This function uses the linear regression method on the log of the rescaled range to estimate the Hurst coefficient.
The Hurst coefficient is related to the long memory parameter d by the formula H = d + 1/2.

# Examples    
```julia
julia> rescaled_range_est(randn(100))
```
"""
function rescaled_range_est(x::Array)
    T = length(x)

    RS = rescaled_range(x)

    Y = log.(RS)

    X = [ones(T - 1, 1) log.(collect(2:T))]

    beta = X \ Y

    return beta[2]
end


"""
    rescaled_range(x::Array)

Computes the rescaled range (R/S) statistic of a time series.

# Arguments
- `x::Array`: The time series.

# Output
- `RS::Array`: The rescaled range statistic.

# Examples    
```julia
julia> rescaled_range(randn(100))
```
"""
function rescaled_range(x::Array)
    T = length(x)

    μ = smean(x)

    Y = x .- μ
    Z = cumsum(Y, dims=1)

    RS = zeros(T - 1, 1)

    for ii = 2:T
        RS[ii-1, 1] = (maximum(Z[1:ii]) - minimum(Z[1:ii])) / sstd(x[1:ii])
    end

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