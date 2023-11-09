"""
# ClassicEstimators

This module contains functions to estimate them Hurst coefficient, closely related to the long memory parameter, of a time series using the rescaled range (R/S) statistic.

## Author
[J. Eduardo Vera-Valdés](https://everval.github.io/)

"""
module ClassicEstimators

export rescaled_range_est, rescaled_range, my_std


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

    X = [ones(T-1,1) log.( collect(2:T) )]

    beta = X\Y

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

    μ = sum(x)/T

    Y = x.-μ
    Z = cumsum(Y,dims=1)

    RS = zeros(T-1,1)

    for ii = 2:T
        RS[ii-1,1] = ( maximum(Z[1:ii]) - minimum(Z[1:ii]) )/ my_std(x[1:ii])
    end

    return RS
end


"""
    my_std(x::Array)

Computes the sample standard deviation of a time series.

# Arguments
- `x::Array`: The time series.

# Output
- `std::Real`: The sample standard deviation.

# Notes
This function divides by T instead of T-1.

# Examples    
```julia
julia> my_std(randn(100))
```
"""
function my_std(x::Array)
    T = length(x)
    μ = sum(x)/T
    return sqrt( sum( (x.-μ).^2 ) / T )
end

end # module