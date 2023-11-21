"""
# Forecasters

This module contains functions to forecast a long memory time series using the fractional differencing method and the CSA method.

## Author
[J. Eduardo Vera-Valdés](https://everval.github.io/)

"""
module Forecasters

export fi_ar_coefs, fi_forecast, csa_forecast, har_forecast


"""
    fi_forecast(x::Array, h::Int, d::Real)

Computes the forecast of a long memory time series using the fractional differencing method.

# Arguments
- `x::Array`: The time series.
- `h::Int`: The number of periods to forecast.
- `d::Real`: The fractional differencing parameter.

# Output
- `xfor::Array`: The forecast of the time series as a column vector. The first T elements are the original time series.

# Notes
Multiple dispatch is used to compute the forecast of a time series with or without confidence bands.

# Examples    
```julia
julia> fi_forecast(figen(100,0.2), 10, 0.2)
```
"""
function fi_forecast(x::Array, h::Int, d::Real)
    T = length(x)

    maxlags = T+h

    vars = fi_ar_coefs(d,maxlags)

    xfor = zeros(maxlags,1)
    xfor[1:T,1] = x

    for ii = T+1:maxlags
        xfor[ii,1] = sum( reverse(xfor[1:(ii-1),1]).*vars[1:ii-1] )
    end

    return xfor

end


"""
    fi_forecast(x::Array, h::Int, d::Real, σ::Real)

Computes the forecast of a long memory time series using the fractional differencing method.

# Arguments
- `x::Array`: The time series.
- `h::Int`: The number of periods to forecast.
- `d::Real`: The fractional differencing parameter.
- `σ::Real`: The standard deviation of the forecast errors.

# Output
- `xfor::Array`: The forecast of the time series as a matrix where the first column is the forecast, the second column is the lower confidence band, and the third column is the upper confidence band. The first T elements are the original time series.

# Notes
Multiple dispatch is used to compute the forecast of a time series with or without confidence bands.

# Examples    
```julia
julia> fi_forecast(figen(100,0.2), 10, 0.2)
```
"""
function fi_forecast(x::Array, h::Int, d::Real, σ::Real)
    T = length(x)

    maxlags = T+h


    xfor = zeros(maxlags,3)
    xfor[:,1] = fi_forecast(x,h,d)

    xfor[T+1:maxlags,2] = xfor[T+1:maxlags,1] - 2*σ*ones(maxlags-T,1)
    xfor[T+1:maxlags,3] = xfor[T+1:maxlags,1] + 2*σ*ones(maxlags-T,1)

    return xfor

end

"""
    fi_ar_coefs(d::Real, maxlags::Int)

Computes the AR coefficients of the fractional differenced process with parameter `d` at lags 1, ..., `maxlags`.

# Arguments
- `d::Real`: The fractional differencing parameter.
- `maxlags::Int`: The number of lags to compute.

# Output
- `ar_coefs::Array`: The AR coefficients of the fractional differenced process with parameter `d` at lags 1, ..., `maxlags`.

# Notes
The AR coefficients are computed using the recursive formula for speed.
The zero lag coefficient is not computed, it is theoretically 1.

# Examples    
```julia
julia> fi_ar_coefs(0.2, 5)
``` 
"""
function fi_ar_coefs(d::Real, maxlags::Int)
    ar_coefs = zeros(maxlags)

    ar_coefs[1] = d
    for ii = 2:maxlags
        ar_coefs[ii] = ar_coefs[ii-1]*(ii-1-d)/(ii)
    end

    return ar_coefs
end


"""
    csa_forecast(x::Array, h::Int, p::Real, q::Real)

Computes the forecast of a long memory time series using the CSA method.

# Arguments
- `x::Array`: The time series.
- `h::Int`: The number of periods to forecast.
- `p::Real`: The parameter p of the CSA process.
- `q::Real`: The parameter q of the CSA process.

# Output
- `xfor::Array`: The forecast of the time series as a column vector. The first T elements are the original time series.

# Notes
Multiple dispatch is used to compute the forecast of a time series with or without confidence bands.

# Examples    
```julia
julia> csa_forecast(csafigen(100,1.4,1.4), 10, 1.4, 1.4)
```
"""
function csa_forecast(x::Array, h::Int, p::Real, q::Real)
    T = length(x)

    maxlags = T+h

    vars = csa_ma_coefs(p,q,maxlags)

    matvar = my_half_toeplitz(vars[1:T,1])

    errs = zeros(maxlags,1)

    errs[1:T,1] = matvar\x

    xfor = zeros(maxlags,1)
    xfor[1:T,1] = x

    for ii = T+1:maxlags
        xfor[ii,1] = sum( reverse(errs[1:(ii-1),1]).*vars[1:ii-1] )
    end

    return xfor

end

"""
    csa_forecast(x::Array, h::Int, p::Real, q::Real, σ::Real)

Computes the forecast of a long memory time series using the CSA method.

# Arguments
- `x::Array`: The time series.
- `h::Int`: The number of periods to forecast.
- `p::Real`: The parameter p of the CSA process.
- `q::Real`: The parameter q of the CSA process.
- `σ::Real`: The standard deviation of the forecast errors.

# Output
- `xfor::Array`: The forecast of the time series as a matrix where the first column is the forecast, the second column is the lower confidence band, and the third column is the upper confidence band. The first T elements are the original time series.

# Notes
Multiple dispatch is used to compute the forecast of a time series with or without confidence bands.

# Examples    
```julia
julia> csa_forecast(csafigen(100,1.4,1.4), 10, 1.4, 1.4, 1.0)
```
"""
function csa_forecast(x::Array, h::Int, p::Real, q::Real, σ::Real)
    T = length(x)

    maxlags = T+h

    xfor = zeros(maxlags,3)
    xfor[:,1] = csa_forecast(x,h,p,q)

    xfor[T+1:maxlags,2] = xfor[T+1:maxlags,1] - 2*σ*ones(maxlags-T,1)
    xfor[T+1:maxlags,3] = xfor[T+1:maxlags,1] + 2*σ*ones(maxlags-T,1)

    return xfor

end


"""
    csa_ma_coefs(p::Real, q::Real, maxlags::Int)

Computes the MA coefficients of the CSA process with parameters `p` and `q` at lags 0, ..., `maxlags-1`.

# Arguments
- `p::Real`: The parameter p of the CSA process.
- `q::Real`: The parameter q of the CSA process.
- `maxlags::Int`: The number of lags to compute.

# Output
- `ma_coefs::Array`: The MA coefficients of the CSA process with parameters `p` and `q` at lags 0, ..., `maxlags-1`.

# Notes
The MA coefficients are computed using the recursive formula for speed.
The zero lag coefficient is included, it is theoretically 1.

# Examples    
```julia
julia> csa_ma_coefs(1.4, 1.4, 5)
``` 
"""
function csa_ma_coefs(p::Real , q::Real, maxlags::Int)
    ma_coefs = zeros(maxlags)

    ma_coefs[1] = 1
    for t in 2:maxlags
        ma_coefs[t] = ma_coefs[t-1] * ((p + t - 2) / (p + t - 2 + q))^(1 / 2)
    end

    return ma_coefs
end


"""
    my_half_toeplitz(coefs::Array)

Constructs a bottom Toeplitz matrix from the given coefficients.

# Arguments
- `coefs::Array`: An array of coefficients.

# Output
- `Toep::Array`: The bottom Toeplitz matrix constructed from the given coefficients.

# Examples    
```julia
julia> my_half_toeplitz([1, 2, 3])
```
"""
function my_half_toeplitz(coefs::Array)
    N = length(coefs)
    Toep = zeros(N, N)

    for ii = 1:N
        for jj = 1:N
            if ii >= jj
                Toep[ii, jj] = coefs[ii-jj+1]
            end
        end
    end

    return Toep
end


"""
    har_forecast(x::Array, h::Int, m::Array=[1,5,22])

Computes the forecast of a time series by fitting and recursevely forecasting the HAR model.

# Arguments
- `x::Array`: The time series.
- `h::Int`: The number of periods to forecast.

# Output
- `xfor::Array`: The forecast of the time series as a column vector. The first T-max(m) elements are the original time series.

# Optional Arguments
- `m::Array`: The lags to include in the HAR model. The default is [1,5,22].

# Examples    
```julia
julia> har_forecast(figen(100,0.2), 10)
```
"""
function har_forecast(x::Array, h::Int, m::Array=[1,5,22])
    T = length(x)
    
    mm = maximum(m)
    n = length(m)

    ## Estimation because the matrix are needed for forecasting
    X = zeros(T-mm+h, n+1)
    X[:,1] = ones(T-mm+h, 1)

    for ii = 1:n
        cm = m[ii]
        aux = zeros(T-mm, 1)
        for jj = 1:cm
            aux = aux + x[mm-jj+1:T-jj,1]
        end
        X[1:T-mm,ii+1] = aux/cm
    end

    Y = zeros(T-mm+h,1)

    Y[1:T-mm,1] = x[mm+1:T,1]

    betas = X[1:T-mm,:]\Y[1:T-mm,1]
    err = Y[1:T-mm,1]-X[1:T-mm,:]*betas

    sigma = (err'*err)/(T-mm-n-1)

    ## Forecasting
    for ii = 1:h
        Y[T-mm+ii,1] = X[T-mm+ii-1,:]'*betas 
        for jj = 1:n
            cm = m[jj]
            X[T-mm+ii,jj+1] = sum(Y[T-mm+ii-cm:T-mm+ii-1,1])/cm
        end
    end

    return Y

end



"""
    har_forecast(x::Array, h::Int, σ::Real, m::Array=[1,5,22])

Computes the forecast of a time series by fitting and recursevely forecasting the HAR model.

# Arguments
- `x::Array`: The time series.
- `h::Int`: The number of periods to forecast.
- `σ::Real`: The standard deviation of the forecast errors.

# Output
- `xfor::Array`: The forecast of the time series as a matrix where the first column is the forecast, the second column is the lower confidence band, and the third column is the upper confidence band. The first T-max(m) elements are the original time series.

# Optional Arguments
- `m::Array`: The lags to include in the HAR model. The default is [1,5,22].

# Examples    
```julia
julia> har_forecast(figen(100,0.2), 10, 1)
```
"""
function har_forecast(x::Array, h::Int, σ::Real, m::Array=[1,5,22])
    T = length(x)

    mm = maximum(m)

    Y = zeros(T-mm+h,3)
    Y[:,1] = har_forecast(x,h,m)
    
    Y[T-mm+1:T-mm+h,2] = Y[T-mm+1:T-mm+h,1] - 2*σ*ones(h,1)
    Y[T-mm+1:T-mm+h,3] = Y[T-mm+1:T-mm+h,1] + 2*σ*ones(h,1)

    return Y

end

end # module Forecasters