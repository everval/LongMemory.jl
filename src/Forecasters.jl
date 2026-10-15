"""
# Forecasters

This module contains functions to forecast a long memory time series using the fractional differencing method and the CSA method.

## Author

[J. Eduardo Vera-Valdés](https://everval.github.io/)

"""
module Forecasters

import ..plots_required

using ..ClassicEstimators: smean

using FFTW

export fi_ar_coefs, fi_ma_coefs, fi_forecast, fi_forecast_plot, csa_forecast, csa_forecast_plot, csa_ma_coefs, har_forecast, har_forecast_plot, gregory_coeffs, hitransform, arfima_forecast, arfima_forecast_plot, arfima_ar_coefs, arfima_ma_coefs

# Plotting functions. Their methods are defined in the package extension LongMemoryPlotsExt, which is loaded with `using Plots`.
function fi_forecast_plot end
function csa_forecast_plot end
function har_forecast_plot end
function arfima_forecast_plot end

fi_forecast_plot(args...; kwargs...) = plots_required(fi_forecast_plot, args)
csa_forecast_plot(args...; kwargs...) = plots_required(csa_forecast_plot, args)
har_forecast_plot(args...; kwargs...) = plots_required(har_forecast_plot, args)
arfima_forecast_plot(args...; kwargs...) = plots_required(arfima_forecast_plot, args)


"""
    fi_forecast(x::AbstractArray, h::Int, d::Real)

Computes the forecast of a long memory time series using the fractional differencing method.

# Arguments
- `x::Array`: The time series.
- `h::Int`: The number of periods to forecast.
- `d::Real`: The fractional differencing parameter.

# Output
- `xfor::Array`: The forecast of the time series as a vector. The first T elements are the original time series.

# Notes
Multiple dispatch is used to compute the forecast of a time series with or without confidence bands.
The series is demeaned with the sample mean before forecasting and the mean is added back to the forecasts.

# Examples
```julia
julia> fi_forecast(fi_gen(100,0.2), 10, 0.2)
```
"""
function fi_forecast(x::AbstractArray, h::Int, d::Real)
    T = length(x)
    μ = smean(x)

    maxlags = T + h

    vars = fi_ar_coefs(d, maxlags)

    xfor = zeros(maxlags)
    xfor[1:T] = x .- μ

    for ii = T+1:maxlags
        xfor[ii] = sum(reverse(xfor[1:(ii-1)]) .* vars[1:ii-1])
    end

    xfor[T+1:maxlags] .+= μ
    xfor[1:T] = x

    return xfor

end


"""
    fi_forecast(x::AbstractArray, h::Int, d::Real, σ::Real)

Computes the forecast of a long memory time series using the fractional differencing method.

# Arguments
- `x::Array`: The time series.
- `h::Int`: The number of periods to forecast.
- `d::Real`: The fractional differencing parameter.
- `σ::Real`: The standard deviation of the innovations.

# Output
- `xfor::Array`: The forecast of the time series as a matrix where the first column is the forecast, the second column is the lower confidence band, and the third column is the upper confidence band. The first T elements are the original time series.

# Notes
Multiple dispatch is used to compute the forecast of a time series with or without confidence bands.
The bands are the forecast ± 2 standard errors of the `k`-step forecast error, `σ (ψ₀² + ... + ψ_{k-1}²)^{1/2}`, where `ψ` are the MA coefficients from `fi_ma_coefs`. Hence, they widen with the horizon.

# Examples
```julia
julia> fi_forecast(fi_gen(100,0.2), 10, 0.2, 1.0)
```
"""
function fi_forecast(x::AbstractArray, h::Int, d::Real, σ::Real)
    T = length(x)

    maxlags = T + h


    xfor = zeros(maxlags, 3)
    xfor[:, 1] = fi_forecast(x, h, d)

    se = σ .* sqrt.(cumsum(fi_ma_coefs(d, h) .^ 2)) # standard errors of the forecast errors

    xfor[T+1:maxlags, 2] = xfor[T+1:maxlags, 1] .- 2 .* se
    xfor[T+1:maxlags, 3] = xfor[T+1:maxlags, 1] .+ 2 .* se

    return xfor

end

"""
    fi_forecast_plot(x::AbstractArray, h::Int, d::Real)

Plots the forecast of a long memory time series using the fractional differencing method.

Requires Plots.jl to be loaded, e.g. with `using LongMemory, Plots`.

# Arguments
- `x::Array`: The time series.
- `h::Int`: The number of periods to forecast.
- `d::Real`: The fractional differencing parameter.

# Output
- `p1::Plot`: The plot of the original time series and the forecast.

# Notes
Multiple dispatch is used to plot the forecast of a time series with or without confidence bands.
Plotting is done by calling the `fi_forecast` function.

# Examples    
```julia
julia> fi_forecast_plot(fi_gen(100,0.2), 10, 0.2)
```
"""
fi_forecast_plot(x::AbstractArray, h::Int, d::Real)

"""
    fi_forecast_plot(x::AbstractArray, h::Int, d::Real, σ::Real)

Plots the forecast of a long memory time series using the fractional differencing method.

Requires Plots.jl to be loaded, e.g. with `using LongMemory, Plots`.

# Arguments
- `x::Array`: The time series.
- `h::Int`: The number of periods to forecast.
- `d::Real`: The fractional differencing parameter.
- `σ::Real`: The standard deviation of the innovations.

# Output
- `p1::Plot`: The plot of the original time series and the forecast with confidence bands.

# Notes
Multiple dispatch is used to plot the forecast of a time series with or without confidence bands.
Plotting is done by calling the `fi_forecast` function.

# Examples    
```julia
julia> fi_forecast_plot(fi_gen(100,0.2), 10, 0.2, 1.0)
```
"""
fi_forecast_plot(x::AbstractArray, h::Int, d::Real, σ::Real)

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
        ar_coefs[ii] = ar_coefs[ii-1] * (ii - 1 - d) / (ii)
    end

    return ar_coefs
end


"""
    fi_ma_coefs(d::Real, maxlags::Int)

Computes the MA coefficients of the fractional differenced process with parameter `d` at lags 0, ..., `maxlags-1`.

# Arguments
- `d::Real`: The fractional differencing parameter.
- `maxlags::Int`: The number of lags to compute.

# Output
- `ma_coefs::Array`: The MA coefficients of the fractional differenced process with parameter `d` at lags 0, ..., `maxlags-1`, i.e. the coefficients of `(1-L)^{-d}`.

# Notes
The MA coefficients are computed using the recursive formula `ψ_j = ψ_{j-1} (j - 1 + d) / j` for speed.
The zero lag coefficient is included, it is 1. This matches the convention of `csa_ma_coefs`.

# Examples
```julia
julia> fi_ma_coefs(0.2, 5)
```
"""
function fi_ma_coefs(d::Real, maxlags::Int)
    ma_coefs = zeros(maxlags)

    ma_coefs[1] = 1
    for ii = 2:maxlags
        ma_coefs[ii] = ma_coefs[ii-1] * (ii - 2 + d) / (ii - 1)
    end

    return ma_coefs
end


"""
    csa_forecast(x::AbstractArray, h::Int, p::Real, q::Real)

Computes the forecast of a long memory time series using the CSA method.

# Arguments
- `x::Array`: The time series.
- `h::Int`: The number of periods to forecast.
- `p::Real`: The parameter p of the CSA process.
- `q::Real`: The parameter q of the CSA process.

# Output
- `xfor::Array`: The forecast of the time series as a vector. The first T elements are the original time series.

# Notes
Multiple dispatch is used to compute the forecast of a time series with or without confidence bands.
The series is demeaned with the sample mean before forecasting and the mean is added back to the forecasts.

# Examples
```julia
julia> csa_forecast(csa_gen(100,1.4,1.4), 10, 1.4, 1.4)
```
"""
function csa_forecast(x::AbstractArray, h::Int, p::Real, q::Real)
    T = length(x)
    μ = smean(x)

    maxlags = T + h

    vars = csa_ma_coefs(p, q, maxlags)

    matvar = my_half_toeplitz(vars[1:T, 1])

    errs = zeros(maxlags)

    errs[1:T] = matvar \ (vec(x) .- μ)

    xfor = zeros(maxlags)
    xfor[1:T] = x

    for ii = T+1:maxlags
        xfor[ii] = μ + sum(reverse(errs[1:ii]) .* vars[1:ii]) # future shocks errs[T+1:ii] are zero
    end

    return xfor

end

"""
    csa_forecast(x::AbstractArray, h::Int, p::Real, q::Real, σ::Real)

Computes the forecast of a long memory time series using the CSA method.

# Arguments
- `x::Array`: The time series.
- `h::Int`: The number of periods to forecast.
- `p::Real`: The parameter p of the CSA process.
- `q::Real`: The parameter q of the CSA process.
- `σ::Real`: The standard deviation of the innovations.

# Output
- `xfor::Array`: The forecast of the time series as a matrix where the first column is the forecast, the second column is the lower confidence band, and the third column is the upper confidence band. The first T elements are the original time series.

# Notes
Multiple dispatch is used to compute the forecast of a time series with or without confidence bands.
The bands are the forecast ± 2 standard errors of the `k`-step forecast error, `σ (ψ₀² + ... + ψ_{k-1}²)^{1/2}`, where `ψ` are the MA coefficients from `csa_ma_coefs`. Hence, they widen with the horizon.

# Examples
```julia
julia> csa_forecast(csa_gen(100,1.4,1.4), 10, 1.4, 1.4, 1.0)
```
"""
function csa_forecast(x::AbstractArray, h::Int, p::Real, q::Real, σ::Real)
    T = length(x)

    maxlags = T + h

    xfor = zeros(maxlags, 3)
    xfor[:, 1] = csa_forecast(x, h, p, q)

    se = σ .* sqrt.(cumsum(csa_ma_coefs(p, q, h) .^ 2)) # standard errors of the forecast errors

    xfor[T+1:maxlags, 2] = xfor[T+1:maxlags, 1] .- 2 .* se
    xfor[T+1:maxlags, 3] = xfor[T+1:maxlags, 1] .+ 2 .* se

    return xfor

end

"""
    csa_forecast_plot(x::AbstractArray, h::Int, p::Real, q::Real)

Plots the forecast of a long memory time series using the CSA method.

Requires Plots.jl to be loaded, e.g. with `using LongMemory, Plots`.

# Arguments
- `x::Array`: The time series.
- `h::Int`: The number of periods to forecast.
- `p::Real`: The parameter p of the CSA process.
- `q::Real`: The parameter q of the CSA process.

# Output
- `p1::Plot`: The plot of the original time series and the forecast.

# Notes
Multiple dispatch is used to plot the forecast of a time series with or without confidence bands.
Plotting is done by calling the `csa_forecast` function.

# Examples    
```julia
julia> csa_forecast_plot(csa_gen(100,1.4,1.4), 10, 1.4, 1.4)
```
"""
csa_forecast_plot(x::AbstractArray, h::Int, p::Real, q::Real)

"""
    csa_forecast_plot(x::AbstractArray, h::Int, p::Real, q::Real, σ::Real)

Plots the forecast of a long memory time series using the CSA method.

Requires Plots.jl to be loaded, e.g. with `using LongMemory, Plots`.

# Arguments
- `x::Array`: The time series.
- `h::Int`: The number of periods to forecast.
- `p::Real`: The parameter p of the CSA process.
- `q::Real`: The parameter q of the CSA process.
- `σ::Real`: The standard deviation of the innovations.

# Output
- `p1::Plot`: The plot of the original time series and the forecast with confidence bands.

# Notes
Multiple dispatch is used to plot the forecast of a time series with or without confidence bands.
Plotting is done by calling the `csa_forecast` function.

# Examples    
```julia
julia> csa_forecast_plot(csa_gen(100,1.4,1.4), 10, 1.4, 1.4, 1.0)
```
"""
csa_forecast_plot(x::AbstractArray, h::Int, p::Real, q::Real, σ::Real)


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
function csa_ma_coefs(p::Real, q::Real, maxlags::Int)
    ma_coefs = zeros(maxlags)

    ma_coefs[1] = 1
    for t in 2:maxlags
        ma_coefs[t] = ma_coefs[t-1] * ((p + t - 2) / (p + t - 2 + q))^(1 / 2)
    end

    return ma_coefs
end


"""
    my_half_toeplitz(coefs::AbstractArray)

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
function my_half_toeplitz(coefs::AbstractArray)
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
    har_forecast(x::AbstractArray, h::Int; flag::Bool = false, m::AbstractArray=[1,5,22])

Computes the forecast of a time series by fitting and recursively forecasting the HAR model.

# Arguments
- `x::Array`: The time series.
- `h::Int`: The number of periods to forecast.

# Output
- `xfor::Array`: The forecast of the time series. The first T-max(m) elements are the original time series.

# Optional Arguments
- `flag::Bool`: If true, the output is a matrix where the first column is the forecast, the second column is the lower confidence band, and the third column is the upper confidence band. The default is false.	
- `m::Array`: The lags to include in the HAR model. The default is [1,5,22].

# Notes
The bands are the forecast ± 2 standard errors of the `k`-step forecast error, computed from the MA representation of the fitted HAR model (an AR(max(m)) model) and the residual standard deviation. Hence, they widen with the horizon. Parameter uncertainty is not taken into account.

# Examples
```julia
julia> har_forecast(fi_gen(100,0.2), 10)
```
"""
function har_forecast(x::AbstractArray, h::Int; flag::Bool=false, m::AbstractArray=[1, 5, 22])
    T = length(x)

    mm = maximum(m)
    n = length(m)

    ## Estimation
    X = zeros(T - mm, n + 1)
    X[:, 1] = ones(T - mm, 1)

    for ii = 1:n
        cm = m[ii]
        aux = zeros(T - mm, 1)
        for jj = 1:cm
            aux = aux + x[mm-jj+1:T-jj, 1]
        end
        X[1:T-mm, ii+1 ] = aux / cm
    end

    Y = x[mm+1:T, 1]

    betas = X \ Y

    err = Y - X * betas

    σ = sqrt((err' * err) / (T - mm - n - 1))

    ## Forecasting: x̂_t = β₀ + Σ_j β_j mean(x_{t-m_j}, ..., x_{t-1}), with forecasts replacing unknown values
    xfor = zeros(T + h, 3)
    xfor[1:T, 1] = x[1:T, 1]
    for ii = T+1:T+h
        xfor[ii, 1] = betas[1]
        for jj = 1:n
            cm = m[jj]
            xfor[ii, 1] += betas[jj+1] * sum(@view xfor[ii-cm:ii-1, 1]) / cm
        end
    end
    ## Standard errors: the HAR model is an AR(mm) model with coefficients a_k = Σ_{j: m_j ≥ k} β_j / m_j
    a = zeros(mm)
    for jj = 1:n
        cm = m[jj]
        a[1:cm] .+= betas[jj+1] / cm
    end

    ψ = zeros(h) # MA coefficients ψ_0, ..., ψ_{h-1}
    ψ[1] = 1
    for ii = 2:h
        for kk = 1:min(ii - 1, mm)
            ψ[ii] += a[kk] * ψ[ii-kk]
        end
    end

    se = σ .* sqrt.(cumsum(ψ .^ 2))

    xfor[T+1:T+h, 2] = xfor[T+1:T+h, 1] .- 2 .* se
    xfor[T+1:T+h, 3] = xfor[T+1:T+h, 1] .+ 2 .* se

    if flag == true
        return xfor
    elseif flag == false
        return xfor[:, 1]
    end

end


"""
    har_forecast_plot(x::AbstractArray, h::Int; flag::Bool = true, m::AbstractArray=[1,5,22])

Plots the forecast of a time series by fitting and recursively forecasting the HAR model.

Requires Plots.jl to be loaded, e.g. with `using LongMemory, Plots`.

# Arguments
- `x::Array`: The time series.
- `h::Int`: The number of periods to forecast.

# Output
- `p1::Plot`: The plot of the original time series and the forecast.

# Optional Arguments
- `flag::Bool`: If true, the output is a matrix where the first column is the forecast, the second column is the lower confidence band, and the third column is the upper confidence band. The default is true.
- `m::Array`: The lags to include in the HAR model. The default is [1,5,22].

# Notes
Multiple dispatch is used to plot the forecast of a time series with or without confidence bands.
Plotting is done by calling the `har_forecast` function.

# Examples    
```julia
julia> har_forecast_plot(fi_gen(100,0.2), 10)
```
"""
har_forecast_plot(x::AbstractArray, h::Int)

"""
    arfima_forecast(x::AbstractArray, h::Int, d::Real, AR::AbstractArray, MA::AbstractArray)

Computes the forecast of a time series using the ARFIMA(p,d,q) model.

# Arguments
- `x::Array`: The time series.
- `h::Int`: The number of periods to forecast.
- `d::Real`: The fractional differencing parameter.
- `AR::Array`: The AR coefficients `φ`, with the AR polynomial `1 - φ₁L - ... - φₚLᵖ`.
- `MA::Array`: The MA coefficients `θ`, with the MA polynomial `1 + θ₁L + ... + θ_qL^q`.

# Output
- `xfor::Array`: The forecast of the time series as a vector. The first T elements are the original time series.

# Notes
Multiple dispatch is used to compute the forecast of a time series with or without confidence bands.
The forecasts use the AR(∞) representation of the model, `φ(L)(1-L)^d / θ(L)`, truncated at the start of the sample, see `arfima_ar_coefs`. With `AR = MA = []` they coincide with `fi_forecast`.
The series is demeaned with the sample mean before forecasting and the mean is added back to the forecasts.
The parametrization is the one used by `arfima_gen` and `arfima_whittle_est`.

# Examples
```julia
julia> arfima_forecast(arfima_gen(100, 0, [0.5], 0.2, []), 10, 0.2, [0.5], [])
```
"""
function arfima_forecast(x::AbstractArray, h::Int, d::Real, AR::AbstractArray, MA::AbstractArray)
    T = length(x)
    μ = smean(x)

    maxlags = T + h

    vars = arfima_ar_coefs(d, AR, MA, maxlags)

    xfor = zeros(maxlags)
    xfor[1:T] = x .- μ

    for ii = T+1:maxlags
        xfor[ii] = sum(reverse(xfor[1:(ii-1)]) .* vars[1:ii-1])
    end

    xfor[T+1:maxlags] .+= μ
    xfor[1:T] = x

    return xfor

end


"""
    arfima_forecast(x::AbstractArray, h::Int, d::Real, AR::AbstractArray, MA::AbstractArray, σ::Real)

Computes the forecast of a time series using the ARFIMA(p,d,q) model, with confidence bands.

# Arguments
- `x::Array`: The time series.
- `h::Int`: The number of periods to forecast.
- `d::Real`: The fractional differencing parameter.
- `AR::Array`: The AR coefficients.
- `MA::Array`: The MA coefficients.
- `σ::Real`: The standard deviation of the innovations.

# Output
- `xfor::Array`: The forecast of the time series as a matrix where the first column is the forecast, the second column is the lower confidence band, and the third column is the upper confidence band. The first T elements are the original time series.

# Notes
Multiple dispatch is used to compute the forecast of a time series with or without confidence bands.
The bands are the forecast ± 2 standard errors of the `k`-step forecast error, `σ (ψ₀² + ... + ψ_{k-1}²)^{1/2}`, where `ψ` are the MA coefficients from `arfima_ma_coefs`.

# Examples
```julia
julia> arfima_forecast(arfima_gen(100, 0, [0.5], 0.2, []), 10, 0.2, [0.5], [], 1.0)
```
"""
function arfima_forecast(x::AbstractArray, h::Int, d::Real, AR::AbstractArray, MA::AbstractArray, σ::Real)
    T = length(x)

    maxlags = T + h

    xfor = zeros(maxlags, 3)
    xfor[:, 1] = arfima_forecast(x, h, d, AR, MA)

    se = σ .* sqrt.(cumsum(arfima_ma_coefs(d, AR, MA, h) .^ 2)) # standard errors of the forecast errors

    xfor[T+1:maxlags, 2] = xfor[T+1:maxlags, 1] .- 2 .* se
    xfor[T+1:maxlags, 3] = xfor[T+1:maxlags, 1] .+ 2 .* se

    return xfor

end


"""
    arfima_forecast_plot(x::AbstractArray, h::Int, d::Real, AR::AbstractArray, MA::AbstractArray)

Plots the forecast of a time series using the ARFIMA(p,d,q) model.

Requires Plots.jl to be loaded, e.g. with `using LongMemory, Plots`.

# Arguments
- `x::Array`: The time series.
- `h::Int`: The number of periods to forecast.
- `d::Real`: The fractional differencing parameter.
- `AR::Array`: The AR coefficients.
- `MA::Array`: The MA coefficients.

# Output
- `p1::Plot`: The plot of the original time series and the forecast.

# Notes
Multiple dispatch is used to plot the forecast of a time series with or without confidence bands.
Plotting is done by calling the `arfima_forecast` function.

# Examples
```julia
julia> arfima_forecast_plot(arfima_gen(100, 0, [0.5], 0.2, []), 10, 0.2, [0.5], [])
```
"""
arfima_forecast_plot(x::AbstractArray, h::Int, d::Real, AR::AbstractArray, MA::AbstractArray)


"""
    arfima_forecast_plot(x::AbstractArray, h::Int, d::Real, AR::AbstractArray, MA::AbstractArray, σ::Real)

Plots the forecast of a time series using the ARFIMA(p,d,q) model, with confidence bands.

Requires Plots.jl to be loaded, e.g. with `using LongMemory, Plots`.

# Arguments
- `x::Array`: The time series.
- `h::Int`: The number of periods to forecast.
- `d::Real`: The fractional differencing parameter.
- `AR::Array`: The AR coefficients.
- `MA::Array`: The MA coefficients.
- `σ::Real`: The standard deviation of the innovations.

# Output
- `p1::Plot`: The plot of the original time series and the forecast with confidence bands.

# Notes
Multiple dispatch is used to plot the forecast of a time series with or without confidence bands.
Plotting is done by calling the `arfima_forecast` function.

# Examples
```julia
julia> arfima_forecast_plot(arfima_gen(100, 0, [0.5], 0.2, []), 10, 0.2, [0.5], [], 1.0)
```
"""
arfima_forecast_plot(x::AbstractArray, h::Int, d::Real, AR::AbstractArray, MA::AbstractArray, σ::Real)


"""
    arfima_ar_coefs(d::Real, AR::AbstractArray, MA::AbstractArray, maxlags::Int)

Computes the AR(∞) coefficients of the ARFIMA(p,d,q) model at lags 1, ..., `maxlags`.

# Arguments
- `d::Real`: The fractional differencing parameter.
- `AR::Array`: The AR coefficients.
- `MA::Array`: The MA coefficients.
- `maxlags::Int`: The number of lags to compute.

# Output
- `ar_coefs::Array`: The coefficients `π_j` in `x_t = Σ_j π_j x_{t-j} + ε_t`, at lags 1, ..., `maxlags`.

# Notes
The coefficients are those of `1 - φ(L)(1-L)^d / θ(L)`. They are computed by multiplying the polynomials and dividing by `θ(L)` recursively.
With `AR = MA = []` they coincide with `fi_ar_coefs`.

# Examples
```julia
julia> arfima_ar_coefs(0.2, [0.5], [0.3], 5)
```
"""
function arfima_ar_coefs(d::Real, AR::AbstractArray, MA::AbstractArray, maxlags::Int)
    p = length(AR)
    q = length(MA)

    # Coefficients of (1-L)^d at lags 0, ..., maxlags
    c = zeros(maxlags + 1)
    c[1] = 1
    for ii = 2:maxlags+1
        c[ii] = c[ii-1] * (ii - 2 - d) / (ii - 1)
    end

    # Multiplying by φ(L)
    a = copy(c)
    for ii = 2:maxlags+1
        for kk = 1:min(ii - 1, p)
            a[ii] -= AR[kk] * c[ii-kk]
        end
    end

    # Dividing by θ(L)
    b = copy(a)
    for ii = 2:maxlags+1
        for kk = 1:min(ii - 1, q)
            b[ii] -= MA[kk] * b[ii-kk]
        end
    end

    ar_coefs = -b[2:maxlags+1]

    return ar_coefs
end


"""
    arfima_ma_coefs(d::Real, AR::AbstractArray, MA::AbstractArray, maxlags::Int)

Computes the MA(∞) coefficients of the ARFIMA(p,d,q) model at lags 0, ..., `maxlags-1`.

# Arguments
- `d::Real`: The fractional differencing parameter.
- `AR::Array`: The AR coefficients.
- `MA::Array`: The MA coefficients.
- `maxlags::Int`: The number of lags to compute.

# Output
- `ma_coefs::Array`: The coefficients of `θ(L)(1-L)^{-d} / φ(L)` at lags 0, ..., `maxlags-1`.

# Notes
The zero lag coefficient is included, it is 1. With `AR = MA = []` they coincide with `fi_ma_coefs`.

# Examples
```julia
julia> arfima_ma_coefs(0.2, [0.5], [0.3], 5)
```
"""
function arfima_ma_coefs(d::Real, AR::AbstractArray, MA::AbstractArray, maxlags::Int)
    p = length(AR)
    q = length(MA)

    ψd = fi_ma_coefs(d, maxlags)

    # Multiplying by θ(L)
    e = copy(ψd)
    for ii = 2:maxlags
        for kk = 1:min(ii - 1, q)
            e[ii] += MA[kk] * ψd[ii-kk]
        end
    end

    # Dividing by φ(L)
    ma_coefs = copy(e)
    for ii = 2:maxlags
        for kk = 1:min(ii - 1, p)
            ma_coefs[ii] += AR[kk] * ma_coefs[ii-kk]
        end
    end

    return ma_coefs
end


"""
    gregory_coeffs(N::Int)

Computes the absolute values of the Gregory coefficients, `|G_n|` for `n = 1, ..., N`, used in the harmonic inverse transformation.

# Arguments
- `N::Int`: The number of coefficients to compute.

# Output
- `coefs::Array`: The Gregory coefficients.

# Notes
The first coefficient is 1/2 and the rest are computed using the recursive formula. The zero lag coefficient is not computed, it is theoretically 1.
The inverse of the harmonically weighted filter is `1 - Σ_{n≥1} |G_n| Lⁿ`.

# Examples
```julia
julia> gregory_coeffs(5)
```
"""
function gregory_coeffs(N::Int)
    coefs = zeros(N)

    coefs[1] = 1/2
    for jj = 1:N
        gs = 0
        for ii = 1:(jj-1)
            gs = gs + coefs[ii] / (jj - ii + 1)
        end
        coefs[jj] = 1/(jj + 1) - gs
    end

    return coefs
end

"""
    hitransform(x::AbstractArray)

Computes the harmonic inverse transformation of a time series. For more details see Hassler and Hosseinkouchack (2020).

# Arguments
- `x::Array`: The time series.

# Output
- `result::Array`: The time series after the harmonic inverse transformation (new array, input not mutated).

# Notes
The harmonic inverse transformation is computed using the FFT.

# Examples
```julia
julia> hitransform(fi_gen(100,0.2))
```
"""
function hitransform(x::AbstractArray)

    T = length(x)-1

    if T < 1
        return float.(vec(x))
    end

    np2 = nextpow(2, 2 * T - 1)

    coefs = gregory_coeffs(T)

    padcoefs = [coefs; zeros(np2 - T, 1)]
    padx = [x[1:T]; zeros(np2 - T, 1)]

    gxs = irfft(rfft(padx) .* rfft(padcoefs), np2)
    gxs = gxs[1:T]

    result = float.(vec(x)) # new array; also avoids InexactError for integer input
    result[2:T+1] = result[2:T+1] .- gxs[1:T]

    return result

end

end # module Forecasters