"""
# ParametricEstimators

This module contains functions for estimating the parameters of the fractional differenced process à la ARFIMA and the CSA process.

## Author

[J. Eduardo Vera-Valdés](https://everval.github.io/)

"""
module ParametricEstimators

using ..LogPeriodEstimators: gph_est, periodogram
using ..ClassicEstimators: smean
using ..GeneratingFunctions: davies_harte

using Optim, LinearAlgebra, SpecialFunctions

export fi_mle_est, csa_mle_est, har_est, fi_var_vals, csa_var_vals, fi_cor_vals, csa_cor_vals, fi_var_matrix, csa_var_matrix, fi_llk, csa_llk, durbin_levinson, fi_exact_gen, csa_exact_gen, arfima_whittle_est, arfima_var_vals, arfima_cor_vals, arfima_var_matrix, arfima_llk, arfima_mle_est, arfima_exact_gen

"""
    my_toeplitz(coefs::AbstractArray)

Constructs a Toeplitz matrix from the given coefficients.

# Arguments
- `coefs::Array`: An array of coefficients.

# Output
- `Toep::Array`: The Toeplitz matrix constructed from the given coefficients.

# Notes
This is an internal helper function; it is not exported.

# Examples
```julia
julia> LongMemory.ParametricEstimators.my_toeplitz([1, 2, 3])
```
"""
function my_toeplitz(coefs::AbstractArray)
    N = length(coefs)
    Toep = zeros(N, N)

    @inbounds for ii in 1:N, jj in 1:ii
        Toep[ii, jj] = coefs[ii-jj+1]
        Toep[jj, ii] = coefs[ii-jj+1]
    end

    return Toep
end


"""
    durbin_levinson(acv::AbstractArray, x::AbstractArray)

Computes the one-step-ahead best linear predictions of a stationary time series `x` and their mean squared errors, given the autocovariances `acv`, using the Durbin-Levinson algorithm.

# Arguments
- `acv::AbstractArray`: The autocovariances at lags 0, 1, ..., `T-1` (at least `length(x)` values).
- `x::AbstractArray`: The data (zero mean).

# Output
- `xhat::Array`: The predictions of `x[t]` given `x[1], ..., x[t-1]`, for t = 1, ..., T (`xhat[1] = 0`).
- `v::Array`: The mean squared errors of the predictions.

# Notes
The algorithm requires `O(T^2)` operations and `O(T)` memory, instead of the `O(T^3)` operations and `O(T^2)` memory needed to factorize the covariance matrix.
The Gaussian log-likelihood follows from the prediction error decomposition: `logdet(V) = sum(log.(v))` and `x'V⁻¹x = sum((x - xhat).^2 ./ v)`, where `V` is the Toeplitz covariance matrix.
See Brockwell and Davis (1991), Time Series: Theory and Methods, Section 5.2.

# Examples
```julia
julia> durbin_levinson(fi_var_vals(100, 0.3), fi_gen(100, 0.3))
```
"""
function durbin_levinson(acv::AbstractArray, x::AbstractArray)
    T = length(x)

    if length(acv) < T
        error("At least length(x) autocovariances are needed")
    end

    xhat = zeros(T)
    v = zeros(T)
    phi = zeros(T)
    phiold = zeros(T)

    v[1] = acv[1]

    @inbounds for n = 1:T-1
        # Partial autocorrelation at lag n
        aux = acv[n+1]
        for jj = 1:n-1
            aux -= phiold[jj] * acv[n-jj+1]
        end
        phinn = aux / v[n]

        # Update of the prediction coefficients
        for jj = 1:n-1
            phi[jj] = phiold[jj] - phinn * phiold[n-jj]
        end
        phi[n] = phinn

        v[n+1] = v[n] * (1 - phinn^2)

        for jj = 1:n
            xhat[n+1] += phi[jj] * x[n+1-jj]
        end

        for jj = 1:n
            phiold[jj] = phi[jj]
        end
    end

    return xhat, v
end


"""
    fi_var_vals(T::Int,d::Real)

Computes the autocovariance function of the fractional differenced process with parameter `d` at lags 0, 1, ..., `T-1`.

# Arguments
- `T::Int`: The number of lags to compute.
- `d::Real`: The fractional differencing parameter.

# Output
- `vars::Array`: The autocovariance function of the fractional differenced process with parameter `d` at lags 0, 1, ..., `T-1`.

# Notes
This function uses the recursive formula for the autocovariance function of the fractional differenced process. 

# Examples    
```julia
julia> fi_var_vals(10, 0.4)
```
"""
function fi_var_vals(T::Int, d::Real)
    vars = zeros(T)
    vars[1] = 1
    for k = 1:(T-1)
        vars[k+1] = (d + k - 1) / (k - d) * vars[k]
    end

    vars .*= (gamma(1 - 2 * d) / (gamma(1 - d))^2)
    return vars
end

"""
    fi_cor_vals(T::Int,d::Real)

Computes the autocorrelation function of the fractional differenced process with parameter `d` at lags 0, 1, ..., `T-1`.

# Arguments
- `T::Int`: The number of lags to compute.
- `d::Real`: The fractional differencing parameter.

# Output
- `vars::Array`: The autocorrelation function of the fractional differenced process with parameter `d` at lags 0, 1, ..., `T-1`.

# Notes
This function uses fi_var_vals() to compute the autocovariance function and then normalizes by the variance (first computed value). 

# Examples    
```julia
julia> fi_cor_vals(10, 0.4)
```
"""
function fi_cor_vals(T::Int, d::Real)
    vars = fi_var_vals(T, d)
    cors = vars ./ vars[1]
    return cors
end


"""
    fi_var_matrix(T::Int, d::Real)

Constructs the autocovariance matrix of the fractional differenced process with parameter `d` at lags 0, 1, ..., `T-1`.

# Arguments
- `T::Int`: The number of lags to compute.
- `d::Real`: The fractional differencing parameter.

# Output
- `V::Array`: The autocovariance matrix of the fractional differenced process with parameter `d` at lags 0, 1, ..., `T-1`.

# Examples    
```julia
julia> fi_var_matrix(10, 0.4)
```
"""
function fi_var_matrix(T::Int, d::Real)
    coefs = fi_var_vals(T, d)
    return my_toeplitz(coefs)
end


"""
    fi_llk(d::Real, x::AbstractArray)

Computes the log-likelihood of the fractional differenced process with parameter `d` given the data `x`.

# Arguments
- `d::Real`: The fractional differencing parameter on the unconstrained (logit) scale; it is mapped to `-1/2 + exp(d)/(1+exp(d)) ∈ (-1/2, 1/2)`.
- `x::Array`: The data.

# Output
- `llk::Real`: The negative concentrated log-likelihood (up to constants, divided by `T`), suitable for minimization.

# Notes
This function computes the concentrated log-likelihood function of the fractional differenced process with parameter `d` given the data `x`.

# Examples    
```julia
julia> fi_llk(0.4, randn(100,1))
```
"""
function fi_llk(d::Real, x::AbstractArray)
    d = -1 / 2 + exp(d) / (1 + exp(d))
    T = length(x)
    xhat, v = durbin_levinson(fi_var_vals(T, d), x) # prediction error decomposition
    llk = 1/2 * (sum(log.(v)) / T + log(sum((vec(x) .- xhat) .^ 2 ./ v)))
    return llk
end


"""
    fi_mle_est(x::AbstractArray)

Computes the maximum likelihood estimate of the fractional differencing parameter and the standard deviation of the fractional differenced process given the data `x`.

# Arguments
- `x::Array`: The data.

# Output
- `d::Real`: The maximum likelihood estimate of the fractional differencing parameter.
- `σ::Real`: The maximum likelihood estimate of the standard deviation of the fractional differenced process.

# Notes
This function uses the `Optim` package to minimize the log-likelihood function.
The likelihood is evaluated with the Durbin-Levinson algorithm, see `durbin_levinson`.

# Examples
```julia
julia> fi_mle_est(randn(100,1))
```
"""
function fi_mle_est(x::AbstractArray)
    x = x .- smean(x)
    d0 = gph_est(x)
    d0 = min(0.49, d0)
    d0 = max(-0.49, d0)
    dini = log((d0 + 1 / 2) / (1 / 2 - d0))

    dmle = optimize(d -> fi_llk(first(d), x), [dini]).minimizer[1]
    dmle = - 1 / 2 + exp(dmle) / (1 + exp(dmle))
    xhat, v = durbin_levinson(fi_var_vals(length(x), dmle), x)
    σ = sqrt(sum((vec(x) .- xhat) .^ 2 ./ v) / length(x))

    return dmle, σ
end


"""
    fi_exact_gen(T::Int, d::Real; μ=0, σ=1)

Generate a stationary fractionally differenced process of length `T` with exactly the theoretical autocovariances, using the Davies-Harte method.

# Arguments
- `T::Int`: The length of the time series.
- `d::Real`: The fractional differencing parameter, `d∈(-1/2,1/2)`.

# Optional arguments
- `μ::Real`: The mean of the time series.
- `σ::Real`: The standard deviation of the innovations.

# Output
- `x::Array`: The time series.

# Notes
In contrast to `fi_gen`, which applies the fractional filter truncated at the start of the sample (a nonstationary, so-called type II, process), this function draws from the stationary process. See `davies_harte` for details.

# Examples
```julia
julia> fi_exact_gen(100, 0.4)
```
"""
function fi_exact_gen(T::Int, d::Real; μ=0, σ=1)
    if abs(d) >= 1 / 2
        error("The stationary process requires d∈(-1/2,1/2)")
    end

    x = μ .+ σ .* davies_harte(fi_var_vals(T, d))

    return x
end


"""
    csa_var_vals(T::Int, p::Real, q::Real)

Computes the autocovariance function of the CSA process with parameters `p` and `q` at lags 0, 1, ..., `T-1`.

# Arguments
- `T::Int`: The number of lags to compute.
- `p::Real`: The first parameter of the CSA process.
- `q::Real`: The second parameter of the CSA process.

# Output
- `acf::Array`: The autocovariance function of the CSA process with parameters `p` and `q` at lags 0, 1, ..., `T-1`.

# Notes
This function uses the recursive formula for the autocovariance function of the CSA process, `γ(k) = B(p + k/2, q - 1) / B(p, q)` for unit innovation variance, which requires `q > 1`.

# Examples    
```julia
julia> csa_var_vals(20, 1.4, 1.6)
```
"""
function csa_var_vals(T::Int, p::Real, q::Real)

    s = Int(ceil(T / 2))

    trin_even = zeros(s)
    trin_odd = zeros(s)

    trin_even[1] = gamma(p) / gamma(p + q - 1)
    trin_odd[1] = gamma(p + 1 / 2) / gamma(p + q - 1 / 2)

    for ii in 1:s-1
        trin_even[ii+1] = (p + ii - 1) / (p + q - 2 + ii) * trin_even[ii]
        trin_odd[ii+1] = (p + ii + 1 / 2 - 1) / (p + q + ii - 1 / 2 - 1) * trin_odd[ii]
    end

    trin = zeros(2 * s)
    trin[1:2:2*s] .= trin_even
    trin[2:2:2*s] .= trin_odd

    trin_c = trin[1:T]

    acv = trin_c .* (gamma(p + q) / (gamma(p) * (q - 1))) # γ(k) = B(p + k/2, q - 1) / B(p, q)

    return acv
end


"""
    csa_cor_vals(T::Int, p::Real, q::Real)

Computes the autocorrelation function of the CSA process with parameters `p` and `q` at lags 0, 1, ..., `T-1`.

# Arguments
- `T::Int`: The number of lags to compute.
- `p::Real`: The first parameter of the CSA process.
- `q::Real`: The second parameter of the CSA process.

# Output
- `acf::Array`: The autocorrelation function of the CSA process with parameters `p` and `q` at lags 0, 1, ..., `T-1`.

# Notes
This function uses csa_var_vals() to compute the autocovariance function and then normalizes by the variance (first computed value). 

# Examples    
```julia
julia> csa_cor_vals(20, 1.4, 1.6)
```
"""
function csa_cor_vals(T::Int, p::Real, q::Real)

    vars = csa_var_vals(T, p, q)
    cors = vars ./ vars[1]
    return cors
end


"""
    csa_var_matrix(T::Int, p::Real, q::Real)

Constructs the autocovariance matrix of the CSA process with parameters`p` and `q` at lags 0, 1, ..., `T-1`.

# Arguments
- `T::Int`: The number of lags to compute.
- `p::Real`: The first parameter of the CSA process.
- `q::Real`: The second parameter of the CSA process.

# Output
- `V::Array`: The autocovariance matrix of the CSA process with parameters `p` and `q` at lags 0, 1, ..., `T-1`.

# Examples    
```julia
julia> csa_var_matrix(10, 1.4, 1.8)
```
"""
function csa_var_matrix(T::Int, p::Real, q::Real)
    coefs = csa_var_vals(T, p, q)
    return my_toeplitz(coefs)
end


"""
    csa_llk(p::Real, q::Real, x::AbstractArray)

Computes the log-likelihood of the CSA process with parameters `p` and `q` given the data `x`.

# Arguments
- `p::Real`: The first parameter of the CSA process on the unconstrained scale; it is mapped to `1 + 2exp(p)/(1+exp(p)) ∈ (1, 3)`.
- `q::Real`: The second parameter of the CSA process on the unconstrained scale; it is mapped to `1 + 2exp(q)/(1+exp(q)) ∈ (1, 3)`.
- `x::Array`: The data.

# Output
- `llk::Real`: The negative concentrated log-likelihood (up to constants, divided by `T`), suitable for minimization.

# Notes
This function computes the concentrated log-likelihood function of the CSA process with parameters `p` and `q` given the data `x`.

# Examples    
```julia
julia> csa_llk(1.4, 1.8, randn(100,1))
```
"""
function csa_llk(p::Real, q::Real, x::AbstractArray)
    p = 1 + 2 * (exp(p) / (1 + exp(p)))
    q = 1 + 2 * (exp(q) / (1 + exp(q)))
    T = length(x)
    xhat, v = durbin_levinson(csa_var_vals(T, p, q), x) # prediction error decomposition
    llk = 1/2 * (sum(log.(v)) / T + log(sum((vec(x) .- xhat) .^ 2 ./ v)))
    return llk
end


"""
    csa_mle_est(x::AbstractArray)

Computes the maximum likelihood estimate of the parameters `p` and `q` of the CSA process and the standard deviation of the CSA process given the data `x`.

# Arguments
- `x::Array`: The data.

# Output
- `p::Real`: The maximum likelihood estimate of the first parameter of the CSA process.
- `q::Real`: The maximum likelihood estimate of the second parameter of the CSA process.
- `σ::Real`: The maximum likelihood estimate of the standard deviation of the CSA process.

# Notes
This function uses the `Optim` package to minimize the log-likelihood function.
The parameters are restricted to `p∈(1,3)` and `q∈(1,3)`. The starting values are `p = 2` and `q = 2(1-d)`, with `d` the log-periodogram estimate, so the estimates are reproducible.
The likelihood is evaluated with the Durbin-Levinson algorithm, see `durbin_levinson`.

# Examples
```julia
julia> csa_mle_est(randn(100,1))
```
"""
function csa_mle_est(x::AbstractArray)
    x = x .- smean(x)
    d0 = gph_est(x)
    d0 = min(0.49, d0)
    d0 = max(-0.49, d0)
    pini = 2.0
    qini = 2 * (1 - d0) # q = 2(1-d)

    pini = log((pini - 1) / (3 - pini))
    qini = log((qini - 1) / (3 - qini))

    res = optimize(pq -> csa_llk(first(pq), last(pq), x), [pini, qini]).minimizer

    pmle = 1 + 2 * (exp(res[1]) / (1 + exp(res[1])))
    qmle = 1 + 2 * (exp(res[2]) / (1 + exp(res[2])))

    xhat, v = durbin_levinson(csa_var_vals(length(x), pmle, qmle), x)
    σ = sqrt(sum((vec(x) .- xhat) .^ 2 ./ v) / length(x))

    return pmle, qmle, σ
end


"""
    csa_exact_gen(T::Int, p::Real, q::Real; μ=0, σ=1)

Generate a time series of length `T` with exactly the autocovariances of the cross-sectional aggregated (CSA) process with parameters `p` and `q`, using the Davies-Harte method.

# Arguments
- `T::Int`: The length of the time series.
- `p::Real`: The first parameter of the CSA process.
- `q::Real`: The second parameter of the CSA process, `q > 1`.

# Optional arguments
- `μ::Real`: The mean of the time series.
- `σ::Real`: The standard deviation of the innovations of the individual AR(1) processes.

# Output
- `x::Array`: The time series.

# Notes
The autocovariances are those of `csa_var_vals`, i.e. of the limit of the cross-sectional aggregation of AR(1) processes with squared parameters drawn from a Beta(p,q) distribution, normalized by the number of units.
An error is thrown if the circulant embedding is not nonnegative definite. See `davies_harte` for details.

# Examples
```julia
julia> csa_exact_gen(100, 1.4, 1.6)
```
"""
function csa_exact_gen(T::Int, p::Real, q::Real; μ=0, σ=1)
    if q <= 1
        error("The CSA process has finite variance only for q > 1")
    end

    x = μ .+ σ .* davies_harte(csa_var_vals(T, p, q))

    return x
end


"""
    arfima_whittle_est(x::AbstractArray, p::Int, q::Int)

Estimates the parameters of the ARFIMA(p,d,q) model by minimizing the Whittle approximation to the Gaussian likelihood in the frequency domain.

# Arguments
- `x::AbstractArray`: The data.
- `p::Int`: The order of the autoregressive polynomial.
- `q::Int`: The order of the moving average polynomial.

# Output
- `d::Real`: The estimate of the fractional differencing parameter.
- `AR::Array`: The estimates of the AR coefficients `φ`, with the AR polynomial `1 - φ₁L - ... - φₚLᵖ`.
- `MA::Array`: The estimates of the MA coefficients `θ`, with the MA polynomial `1 + θ₁L + ... + θ_qL^q`.
- `σ::Real`: The estimate of the standard deviation of the innovations.

# Notes
The model is `φ(L)(1-L)^d x_t = θ(L)ε_t`, the same parametrization used by `arfima_gen`.
The objective is the concentrated Whittle function `log(mean(I_j/g_j)) + mean(log(g_j))` over the Fourier frequencies `2πj/T`, `j = 1, ..., floor((T-1)/2)`, where `g` is the spectral density up to the factor `σ²/(2π)`.
Stationarity and invertibility are imposed by optimizing over `d∈(-1/2,1/2)` and over the partial autocorrelations of the AR and MA polynomials (Monahan, 1984). Hence, nonstationary data should be differenced before estimation.
The estimator does not depend on the mean of `x`. It uses the `Optim` package, started at the log-periodogram estimate of `d` and zero ARMA coefficients.

# References
Fox, R., and Taqqu, M. S. (1986). Large-sample properties of parameter estimates for strongly dependent stationary Gaussian time series. The Annals of Statistics, 14(2), 517-532.

# Examples
```julia
julia> arfima_whittle_est(arfima_gen(500, 0, [0.5], 0.3, []), 1, 0)
```
"""
function arfima_whittle_est(x::AbstractArray, p::Int, q::Int)
    T = length(x)

    I_w0, w0 = periodogram(x)

    nf = floor(Int, (T - 1) / 2) # excludes the zero and Nyquist frequencies
    I_w = I_w0[2:nf+1]
    w = w0[2:nf+1]

    d0 = gph_est(x)
    d0 = min(0.49, d0)
    d0 = max(-0.49, d0)
    θini = [log((d0 + 1 / 2) / (1 / 2 - d0)); zeros(p + q)]

    res = optimize(θ -> arfima_whittle_llk(θ, I_w, w, p, q), θini).minimizer

    d, AR, MA = arfima_params(res, p, q)

    g = arfima_spec(w, d, AR, MA)
    σ = sqrt(smean(I_w ./ g))

    return d, AR, MA, σ
end


## Internal functions for the ARFIMA Whittle estimator
"""
    arfima_whittle_llk(θ::AbstractArray, I_w::AbstractArray, w::AbstractArray, p::Int, q::Int)

Computes the concentrated Whittle objective of the ARFIMA(p,d,q) model for the unconstrained parameters `θ`. Internal function used by `arfima_whittle_est`.
"""
function arfima_whittle_llk(θ::AbstractArray, I_w::AbstractArray, w::AbstractArray, p::Int, q::Int)
    d, AR, MA = arfima_params(θ, p, q)
    g = arfima_spec(w, d, AR, MA)

    Q = log(smean(I_w ./ g)) + smean(log.(g))

    return Q
end

"""
    arfima_params(θ::AbstractArray, p::Int, q::Int)

Maps the unconstrained parameters `θ` to `d∈(-1/2,1/2)`, stationary AR coefficients and invertible MA coefficients. Internal function used by `arfima_whittle_est`.
"""
function arfima_params(θ::AbstractArray, p::Int, q::Int)
    d = -1 / 2 + exp(θ[1]) / (1 + exp(θ[1]))
    AR = pacf_to_coefs(tanh.(θ[2:p+1]))
    MA = -pacf_to_coefs(tanh.(θ[p+2:p+q+1])) # 1 + Σθ_k z^k = 1 - Σψ_k z^k

    return d, AR, MA
end

"""
    pacf_to_coefs(r::AbstractArray)

Maps partial autocorrelations `r∈(-1,1)` to the coefficients `φ` of a stationary AR polynomial `1 - φ₁z - ... - φₚzᵖ` using the Durbin-Levinson recursion. Internal function used by `arfima_whittle_est`.
"""
function pacf_to_coefs(r::AbstractArray)
    p = length(r)
    φ = zeros(p)

    for k = 1:p
        aux = φ[1:k-1]
        for jj = 1:k-1
            φ[jj] = aux[jj] - r[k] * aux[k-jj]
        end
        φ[k] = r[k]
    end

    return φ
end

"""
    coefs_to_pacf(φ::AbstractArray)

Maps the coefficients `φ` of a stationary AR polynomial `1 - φ₁z - ... - φₚzᵖ` to its partial autocorrelations, the inverse of `pacf_to_coefs`. Internal function used by `arfima_mle_est`.
"""
function coefs_to_pacf(φ::AbstractArray)
    p = length(φ)
    r = zeros(p)
    aux = float.(collect(φ))

    for k = p:-1:1
        r[k] = aux[k]
        aux = [(aux[jj] + r[k] * aux[k-jj]) / (1 - r[k]^2) for jj = 1:k-1]
    end

    return r
end

"""
    arfima_spec(w::AbstractArray, d::Real, AR::AbstractArray, MA::AbstractArray)

Computes the spectral density of the ARFIMA(p,d,q) model at the frequencies `w`, up to the factor `σ²/(2π)`. Internal function used by `arfima_whittle_est`.
"""
function arfima_spec(w::AbstractArray, d::Real, AR::AbstractArray, MA::AbstractArray)
    z = exp.(-im .* w)

    φz = ones(ComplexF64, length(w))
    for k in eachindex(AR)
        φz .-= AR[k] .* z .^ k
    end

    θz = ones(ComplexF64, length(w))
    for k in eachindex(MA)
        θz .+= MA[k] .* z .^ k
    end

    g = abs2.(θz) ./ abs2.(φz) .* (2 .* sin.(w ./ 2)) .^ (-2 * d)

    return g
end


"""
    arfima_var_vals(T::Int, d::Real, AR::AbstractArray, MA::AbstractArray; σ=1)

Computes the autocovariance function of the stationary ARFIMA(p,d,q) process at lags 0, 1, ..., `T-1`.

# Arguments
- `T::Int`: The number of lags to compute.
- `d::Real`: The fractional differencing parameter, `d∈(-1/2,1/2)`.
- `AR::AbstractArray`: The AR coefficients `φ`, with the AR polynomial `1 - φ₁L - ... - φₚLᵖ`, which must be stationary.
- `MA::AbstractArray`: The MA coefficients `θ`, with the MA polynomial `1 + θ₁L + ... + θ_qL^q`.

# Optional arguments
- `σ::Real`: The standard deviation of the innovations.

# Output
- `acv::Array`: The autocovariance function of the ARFIMA(p,d,q) process at lags 0, 1, ..., `T-1`.

# Notes
The model is `φ(L)(1-L)^d x_t = θ(L)ε_t`, the same parametrization used by `arfima_gen`, `arfima_whittle_est` and `arfima_forecast`.
The autocovariances are computed in two steps. First, the autocovariances of the fractionally integrated MA process `u_t = θ(L)(1-L)^{-d}ε_t` are a finite combination of those of `fi_var_vals`. Second, the autocovariances of `x_t = φ(L)^{-1}u_t` are `γ_x(h) = Σ_k g(k) γ_u(h - k)`, where `g` are the autocovariances of the AR(p) process, which decay geometrically. The sum is truncated when the AR autocovariances are below `1e-15` relative to their variance, so the result is exact up to rounding error.
This is equivalent to the hypergeometric representation of Sowell (1992), but it does not require the roots of the AR polynomial to be distinct.
The computing time grows as the roots of the AR polynomial approach the unit circle.

# References
Sowell, F. (1992). Maximum likelihood estimation of stationary univariate fractionally integrated time series models. Journal of Econometrics, 53(1-3), 165-188.

# Examples
```julia
julia> arfima_var_vals(10, 0.3, [0.5], [0.2])
```
"""
function arfima_var_vals(T::Int, d::Real, AR::AbstractArray, MA::AbstractArray; σ=1)
    if abs(d) >= 1 / 2
        error("The stationary process requires d∈(-1/2,1/2)")
    end

    p = length(AR)
    q = length(MA)

    # Truncation lag of the AR autocovariances
    K = 0
    if p > 0
        C = zeros(p, p) # companion matrix
        C[1, :] .= AR
        for ii = 2:p
            C[ii, ii-1] = 1
        end
        r = maximum(abs.(eigvals(C)))

        if r >= 1
            error("The AR polynomial must be stationary")
        end

        while r^K * (K + 1)^(p - 1) > 1e-15
            K += 1
        end
    end

    # Autocovariances of u_t = θ(L)(1-L)^{-d}ε_t
    γfi = fi_var_vals(T + K + q, d)

    θ = [1; MA]
    ψ = zeros(q + 1) # autocovariances of the MA polynomial
    for l = 0:q, s = 0:q-l
        ψ[l+1] += θ[s+1] * θ[s+l+1]
    end

    γu = zeros(T + K)
    @inbounds for h = 0:T+K-1
        γu[h+1] = ψ[1] * γfi[h+1]
        for l = 1:q
            γu[h+1] += ψ[l+1] * (γfi[abs(h - l)+1] + γfi[h+l+1])
        end
    end

    if p == 0
        acv = γu[1:T]
    else
        # Autocovariances of x_t = φ(L)^{-1}u_t
        g = ar_var_vals(K + 1, AR)

        acv = zeros(T)
        @inbounds for h = 0:T-1
            acv[h+1] = g[1] * γu[h+1]
            for k = 1:K
                acv[h+1] += g[k+1] * (γu[abs(h - k)+1] + γu[h+k+1])
            end
        end
    end

    acv .*= σ^2

    return acv
end


"""
    ar_var_vals(T::Int, AR::AbstractArray)

Computes the autocovariance function of the stationary AR(p) process with unit innovation variance at lags 0, 1, ..., `T-1`, solving the Yule-Walker equations. Internal function used by `arfima_var_vals`.
"""
function ar_var_vals(T::Int, AR::AbstractArray)
    p = length(AR)

    # Yule-Walker equations: g(k) - Σ_j φ_j g(|k-j|) = 1{k = 0}, k = 0, ..., p
    A = zeros(p + 1, p + 1)
    for k = 0:p
        A[k+1, k+1] += 1
        for jj = 1:p
            A[k+1, abs(k - jj)+1] -= AR[jj]
        end
    end
    b = zeros(p + 1)
    b[1] = 1

    g0 = A \ b

    g = zeros(max(T, p + 1))
    g[1:p+1] .= g0
    for k = p+1:T-1
        for jj = 1:p
            g[k+1] += AR[jj] * g[k-jj+1]
        end
    end

    return g[1:T]
end


"""
    arfima_cor_vals(T::Int, d::Real, AR::AbstractArray, MA::AbstractArray)

Computes the autocorrelation function of the stationary ARFIMA(p,d,q) process at lags 0, 1, ..., `T-1`.

# Arguments
- `T::Int`: The number of lags to compute.
- `d::Real`: The fractional differencing parameter, `d∈(-1/2,1/2)`.
- `AR::AbstractArray`: The AR coefficients, with the AR polynomial `1 - φ₁L - ... - φₚLᵖ`.
- `MA::AbstractArray`: The MA coefficients, with the MA polynomial `1 + θ₁L + ... + θ_qL^q`.

# Output
- `cors::Array`: The autocorrelation function of the ARFIMA(p,d,q) process at lags 0, 1, ..., `T-1`.

# Notes
This function uses arfima_var_vals() to compute the autocovariance function and then normalizes by the variance (first computed value).

# Examples
```julia
julia> arfima_cor_vals(10, 0.3, [0.5], [0.2])
```
"""
function arfima_cor_vals(T::Int, d::Real, AR::AbstractArray, MA::AbstractArray)
    vars = arfima_var_vals(T, d, AR, MA)
    cors = vars ./ vars[1]
    return cors
end


"""
    arfima_var_matrix(T::Int, d::Real, AR::AbstractArray, MA::AbstractArray)

Constructs the autocovariance matrix of the stationary ARFIMA(p,d,q) process at lags 0, 1, ..., `T-1`.

# Arguments
- `T::Int`: The number of lags to compute.
- `d::Real`: The fractional differencing parameter, `d∈(-1/2,1/2)`.
- `AR::AbstractArray`: The AR coefficients, with the AR polynomial `1 - φ₁L - ... - φₚLᵖ`.
- `MA::AbstractArray`: The MA coefficients, with the MA polynomial `1 + θ₁L + ... + θ_qL^q`.

# Output
- `V::Array`: The autocovariance matrix of the ARFIMA(p,d,q) process at lags 0, 1, ..., `T-1`.

# Examples
```julia
julia> arfima_var_matrix(10, 0.3, [0.5], [0.2])
```
"""
function arfima_var_matrix(T::Int, d::Real, AR::AbstractArray, MA::AbstractArray)
    coefs = arfima_var_vals(T, d, AR, MA)
    return my_toeplitz(coefs)
end


"""
    arfima_llk(θ::AbstractArray, x::AbstractArray, p::Int, q::Int)

Computes the log-likelihood of the ARFIMA(p,d,q) process given the data `x`.

# Arguments
- `θ::AbstractArray`: The `1 + p + q` parameters on the unconstrained scale: `d` on the logit scale, mapped to `-1/2 + exp(θ₁)/(1+exp(θ₁)) ∈ (-1/2, 1/2)`, followed by the partial autocorrelations of the AR and MA polynomials on the `atanh` scale (Monahan, 1984).
- `x::AbstractArray`: The data (zero mean).
- `p::Int`: The order of the autoregressive polynomial.
- `q::Int`: The order of the moving average polynomial.

# Output
- `llk::Real`: The negative concentrated log-likelihood (up to constants, divided by `T`), suitable for minimization.

# Notes
This function computes the exact concentrated log-likelihood function of the ARFIMA(p,d,q) process, using the autocovariances of `arfima_var_vals` and the Durbin-Levinson algorithm, see `durbin_levinson`.
The parametrization imposes stationarity and invertibility.

# Examples
```julia
julia> arfima_llk([0.4, 0.2, -0.1], randn(100), 1, 1)
```
"""
function arfima_llk(θ::AbstractArray, x::AbstractArray, p::Int, q::Int)
    d, AR, MA = arfima_params(θ, p, q)
    T = length(x)
    xhat, v = durbin_levinson(arfima_var_vals(T, d, AR, MA), x) # prediction error decomposition
    llk = 1/2 * (sum(log.(v)) / T + log(sum((vec(x) .- xhat) .^ 2 ./ v)))
    return llk
end


"""
    arfima_mle_est(x::AbstractArray, p::Int, q::Int)

Computes the exact maximum likelihood estimates of the parameters of the ARFIMA(p,d,q) model given the data `x`.

# Arguments
- `x::AbstractArray`: The data.
- `p::Int`: The order of the autoregressive polynomial.
- `q::Int`: The order of the moving average polynomial.

# Output
- `d::Real`: The maximum likelihood estimate of the fractional differencing parameter.
- `AR::Array`: The maximum likelihood estimates of the AR coefficients `φ`, with the AR polynomial `1 - φ₁L - ... - φₚLᵖ`.
- `MA::Array`: The maximum likelihood estimates of the MA coefficients `θ`, with the MA polynomial `1 + θ₁L + ... + θ_qL^q`.
- `σ::Real`: The maximum likelihood estimate of the standard deviation of the innovations.

# Notes
The model is `φ(L)(1-L)^d x_t = θ(L)ε_t`, the same parametrization used by `arfima_gen` and `arfima_whittle_est`.
The series is demeaned with the sample mean, and the Gaussian likelihood is evaluated exactly with the autocovariances of `arfima_var_vals` and the Durbin-Levinson algorithm (Sowell, 1992).
Stationarity and invertibility are imposed by optimizing over `d∈(-1/2,1/2)` and over the partial autocorrelations of the AR and MA polynomials (Monahan, 1984). Hence, nonstationary data should be differenced before estimation.
This function uses the `Optim` package to minimize the log-likelihood function, started at the Whittle estimates of `arfima_whittle_est`.
With `p = q = 0` it estimates the same model as `fi_mle_est`.
As usual for the exact maximum likelihood estimator with an estimated mean, `d` is biased downward in small samples, more so with AR terms that compete with `d` for the low-frequency dynamics (in simulations of an ARFIMA(1,0.3,0) with `φ₁ = 0.5`, the mean estimate of `d` is about 0.19 for `T = 500` and 0.27 for `T = 2000`). In such cases the likelihood can also have several local maxima. The Whittle estimator `arfima_whittle_est` does not depend on the mean and has a smaller bias.

# References
Sowell, F. (1992). Maximum likelihood estimation of stationary univariate fractionally integrated time series models. Journal of Econometrics, 53(1-3), 165-188.

# Examples
```julia
julia> arfima_mle_est(arfima_gen(500, 0, [0.5], 0.3, []), 1, 0)
```
"""
function arfima_mle_est(x::AbstractArray, p::Int, q::Int)
    x = x .- smean(x)

    d0, AR0, MA0, = arfima_whittle_est(x, p, q)
    d0 = min(0.49, d0)
    d0 = max(-0.49, d0)
    r0 = clamp.([coefs_to_pacf(AR0); coefs_to_pacf(-MA0)], -0.99, 0.99)
    θini = [log((d0 + 1 / 2) / (1 / 2 - d0)); atanh.(r0)]

    res = optimize(θ -> arfima_llk(θ, x, p, q), θini).minimizer

    d, AR, MA = arfima_params(res, p, q)

    xhat, v = durbin_levinson(arfima_var_vals(length(x), d, AR, MA), x)
    σ = sqrt(sum((vec(x) .- xhat) .^ 2 ./ v) / length(x))

    return d, AR, MA, σ
end


"""
    arfima_exact_gen(T::Int, d::Real, AR::AbstractArray, MA::AbstractArray; μ=0, σ=1)

Generate a stationary ARFIMA(p,d,q) process of length `T` with exactly the theoretical autocovariances, using the Davies-Harte method.

# Arguments
- `T::Int`: The length of the time series.
- `d::Real`: The fractional differencing parameter, `d∈(-1/2,1/2)`.
- `AR::AbstractArray`: The AR coefficients, with the AR polynomial `1 - φ₁L - ... - φₚLᵖ`, which must be stationary.
- `MA::AbstractArray`: The MA coefficients, with the MA polynomial `1 + θ₁L + ... + θ_qL^q`.

# Optional arguments
- `μ::Real`: The mean of the time series.
- `σ::Real`: The standard deviation of the innovations.

# Output
- `x::Array`: The time series.

# Notes
In contrast to `arfima_gen`, which applies the fractional filter truncated at the start of the sample (a nonstationary, so-called type II, process), this function draws from the stationary process, with the autocovariances of `arfima_var_vals`.
Note the order of the arguments, `(d, AR, MA)` as in `arfima_var_vals` and `arfima_forecast`, while `arfima_gen` takes `(μ, AR, d, MA)`.
An error is thrown if the circulant embedding is not nonnegative definite. See `davies_harte` for details.

# Examples
```julia
julia> arfima_exact_gen(100, 0.3, [0.5], [0.2])
```
"""
function arfima_exact_gen(T::Int, d::Real, AR::AbstractArray, MA::AbstractArray; μ=0, σ=1)
    x = μ .+ σ .* davies_harte(arfima_var_vals(T, d, AR, MA))

    return x
end


"""
    har_est(x::AbstractArray; m::AbstractArray = [1 , 5 , 22])

Estimates the parameters of the Heterogenous Autoregressive (HAR) model given the data `x`. See [Corsi (2009)](https://academic.oup.com/jfec/article/7/2/174/856522).

# Arguments
- `x::Array`: The data.

# Optional arguments
- `m::Array`: An array with the lags to use in the estimation. By default, the lags are 1, 5, and 22; as suggested by the original paper.

# Output
- `β::Array`: The estimated parameters of the HAR model.
- `σ::Real`: The estimated standard deviation of the HAR model.

# Examples    
```julia
julia> har_est(randn(100,1))
```
"""
function har_est(x::AbstractArray; m::AbstractArray=[1, 5, 22])
    T = length(x)
    n = length(m)
    m = sort(m)

    mm = maximum(m)
    
    X = zeros(T - mm, n + 1)
    X[:, 1] = ones(T - mm, 1)

    for ii = 1:n
        cm = m[ii]
        aux = zeros(T - mm)
        for jj = 1:cm
            aux .+= @view x[mm-jj+1:T-jj, 1]
        end
        X[:, ii+1] = aux / cm
    end

    Y = x[mm+1:T, 1]

    β = X \ Y
    err = Y - X * β

    σ = sqrt( (err' * err) / (T - mm - n - 1) )

    return β, σ

end


end # module ParametricEstimators
