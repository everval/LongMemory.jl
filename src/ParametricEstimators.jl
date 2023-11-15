"""
# ParametricEstimators

This module contains functions for estimating the parameters of the fractional differenced process à la ARFIMA and the CSA process.

## Author
[J. Eduardo Vera-Valdés](https://everval.github.io/)

"""
module ParametricEstimators

include("LogPeriodEstimators.jl")
import .LogPeriodEstimators: gph_est

using Optim, LinearAlgebra, SpecialFunctions

export fimle_est, csamle_est, har_est

"""
    my_toeplitz(coefs::Array)

Constructs a Toeplitz matrix from the given coefficients.

# Arguments
- `coefs::Array`: An array of coefficients.

# Output
- `Toep::Array`: The Toeplitz matrix constructed from the given coefficients.

# Examples    
```julia
julia> my_toeplitz([1, 2, 3])
```
"""
function my_toeplitz(coefs::Array)
    N = length(coefs)
    Toep = zeros(N, N)

    for ii = 1:N
        for jj = 1:N
            if ii >= jj
                Toep[ii, jj] = coefs[ii-jj+1]
            else
                Toep[ii, jj] = coefs[jj-ii+1]
            end
        end
    end

    return Toep
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

    vars = vars.*( gamma(1-2*d)/(gamma(1-d))^2 )
    return vars
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
    fi_llk(x::Array,d::Real)

Computes the log-likelihood of the fractional differenced process with parameter `d` given the data `x`.

# Arguments
- `x::Array`: The data.
- `d::Real`: The fractional differencing parameter.

# Output
- `llk::Real`: The log-likelihood of the fractional differenced process with parameter `d` given the data `x`.

# Notes
This function computes the concentrated log-likelihood function of the fractional differenced process with parameter `d` given the data `x`.
The function is inspired by the `arfima.Estimate()` function in Ox; see [Doornik (1999)](http://fmwww.bc.edu/ec-p/software/ox/Ox.arfima.v2.1.pdf).	

# Examples    
```julia
julia> fi_llk(randn(100,1), 0.4)
```
"""
function fi_llk(d::Real, x::Array)
    d = -1 / 2 + exp(d) / (1 + exp(d))
    T = length(x)
    V = fi_var_matrix(T, d)
    llk = logdet(V) / T + log((x'/V*x)[1, 1] / T)
    # σ2 = (x'/V*x)[1, 1] /T
    return llk
end


"""
    fimle_est(x::Array)

Computes the maximum likelihood estimate of the fractional differencing parameter and the variance of the fractional differenced process given the data `x`.

# Arguments
- `x::Array`: The data.

# Output
- `d::Real`: The maximum likelihood estimate of the fractional differencing parameter.
- `σ2::Real`: The maximum likelihood estimate of the variance of the fractional differenced process.

# Notes
This function uses the `Optim` package to minimize the log-likelihood function.
The function is inspired by the `arfima.Estimate()` function in Ox; see [Doornik (1999)](http://fmwww.bc.edu/ec-p/software/ox/Ox.arfima.v2.1.pdf).	

# Examples    
```julia
julia> fimle_est(randn(100,1))
```
"""
function fimle_est(x::Array)
    d0 = gph_est(x)
    d0 = min(0.49, d0)
    d0 = max(-0.49, d0)
    dini = log((d0 + 1 / 2) / (1 / 2 - d0))

    dmle = optimize(d -> fi_llk(first(d), x), [dini]).minimizer[1]
    dmle = -1 / 2 + exp(dmle) / (1 + exp(dmle))
    V = fi_var_matrix(length(x), dmle)
    σ2 = (x'/V*x)[1, 1] / length(x)

    return dmle, σ2
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
This function uses the recursive formula for the autocovariance function of the CSA process.

# Examples    
```julia
julia> csa_var_vals(20, 0.4, 0.6)
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

    acv = trin_c .* (gamma(p + q) / (gamma(p) * (q - 1))) * (beta(p, q))

    return acv
end

"""
    csa_var_matrix(T::Int, d::Real)

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
    csa_llk(p::Real, q::Real, x::Array)

Computes the log-likelihood of the CSA process with parameters `p` and `q` given the data `x`.

# Arguments
- `p::Real`: The first parameter of the CSA process.
- `q::Real`: The second parameter of the CSA process.
- `x::Array`: The data.

# Output
- `llk::Real`: The log-likelihood of the CSA process with parameters `p` and `q` given the data `x`.

# Notes
This function computes the concentrated log-likelihood function of the CSA process with parameters `p` and `q` given the data `x`.

# Examples    
```julia
julia> csa_llk(1.4, 1.8, randn(100,1))
```
"""
function csa_llk(p::Real, q::Real, x::Array)
    p = 1 + 2 * ( exp(p) / (1+exp(p)) )
    q = 1 + 2 * ( exp(q) / (1+exp(q)) )
    T = length(x)
    V = csa_var_matrix(T, p, q)
    llk = logdet(V) / T + log((x'/V*x)[1, 1] / T)
    return llk
end


"""
    csa_mle_est(x::Array)

Computes the maximum likelihood estimate of the parameters `p` and `q` of the CSA process and the variance of the CSA process given the data `x`.

# Arguments
- `x::Array`: The data.

# Output
- `p::Real`: The maximum likelihood estimate of the first parameter of the CSA process.
- `q::Real`: The maximum likelihood estimate of the second parameter of the CSA process.
- `σ2::Real`: The maximum likelihood estimate of the variance of the CSA process.

# Notes
This function uses the `Optim` package to minimize the log-likelihood function.

# Examples    
```julia
julia> csa_mle_est(randn(100,1))
```
"""
function csamle_est(x::Array)
    pini = 1+rand()
    qini = 1+rand()

    pini = log((pini-1)/(3-pini))
    qini = log((qini-1)/(3-qini))

    res = optimize(pq -> csa_llk(first(pq), last(pq), x), [pini, qini]).minimizer

    pmle = 1 + 2 * ( exp(res[1]) / (1+exp(res[1])) );
    qmle = 1 + 2 * ( exp(res[2]) / (1+exp(res[2])) );

    V = csa_var_matrix(length(x), pmle, qmle)
    σ2 = (x'/V*x)[1, 1] / length(x)

    return pmle, qmle, σ2
end


"""
    har_est(x::Array; m::Array = [1 , 5 , 22])

Estimates the parameters of the Heterogenous Autoregressive (HAR) model given the data `x`. See [Corsi (2009)](https://academic.oup.com/jfec/article/7/2/174/856522).

# Arguments
- `x::Array`: The data.

# Optional arguments
- `m::Array`: An array with the lags to use in the estimation. By default, the lags are 1, 5, and 22; as suggested by the original paper.

# Output
- `betas::Array`: The estimated parameters of the HAR model.
- `sigma::Real`: The estimated variance of the HAR model.

# Examples    
```julia
julia> har_est(randn(100,1))
```
"""
function har_est(x::Array; m::Array=[1,5,22])
    T = length(x)
    n = length(m)
    sort!(m)

    mm = maximum(m)

    if mm != m[end]
        error("The maximum lag must be the last value in the array.")
    end

    X = zeros(T-mm, n+1)
    X[:,1] = ones(T-mm, 1)

    for ii = 1:n
        cm = m[ii]
        aux = zeros(T-mm, 1)
        for jj = 1:cm
            aux = aux + x[mm-jj+1:T-jj,1]
        end
        X[:,ii+1] = aux/cm
    end

    Y = x[mm+1:T,1]

    betas = X\Y
    err = Y-X*betas

    sigma = (err'*err)/(T-mm-n-1)

    return betas, sigma

end


end # module ParametricEstimators
