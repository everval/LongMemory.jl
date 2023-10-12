using Optim, LinearAlgebra

export fi_mle_est

"""
    fi_var_vals(T::Int,d::Real)

Computes the autocovariance function of the fractional differenced process with parameter `d` at lags 0, 1, ..., `T-1`.

# Arguments
- `T::Int`: The number of lags to compute.
- `d::Real`: The fractional differencing parameter.

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

    return vars

end


"""
    my_toeplitz(coefs::Array)

Constructs a Toeplitz matrix from the given coefficients.

# Arguments
- `coefs::Array`: An array of coefficients.

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
    fi_var_matrix(T::Int, d::Real)

Constructs the autocovariance matrix of the fractional differenced process with parameter `d` at lags 0, 1, ..., `T-1`.

# Arguments
- `T::Int`: The number of lags to compute.
- `d::Real`: The fractional differencing parameter.

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

# Notes
This function computes the concentrated log-likelihood function of the fractional differenced process with parameter `d` given the data `x`.

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

# Notes
This function uses the `Optim` package to optimize the log-likelihood function.

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
    dmle = -1/2 + exp(dmle) / (1 + exp(dmle))
    V = fi_var_matrix(length(x), dmle)
    σ2 = (x'/V*x)[1, 1] / length(x)

    return dmle, σ2

end
