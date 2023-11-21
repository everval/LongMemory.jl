"""
# GeneratingFunctions

This module contains functions to generate time series with long memory.


## Author
- [J. Eduardo Vera-Valdés](https://everval.github.io/)

"""
module GeneratingFunctions

using FFTW, Distributions

export fracdiff, csadiff, csa_gen, edm_gen, fi, fi_gen, arfi_gen, arfima_gen


"""
    fracdiff(x,d)

Compute the fractional difference of a time series `x` with fractional order `d∈(-1/2,1/2)`.

# Arguments
- `x::Vector`: time series
- `d::Float64`: fractional difference parameter

# Output
- `dx::Vector`: fractional difference of `x`

# Notes
The function uses the fast Fourier transform to compute the convolution of the time series with the fractional difference filter. 
See [Jensen and Nielsen (2014)](https://onlinelibrary.wiley.com/doi/10.1111/jtsa.12074) for details.
We use autoregressive formulas to efficiently compute the coefficients of the fractional difference filter.

# Examples
```julia-repl
julia> fracdiff(randn(100,1),0.4)
```
"""
function fracdiff(x::Array, d::Float64)
    T = length(x)

    np2 = nextpow(2, 2 * T - 1)

    coefs = zeros(T)
    coefs[1] = 1
    for t in 1:T-1
        coefs[t+1] = coefs[t] * ((t - d - 1) / t)
    end

    padcoefs = [coefs; zeros(np2 - T, 1)]
    padx = [x; zeros(np2 - T, 1)]

    dx = irfft(rfft(padx) .* rfft(padcoefs), np2)
    dx = dx[1:T]

    return dx
end


"""
    fracdiff(x,d::Int)

Compute the first or null difference of a time series `x`.
Multiple dispatch is used to return the same input or call the function `diff` from the Julia standard library if `d=1` or `d=0`, respectively.

# Arguments
- `x::Vector`: time series
- `d::Int64`: difference parameter

# Output
- `dx::Vector`: first or null difference of `x`

# Examples
```julia-repl
julia> fracdiff(randn(100,1),1)
```
"""
function fracdiff(x::Array, d::Int)
    T = length(x)

    if d == 0
        dx = x
    elseif d == 1
        dx = diff(x, dims=1)
    end

    return dx
end


"""
    csadiff(x,p,q)

Generate long memory by using the moving average representation of the cross-sectional aggregated process using the fast Fourier algorithm. See [Vera-Valdes(2021)](https://www.mdpi.com/2225-1146/9/4/39) for details.

# Arguments
- `x::Vector`: time series
- `p::Float64`: first parameter of the cross-sectional aggregated process
- `q::Float64`: second parameter of the cross-sectional aggregated process, which is related to the fractional difference parameter `d` by `q = 2(1-d)`

# Output
- `dx::Vector`: time series with long memory

# Notes
`q` determines the long memory parameter of the cross-sectional aggregated process. The relation `q = 2(1-d)` holds, where `d` is the fractional difference parameter.
We use autoregressive formulas to efficiently compute the coefficients of the moving average representation of the cross-sectional aggregated process. 

# Examples
```julia-repl
julia> csadiff(randn(100,1),1.2,1.4)
```
"""
function csadiff(x::Array, p, q)
    T = length(x)

    np2 = nextpow(2, 2 * T - 1)

    coefs = zeros(T)
    coefs[1] = 1
    for t in 2:T
        coefs[t] = coefs[t-1] * ((p + t - 2) / (p + t - 2 + q))^(1 / 2)
    end

    padcoefs = [coefs; zeros(np2 - T, 1)]
    padx = [x; zeros(np2 - T, 1)]
    dx = irfft(rfft(padx) .* rfft(padcoefs), np2)
    dx = dx[1:T]

    return dx
end


"""
    csa_gen(T::Int,p,q;μ=0,σ=1)

Generate a time series with long memory parameter `q` and length `T` using the cross-sectional aggregated process. 
See [Vera-Valdes(2021)](https://www.mdpi.com/2225-1146/9/4/39) for details.

# Arguments
- `T::Int`: length of the time series
- `p::Float64`: first parameter of the cross-sectional aggregated process
- `q::Float64`: second parameter of the cross-sectional aggregated process, which is related to the fractional difference parameter `d` by `q = 2(1-d)`

# Optional arguments
- `μ::Float64`: mean of the time series
- `σ::Float64`: standard deviation of the time series

# Output
- `x::Vector`: time series with long memory

# Examples
```julia-repl
julia> csa_gen(100,1.2,1.4)
```
"""
function csa_gen(T::Int, p, q; μ=0, σ=1)
    x = csadiff(rand(Normal(μ, σ), T), p, q)

    return x
end


"""
    csa_gen(T::Int,N::Int,p,q;t=0.01;μ=0,σ=1)

Generate a time series with long memory parameter `q` and length `T` using the cross-sectional aggregation of 'N' AR(1) processes à la Granger (1980).

# Arguments
- `T::Int`: length of the time series
- `N::Int`: number of AR(1) processes
- `p::Float64`: first parameter of the cross-sectional aggregated process
- `q::Float64`: second parameter of the cross-sectional aggregated process, which is related to the fractional difference parameter `d` by `q = 2(1-d)`

# Optional arguments
- `t::Float64`: taper length
- `μ::Float64`: mean of the time series
- `σ::Float64`: standard deviation of the time series

# Notes
Multiple dispatch is used to generate the finite sample process if 'N' is included in the arguments.

# Output
- `x::Vector`: time series with long memory

# Examples
```julia-repl
julia> csa_gen(100,100,1.2,1.4)
```
"""
function csa_gen(T::Int, N::Int, p, q; t=0.01, μ=0, σ=1)

    params = sqrt.(rand(Beta(p, q), N)) #generate AR parameters from a Beta distribution

    t = max(Int(round(T * t)), 10) #taper length

    X = zeros(T + t, N)

    errors = rand(Normal(μ, σ), T + t, N)

    for i = 1:N
        for j = 2:T+t
            X[j, i] = params[i] * X[j-1, i] + errors[j, i]
        end
    end

    x = sum(X[t+1:T+t, :], dims=2)

    return x

end


"""
    fi(T,d;μ=0,σ=1)

Generate a time series with long memory parameter `d` and length `T` using the fractional difference filter.

# Arguments
- `T::Int`: length of the time series
- `d::Float64`: fractional difference parameter

# Optional arguments
- `μ::Float64`: mean of the time series
- `σ::Float64`: standard deviation of the time series

# Output
- `x::Vector`: time series with long memory

# Notes     
Multiple dispatch is used for generation: If `d` is an integer, the function returns a time series with first or null difference.
See `fracdiff` for details.

# Examples
```julia-repl
julia> fi(100,0.4)
```
"""
function fi(T::Int, d; μ=0, σ=1)
    x = fracdiff(rand(Normal(μ, σ), T), -d)

    return x
end


"""
    fi_gen(T,d;μ=0,σ=1)

Generate a time series with long memory parameter `d` and length `T` using the fractional difference filter.

# Arguments
- `T::Int`: length of the time series
- `d::Float64`: fractional difference parameter

# Optional arguments
- `μ::Float64`: mean of the time series
- `σ::Float64`: standard deviation of the time series

# Output
- `x::Vector`: time series with long memory

# Notes     
Multiple dispatch is used for generation: If `d` is an integer, the function returns a time series with first or null difference.
See `fracdiff` for details.

# Examples
```julia-repl
julia> fi_gen(100,0.4)
```
"""
function fi_gen(T::Int, d; μ=0, σ=1)
    x = fracdiff(rand(Normal(μ, σ), T), -d)

    return x
end


"""
    arfima_gen(T::Int, μ::Real, AR::Array, d::Real, MA::Array; σ=1)

Generate a time series with long memory parameter `d` and length `T` using the ARFIMA(p,d,q) model.

# Arguments
- `T::Int`: length of the time series
- `μ::Float64`: mean of the time series
- `AR::Array`: AR coefficients
- `d::Float64`: fractional difference parameter
- `MA::Array`: MA coefficients

# Optional arguments
- `σ::Float64`: standard deviation of the time series

# Output
- `x::Vector`: time series with long memory

# Notes
The code is inspired by the function `dgp_arfima.m` by [Carlos Vladimir Rodríguez Caballero (2023)](https://www.mathworks.com/matlabcentral/fileexchange/53301-arfima-p-d-q)

# Examples
```julia-repl
julia> arfima_gen(100, 0, [0.2; -0.5], 0.4, [-0.3; 0.1]])
```
"""
function arfima_gen(T::Int, μ::Real, AR::Array, d::Real, MA::Array; σ=1)
    p = length(AR)
    q = length(MA)

    if q == 0
        u = rand(Normal(0, σ), T + p, 1)
    else
        uM = zeros(T, q + 1)
        u = rand(Normal(0, σ), T + p + q, 1)
        MAb = [1; MA]

        for ii in 1:q+1
            for jj in 1:ii
                uM[:, ii] .+= MAb[jj] * u[q+1-jj+1:T+q-jj+1]
            end
        end
        u = sum(uM, dims=2)
    end

    fiu = fracdiff(u, -d)  #fracdiff with multiple dispatch

    if p == 0
        x = fiu
    else
        xb = zeros(T + p)
        for ii in p+1:T+p
            for jj in 1:p
                xb[ii] += AR[jj] * xb[ii-jj]
            end
            xb[ii] = μ + xb[ii] + fiu[ii-p]
        end
        x = xb[p+1:T+p]
    end
    return x
end


"""
    arfi_gen(T::Int, μ::Real, AR::Array, d::Real; σ=1)

Generate a time series with long memory parameter `d` and length `T` using the ARFIMA(p,d,0) model.

# Arguments
- `T::Int`: length of the time series
- `μ::Float64`: mean of the time series
- `AR::Array`: AR coefficients
- `d::Float64`: fractional difference parameter

# Optional arguments
- `σ::Float64`: standard deviation of the time series

# Output
- `x::Vector`: time series with long memory

# Notes
The code is inspired by the function `dgp_arfima.m` by [Carlos Vladimir Rodríguez Caballero (2023)](https://www.mathworks.com/matlabcentral/fileexchange/53301-arfima-p-d-q)

# Examples
```julia-repl
julia> arfi_gen(100, 0, [0.2; -0.5], 0.4])
```
"""
function arfi_gen(T::Int, μ::Real, AR, d::Real; σ=1)
    p = length(AR)

    u = rand(Normal(0, σ), T + p, 1)

    fiu = fracdiff(u, -d)  #fracdiff with multiple dispatch

    xb = zeros(T + p)
    for ii in p+1:T+p
        for jj in 1:p
            xb[ii] += AR[jj] * xb[ii-jj]
        end
        xb[ii] = μ + xb[ii] + fiu[ii-p]
    end
    x = xb[p+1:T+p]

    return x
end


"""
    edm_gen(T::Int,d; t=0.5, μ=0, σ=1)

Generate a time series with long memory parameter `d` and length `T` using the error duration model à la Parke (1999).

# Arguments
- `T::Int`: length of the time series
- `d::Float64`: long memory parameter

# Optional arguments
- `t::Float64`: taper length
- `μ::Float64`: mean of the time series
- `σ::Float64`: standard deviation of the time series

# Output
- `x::Vector`: time series with long memory

# Notes
The taper length `t` is the proportion of the time series that is pre-sampled to avoid the initial bias of the error duration model.

# Examples
```julia-repl
julia> edm_gen(100,0.4)
```
"""
function edm_gen(T::Int, d; t=0.5, μ=0, σ=1)

    t = max(Int(round(T * t)), 200) # Pre-sample

    p = fi_survival_probs(T + t, d) # Survival probabilities

    u = rand(T + t, 1)

    tiempos = searchsortedlast.(Ref(p[:, 1]), u, rev=true) # Random times

    x = zeros(T + t, 1) # Pre-allocate
    err = rand(Normal(μ, σ), T + t) # Error term

    for ii = 1:(T+t)
        uxu = zeros(T + t, 1) # Pre-allocate
        uxu[ii:min((ii + tiempos[ii] - 1), T + t), 1] .= err[ii] # Error term surviving
        x = x + uxu
    end

    return x

end


""""
    fi_survival_probs(N::Int,d)

Generate the survival probabilities of the error duration model à la Parke (1999).

# Arguments
- `N::Int`: length of the time series
- `d::Float64`: fractional difference parameter

# Output
- `p::Vector`: survival probabilities

# Notes
The survival probabilities are computed using the recursive formula `p_{t+1} = p_t * (t + d - 1) / (t + 1 - d)` to avoid numerical overflow.

# Examples
```julia-repl
julia> fi_survival_probs(100,0.4)
```
"""
function fi_survival_probs(N::Int, d)

    p = zeros(N, 1)
    p[1, 1] = 1
    for ii = 1:N-1
        p[ii+1, 1] = p[ii, 1] * (ii - 1 + d) / (ii + 1 - d)
    end

    return p

end

end #module
