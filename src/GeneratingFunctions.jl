using FFTW, SpecialFunctions, Distributions

export fracdiff, csadiff, csagen, fi


"""
    fracdiff(x,d)

Compute the fractional difference of a time series `x` with fractional order `d∈(-1/2,1/2)`.

The function uses the fast Fourier transform to compute the convolution of the time series with the
fractional difference filter. See [Jensen and Nielsen (2014)](https://onlinelibrary.wiley.com/doi/10.1111/jtsa.12074) for details.

# Arguments
- `x::Vector`: time series
- `d::Float64`: fractional difference parameter

# Examples
```julia-repl
julia> fracdiff(randn(100,1),0.4)
```
"""
function fracdiff(x::Array, d::Float64)
    T = length(x)

    np2 = nextpow(2, 2 * T - 1)
    k = 1:(T-1)
    b = [1; cumprod((k .- d .- 1) ./ k)]
    padb = [b; zeros(np2 - T, 1)]
    padx = [x; zeros(np2 - T, 1)]
    dx = irfft(rfft(padx) .* rfft(padb), np2)
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

# Notes
`q` determines the long memory parameter of the cross-sectional aggregated process. The relation `q = 2(1-d)` holds, where `d` is the fractional difference parameter.	

# Examples
```julia-repl
julia> csadiff(randn(100,1),1.2,1.4)
```
"""
function csadiff(x::Array, p, q)
    T = length(x)

    np2 = nextpow(2, 2 * T - 1)
    coefs = (beta.(p .+ (0:T-1), q) ./ beta(p, q)) .^ (1 / 2)
    padcoefs = [coefs; zeros(np2 - T, 1)]
    padx = [x; zeros(np2 - T, 1)]
    dx = irfft(rfft(padx) .* rfft(padcoefs), np2)
    dx = dx[1:T]

    return dx
end


"""
    csagen(T::Int,p,q;μ=0,σ=1)

Generate a time series with long memory parameter `q` and length `T` using the cross-sectional aggregated process. See [Vera-Valdes(2021)](https://www.mdpi.com/2225-1146/9/4/39) for details.

# Arguments
- `T::Int`: length of the time series
- `p::Float64`: first parameter of the cross-sectional aggregated process
- `q::Float64`: second parameter of the cross-sectional aggregated process, which is related to the fractional difference parameter `d` by `q = 2(1-d)`

# Optional arguments
- `μ::Float64`: mean of the time series
- `σ::Float64`: standard deviation of the time series

# Output
- `x::Vector`: time series

# Examples
```julia-repl
julia> csagen(100,1.2,1.4)
```
"""
function csagen(T::Int, p, q; μ=0, σ=1)
    x = csadiff(rand(Normal(μ, σ), T), p, q)

    return x
end


"""
    csagen(T::Int,N::Int,p,q;t=0.01;μ=0,σ=1)

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
- `x::Vector`: time series

# Examples
```julia-repl
julia> csagen(100,100,1.2,1.4)
```
"""
function csagen(T::Int, N::Int, p, q; t=0.01, μ=0, σ=1)

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
- `x::Vector`: time series

# Notes     
Multiple dispatch is used for generation: If `d` is an integer, the function returns a time series with first or null difference.
See `fracdiff` for details.

# Examples
```julia-repl
julia> fi(100,0.4)
```
"""
function fi(T::Int, d; μ=0, σ=1)
    x = fracdiff(rand(Normal(μ, σ), T), d)

    return x
end



function edmgen(T::Int, d; t=0.05)

    first = Int(round(T * t))

    p = fi_survival_probs(T, d)

    x = 0

    return x

end


""""
    fi_survival_probs(N::Int,d)

Generate the survival probabilities of the fractional difference filter à la Parke (1999).

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