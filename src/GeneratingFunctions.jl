using FFTW, SpecialFunctions

export fracdiff, csadiff

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
function fracdiff(x::Array,d)
    T = length(x)

    np2 = nextpow(2,2*T-1)
    k = 1:(T-1)
    b = [1;cumprod((k.-d.-1)./k)]
    padb = [b;zeros(np2-T,1)]
    padx = [x;zeros(np2-T,1)]
    dx = irfft(rfft(padx).*rfft(padb),np2)
    dx = dx[1:T]

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
function csadiff(x::Array,p,q)
    T = length(x)

    np2 = nextpow(2,2*T-1)
    coefs = ( beta.(p.+(0:T-1),q) ./ beta(p,q) ).^(1/2);
    padcoefs = [coefs;zeros(np2-T,1)]
    padx = [x;zeros(np2-T,1)]
    dx = irfft(rfft(padx).*rfft(padcoefs),np2)
    dx = dx[1:T]

    return dx
end
