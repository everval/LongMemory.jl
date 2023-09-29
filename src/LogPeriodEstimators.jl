using FFTW

export periodogram, gph_est

"""
    periodogram(x)

Compute the periodogram of a time series `x` using the fast Fourier transform.

# Arguments
- `x::Vector`: time series

# Examples
```julia-repl
julia> periodogram(randn(100,1))
```
"""
function periodogram(x::Array)
    T = length(x)

    I_w = abs.(rfft(x)) .^ 2 ./ T # periodogram
    w = 2pi * (0:T-1) ./ T  # Fourier frequencies

    # int rounds to nearest integer. We want to round up or take 1/2 + 1 to
    # make sure we get the whole interval from [0, pi]
    ind = iseven(T) ? round(Int, T / 2 + 1) : ceil(Int, T / 2)
    I_w, w = I_w[1:ind], w[1:ind]
    return I_w, w
end

"""
    gph_est(x;m=0.5,l=0,br=0)

Estimate the long memory parameter of a time series `x` using the log-periodogram estimator. See [Geweke and Porter-Hudak (1983)](https://onlinelibrary.wiley.com/doi/10.1111/j.1467-9892.1983.tb00371.x) and [Andrews and Guggenberger (2003)](https://www.jstor.org/stable/3082070) for details.

# Arguments
- `x::Vector`: time series
- `m∈(0,1)::Float64`: taper final
- `l∈(0,1)::Float64`: taper initial
- `br::Int64`: number of bias reduction terms

# Notes
The function considers the periodogram of the time series `x` for frequencies in the interval `[T^l,T^m]`. The zero frequency is always excluded.
The default values of `m` and `l` are 0.5 and 0, respectively.
The condition `m < l` must hold.

The default value of `br` is 0 which returns the original GPH log-periodogram estimator.

# Examples
```julia-repl
julia> gph_est(randn(100,1))
```
"""
function gph_est(x::Array;m=0.5,l=0,br=0::Int)
    T = length(x)

    if m < l
       error("Taper initial is greater than final")
    end

    last = Int(round(T^m))
    first = Int(max(round(T^l),2))

    I_w, w = periodogram(x)

    I_w = I_w[first:last]
    w = w[first:last]

    Y = log.(I_w)
    X = [-2*log.(w) w.^[0:2:(2*br)...]']  
    β = X\Y

    return β[1]
end
