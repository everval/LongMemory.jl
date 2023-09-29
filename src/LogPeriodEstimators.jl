using FFTW

export periodogram

function periodogram(x)
    T,c = size(x)

    if c > 1
        error("x must be a vector")
    end

    I_w = abs.(rfft(x)) .^ 2 ./ T # periodogram
    w = 2pi * (0:T-1) ./ T  # Fourier frequencies

    # int rounds to nearest integer. We want to round up or take 1/2 + 1 to
    # make sure we get the whole interval from [0, pi]
    ind = iseven(T) ? round(Int, T / 2 + 1) : ceil(Int, T / 2)
    w, I_w = w[1:ind], I_w[1:ind]
    return w, I_w
end

function gph_est(x)
    T,c = size(x)

    if c > 1
        error("x must be a vector")
    end

    I_w = abs.(rfft(x)) .^ 2 ./ T # periodogram
    w = 2pi * (0:T-1) ./ T  # Fourier frequencies
    
    return 0
end