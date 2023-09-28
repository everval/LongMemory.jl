using FFTW

export fracdiff

function fracdiff(x,d)
    T,c = size(x)

    if c > 1
        error("x must be a column vector")
    end

    np2 = nextpow(2,2*T-1)
    k = 1:(T-1)
    b = [1;cumprod((k.-d.-1)./k)]
    padb = [b;zeros(np2-T,1)]
    padx = [x;zeros(np2-T,1)]
    dx = irfft(rfft(padx).*rfft(padb),np2)
    dx = dx[1:T]

    return dx
end