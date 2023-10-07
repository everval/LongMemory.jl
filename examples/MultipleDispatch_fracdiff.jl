using FFTW
# File to benchmark the performance of the fracdiff function using multiple dispatch.

# Float64 version
function fracdiff(x::Array,d::Float64)
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

# Int64 version
function fracdiff(x::Array, d::Int)
    T = length(x)

    if d == 0
        dx = x
    elseif d == 1
        dx = diff(x, dims=1)
    end

    return dx
end

# Parameters
R = 10^4; T = 10^3;

# Time using the Float64 version
@time begin
	for i = 1:R
		fracdiff(randn(T),0.0)
	end
end
		
# Time using the Int64 version
@time begin
	for i = 1:R
		fracdiff(randn(T),0)
	end
end