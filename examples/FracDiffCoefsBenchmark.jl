d = 0.2; T = 10^4

# Recursive
function fi_cum(d, T::Int)
    k = 1:(T-1)
    coefs = [1; cumprod((k .- d .- 1) ./ k)]

    return coefs
end

function fi_rec(d,T::Int)
    
    coefs = zeros(T)
    coefs[1] = 1
    for t in 1:T-1
        coefs[t+1] = coefs[t] * ( (t - d - 1) / t )
    end

    return coefs
end


R = 10^4

@time begin
    for r in 1:R
        coefsr = fi_rec(d, T)
    end
end

@time begin
    for r in 1:R
        coefsc = fi_cum(d, T)
    end
end

fi_cum(d,T) ≈  fi_rec(d,T)