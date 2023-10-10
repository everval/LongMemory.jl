using SpecialFunctions

p = 1.2; q = 1.8; T = 10^3
# To check speed of betas against recursive
function beta_non(p,q,T::Int)
    coefs = (beta.(p .+ (0:T-1), q) ./ beta(p, q)) .^ (1 / 2)
end

# Recursive
function beta_rec(p, q, T::Int)
    coefs = zeros(T)
    coefs[1] = 1
    for t in 2:T
        coefs[t] = coefs[t-1] * ( (p + t - 2) / (p + t - 2 + q) )^(1/2)
    end
    coefs
end

function beta_cum(p,q,T::Int)
       
    coefs = [1; cumprod(  ((p .+ (2:T) .- 2) ./ (p .+ (2:T) .- 2 .+ q) ).^(1/2) ) ]

    return coefs
end


R = 10^4

@time begin
    for r in 1:R
        coefs = beta_non(p,q,T)
    end
end

@time begin
    for r in 1:R
        coefsr = beta_rec(p, q, T)
    end
end

@time begin
    for r in 1:R
        coefsc = beta_cum(p, q, T)
    end
end

T=10;[beta_non(p,q,T) beta_rec(p, q, T) beta_cum(p, q, T)]