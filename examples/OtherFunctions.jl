using SpecialFunctions

function fi_var_vals_re(T::Int, d::Real)
    k = 0:(T-1)

    vars = (gamma.(k .+ d) .* gamma(1 - d)) ./ (gamma.(k .- d .+ 1) .* gamma(d))

    return vars

end

function fi_var_vals(T::Int, d::Real)

    vars = zeros(T)
    vars[1] = 1
    for k = 1:(T-1)
        vars[k+1] = (d + k - 1) / (k - d) * vars[k]
    end

    return vars

end


T = 10^5
R = 10^3
d = 0.4
@time begin
    for ii = 1:R
        fi_var_vals_re(T, 0.4)
    end
end

@time begin
    for ii = 1:R
        fi_var_vals(T, 0.4)
    end
end

fi_var_vals_re(10, d) ≈ fi_var_vals(10, d)


N = 4
coefs = fi_var_vals(N, d)

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

Gam = my_toeplitz(coefs)


d = 400;theta = -0.5+(exp(d))/(exp(d))