using Random, Plots, FFTW

function arfimagen(T::Int, μ::Real, AR::Array, d::Real, MA::Array; σ=1)
    p = length(AR)
    q = length(MA)

    if q == 0
        u = rand(Normal(0, σ), T + p, 1)
    else
        uM = zeros(T, q + 1)
        u = rand(Normal(0, σ), T + p + q, 1)
        MAb = [1; MA]

        for j in 1:q+1
            for i in 1:j
                uM[:, j] .+= MAb[i] * u[q+1-i+1:T+q-i+1]
            end
        end
        u = sum(uM, dims=2)
    end

    fiu = fracdiff(u, -d)  #fracdiff with multiple dispatch

    if p == 0
        y = fiu
    else
        yb = zeros(T + p)
        for j in p+1:T+p
            for i in 1:p
                yb[j] += AR[i] * yb[j-i]
            end
            yb[j] = μ + yb[j] + fiu[j-p]
        end
        y = yb[p+1:T+p]
    end
    return y
end

function arfigen(T::Int, μ::Real, AR::Array, d::Real; σ=1)
    p = length(AR)

    u = rand(Normal(0, σ), T + p, 1)

    fiu = fracdiff(u, -d)  #fracdiff with multiple dispatch

    yb = zeros(T + p)
    for ii in p+1:T+p
        for jj in 1:p
            yb[ii] += AR[jj] * yb[ii-jj]
        end
        yb[ii] = μ + yb[ii] + fiu[ii-p]
    end
    y = yb[p+1:T+p]

    return y
end




gph_est((arfigen(1000, 0, [], 0.4)))