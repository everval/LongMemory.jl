using Random, Plots, FFTW

function dgp_arfima(c, ARo, MAo, T, SD, d, F)
    p = length(ARo)
    q = length(MAo)
    
    # Simulation of the process
    if SD == 0
        SD = 1
    end
    
    if q == 0
        u = randn(T + p)
    else
        uM = zeros(T, q + 1)
        u = randn(T + q)
        MAb = [1; MAo]
        
        for j in 1:q + 1
            for i in 1:j
                uM[:, j] .+= MAb[i] * u[q + 1 - i + 1:T + q - i + 1]
            end
        end
        u = sum(uM, dims=2)
    end
    
    if d == 0
        fu = u
    else
        fu = fracdiff(u, -d)  # You should define the fastfrac function separately
    end
    
    if p == 0
        y = fu
    else
        yb = zeros(T + p)
        for j in p + 1:T + p
            for i in 1:p
                yb[j] += ARo[i] * yb[j - i]
            end
            yb[j] = c + yb[j] + fu[j - p]
        end
        y = yb[p + 1:T + p]
    end
    return y
end

gph_est((dgp_arfima(0,[],[],10000,1,-0.4)))