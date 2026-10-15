"""
# LongMemoryPlotsExt

Package extension with the plotting functions of LongMemory.jl. It is loaded automatically when both LongMemory and Plots are loaded, e.g. with `using LongMemory, Plots`.

## Author

[J. Eduardo Vera-Valdés](https://everval.github.io/)

"""
module LongMemoryPlotsExt

using LongMemory, Plots

import LongMemory: autocorrelation_plot, log_variance_plot, rescaled_range_plot, dfa_plot, periodogram_plot, fi_forecast_plot, csa_forecast_plot, har_forecast_plot, arfima_forecast_plot, NileDataPlot, NHTempDataPlot, LMPlot


## Plots of the ClassicEstimators module

function autocorrelation_plot(x::AbstractArray, k::Int=30)
    acv = autocovariance(x, k)
    acf = acv ./ acv[1]

    p1 = plot(0:k-1, acf, line=:stem, xlabel="Lags", ylabel="Autocorrelation function", legend=false, marker=:circle)

    return p1
end


function log_variance_plot(x::AbstractArray; m::Int=100, slope::Bool=false, slope2::Bool=false)
    T = length(x)
    if m >= T
        error("m must be less than T")
    end
    k = floor(Int, T / 2) - 1
    step = round(Int, T / m)

    nt = collect(1:step:k)
    cr = length(nt)

    if cr - 1 < 3
        error("At least 3 block sizes are needed for the regression; increase the length of the series or m")
    end

    Y = zeros(cr - 1, 1)
    for ii = 2:cr
        Y[ii-1, 1] = sstdk(x, nt[ii] )
    end

    X = [ones(cr-1,1) log.(nt[2:end])]
    Z = log.(Y)
    beta = X \ Z

    p1 = plot(X[:,2], Z, line=:scatter, label="", xlabel="log-sampling", ylabel="log-variance")
    if slope == true
        oldylims = ylims(p1)
        plot!(X[:,2], X * beta, line=:dash, label=string("Slope = ", round(beta[2], digits=4)), linewidth=3)
        ylims!(p1, oldylims)
    end
    if slope2 == true
        oldylims = ylims(p1)
        plot!(X[:,2], X *  [beta[1]; -1], line=:dashdot, label="Slope = -1", linewidth=3)
        ylims!(p1, oldylims)
    end
    return p1

end


function rescaled_range_plot(x::AbstractArray; k::Int = 100, slope = false, slope2 = false)
    T = length(x)

    if k >= T
        error("k must be less than T")
    end

    RS = rescaled_range(x; k=k)

    step = round(Int, T / k)

    position = collect(1:step:T)
    popfirst!(position)

    Y = log.(RS)

    X = [ones(length(Y), 1) log.(position)]

    beta = X \ Y

    p1 = plot(X[:,2], Y, line=:scatter, label="", xlabel="log-sampling", ylabel="log-rescaled range")
    if slope == true
        oldylims = ylims(p1)
        plot!(X[:,2], X * beta, line=:dash, label=string("Slope = ", round(beta[2], digits=4)), linewidth=3)
        ylims!(p1, oldylims)
    end
    if slope2 == true
        oldylims = ylims(p1)
        plot!(X[:,2], X *  [beta[1]; 0.5], line=:dashdot, label="Slope = 0.500", linewidth=3)
        ylims!(p1, oldylims)
    end
    return p1

end


function dfa_plot(x::AbstractArray; k::Int = 20, slope::Bool = false, slope2::Bool = false)
    nt, F = dfa_fluctuation(x; k=k)

    X = [ones(length(nt), 1) log.(nt)]
    Y = log.(F)

    beta = X \ Y

    p1 = plot(X[:,2], Y, line=:scatter, label="", xlabel="log-box size", ylabel="log-fluctuation")
    if slope == true
        oldylims = ylims(p1)
        plot!(X[:,2], X * beta, line=:dash, label=string("Slope = ", round(beta[2], digits=4)), linewidth=3)
        ylims!(p1, oldylims)
    end
    if slope2 == true
        oldylims = ylims(p1)
        plot!(X[:,2], X * [beta[1]; 0.5], line=:dashdot, label="Slope = 0.500", linewidth=3)
        ylims!(p1, oldylims)
    end
    return p1
end



## Plots of the LogPeriodEstimators module

function periodogram_plot(x::AbstractArray; slope::Bool=false)
    T = length(x)

    I_wc, wc = periodogram(x)

    I_w = I_wc[2:end]
    w = wc[2:end]

    if slope == false
        p1 = plot(w, I_w, xlabel="log-frequency", ylabel="log-periodogram", label="", xaxis=:log, yaxis=:log, line=:scatter)
    end

    if slope == true
        p1 = plot(log.(w), log.(I_w), xlabel="log-frequency", ylabel="log-periodogram", label="", line=:scatter)
        oldylims = ylims(p1)
        X = [ones(length(w)) log.(w)]
        beta = X \ log.(I_w)
        plot!(log.(w), X*beta, line=:dash, label=string("Slope = ", round(beta[2],digits=4)), linewidth = 3)
        ylims!(p1, oldylims)
    end

    return p1
end



## Plots of the Forecasters module

function fi_forecast_plot(x::AbstractArray, h::Int, d::Real)
    x = x .- smean(x)
    T = length(x)

    xfor = fi_forecast(x, h, d)

    p1 = plot(xfor[1:T], xlabel="Time", ylabel="Time series", label="Original time series", line=:solid)
    plot!(T+1:T+h, xfor[T+1:T+h, 1], label="Forecast", color=:red, line=:dash)

    return p1
end


function fi_forecast_plot(x::AbstractArray, h::Int, d::Real, σ::Real)
    x = x .- smean(x)
    T = length(x)

    xfor = fi_forecast(x, h, d, σ)

    p1 = plot(xfor[1:T, 1], xlabel="Time", ylabel="Time series", label="Original time series", line=:solid)
    plot!(T+1:T+h, xfor[T+1:T+h, 1], label="Forecast", color=:red, line=:dash)
    plot!(T+1:T+h, xfor[T+1:T+h, 2], fillrange=xfor[T+1:T+h, 3], label="Confidence band", fillalpha=0.25, color=:black)
    plot!(T+1:T+h, xfor[T+1:T+h, 3], label="", color=:black, line=:solid)

    return p1
end


function csa_forecast_plot(x::AbstractArray, h::Int, p::Real, q::Real)
    x = x .- smean(x)
    T = length(x)

    xfor = csa_forecast(x, h, p, q)

    p1 = plot(xfor[1:T], xlabel="Time", ylabel="Time series", label="Original time series", line=:solid)
    plot!(T+1:T+h, xfor[T+1:T+h, 1], label="Forecast", color=:red, line=:dash)

    return p1
end


function csa_forecast_plot(x::AbstractArray, h::Int, p::Real, q::Real, σ::Real)
    x = x .- smean(x)
    T = length(x)

    xfor = csa_forecast(x, h, p, q, σ)

    p1 = plot(xfor[1:T, 1], xlabel="Time", ylabel="Time series", label="Original time series", line=:solid)
    plot!(T+1:T+h, xfor[T+1:T+h, 1], label="Forecast", color=:red, line=:dash)
    plot!(T+1:T+h, xfor[T+1:T+h, 2], fillrange=xfor[T+1:T+h, 3], label="Confidence band", fillalpha=0.25, color=:black)
    plot!(T+1:T+h, xfor[T+1:T+h, 3], label="", color=:black, line=:solid)

    return p1
end


function har_forecast_plot(x::AbstractArray, h::Int; flag::Bool=true, m::AbstractArray=[1, 5, 22])
    x = x .- smean(x)
    T = length(x)

    xfor = har_forecast(x, h; flag=flag, m=m)

    p1 = plot(xfor[1:T], xlabel="Time", ylabel="Time series", label="Original time series", line=:solid)
    plot!(T+1:T+h, xfor[T+1:T+h, 1], label="Forecast", color=:red, line=:dash)

    if flag == true
        plot!(T+1:T+h, xfor[T+1:T+h, 2], fillrange=xfor[T+1:T+h, 3], label="Confidence band", fillalpha=0.25, color=:black)
        plot!(T+1:T+h, xfor[T+1:T+h, 3], label="", color=:black, line=:solid)
    end

    return p1
end


function arfima_forecast_plot(x::AbstractArray, h::Int, d::Real, AR::AbstractArray, MA::AbstractArray)
    T = length(x)

    xfor = arfima_forecast(x, h, d, AR, MA)

    p1 = plot(xfor[1:T], xlabel="Time", ylabel="Time series", label="Original time series", line=:solid)
    plot!(T+1:T+h, xfor[T+1:T+h, 1], label="Forecast", color=:red, line=:dash)

    return p1
end


function arfima_forecast_plot(x::AbstractArray, h::Int, d::Real, AR::AbstractArray, MA::AbstractArray, σ::Real)
    T = length(x)

    xfor = arfima_forecast(x, h, d, AR, MA, σ)

    p1 = plot(xfor[1:T, 1], xlabel="Time", ylabel="Time series", label="Original time series", line=:solid)
    plot!(T+1:T+h, xfor[T+1:T+h, 1], label="Forecast", color=:red, line=:dash)
    plot!(T+1:T+h, xfor[T+1:T+h, 2], fillrange=xfor[T+1:T+h, 3], label="Confidence band", fillalpha=0.25, color=:black)
    plot!(T+1:T+h, xfor[T+1:T+h, 3], label="", color=:black, line=:solid)

    return p1
end



## Plots of the DataExamples module

function NileDataPlot()
    data = NileData()
    l = @layout [a b; c d]

    p1 = plot(data.Year, data.NileMin, ylabel="Nile River Min.", xlabel="Year", legend=false)
    p2 = autocorrelation_plot(data.NileMin,50)
    ylabel!(p2, "Autocorrelation")
    p3 = periodogram_plot(data.NileMin)
    p4 = log_variance_plot(data.NileMin;m=200)

    pp = plot(p1, p2, p3, p4, layout=l)

    return pp
end


function NHTempDataPlot()
    data = NHTempData()
    l = @layout [a b; c d]

    p1 = plot(data.NHTemp, ylabel="North Hem. Temp.", xlabel="Months", legend=false)
    p2 = autocorrelation_plot(data.NHTemp,50)
    ylabel!(p2, "Autocorrelation") # This is how you add a
    p3 = periodogram_plot(data.NHTemp)
    p4 = log_variance_plot(data.NHTemp;m=500)

    pp = plot(p1, p2, p3, p4, layout=l)

    return pp
end


function LMPlot(x::AbstractArray;lags::Int=50,name::String="Time series")
    l = @layout [a b; c d]

    p1 = plot(x, ylabel=name, xlabel="Time", legend=false)
    p2 = autocorrelation_plot(x,lags)
    ylabel!(p2, "Autocorrelation")
    p3 = periodogram_plot(x)
    p4 = log_variance_plot(x;m=floor(Int,length(x)/2)-1)

    pp = plot(p1, p2, p3, p4, layout=l)

    return pp
end


end # module LongMemoryPlotsExt
