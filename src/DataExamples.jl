"""
# DataExamples

This package contains a collection of data sets that are used in the LongMemory.jl package.

## Author

[J. Eduardo Vera-Valdés](https://everval.github.io/)
"""
module DataExamples

using DataFrames, CSV, Artifacts, Plots

include("ClassicEstimators.jl")
import .ClassicEstimators: autocorrelation_plot, variance_plot

include("LogPeriodEstimators.jl")
import .LogPeriodEstimators: periodogram_plot

export NileData, NHTempData, NileDataPlot, NHTempDataPlot, LMPlot


"""
    NileData()

Returns the annual minimum level of the Nile River at Aswan, Egypt, from 622 to 1284.

# Output
- `data::DataFrame`: The Nile River data.

# References
- [Hurst, H. E. (1951). Long-term storage capacity of reservoirs. Transactions of the American Society of Civil Engineers, 116, 770-799.](https://ascelibrary.org/doi/abs/10.1061/TACEAT.0006518)

# Examples
```julia
julia> NileData()
```
"""
function NileData()
    return CSV.read(joinpath(artifact"Data", "NileRiverMin.csv"), DataFrame)
end


"""
    NHTempData()

Returns the monthly mean temperature of the Northern Hemisphere from January 1850 to September 2023.

# Output
- `data::DataFrame`: The Northern Hemisphere temperature data.

# References
- [Morice et al. (2021). An Updated Assessment of Near-Surface Temperature Change From 1850: The HadCRUT5 Data Set. Journal of Geophysical Research: Atmospheres, 126, e2019JD032361. https://doi.org/10.1029/2019JD032361](https://doi.org/10.1029/2019JD032361)
"""
function NHTempData()
    return CSV.read(joinpath(artifact"Data", "NorthernHemTemp.csv"), DataFrame)
end


"""
    NileDataPlot()

Returns a plot of the yearly Nile River minimum data from 622 to 1284.

# Output
- `pp::Plot`: The plot of the Nile River data.
- Four plots are returned in a 2x2 grid: 
    - The first plot is the time series of the Nile River data. 
    - The second plot is the autocorrelation function of the Nile River data using the `autocorrelation_plot` function. 
    - The third plot is the periodogram of the Nile River data using the `periodogram_plot` function. 
    - The fourth plot is the variance plot of the Nile River data using the `variance_plot` function.

# Notes
This function uses the `NileData()` function to load the data.

# Examples
```julia
julia> NileDataPlot()
```
"""
function NileDataPlot()
    data = NileData()
    l = @layout [a b; c d]

    p1 = plot(data.Year, data.NileMin, ylabel="Nile River Min.", xlabel="Year", legend=false)
    p2 = autocorrelation_plot(data.NileMin,50)
    ylabel!(p2, "Autocorrelation")
    p3 = periodogram_plot(data.NileMin)
    p4 = variance_plot(data.NileMin)[1]

    pp = plot(p1, p2, p3, p4, layout=l)
    display(p1)
    return pp
end


"""
    NHTempDataPlot()

Returns a plot of the monthly Northern Hemisphere temperature data from January 1850 to September 2023.

# Output
- `pp::Plot`: The plot of the Nile River data.
- Four plots are returned in a 2x2 grid: 
    - The first plot is the time series of the Northern Hemisphere temperature data. 
    - The second plot is the autocorrelation function of the Northern Hemisphere temperature data using the `autocorrelation_plot` function. 
    - The third plot is the periodogram of the Northern Hemisphere temperature data using the `periodogram_plot` function. 
    - The fourth plot is the variance plot of the Northern Hemisphere temperature data using the `variance_plot` function.

# Notes
This function uses the `NHTempData()` function to load the data.

# Examples
```julia
julia> NHTempDataPlot()
```
"""
function NHTempDataPlot()
    data = NHTempData()
    l = @layout [a b; c d]

    p1 = plot(data.NHTemp, ylabel="North Hem. Temp.", xlabel="Months", legend=false)
    p2 = autocorrelation_plot(data.NHTemp,50)
    ylabel!(p2, "Autocorrelation") # This is how you add a
    p3 = periodogram_plot(data.NHTemp)
    p4 = variance_plot(data.NHTemp)[1]

    pp = plot(p1, p2, p3, p4, layout=l)
    display(p1)
    return pp
end

"""
    LMPlot(x::Array;lags::Int=50)

Returns a plot of the time series `x`.

# Arguments
- `x::Array`: The time series.

# Optional arguments
- `lags::Int`: The number of lags to use in the autocorrelation function.

# Output
- `pp::Plot`: The plot of the time series `x`.
- Four plots are returned in a 2x2 grid: 
    - The first plot is the time series of `x`. 
    - The second plot is the autocorrelation function of `x` using the `autocorrelation_plot` function. 
    - The third plot is the periodogram of `x` using the `periodogram_plot` function. 
    - The fourth plot is the variance plot of `x` using the `variance_plot` function.

# Examples
```julia
julia> LMPlot(randn(100))
```
"""
function LMPlot(x::Array;lags::Int=50,name::String="Time series")
    l = @layout [a b; c d]

    p1 = plot(x, ylabel=name, xlabel="Time", legend=false)
    p2 = autocorrelation_plot(x,lags)
    ylabel!(p2, "Autocorrelation")
    p3 = periodogram_plot(x)
    p4 = variance_plot(x)[1]

    pp = plot(p1, p2, p3, p4, layout=l)
    display(p1)
    return pp
end


end # module DataExamples
