"""
# DataExamples

This package contains a collection of data sets that are used in the LongMemory.jl package.

## Author

[J. Eduardo Vera-Valdés](https://everval.github.io/)
"""
module DataExamples

import ..plots_required

using DelimitedFiles, Artifacts

export NileData, NHTempData, NileDataPlot, NHTempDataPlot, LMPlot

# Plotting functions. Their methods are defined in the package extension LongMemoryPlotsExt, which is loaded with `using Plots`.
function NileDataPlot end
function NHTempDataPlot end
function LMPlot end

NileDataPlot(args...; kwargs...) = plots_required(NileDataPlot, args)
NHTempDataPlot(args...; kwargs...) = plots_required(NHTempDataPlot, args)
LMPlot(args...; kwargs...) = plots_required(LMPlot, args)


"""
    NileData()

Returns the annual minimum level of the Nile River at Aswan, Egypt, from 622 to 1284.

# Output
- `data::NamedTuple`: The Nile River data, with columns `Year` and `NileMin`.

# Notes
The data is returned as a named tuple of vectors, e.g. `NileData().NileMin`. It is a Tables.jl table, so `DataFrame(NileData())` returns a `DataFrame` if DataFrames.jl is loaded.

# References
- [Hurst, H. E. (1951). Long-term storage capacity of reservoirs. Transactions of the American Society of Civil Engineers, 116, 770-799.](https://ascelibrary.org/doi/abs/10.1061/TACEAT.0006518)

# Examples
```julia
julia> NileData()
```
"""
function NileData()
    data = readdlm(joinpath(artifact"Data", "NileRiverMin.csv"), ',', Float64; header=true)[1]
    return (Year = Int.(data[:, 1]), NileMin = data[:, 2])
end


"""
    NHTempData()

Returns the monthly mean temperature of the Northern Hemisphere from January 1850 to September 2023.

# Output
- `data::NamedTuple`: The Northern Hemisphere temperature data, with column `NHTemp`.

# Notes
The data is returned as a named tuple of vectors, e.g. `NHTempData().NHTemp`. It is a Tables.jl table, so `DataFrame(NHTempData())` returns a `DataFrame` if DataFrames.jl is loaded.

# Examples
```julia
julia> NHTempData()
```

# References
- [Morice et al. (2021). An Updated Assessment of Near-Surface Temperature Change From 1850: The HadCRUT5 Data Set. Journal of Geophysical Research: Atmospheres, 126, e2019JD032361. https://doi.org/10.1029/2019JD032361](https://doi.org/10.1029/2019JD032361)
"""
function NHTempData()
    data = readdlm(joinpath(artifact"Data", "NorthernHemTemp.csv"), ',', Float64; header=true)[1]
    return (NHTemp = data[:, 1],)
end


"""
    NileDataPlot()

Returns a plot of the yearly Nile River minimum data from 622 to 1284.

Requires Plots.jl to be loaded, e.g. with `using LongMemory, Plots`.

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
NileDataPlot()


"""
    NHTempDataPlot()

Returns a plot of the monthly Northern Hemisphere temperature data from January 1850 to September 2023.

Requires Plots.jl to be loaded, e.g. with `using LongMemory, Plots`.

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
NHTempDataPlot()

"""
    LMPlot(x::AbstractArray;lags::Int=50,name::String="Time series")

Returns a plot of the time series `x`.

Requires Plots.jl to be loaded, e.g. with `using LongMemory, Plots`.

# Arguments
- `x::Array`: The time series.

# Optional arguments
- `lags::Int`: The number of lags to use in the autocorrelation function.
- `name::String`: The name of the time series to use in the plot.

# Output
- `pp::Plot`: Four plots are returned in a 2x2 grid: 
    - The first plot is the time series of `x`. 
    - The second plot is the autocorrelation function of `x` using the `autocorrelation_plot` function. 
    - The third plot is the periodogram of `x` using the `periodogram_plot` function. 
    - The fourth plot is the variance plot of `x` using the `variance_plot` function.

# Examples
```julia
julia> LMPlot(randn(100))
```
"""
LMPlot(x::AbstractArray)


end # module DataExamples
