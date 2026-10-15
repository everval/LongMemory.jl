```@meta
CurrentModule = LongMemory
```

# Data Available

The package includes the following data sets:

- Nile River Minima (622-1284). Available as `NileRiverMin.csv`. The data set contains two columns: `Year`and `NileMin`.

- Northern Hemisphere monthly temperature anomalies. Available as `NorthernHemTemp.csv` in a long time series format starting from January 1850.

## Loading the Data

The data are stored in `CSV` files and are returned as named tuples of vectors, so no additional packages are needed.

You can load the data sets with the following commands:

```julia
using LongMemory
NileMin = NileData()          # NileMin.Year, NileMin.NileMin
NorthernHemTemp = NHTempData() # NorthernHemTemp.NHTemp
```

The named tuples are [Tables.jl](https://tables.juliadata.org/stable/) tables, so they can be converted to a `DataFrame` with [DataFrames.jl](https://dataframes.juliadata.org/stable/):

```julia
using DataFrames
DataFrame(NileData())
```

Moreover, you can generate plots of the data and associated long memory processes with the following commands, which require [Plots.jl](https://docs.juliaplots.org/stable/):

```julia
using Plots
NileDataPlot()
NHTempDataPlot()
```

# Sources

- Northern Hemisphere Temperature Anomaly: [GISTEMP](https://data.giss.nasa.gov/gistemp/). [Morice et al. (2021)](https://agupubs.onlinelibrary.wiley.com/action/showCitFormats?doi=10.1029%2F2019JD032361)

Documentation for [LongMemory.jl](https://github.com/everval/LongMemory.jl).
