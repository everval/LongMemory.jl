```@meta
CurrentModule = LongMemory
```

# Data Available

The package includes the following data sets:

- Nile River Minima (622-1284). Available as `NileRiverMin.csv`. The data set contains two columns: `Year`and `NileMin`.

- Northern Hemisphere monthly temperature anomalies. Available as `NorthernHemTemp.csv` in a long time series format starting from January 1850.


## Loading the Data

The data are stored in `CSV` files.

Hence, loading the data requires the [CSV.jl](https://csv.juliadata.org/stable/) and [DataFrames.jl](https://dataframes.juliadata.org/stable/) packages.

You can load the data sets with the following commands:

```julia
using LongMemory
using CSV, DataFrames
NileMin = CSV.read(joinpath(artifact"Data","NileRiverMin.csv"),DataFrame)
NorthernHemTemp = CSV.read(joinpath(artifact"Data","NorthernHemTemp.csv"),DataFrame)
```

# Sources

- Northern Hemisphere Temperature Anomaly: [GISTEMP](https://data.giss.nasa.gov/gistemp/). [Morice et al. (2021)](https://agupubs.onlinelibrary.wiley.com/action/showCitFormats?doi=10.1029%2F2019JD032361)

Documentation for [LongMemory.jl](https://github.com/everval/LongMemory.jl).
