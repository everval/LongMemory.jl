```@meta
CurrentModule = LongMemory
```

# Data Available

The package includes the following data sets:
- Nile River Minima (622-1284). Available as `NileRiverMin.csv`.


The Nile River Minima contains the following columns:
- `Year`: Year of the observation.
- `NileMin`: Value of the observation.

## Loading the Data

The data are stored in `CSV` files as `DataFrame`.

Hence, loading the data requires the [CSV.jl](https://csv.juliadata.org/stable/) and [DataFrames.jl](https://dataframes.juliadata.org/stable/) packages.

You can load the data sets with the following commands:
```julia
using LongMemory
using CSV, DataFrames
NileMin = CSV.read("examples/NileRiverMin.csv",DataFrame)
```

Documentation for [LongMemory.jl](https://github.com/everval/LongMemory.jl).
