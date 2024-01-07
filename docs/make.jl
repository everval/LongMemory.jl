using LongMemory
using Documenter

DocMeta.setdocmeta!(LongMemory, :DocTestSetup, :(using LongMemory); recursive=true)

makedocs(
    #modules=[LongMemory, LongMemory.GeneratingFunctions],
    authors="J. Eduardo Vera-Valdés",
    repo="https://github.com/everval/LongMemory.jl/blob/{commit}{path}#{line}",
    sitename="LongMemory.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://everval.github.io/LongMemory.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Plotting" => "plotting.md",
        "Generating Functions" => "generating.md",
        "Log-Periodogram Estimation" => "logperiod.md",
        "Parametric Estimation" => "parametric.md",
        "Classic Estimators" => "classicest.md",
        "Forecasting" => "forecasting.md",
        "Data Available" => "data.md",
        "Illustrative Examples" => "examples.md",
        "List of Functions" => "functionlist.md",
    ],
)

deploydocs(;
    repo="github.com/everval/LongMemory.jl",
    devbranch="master",
)
