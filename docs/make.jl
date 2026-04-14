using LongMemory
using Documenter

DocMeta.setdocmeta!(LongMemory, :DocTestSetup, :(using LongMemory); recursive=true)

 Changelog.generate(
    Changelog.Documenter(),                 # output type
    joinpath(@__DIR__, "../CHANGELOG.md"),  # input file
    joinpath(@__DIR__, "src/CHANGELOG.md"); # output file
    repo = "everval/LongMemory.jl",        # default repository for links
)

makedocs(
    #modules=[LongMemory, LongMemory.GeneratingFunctions],
    authors="J. Eduardo Vera-Valdés",
    repo="https://github.com/everval/LongMemory.jl/blob/{commit}{path}#{line}",
    sitename="LongMemory.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://everval.github.io/LongMemory.jl",
        repolink="https://everval.github.io/LongMemory.jl",
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
        "Structural Changes" => "structchange.md",
        "Data Available" => "data.md",
        "Illustrative Examples" => "examples.md",
        "List of Functions" => "functionlist.md",
    ],
)

deploydocs(;
    repo="github.com/everval/LongMemory.jl",
    devbranch="master",
)
