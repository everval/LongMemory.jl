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
        "Home" => "about.md",
        "Generating Functions" => "generating.md",
        "Log-Periodogram Estimation" => "logperiod.md",
        "Parameter Estimation" => "parametric.md",
        "Classic Estimation" => "classicest.md",
        "Data Available" => "data.md",
        "Index" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/everval/LongMemory.jl",
    devbranch="master",
)
