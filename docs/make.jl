using GeneralizedSDistributions
using Documenter

makedocs(;
    modules=[GeneralizedSDistributions],
    authors="Alex Knudson <aknudson@nevada.unr.edu>",
    repo="https://github.com/adknudson/GeneralizedSDistributions.jl/blob/{commit}{path}#L{line}",
    sitename="GeneralizedSDistributions.jl",
    doctest=false,
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://adknudson.github.io/GeneralizedSDistributions.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/adknudson/GeneralizedSDistributions.jl",
    devbranch="main",
    versions = ["stable" => "v^", "v#.#", "dev" => "main"]
)
