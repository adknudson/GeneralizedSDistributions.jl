using Documenter, GeneralizedSDistributions

makedocs(;
    modules=[GeneralizedSDistributions],
    authors="Alex Knudson <aknudson@nevada.unr.edu>",
    sitename="GeneralizedSDistributions.jl",
    doctest=false,
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/adknudson/GeneralizedSDistributions.jl",
    devbranch="main",
    versions = ["stable" => "v^", "v#.#", "dev" => "main"]
)
