using GSDistributions
using Documenter

makedocs(;
    modules=[GSDistributions],
    authors="Alex Knudson <aknudson@nevada.unr.edu>",
    repo="https://github.com/adknudson/GSDistributions.jl/blob/{commit}{path}#L{line}",
    sitename="GSDistributions.jl",
    doctest=false,
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://adknudson.github.io/GSDistributions.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/adknudson/GSDistributions.jl",
)
