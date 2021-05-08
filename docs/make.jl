using Documenter
using SimpleANOVA

DocMeta.setdocmeta!(
    SimpleANOVA,
    :DocTestSetup,
    :(using SimpleANOVA);
    recursive=true
)

makedocs(
    sitename = "SimpleANOVA.jl",
    pages = [
        "SimpleANOVA" => "index.md",
        "API" => "api.md"
    ],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    modules = [SimpleANOVA],
    strict = true,
    checkdocs = :exports
)

deploydocs(repo="github.com/BioTurboNick/SimpleANOVA.jl")
