using FiniteSizeScaling
using Documenter

DocMeta.setdocmeta!(FiniteSizeScaling, :DocTestSetup, :(using FiniteSizeScaling); recursive=true)

makedocs(;
    modules=[FiniteSizeScaling],
    authors="Owen Bradley",
    repo="https://github.com/owenpb/FiniteSizeScaling.jl/blob/{commit}{path}#{line}",
    sitename="FiniteSizeScaling.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://owenpb.github.io/FiniteSizeScaling.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/owenpb/FiniteSizeScaling.jl",
    devbranch="main",
)
