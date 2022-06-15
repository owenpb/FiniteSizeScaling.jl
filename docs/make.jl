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
        assets=String[]
    ),
    pages=[
        "Home" => "index.md",
        "Methods" => "methods.md",
        "Example: Holstein model" => "example_page.md",
        "Demo 1: One-parameter scaling" => "demo_1.md",
        "Demo 2: Two-parameter scaling" => "demo_2.md"
    ]
)

deploydocs(;
    repo="github.com/owenpb/FiniteSizeScaling.jl",
    devbranch="main"
)
