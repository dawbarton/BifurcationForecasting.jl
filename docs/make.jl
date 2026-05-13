# docs/make.jl
#
# Build script for BifurcationForecasting.jl documentation using Documenter.jl.
#
# Run from the package root:
#   julia --project=docs docs/make.jl
#
# The generated HTML is written to docs/build/.

using Documenter
using BifurcationForecasting

DocMeta.setdocmeta!(
    BifurcationForecasting,
    :DocTestSetup,
    :(using BifurcationForecasting);
    recursive = true
)

makedocs(;
    modules  = [BifurcationForecasting],
    sitename = "BifurcationForecasting.jl",
    authors  = "David A. W. Barton",
    format   = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
    ),
    remotes  = nothing,
    pages = [
        "Home"     => "index.md",
        "Workflow" => "workflow.md",
        "API"      => "api.md",
    ],
    checkdocs = :exports,
    warnonly  = false,
)
