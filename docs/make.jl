using StrongFieldDynamics
using Documenter
using DocumenterVitepress

DocMeta.setdocmeta!(StrongFieldDynamics, :DocTestSetup, :(using StrongFieldDynamics); recursive=true)

makedocs(;
    modules=[StrongFieldDynamics],
    authors="Aloka Kumar Sahoo <aloka_s@ph.iitr.ac.in>",
    sitename="StrongFieldDynamics.jl",
    format = DocumenterVitepress.MarkdownVitepress(
        repo="https://AlokaSahoo.github.io/StrongFieldDynamics.jl",
    ),
    pages=[
        "Home" => "index.md",
    ],
)

DocumenterVitepress.deploydocs(;
    repo = "github.com/AlokaSahoo/StrongFieldDynamics.jl",
    target = joinpath(@__DIR__, "build"),
    devbranch = "main",
    branch = "gh-pages",
    push_preview = true,
)
