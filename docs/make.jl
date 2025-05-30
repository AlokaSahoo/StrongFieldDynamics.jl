using StrongFieldDynamics
using Documenter

DocMeta.setdocmeta!(StrongFieldDynamics, :DocTestSetup, :(using StrongFieldDynamics); recursive=true)

makedocs(;
    modules=[StrongFieldDynamics],
    authors="Aloka Kumar Sahoo <aloka_s@ph.iitr.ac.in>",
    sitename="StrongFieldDynamics.jl",
    format=Documenter.HTML(;
        canonical="https://AlokaSahoo.github.io/StrongFieldDynamics.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/AlokaSahoo/StrongFieldDynamics.jl",
    devbranch="main",
)