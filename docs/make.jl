using Pkg
Pkg.activate(@__DIR__)                     
Pkg.develop(PackageSpec(path=joinpath(@__DIR__, "..")))
Pkg.instantiate()

using AttainableRegionTutorial, Documenter

DocMeta.setdocmeta!(AttainableRegionTutorial, :DocTestSetup, :(using AttainableRegionTutorial); recursive=true)

makedocs(;
    modules=[AttainableRegionTutorial],
    authors="Thomas Vetter <vettert85@gmail.com> and contributors",
    sitename="AttainableRegionTutorial.jl",
    format=Documenter.HTML(;
        canonical="https://github.com/thomvet/AttainableRegionTutorial.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/thomvet/AttainableRegionTutorial.jl",
    devbranch="main",
)
