using Documenter, GXBeam

makedocs(;
    modules = [GXBeam],
    pages = [
        "Home" => "index.md",
        "Getting Started" => "guide.md",
        "Examples" => "examples.md",
        "Using GXBeam with DifferentialEquations.jl" => "diffeq.md",
        "Library" => "library.md"
    ],
    sitename = "GXBeam.jl",
    authors = "Taylor McDonnell <taylormcd@byu.edu>",
)

deploydocs(
    repo = "github.com/byuflowlab/GXBeam.jl.git",
)
