using Documenter, Literate, GXBeam

const is_ci = haskey(ENV, "GITHUB_ACTIONS")

# Pre-install matplotlib
import Plots; Plots.pyplot()

# Generate examples
include("generate.jl")

# Build documentation
makedocs(;
    modules = [GXBeam],
    pages = [
        "Home" => "index.md",
        "Getting Started" => joinpath("examples", "guide.md"),
        "Examples" => [
            joinpath("examples", "cantilever.md"),
            joinpath("examples", "overdetermined.md"),
            joinpath("examples", "tipforce.md"),
            joinpath("examples", "tipmoment.md"),
            joinpath("examples", "curved.md"),
            joinpath("examples", "rotating.md"),
            joinpath("examples", "wind-turbine-blade.md"),
            joinpath("examples", "static-joined-wing.md"),
            joinpath("examples", "dynamic-joined-wing.md"),
            joinpath("examples", "vertical-axis-wind-turbine.md"),
        ],
        "Using GXBeam with DifferentialEquations.jl" => joinpath("examples", "diffeq.md"),
        "API Reference" => joinpath("reference", "reference.md"),
    ],
    sitename = "GXBeam.jl",
    authors = "Taylor McDonnell <taylormcd@byu.edu>",
)

deploydocs(
    repo = "github.com/byuflowlab/GXBeam.jl.git",
)
