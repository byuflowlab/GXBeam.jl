using Documenter, Literate, GXBeam

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
        "Section Properties and Strain Recovery" => joinpath("examples", "section.md"),
        "Sensitivity Analysis" => joinpath("examples", "sensitivities.md"),
        "DifferentialEquations" => joinpath("examples", "diffeq.md"),
        "Examples" => [
            joinpath("examples", "cantilever.md"),
            joinpath("examples", "overdetermined.md"),
            joinpath("examples", "tipforce.md"),
            joinpath("examples", "tipmoment.md"),
            joinpath("examples", "curved.md"),
            joinpath("examples", "rotating.md"),
            joinpath("examples", "excited.md"),
            joinpath("examples", "wind-turbine-blade.md"),
            joinpath("examples", "static-joined-wing.md"),
            joinpath("examples", "dynamic-joined-wing.md"),
            joinpath("examples", "vertical-axis-wind-turbine.md"),
        ],
        "API Reference" => [
            joinpath("reference", "public.md"),
            joinpath("reference", "private.md"),
        ],
    ],
    format = Documenter.HTMLWriter.HTML(
        size_threshold_warn = 307200,  # 300 KiB
        size_threshold = 409600 # 400 KiB
    ),
    sitename = "GXBeam.jl",
    authors = "Taylor McDonnell <taylormcd@byu.edu>",
)

deploydocs(;
    repo = "github.com/byuflowlab/GXBeam.jl.git",
    devbranch = "master"
)
