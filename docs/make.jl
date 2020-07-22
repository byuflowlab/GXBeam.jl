using Documenter, GEBT

makedocs(;
    modules = [GEBT],
    format = Documenter.HTML(),
    pages = [
        "Home" => "index.md",
        "Examples" => [
            "Linear" => "linear.md",
            "Static" => "static.md",
            "Steady State" => "steady_state.md",
            "Stability" => "stability.md",
            "Time Marching" => "time_marching.md"
        ],
        "Library" => "library.md"
    ],
    sitename = "GEBT.jl",
    authors = "Taylor McDonnell <taylormcd@byu.edu>",
)

deploydocs(
    repo="github.com/byuflowlab/GEBT.jl.git",
)
