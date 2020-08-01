using Documenter, GEBT

makedocs(;
    modules = [GEBT],
    pages = [
        "Home" => "index.md",
        "Examples" => "examples.md",
        "Library" => "library.md"
    ],
    sitename = "GEBT.jl",
    authors = "Taylor McDonnell <taylormcd@byu.edu>",
)

deploydocs(
    deps = Deps.pip("python-markdown-math"),
    repo = "github.com/byuflowlab/GEBT.jl.git",
)
