using GXBeam, DifferentialEquations, LinearAlgebra, ForwardDiff, Random, Test

# input/output tests
include("io.jl")

# jacobian tests
include("jacobians.jl")

# cross-section tests
include("section.jl")

# examples
include(joinpath("examples", "cantilever.jl"))
include(joinpath("examples", "overdetermined.jl"))
include(joinpath("examples", "tipforce.jl"))
include(joinpath("examples", "tipmoment.jl"))
include(joinpath("examples", "curved.jl"))
include(joinpath("examples", "rotating.jl"))
include(joinpath("examples", "wind-turbine-blade.jl"))
include(joinpath("examples", "static-joined-wing.jl"))
include(joinpath("examples", "dynamic-joined-wing.jl"))

# interfaces
include(joinpath("interfaces", "diffeq.jl"))
include(joinpath("interfaces", "forwarddiff.jl"))

# issues
include(joinpath("issues", "zeros.jl"))
include(joinpath("issues", "pointmass.jl"))




