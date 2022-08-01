using GXBeam, DifferentialEquations, Elliptic, LinearAlgebra, ForwardDiff, Random, Test

# input/output tests
include("io.jl")

# jacobian tests
include("jacobians.jl")

# cross-section tests
include("section.jl")

# examples
include("examples/cantilever.jl")
include("examples/overdetermined.jl")
include("examples/tipforce.jl")
include("examples/tipmoment.jl")
include("examples/curved.jl")
include("examples/rotating.jl")
# include("examples/wind-turbine-blade.jl")
include("examples/static-joined-wing.jl")
# include("examples/dynamic-joined-wing.jl")

# interfaces
# include("interfaces/diffeq.jl")
include("interfaces/forwarddiff.jl")

# issues
include("issues/zeros.jl")
include("issues/pointmass.jl")




