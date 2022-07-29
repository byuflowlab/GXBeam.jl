using GXBeam
using LinearAlgebra
using StaticArrays
using DifferentialEquations
using Test
import Elliptic
using ForwardDiff
using Random

const RNG = MersenneTwister(1234)

# internals
include("internals/io.jl")
include("internals/rotations.jl")
include("internals/jacobians.jl")

# examples
include("examples/cantilever.jl")
include("examples/overdetermined.jl")
include("examples/tipforce.jl")
include("examples/tipmoment.jl")
include("examples/curved.jl")
include("examples/rotating.jl")
include("examples/wind-turbine-blade.jl")
include("examples/static-joined-wing.jl")
include("examples/dynamic-joined-wing.jl")

# interfaces
include("interfaces/diffeq.jl")
include("interfaces/forwarddiff.jl")

# issues
include("issues/zeros.jl")
include("issues/pointmass.jl")




include("sectiontests.jl")