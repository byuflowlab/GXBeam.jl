module GEBT

using Arpack
using FLOWMath
using LinearAlgebra
using LinearMaps
using NLsolve
using SparseArrays
using StaticArrays
using OffsetArrays
using WriteVTK

export curve_length, discretize_beam
export System, Assembly, PrescribedConditions, DistributedLoads, TimeFunction
export static_analysis, steady_state_analysis, eigenvalue_analysis, time_domain_analysis
export static_analysis!, steady_state_analysis!, eigenvalue_analysis!, time_domain_analysis!
export AssemblyState, left_eigenvectors, wiener_milenkovic, write_vtk

# Constant used for scaling forces/moments
# this is needed because sparse arrays don't have scaling methods defined yet
const FORCE_SCALING = 1e3

# Global constants used for distributed-load integration
N_GAUSS = 6
X_GAUSS = SVector(
    -0.9324695142031521,
    -0.6612093864662645,
    -0.2386191860831969,
     0.2386191860831969,
     0.6612093864662645,
     0.9324695142031521
     )
W_GAUSS = SVector(
    0.17132449237917044,
    0.36076157304813855,
    0.46791393457269115,
    0.46791393457269109,
    0.36076157304813855,
    0.17132449237917044
    )
W1 = ((1 .- X_GAUSS) ./ 4 .* W_GAUSS)
W2 = ((1 .+ X_GAUSS) ./ 4 .* W_GAUSS)

include("math.jl")
include("element.jl")
include("point.jl")
include("loads.jl")
include("assembly.jl")
include("system.jl")
include("analyses.jl")
include("postprocess.jl")

end
