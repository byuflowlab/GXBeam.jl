module GXBeam

using ArnoldiMethod
using FLOWMath
using LinearAlgebra
using LinearMaps
using NLsolve
using SparseArrays
using StaticArrays
using WriteVTK

export curve_length, discretize_beam
export System, Assembly, PrescribedConditions, DistributedLoads
export reset_state!
export static_analysis, steady_state_analysis, eigenvalue_analysis, time_domain_analysis
export static_analysis!, steady_state_analysis!, eigenvalue_analysis!, time_domain_analysis!
export AssemblyState, left_eigenvectors, wiener_milenkovic, write_vtk, correlate_eigenmodes

# Constant used for scaling forces/moments
# this is needed because sparse arrays don't seem to be able to handle
# mismatched units very well
const FORCE_SCALING = nextpow(2, 1e3)

const GAUSS_NODES = SVector(-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526)
const GAUSS_WEIGHTS = SVector(0.34785484513745385, 0.6521451548625462, 0.6521451548625462, 0.34785484513745385)

include("math.jl")
include("element.jl")
include("point.jl")
include("loads.jl")
include("assembly.jl")
include("system.jl")
include("analyses.jl")
include("postprocess.jl")

end
