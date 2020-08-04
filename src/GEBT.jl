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
export System, Assembly, PrescribedConditions, DistributedLoads
export static_analysis, steady_state_analysis, eigenvalue_analysis, time_domain_analysis
export static_analysis!, steady_state_analysis!, eigenvalue_analysis!, time_domain_analysis!
export AssemblyState, left_eigenvectors, wiener_milenkovic, write_vtk, correlate_eigenmodes

# Constant used for scaling forces/moments
# this is needed because sparse arrays don't seem to be able to handle
# mismatched units very well
const FORCE_SCALING = 1e3

include("math.jl")
include("element.jl")
include("point.jl")
include("loads.jl")
include("assembly.jl")
include("system.jl")
include("analyses.jl")
include("postprocess.jl")

end
