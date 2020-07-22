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
export left_eigenvectors
export AssemblyState, write_vtk

# Constant used for scaling forces/moments
# this is needed because sparse arrays don't have scaling methods defined yet
const FORCE_SCALING = 1e3

# Global constants used for distributed-load integration
const N_GAUSS = 6
const X_GAUSS = @SVector [-cos((2*i-1)/(2*N_GAUSS)*pi) for i = 1:N_GAUSS]
const W_GAUSS = @SVector [pi/N_GAUSS for i = 1:N_GAUSS]
const W1 = (1 .- X_GAUSS) ./ 4 .* W_GAUSS
const W2 = (1 .+ X_GAUSS) ./ 4 .* W_GAUSS

include("math.jl")
include("element.jl")
include("point.jl")
include("loads.jl")
include("assembly.jl")
include("system.jl")
include("analyses.jl")
include("postprocess.jl")

end
