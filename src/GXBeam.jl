module GXBeam

using ArnoldiMethod
using LinearAlgebra
using LinearMaps
using NLsolve
using SparseArrays
using StaticArrays
using WriteVTK
using DifferentialEquations

export curve_length, discretize_beam
export System, Assembly, PrescribedConditions, DistributedLoads
export combine_loads
export system_state, reset_state!, set_state!
export set_element_deflections!, set_element_rotations!
export set_element_forces!, set_element_moments!
export set_element_linear_momenta!, set_element_angular_momenta!
export set_point_deflections!, set_point_rotations!
export set_point_forces!, set_point_moments!
export static_analysis, static_analysis!
export steady_state_analysis, steady_state_analysis!
export eigenvalue_analysis, eigenvalue_analysis!
export initial_condition_analysis, initial_condition_analysis!
export time_domain_analysis, time_domain_analysis!
export AssemblyState
export PointState, extract_point_state, extract_point_states, extract_point_states!
export ElementState, extract_element_state, extract_element_states, extract_element_states!
export wiener_milenkovic, angular_velocities
export rotate, rotate!
export translate, translate!
export cross_section_velocities, cross_section_velocities!
export deform_cross_section, deform_cross_section!


export write_vtk
export left_eigenvectors, correlate_eigenmodes

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

include("interfaces/differentialequations.jl")

end
