module GXBeam

using ArnoldiMethod
using ForwardDiff
using LinearAlgebra
using LinearMaps
using NLsolve
import Roots
using SparseArrays
using StaticArrays
using WriteVTK
using SciMLBase
using UnPack

export AbstractSystem, StaticSystem, DynamicSystem, ExpandedSystem 
export Assembly, PrescribedConditions, DistributedLoads, PointMass
export AssemblyState, PointState, ElementState

export reset_state!, copy_state!, set_state!
export set_linear_deflection!, set_angular_deflection!
export set_external_forces!, set_external_moments!
export set_linear_velocity!, set_angular_velocity!
export set_internal_forces!, set_internal_moments!
export set_start_forces!, set_start_moments!
export set_end_forces!, set_end_moments!
export set_point_linear_velocity!, set_point_angular_velocity!
export set_element_linear_velocity!, set_element_angular_velocity!

export static_analysis, static_analysis!
export steady_state_analysis, steady_state_analysis!
export linearize!, solve_eigensystem, left_eigenvectors, correlate_eigenmodes
export eigenvalue_analysis, eigenvalue_analysis!
export initial_condition_analysis, initial_condition_analysis!
export time_domain_analysis, time_domain_analysis!
export AssemblyState
export PointState, extract_point_state, extract_point_states, extract_point_states!
export ElementState, extract_element_state, extract_element_states, extract_element_states!

export curve_length, discretize_beam
export combine_loads, combine_masses

export wiener_milenkovic
export rotate, rotate!
export translate, translate!
export cross_section_velocities, cross_section_velocities!
export deform_cross_section, deform_cross_section!

export write_vtk

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

include("interfaces/diffeq.jl")
include("interfaces/forwarddiff.jl")

end
