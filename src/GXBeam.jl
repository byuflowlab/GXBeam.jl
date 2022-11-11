module GXBeam

using ArnoldiMethod
using ForwardDiff
using ChainRulesCore
using FillArrays
using FLOWMath
using ImplicitAD
using LinearAlgebra
using LinearMaps
using NLsolve
import Roots
using ReverseDiff
using SparseArrays
using StaticArrays
using SciMLBase
using UnPack
using WriteVTK

export AbstractSystem, StaticSystem, DynamicSystem, ExpandedSystem 

export static_analysis, static_analysis!
export steady_state_analysis, steady_state_analysis!
export linearize!, solve_eigensystem, left_eigenvectors, correlate_eigenmodes
export eigenvalue_analysis, eigenvalue_analysis!
export initial_condition_analysis, initial_condition_analysis!
export time_domain_analysis, time_domain_analysis!

export reset_state!, copy_state!, set_state!

export set_linear_displacement!, set_angular_displacement!
export set_external_forces!, set_external_moments!
export set_linear_velocity!, set_angular_velocity!
export set_internal_forces!, set_internal_moments!
export set_start_forces!, set_start_moments!
export set_end_forces!, set_end_moments!
export set_point_linear_velocity!, set_point_angular_velocity!
export set_element_linear_velocity!, set_element_angular_velocity!

export Assembly
export curve_length, discretize_beam

export PrescribedConditions, DistributedLoads, PointMass
export combine_loads, combine_masses

export AssemblyState
export PointState, extract_point_state, extract_point_states, extract_point_states!
export ElementState, extract_element_state, extract_element_states, extract_element_states!
export body_accelerations

export wiener_milenkovic
export rotate, rotate!
export translate, translate!
export deform_cross_section, deform_cross_section!

export write_vtk

# from section
export Material, Node, MeshElement
export initialize_cache, compliance_matrix, mass_matrix, plotmesh

# from afmesh
export Layer
export afmesh

const GAUSS_NODES = SVector(-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526)
const GAUSS_WEIGHTS = SVector(0.34785484513745385, 0.6521451548625462, 0.6521451548625462, 0.34785484513745385)

# common math functions
include("math.jl")

# section properties
include("section.jl")

# airfoil meshing
include("afmesh.jl")

# assembly creation
include("assembly.jl")

# prescribed conditions, distributed loads, and point masses
include("loads.jl")

# system storage and pointers 
include("system.jl")

# state variable input and output methods
include("input.jl")
include("output.jl")

# point residuals and jacobians
include("point.jl")

# element residuals and jacobians
include("element.jl")

# system analyses
include("analyses.jl")

# DifferentialEquations Interface
include("interfaces/diffeq.jl")

# ReverseDiff Overloads
include("interfaces/reversediff.jl")

# WriteVTK Interface
include("interfaces/writevtk.jl")

end
