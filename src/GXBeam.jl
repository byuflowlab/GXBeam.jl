module GXBeam

using ArnoldiMethod
using ForwardDiff
using FillArrays
using LinearAlgebra
using LinearMaps
using NLsolve
import Roots
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

export set_body_linear_displacement!, set_body_angular_displacement!
export set_body_linear_velocity!, set_body_angular_velocity!
export set_body_linear_acceleration!, set_body_angular_acceleration!
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
export BodyState, extract_body_state
export PointState, extract_point_state, extract_point_states, extract_point_states!
export ElementState, extract_element_state, extract_element_states, extract_element_states!

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

# airfoil meshing
include("afmesh.jl")

# section properties
include("section.jl")

# assembly creation
include("assembly.jl")

# prescribed conditions, distributed loads, and point masses
include(joinpath("loads", "prescribed.jl"))
include(joinpath("loads", "distributed.jl"))
include(joinpath("loads", "pointmass.jl"))

# system storage and pointers 
include("system.jl")

# input and output methods
include("preprocess.jl")
include("postprocess.jl")

# body residuals and jacobians
include("body.jl")

# point residuals
include(joinpath("point", "residual", "static.jl"))
include(joinpath("point", "residual", "steady.jl"))
include(joinpath("point", "residual", "initial.jl"))
include(joinpath("point", "residual", "newmark.jl"))
include(joinpath("point", "residual", "dynamic.jl"))
include(joinpath("point", "residual", "expanded.jl"))

# point jacobians
include(joinpath("point", "jacobian", "static.jl"))
include(joinpath("point", "jacobian", "steady.jl"))
include(joinpath("point", "jacobian", "initial.jl"))
include(joinpath("point", "jacobian", "newmark.jl"))
include(joinpath("point", "jacobian", "dynamic.jl"))
include(joinpath("point", "jacobian", "expanded.jl"))

# element residuals
include(joinpath("element", "residual", "static.jl"))
include(joinpath("element", "residual", "steady.jl"))
include(joinpath("element", "residual", "initial.jl"))
include(joinpath("element", "residual", "newmark.jl"))
include(joinpath("element", "residual", "dynamic.jl"))
include(joinpath("element", "residual", "expanded.jl"))

# element jacobians
include(joinpath("element", "jacobian", "static.jl"))
include(joinpath("element", "jacobian", "steady.jl"))
include(joinpath("element", "jacobian", "initial.jl"))
include(joinpath("element", "jacobian", "newmark.jl"))
include(joinpath("element", "jacobian", "dynamic.jl"))
include(joinpath("element", "jacobian", "expanded.jl"))

# system analyses
include(joinpath("analyses", "static.jl"))
include(joinpath("analyses", "steady.jl"))
include(joinpath("analyses", "initial.jl"))
include(joinpath("analyses", "newmark.jl"))
include(joinpath("analyses", "eigen.jl"))

# DifferentialEquations Interface
include(joinpath("interfaces", "diffeq.jl"))

# ForwardDiff Interface
include(joinpath("interfaces", "forwarddiff.jl"))

# WriteVTK Interface
include(joinpath("interfaces", "writevtk.jl"))

end
