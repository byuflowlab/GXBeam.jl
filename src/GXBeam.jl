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

const GAUSS_NODES = SVector(-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526)
const GAUSS_WEIGHTS = SVector(0.34785484513745385, 0.6521451548625462, 0.6521451548625462, 0.34785484513745385)

# common math functions
include("math.jl")

# body residuals and jacobians
include("body.jl")

# assembly creation
include("assembly.jl")

# prescribed conditions, distributed loads, and point masses
include("loads/prescribed.jl")
include("loads/distributed.jl")
include("loads/pointmass.jl")

# system storage and pointers 
include("system.jl")

# input and output methods
include("preprocess.jl")
include("postprocess.jl")

# point residuals
include("point/residual/static.jl")
include("point/residual/steady.jl")
include("point/residual/initial.jl")
include("point/residual/newmark.jl")
include("point/residual/dynamic.jl")
include("point/residual/expanded.jl")

# point jacobians
include("point/jacobian/static.jl")
include("point/jacobian/steady.jl")
include("point/jacobian/initial.jl")
include("point/jacobian/newmark.jl")
include("point/jacobian/dynamic.jl")
include("point/jacobian/expanded.jl")

# element residuals
include("element/residual/static.jl")
include("element/residual/steady.jl")
include("element/residual/initial.jl")
include("element/residual/newmark.jl")
include("element/residual/dynamic.jl")
include("element/residual/expanded.jl")

# element jacobians
include("element/jacobian/static.jl")
include("element/jacobian/steady.jl")
include("element/jacobian/initial.jl")
include("element/jacobian/newmark.jl")
include("element/jacobian/dynamic.jl")
include("element/jacobian/expanded.jl")

# system analyses
include("analyses/static.jl")
include("analyses/steady.jl")
include("analyses/initial.jl")
include("analyses/newmark.jl")
include("analyses/eigen.jl")

# DifferentialEquations Interface
include("interfaces/diffeq.jl")

# ForwardDiff Interface
include("interfaces/forwarddiff.jl")

# WriteVTK Interface
include("interfaces/writevtk.jl")

end
