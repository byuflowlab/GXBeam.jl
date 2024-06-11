module GXBeamCS

import ForwardDiff
using UnPack
using SparseArrays
using StaticArrays
using FLOWMath
using LinearAlgebra

# from afmesh
export Layer
export afmesh

# from fem
export Material, Node, MeshElement
export initialize_cache, compliance_matrix, mass_matrix, plotmesh
export strain_recovery, plotsoln, tsai_wu

# from clt
export MaterialPlane, Lamina, BeamSection
export beamstiffness

# section properties using finite element based approach
include("fem.jl")

# airfoil meshing
include("afmesh.jl")

# section properties using classical laminate theory
include("clt.jl")


end