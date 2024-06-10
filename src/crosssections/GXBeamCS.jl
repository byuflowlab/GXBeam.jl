module GXBeamCS

import ForwardDiff
using UnPack
using SparseArrays
using StaticArrays

# from fem
export Material, Node, MeshElement
export initialize_cache, compliance_matrix, mass_matrix, plotmesh
export strain_recovery, plotsoln, tsai_wu

# from afmesh
export Layer
export afmesh

# section properties using finite element based approach
include("fem.jl")

# airfoil meshing
include("afmesh.jl")

end