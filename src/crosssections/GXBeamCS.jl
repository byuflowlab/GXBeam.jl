module GXBeamCS

import ForwardDiff
using UnPack
using SparseArrays
using StaticArrays

# from section
export Material, Node, MeshElement
export initialize_cache, compliance_matrix, mass_matrix, plotmesh
export strain_recovery, plotsoln, tsai_wu

# from afmesh
export Layer
export afmesh

# section properties
include("section.jl")

# airfoil meshing
include("afmesh.jl")

end