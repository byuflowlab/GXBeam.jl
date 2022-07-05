using LinearAlgebra: I, Symmetric, det, factorize
using SparseArrays: spzeros, sparse
using StaticArrays

"""
    Material(E1, E2, E2, G12, G13, G23, nu12, nu13, nu23, rho)

General orthotropic material properties. 
The "1" direction is along the beam axis for a fiber orientation of zero (theta=0).
"2" corresponds to the local x-direction and "3" the local y at theta=0. (see documentation for figures).  
If the two in-plane directions have the same stiffness properties then only one 
needs to be specified (e.g., compatible with the plane-stress assumption of CLT-based tools):
`Material(E1, E2, G12, nu12, rho)`

**Arguments**
- `E::float`: Young's modulus along 1st, 2nd and 3rd axes.
- `G::float`: shear moduli
- `nu::float`: Poisson's ratio.  nu_ij E_j = nu_ji E_i
- `rho::float`: density
"""
struct Material{TF}
    E1::TF
    E2::TF
    E3::TF
    G12::TF
    G13::TF
    G23::TF
    nu12::TF
    nu13::TF
    nu23::TF
    rho::TF
end

Material(E1, E2, G12, nu12, rho) = Material(E1, E2, E2, G12, G12, G12, nu12, nu12, nu12, rho)

"""
    Node(x, y, number)

A node in the finite element mesh at location x, y.  If assembled in a vector, the vector index corresponds to the node number.

**Arguments**
- `x::float`: x location of node in global coordinate system
- `y::float`: y location of node in global coordinate system
"""
struct Node{TF}
    x::TF
    y::TF
end

"""
    Element(nodenum, material, theta)

An element in the mesh, consisting of four ordered nodes, a material, and a fiber orientation.

**Arguments**
- `nodenum::Vector{integer}`: a vector of four node numbers corresponding the the four nodes defining this element (vector indices of the nodes). 
    Node order should be counterclockwise starting from the bottom left node using the local coordinate sytem (see figure).
- `material::Material`: material properties of this element
- `theta::float`: fiber orientation
"""
struct Element{VI, TF}
    nodenum::VI
    material::Material
    theta::TF
end

"""
internal cache so allocations happen only once upfront
"""
struct Cache{TM, TSM, TFM, TV}  # matrix, sparse matrix, float matrix, vector
    Q::TM
    Ttheta::TSM
    Tbeta::TSM
    Z::TSM
    S::TFM
    N::TSM
    SZ::TSM
    SN::TSM
    Bksi::TSM
    Beta::TSM
    dNM_dksi::TSM
    dNM_deta::TSM
    BN::TSM
    Ae::TM
    Re::TM
    Ee::TM
    Ce::TM
    Le::TM
    Me::TM
    idx::TV
    A::TM
    R::TM
    E::TSM
    C::TSM
    L::TM
    M::TSM
end

"""
create chace.  set sizes of static matrices, and set sparsity patterns for those that are fixed.
"""
function initializecache(nodes, elements, etype=Float64)

    # create cache
    Q = zeros(etype, 6, 6)
    
    Ttheta = zeros(etype, 6, 6)
    Ttheta[1, 1] = 1.0; Ttheta[1, 4] = 1.0; Ttheta[1, 6] = 1.0
    Ttheta[2, 2] = 1.0
    Ttheta[3, 3] = 1.0; Ttheta[3, 5] = 1.0
    Ttheta[4, 1] = 1.0; Ttheta[4, 4] = 1.0; Ttheta[4, 6] = 1.0
    Ttheta[5, 3] = 1.0; Ttheta[5, 5] = 1.0
    Ttheta[6, 1] = 1.0; Ttheta[6, 4] = 1.0; Ttheta[6, 6] = 1.0
    Ttheta = sparse(Ttheta)
    
    Tbeta = zeros(etype, 6, 6)
    Tbeta[1, 1] = 1.0; Tbeta[1, 2] = 1.0; Tbeta[1, 3] = 1.0
    Tbeta[2, 1] = 1.0; Tbeta[2, 2] = 1.0; Tbeta[2, 3] = 1.0
    Tbeta[3, 1] = 1.0; Tbeta[3, 2] = 1.0; Tbeta[3, 3] = 1.0
    Tbeta[4, 4] = 1.0; Tbeta[4, 5] = 1.0
    Tbeta[5, 4] = 1.0; Tbeta[5, 5] = 1.0
    Tbeta[6, 6] = 1.0
    Tbeta = sparse(Tbeta)
    
    Z = zeros(etype, 3, 6)  # [I zeros(etype, 3, 3)]
    Z[1, 1] = 1.0
    Z[2, 2] = 1.0
    Z[3, 3] = 1.0
    Z[1, 6] = 1.0
    Z[2, 6] = 1.0
    Z[3, 4] = 1.0
    Z[3, 5] = 1.0
    Z = sparse(Z)

    S = [zeros(3, 3); I]
    S = sparse(S)

    N = zeros(etype, 3, 12)
    N[1, 1] = 1.0
    N[2, 2] = 1.0
    N[3, 3] = 1.0
    N[1, 4] = 1.0
    N[2, 5] = 1.0
    N[3, 6] = 1.0
    N[1, 7] = 1.0
    N[2, 8] = 1.0
    N[3, 9] = 1.0
    N[1, 10] = 1.0
    N[2, 11] = 1.0
    N[3, 12] = 1.0
    N = sparse(N)

    SZ = spzeros(etype, 6, 6)
    SN = spzeros(etype, 6, 12)
    
    Bksi = zeros(etype, 6, 3)
    Bksi[1, 1] = 1.0
    Bksi[2, 2] = 1.0
    Bksi[3, 1] = 1.0
    Bksi[3, 2] = 1.0
    Bksi[4, 3] = 1.0
    Bksi[5, 3] = 1.0
    Bksi = sparse(Bksi)
    
    Beta = zeros(etype, 6, 3)
    Beta[1, 1] = 1.0
    Beta[2, 2] = 1.0
    Beta[3, 1] = 1.0
    Beta[3, 2] = 1.0
    Beta[4, 3] = 1.0
    Beta[5, 3] = 1.0
    Beta = sparse(Beta)
    
    dNM_dksi = zeros(etype, 3, 12)
    dNM_dksi[1, 1] = 1.0
    dNM_dksi[2, 2] = 1.0
    dNM_dksi[3, 3] = 1.0
    dNM_dksi[1, 4] = 1.0
    dNM_dksi[2, 5] = 1.0
    dNM_dksi[3, 6] = 1.0
    dNM_dksi[1, 7] = 1.0
    dNM_dksi[2, 8] = 1.0
    dNM_dksi[3, 9] = 1.0
    dNM_dksi[1, 10] = 1.0
    dNM_dksi[2, 11] = 1.0
    dNM_dksi[3, 12] = 1.0
    dNM_dksi = sparse(dNM_dksi)
    
    dNM_deta = zeros(etype, 3, 12)
    dNM_deta[1, 1] = 1.0
    dNM_deta[2, 2] = 1.0
    dNM_deta[3, 3] = 1.0
    dNM_deta[1, 4] = 1.0
    dNM_deta[2, 5] = 1.0
    dNM_deta[3, 6] = 1.0
    dNM_deta[1, 7] = 1.0
    dNM_deta[2, 8] = 1.0
    dNM_deta[3, 9] = 1.0
    dNM_deta[1, 10] = 1.0
    dNM_deta[2, 11] = 1.0
    dNM_deta[3, 12] = 1.0
    dNM_deta = sparse(dNM_deta)
    
    BN = spzeros(etype, 6, 12)

    Ae = zeros(etype, 6, 6)
    Re = zeros(etype, 12, 6)
    Ee = zeros(etype, 12, 12)
    Ce = zeros(etype, 12, 12)
    Le = zeros(etype, 12, 6)
    Me = zeros(etype, 12, 12)

    idx = zeros(Int64, 12)

    # big system matrices
    ne = length(elements) # number of elements
    nn = length(nodes)  # number of nodes
    ndof = 3 * nn  # 3 displacement dof per node

    A = zeros(etype, 6, 6)  # 6 x 6
    R = zeros(etype, ndof, 6)  #  nn*3 x 6
    E = spzeros(etype, ndof, ndof)  # nn*3 x nn*3
    C = spzeros(etype, ndof, ndof)  # nn*3 x nn*3
    L = zeros(etype, ndof, 6)  # nn*3 x 6
    M = spzeros(etype, ndof, ndof)  # nn*3 x nn*3

    # initialize sparsity pattern
    one12 = ones(etype, 12, 12)
    @views for i = 1:ne
        nodenum = elements[i].nodenum
        idx = node2idx(nodenum)

        E[idx, idx] .= one12
        C[idx, idx] .= one12
        M[idx, idx] .= one12
    end
    # reset to zeros
    E .= 0.0
    C .= 0.0
    M .= 0.0

    cache = Cache(Q, Ttheta, Tbeta, Z, S, N, SZ, SN, Bksi, Beta, dNM_dksi, dNM_deta, BN, Ae, Re, Ee, Ce, Le, Me, idx, A, R, E, C, L, M)

    return cache
end

"""
Constituitive matrix of this material using the internal ordering.
"""
function stiffness!(material, cache) 
    E1 = material.E1; E2 = material.E2; E3 = material.E3
    nu12 = material.nu12; nu13 = material.nu13; nu23 = material.nu23
    G12 = material.G12; G13 = material.G13; G23 = material.G23

    nu21 = nu12*E2/E1
    nu31 = nu13*E3/E1
    nu32 = nu23*E3/E2
    delta = 1.0 / (1 - nu12*nu21 - nu23*nu32 - nu13*nu31 - 2*nu21*nu32*nu13)

    cache.Q .= 0.0  # reset
    cache.Q[6, 6] = E1*(1 - nu23*nu32)*delta
    cache.Q[1, 1] = E2*(1 - nu13*nu31)*delta
    cache.Q[2, 2] = E3*(1 - nu12*nu21)*delta
    cache.Q[1, 6] = E1*(nu21 + nu31*nu23)*delta
    cache.Q[6, 1] = E1*(nu21 + nu31*nu23)*delta
    cache.Q[2, 6] = E1*(nu31 + nu21*nu32)*delta
    cache.Q[6, 2] = E1*(nu31 + nu21*nu32)*delta
    cache.Q[1, 2] = E2*(nu32 + nu12*nu31)*delta
    cache.Q[2, 1] = E2*(nu32 + nu12*nu31)*delta
    cache.Q[4, 4] = G12
    cache.Q[5, 5] = G13
    cache.Q[3, 3] = G23

    
    
    return nothing
end

"""
Rotate constituitive matrix by ply angle
"""
function rotate_to_ply!(theta, cache)
    c = cos(theta)
    s = sin(theta)

    cache.Ttheta[1, 1] = c^2
    cache.Ttheta[1, 4] = 2*s*c
    cache.Ttheta[1, 6] = s^2
    cache.Ttheta[2, 2] = 1.0
    cache.Ttheta[3, 3] = c
    cache.Ttheta[3, 5] = s
    cache.Ttheta[4, 1] = -s*c
    cache.Ttheta[4, 4] = c^2 - s^2
    cache.Ttheta[4, 6] = s*c
    cache.Ttheta[5, 3] = -s
    cache.Ttheta[5, 5] = c
    cache.Ttheta[6, 1] = s^2
    cache.Ttheta[6, 4] = -2*s*c
    cache.Ttheta[6, 6] = c^2

    cache.Q .= cache.Ttheta * Symmetric(cache.Q) * cache.Ttheta'
    return nothing
end

"""
Rotate constituitive matrix by element orientation where `c = cos(beta)` and `s = sin(beta)`
"""
function rotate_to_element!(c, s, cache)  # c = cos(beta), s = sin(beta)

    cache.Tbeta[1, 1] = c^2
    cache.Tbeta[1, 2] = s^2
    cache.Tbeta[1, 3] = -2*s*c
    cache.Tbeta[2, 1] = s^2
    cache.Tbeta[2, 2] = c^2
    cache.Tbeta[2, 3] = 2*s*c
    cache.Tbeta[3, 1] = s*c
    cache.Tbeta[3, 2] = -s*c
    cache.Tbeta[3, 3] = c^2-s^2
    cache.Tbeta[4, 4] = c
    cache.Tbeta[4, 5] = -s
    cache.Tbeta[5, 4] = s
    cache.Tbeta[5, 5] = c
    cache.Tbeta[6, 6] = 1.0

    cache.Q .= cache.Tbeta * Symmetric(cache.Q) * cache.Tbeta'
    return nothing
end


"""
Get element constituitive matrix accounting for fiber orientation and element orientation
"""
function elementQ!(material, theta, cbeta, sbeta, cache)
    stiffness!(material, cache)
    rotate_to_ply!(theta, cache)
    rotate_to_element!(cbeta, sbeta, cache)

    return nothing
end

"""
Compute the integrand for a single element with a given ksi, eta.
"""
function addelementintegrand!(ksi, eta, element, nodes, cache)

    # shape functions
    N = zeros(4)
    N[1] = 0.25*(1 - ksi)*(1 - eta)
    N[2] = 0.25*(1 + ksi)*(1 - eta)
    N[3] = 0.25*(1 + ksi)*(1 + eta)
    N[4] = 0.25*(1 - ksi)*(1 + eta)
    
    # x, y position
    x = 0.0
    y = 0.0
    for i = 1:4
        x += N[i]*nodes[i].x
        y += N[i]*nodes[i].y
    end

    # orientation (beta)
    xl = 0.5*(nodes[1].x + nodes[4].x)
    yl = 0.5*(nodes[1].y + nodes[4].y)
    xr = 0.5*(nodes[2].x + nodes[3].x)
    yr = 0.5*(nodes[2].y + nodes[3].y)
    dx = xr - xl
    dy = yr - yl
    ds = sqrt(dx^2 + dy^2)
    cbeta = dx/ds
    sbeta = dy/ds
    
    # basic matrices
    # Z = [I [0.0 0 -y; 0 0 x; y -x 0]]  # translation and cross product
    cache.Z[1, 6] = -y
    cache.Z[2, 6] = x
    cache.Z[3, 4] = y
    cache.Z[3, 5] = -x
    # S = [zeros(3, 3); I]
    elementQ!(element.material, element.theta, cbeta, sbeta, cache)
    cache.N[1, 1] = N[1]
    cache.N[2, 2] = N[1]
    cache.N[3, 3] = N[1]
    cache.N[1, 4] = N[2]
    cache.N[2, 5] = N[2]
    cache.N[3, 6] = N[2]
    cache.N[1, 7] = N[3]
    cache.N[2, 8] = N[3]
    cache.N[3, 9] = N[3]
    cache.N[1, 10] = N[4]
    cache.N[2, 11] = N[4]
    cache.N[3, 12] = N[4]
    cache.SZ .= cache.S * cache.Z
    cache.SN .= cache.S * cache.N

    # derivatives of shape functions
    dN_dksi = zeros(4)
    dN_dksi[1] = -0.25*(1 - eta)
    dN_dksi[2] = 0.25*(1 - eta)
    dN_dksi[3] = 0.25*(1 + eta)
    dN_dksi[4] = -0.25*(1 + eta)

    dN_deta = zeros(4)
    dN_deta[1] = -0.25*(1 - ksi)
    dN_deta[2] = -0.25*(1 + ksi)
    dN_deta[3] = 0.25*(1 + ksi)
    dN_deta[4] = 0.25*(1 - ksi)

    # Jacobian
    dx_dksi = 0.0
    dx_deta = 0.0
    dy_dksi = 0.0
    dy_deta = 0.0
    for i = 1:4
        dx_dksi += dN_dksi[i]*nodes[i].x
        dx_deta += dN_deta[i]*nodes[i].x
        dy_dksi += dN_dksi[i]*nodes[i].y
        dy_deta += dN_deta[i]*nodes[i].y
    end

    J = [dx_dksi dy_dksi;
         dx_deta dy_deta]
    detJ = det(J)
    Jinv = [dy_deta -dy_dksi;
           -dx_deta dx_dksi] / detJ

    # BN matrix
    cache.Bksi[1, 1] = Jinv[1, 1]
    cache.Bksi[2, 2] = Jinv[2, 1]
    cache.Bksi[3, 1] = Jinv[2, 1]
    cache.Bksi[3, 2] = Jinv[1, 1]
    cache.Bksi[4, 3] = Jinv[1, 1]
    cache.Bksi[5, 3] = Jinv[2, 1]

    cache.Beta[1, 1] = Jinv[1, 2]
    cache.Beta[2, 2] = Jinv[2, 2]
    cache.Beta[3, 1] = Jinv[2, 2]
    cache.Beta[3, 2] = Jinv[1, 2]
    cache.Beta[4, 3] = Jinv[1, 2]
    cache.Beta[5, 3] = Jinv[2, 2]

    # dNM_dksi = [dN_dksi[1]*Matrix(1.0I, 3, 3) dN_dksi[2]*I dN_dksi[3]*I dN_dksi[4]*I]
    cache.dNM_dksi[1, 1] = dN_dksi[1]
    cache.dNM_dksi[2, 2] = dN_dksi[1]
    cache.dNM_dksi[3, 3] = dN_dksi[1]
    cache.dNM_dksi[1, 4] = dN_dksi[2]
    cache.dNM_dksi[2, 5] = dN_dksi[2]
    cache.dNM_dksi[3, 6] = dN_dksi[2]
    cache.dNM_dksi[1, 7] = dN_dksi[3]
    cache.dNM_dksi[2, 8] = dN_dksi[3]
    cache.dNM_dksi[3, 9] = dN_dksi[3]
    cache.dNM_dksi[1, 10] = dN_dksi[4]
    cache.dNM_dksi[2, 11] = dN_dksi[4]
    cache.dNM_dksi[3, 12] = dN_dksi[4]
    # dNM_deta = [dN_deta[1]*Matrix(1.0I, 3, 3) dN_deta[2]*I dN_deta[3]*I dN_deta[4]*I]
    cache.dNM_deta[1, 1] = dN_deta[1]
    cache.dNM_deta[2, 2] = dN_deta[1]
    cache.dNM_deta[3, 3] = dN_deta[1]
    cache.dNM_deta[1, 4] = dN_deta[2]
    cache.dNM_deta[2, 5] = dN_deta[2]
    cache.dNM_deta[3, 6] = dN_deta[2]
    cache.dNM_deta[1, 7] = dN_deta[3]
    cache.dNM_deta[2, 8] = dN_deta[3]
    cache.dNM_deta[3, 9] = dN_deta[3]
    cache.dNM_deta[1, 10] = dN_deta[4]
    cache.dNM_deta[2, 11] = dN_deta[4]
    cache.dNM_deta[3, 12] = dN_deta[4]

    cache.BN .= sparse(cache.Bksi) * sparse(cache.dNM_dksi) + sparse(cache.Beta) * sparse(cache.dNM_deta)
    
    # integrands
    cache.Ae .+= cache.SZ' * cache.Q * cache.SZ * detJ
    cache.Re .+= cache.BN' * cache.Q * cache.SZ * detJ
    cache.Ee .+= cache.BN' * cache.Q * cache.BN * detJ
    cache.Ce .+= cache.BN' * cache.Q * cache.SN * detJ  # use Giavottoa definition
    cache.Le .+= cache.SN' * cache.Q * cache.SZ * detJ
    cache.Me .+= cache.SN' * cache.Q * cache.SN * detJ

    return nothing
end

"""
2-point Gauss quadrature of one element
"""
function submatrix!(element, nodes, cache)

    # initialize
    cache.Ae .= 0.0
    cache.Re .= 0.0
    cache.Ee .= 0.0
    cache.Ce .= 0.0
    cache.Le .= 0.0
    cache.Me .= 0.0

    # 2-point Gauss quadrature in both dimensions
    ksi = [1.0/sqrt(3), 1.0/sqrt(3), -1.0/sqrt(3), -1.0/sqrt(3)]
    eta = [1.0/sqrt(3), -1.0/sqrt(3), 1.0/sqrt(3), -1.0/sqrt(3)]
    
    for i = 1:4
        addelementintegrand!(ksi[i], eta[i], element, nodes, cache)
    end
    
    return nothing
end

"""
Convenience function to map node numbers to locations in global matrix.
"""
function node2idx!(nodenums, cache)
    nn = 4
    # idx = Vector{Int64}(undef, nn*3)
    for i = 1:nn
        cache.idx[((i-1)*3+1):i*3] = ((nodenums[i]-1)*3+1):nodenums[i]*3
    end
    return nothing
end

function node2idx(nodenums)
    nn = 4
    idx = Vector{Int64}(undef, nn*3)
    for i = 1:nn
        idx[((i-1)*3+1):i*3] = ((nodenums[i]-1)*3+1):nodenums[i]*3
    end
    return idx
end




"""
    compliance(nodes, elements)

Compute compliance matrix given the finite element mesh described by nodes and elements.

**Arguments**
- `nodes::Vector{Node}`: all the nodes in the mesh
- `elements::Vector{Element}`: all the elements in the mesh

**Returns**
- `S::Matrix`: compliance matrix
- `sc::Vector{float}`: x, y location of shear center
    (location where a transverse/shear force will not produce any torsion, i.e., beam will not twist)
- `tc::Vector{float}`: x, y location of tension center, aka elastic center, aka centroid 
    (location where an axial force will not produce any bending, i.e., beam will remain straight)
"""
function compliance(nodes, elements, cache=initializecache(nodes, elements))

    # initialize
    ne = length(elements) # number of elements
    nn = length(nodes)  # number of nodes
    ndof = 3 * nn  # 3 displacement dof per node

    # place element matrices in global matrices (scatter)
    @views for i = 1:ne
        nodenum = elements[i].nodenum
        submatrix!(elements[i], nodes[nodenum], cache)
        node2idx!(nodenum, cache)

        cache.A .+= cache.Ae
        cache.R[cache.idx, :] .+= cache.Re
        cache.E[cache.idx, cache.idx] .+= cache.Ee
        cache.C[cache.idx, cache.idx] .+= cache.Ce
        cache.L[cache.idx, :] .+= cache.Le
        cache.M[cache.idx, cache.idx] .+= cache.Me
    end
    A = cache.A
    R = cache.R
    E = cache.E
    C = cache.C
    L = cache.L
    M = cache.M

    # assemble displacement constraint matrix
    DT = spzeros(6, ndof)
    for i = 1:nn
        s = 3*(i-1)
        DT[1, s+1] = 1.0
        DT[2, s+2] = 1.0
        DT[3, s+3] = 1.0
        DT[4, s+3] = nodes[i].y
        DT[5, s+3] = -nodes[i].x
        DT[6, s+1] = -nodes[i].y
        DT[6, s+2] = nodes[i].x
    end
    D = sparse(DT')

    # Tr matrix
    Tr = spzeros(6, 6)
    Tr[1, 5] = -1.0
    Tr[2, 4] = 1.0

    # solve first linear system
    AM = [E R D;
        R' A zeros(6, 6);
        DT zeros(6, 12)]
    AM = factorize(Symmetric(sparse(AM)))
    
    B2 = sparse([zeros(ndof, 6); Tr'; zeros(6, 6)])
    X2 = zeros(ndof+12, 6)
    for i = 1:6
        X2[:, i] = AM\B2[:, i]
    end
    dX = X2[1:ndof, :]
    dY = X2[ndof+1:ndof+6, :]

    # solve second linear system
    Bsub1 = [C'-C  L;  # NOTE: error in manual should be C' - C
            -L' zeros(6, 6);  # NOTE: sign fix, error in BECAS documentation
            zeros(6, ndof+6)]
    Bsub2 = [zeros(ndof, 6); I; zeros(6, 6)]
    B1 = sparse(Bsub1*X2[1:end-6, :] + Bsub2)
    X1 = zeros(ndof+12, 6)
    for i = 1:6
        X1[:, i] = AM\B1[:, i]
    end
    X = X1[1:ndof, :]
    Y = X1[ndof+1:ndof+6, :]

    # compliance matrix
    XY = [X; dX; Y]
    S = XY'*[E C R; C' M L; R' L' A]*XY

    # S += 2*X'*L*dY  # NOTE: add additional term from 2015 paper, though I don't see how that follows...

    # # stiffness matrix
    # K = inv(Symmetric(S))

    xs = -S[6, 2]/S[6, 6]
    ys = S[6, 1]/S[6, 6]
    xt = (S[4, 4]*S[5, 3] - S[4, 5]*S[4, 3])/(S[4, 4]*S[5, 5] - S[4, 5]^2)
    yt = (-S[4, 3]*S[5, 5] + S[4, 5]*S[5, 3])/(S[4, 4]*S[5, 5] - S[4, 5]^2) 
    sc = [xs, ys]
    tc = [xt, yt]

    # K11 = [E R D; R' A zeros(6, 6); DT zeros(6, 12)]
    # K12 = [C'-C  -L zeros(ndof, 6);
    #         L' zeros(6, 12);
    #         zeros(6, ndof+12)]
    # W = [K11 K12; zeros(size(K11)) K11] \ [Bsub2; B2]
    # G11 = [E R zeros(ndof, 6); R' A zeros(6, 6); zeros(6, ndof+12)]
    # G12 = [C L zeros(ndof, 6); L' zeros(6, 12); zeros(6, ndof+12)]
    # G22 = [M zeros(ndof, 12); zeros(12, ndof+12)]
    # G = [G11 G12; G12' G22]
    # Fs = W'*G*W

    return S, sc, tc
end

# function centers(S)  #, z, L)
#     # assume free end
#     # xs = -(S[6, 2] + S[6, 4]*(L - z))/S[6, 6]
#     # ys = (S[6, 1] + S[6, 5]*(L - z))/S[6, 6]
#     xs = -S[6, 2]/S[6, 6]
#     ys = S[6, 1]/S[6, 6]
#     xt = (S[4, 4]*S[5, 3] - S[4, 5]*S[4, 3])/(S[4, 4]*S[5, 5] - S[4, 5]^2)
#     yt = (-S[4, 3]*S[5, 5] + S[4, 5]*S[5, 3])/(S[4, 4]*S[5, 5] - S[4, 5]^2) 

#     return xs, ys, xt, yt
# end

"""
Reorder stiffness or compliance matrix from internal order to GXBeam order
"""
function reorder(K)  # reorder to GXBeam format
    idx = [3, 1, 2, 6, 4, 5]
    return K[idx, idx]
end


function plotmesh(nodes, elements; plotnumbers=false)
    ne = length(elements)
    
    figure()
    for i = 1:ne
        node = nodes[elements[i].nodenum]
        for i = 1:4
            iplus = i+1
            if iplus == 5
                iplus = 1 
            end
            plot([node[i].x, node[iplus].x], [node[i].y, node[iplus].y])
        end
        if plotnumbers
            barx = sum([n.x/4 for n in node])
            bary = sum([n.y/4 for n in node])
            text(barx, bary, string(i), color="r")
        end
    end
    if plotnumbers
        nn = length(nodes)
        for i = 1:nn
            text(nodes[i].x, nodes[i].y, string(i))
        end
    end
end



# --------- Unit Tests --------
using Test

E1 = rand()
E2 = rand()
E3 = rand()
G12 = rand()
G13 = rand()
G23 = rand()
nu12 = rand()
nu23 = rand()
nu13 = rand()
rho = 1.0
mat = Material(E1, E2, E3, G12, G13, G23, nu12, nu13, nu23, rho)
cache = initializecache([Node(0.0, 0.0), Node(0.0, 0.0), Node(0.0, 0.0), Node(0.0, 0.0)], [Element([1, 2, 3, 4], mat, 0.0)])
stiffness!(mat, cache)
Q1 = cache.Q

nu21 = nu12*E2/E1
nu31 = nu13*E3/E1
nu32 = nu23*E3/E2
C2 = [1/E2 -nu23/E2 0.0 0 0 -nu12/E1;
       -nu23/E2 1/E3 0 0 0 -nu13/E1;
       0 0 1/G23 0 0 0;
       0 0 0 1/G12 0 0
       0 0 0 0 1/G13 0;
       -nu12/E1 -nu13/E1 0 0 0 1/E1]
Q2 = inv(C2)


@test all(isapprox.(Q1, Q2, atol=1e-6))


# ---- Case 1: Square cross section of isotropic material - S1
iso1 = Material(100.0, 100.0, 100.0, 41.667, 41.667, 41.667, 0.2, 0.2, 0.2, 1.0)
x = range(-0.05, 0.05, length=11)
y = range(-0.05, 0.05, length=11)

nodes = Vector{Node}(undef, 11*11)
elements = Vector{Element}(undef, 10*10)

let
m = 1
for i = 1:11
    for j = 1:11
        nodes[m] = Node(x[i], y[j])
        m += 1
    end
end

m = 1
for i = 1:10
    for j = 1:10
        elements[m] = Element([11*(i-1)+j, 11*(i)+j, 11*(i)+j+1, 11*(i-1)+j+1], iso1, 0.0)
        m += 1
    end
end

end

using PyPlot
close("all"); pygui(true);

# plotmesh(nodes, elements, plotnumbers=true)

S, sc, tc = compliance(nodes, elements)
K = inv(S)

@test isapprox(K[1, 1], 3.4899e-1, atol=0.0001e-1)
@test isapprox(K[2, 2], 3.4899e-1, atol=0.0001e-1)
@test isapprox(K[3, 3], 1.0, atol=1e-4)
@test isapprox(K[4, 4], 8.3384e-4, atol=0.0001e-4)
@test isapprox(K[5, 5], 8.3384e-4, atol=0.0001e-4)
@test isapprox(K[6, 6], 5.9084e-4, atol=0.0001e-4)

@test isapprox(sc[1], 0.0, atol=1e-8)
@test isapprox(sc[2], 0.0, atol=1e-8)
@test isapprox(tc[1], 0.0, atol=1e-8)
@test isapprox(tc[2], 0.0, atol=1e-8)

# --------- Case 2 --------
alpha = 1e1
iso2 = Material(100.0/alpha, 100.0/alpha, 41.667/alpha, 0.2, 1.0)  # note error in user guide for nu

let
m = 1
for i = 1:10
    for j = 1:10
        if i <= 5
            elements[m] = Element([11*(i-1)+j, 11*(i)+j, 11*(i)+j+1, 11*(i-1)+j+1], iso1, 0.0)
        else
            elements[m] = Element([11*(i-1)+j, 11*(i)+j, 11*(i)+j+1, 11*(i-1)+j+1], iso2, 0.0)
        end
        m += 1
    end
end
end


S, sc, tc = compliance(nodes, elements)
K = inv(S)

@test isapprox(K[1, 1], 1.28e-1, atol=0.01e-1)
@test isapprox(K[2, 2], 1.92e-1, atol=0.01e-1)
@test isapprox(K[3, 3], 5.50e-1, atol=0.01e-1)
@test isapprox(K[4, 4], 4.59e-4, atol=0.01e-4)
@test isapprox(K[5, 5], 4.59e-4, atol=0.01e-4)
@test isapprox(K[6, 6], 2.77e-4, atol=0.01e-4)
@test isapprox(K[2, 6], -3.93e-3, atol=0.01e-3)
@test isapprox(K[3, 5], 1.13e-2, atol=0.01e-2)

# ------ case 3 -------
ortho = Material(480.0, 120.0, 120.0, 60.0, 50.0, 60.0, 0.19, 0.26, 0.19, 1.0)

theta = 0.0
let
    m = 1
    for i = 1:10
        for j = 1:10
            elements[m] = Element([11*(i-1)+j, 11*(i)+j, 11*(i)+j+1, 11*(i-1)+j+1], ortho, theta)
            m += 1
        end
    end
end

S, sc, tc = compliance(nodes, elements)
K = inv(S)
@test isapprox(K[1, 1], 5.039E-01, atol=0.001e-1)
@test isapprox(K[2, 2], 4.201E-01, atol=0.001e-1)
@test isapprox(K[3, 3], 4.800E+00, atol=0.001e0)
@test isapprox(K[4, 4], 4.001E-03, atol=0.001e-3)
@test isapprox(K[5, 5], 4.001E-03, atol=0.001e-3)
@test isapprox(K[6, 6], 7.737E-04, atol=0.001e-4)

theta = 22.5*pi/180
let
    m = 1
    for i = 1:10
        for j = 1:10
            elements[m] = Element([11*(i-1)+j, 11*(i)+j, 11*(i)+j+1, 11*(i-1)+j+1], ortho, theta)
            m += 1
        end
    end
end

S, sc, tc = compliance(nodes, elements)
K = inv(S)
@test isapprox(K[1, 1], 7.598E-01, atol=0.001e-1)
@test isapprox(K[2, 2], 4.129E-01, atol=0.001e-1)
@test isapprox(K[3, 3], 3.435E+00, atol=0.001e0)
@test isapprox(K[4, 4], 2.489E-03, atol=0.001e-3)
@test isapprox(K[5, 5], 2.274E-03, atol=0.001e-3)
@test isapprox(K[6, 6], 9.499E-04, atol=0.001e-4)
@test isapprox(K[1, 3], 7.387E-01, atol=0.001e-1)
@test isapprox(K[4, 6], -4.613E-04, atol=0.001e-4)

theta = 90*pi/180
let
    m = 1
    for i = 1:10
        for j = 1:10
            elements[m] = Element([11*(i-1)+j, 11*(i)+j, 11*(i)+j+1, 11*(i-1)+j+1], ortho, theta)
            m += 1
        end
    end
end

S, sc, tc = compliance(nodes, elements)
K = inv(S)
@test isapprox(K[1, 1], 5.0202E-01, atol=0.0001e-1)
@test isapprox(K[2, 2], 5.0406E-01, atol=0.0001e-1)
@test isapprox(K[3, 3], 1.2, atol=0.001e0)
@test isapprox(K[4, 4], 1.0004E-03, atol=0.0001e-3)
@test isapprox(K[5, 5], 1.0002E-03, atol=0.0001e-3)
@test isapprox(K[6, 6], 8.5081E-04, atol=0.0001e-4)


# --------- cylinder -------
R = 0.1
t = 0.01
nr = 4
nt = 120
r = range(R - t, R, length=nr)
theta = range(0, 2*pi, length=nt)

nodes = Vector{Node}(undef, nr*(nt-1))
elements = Vector{Element}(undef, (nr-1)*(nt-1))
let 
m = 1
for i = 1:nt-1
    for j = 1:nr
        nodes[m] = Node(r[j]*cos(theta[i]), r[j]*sin(theta[i]))
        m += 1
    end
end
end

thetam = 0.5* (theta[1:end-1] + theta[2:end])
beta = thetam .- pi/2

let
n = 1
for i = 1:nt-1
    for j = 1:nr-1
        if i == nt-1
            ip = 0
        else
            ip = i
        end

        elements[n] = Element([nr*ip+j, nr*(i-1)+j, nr*(i-1)+j+1, nr*ip+j+1], iso1, 0.0)
        n += 1
    end
end
end

# plotmesh(nodes, elements)

S, sc, tc = compliance(nodes, elements)
K = inv(S)

@test isapprox(K[1, 1], 1.249E-01, atol=0.001e-1/2)
@test isapprox(K[2, 2], 1.249E-01, atol=0.001e-1/2)
@test isapprox(K[3, 3], 5.965E-01, atol=0.0015e-1)
@test isapprox(K[4, 4], 2.697E-03, atol=0.002e-3)  # mesh discretization is not the same
@test isapprox(K[5, 5], 2.697E-03, atol=0.002e-3)
@test isapprox(K[6, 6], 2.248E-03, atol=0.001e-3)


# ---- C2 ----

R = 0.1
t = 0.01
nr = 4
nt = 60
r = range(R - t, R, length=nr)
theta = range(pi/2, 3*pi/2, length=nt)

nodes = Vector{Node}(undef, nr*nt)
elements = Vector{Element}(undef, (nr-1)*(nt-1))
let 
m = 1
for i = 1:nt
    for j = 1:nr
        nodes[m] = Node(r[j]*cos(theta[i]), r[j]*sin(theta[i]))
        m += 1
    end
end
end


let
n = 1
for i = 1:nt-1
    for j = 1:nr-1
        elements[n] = Element([nr*i+j, nr*(i-1)+j, nr*(i-1)+j+1, nr*i+j+1], iso1, 0.0)
        n += 1
    end
end
end

# plotmesh(nodes, elements)

S, sc, tc = compliance(nodes, elements)
K = inv(S)

@test isapprox(K[1, 1], 4.964E-02, atol=0.002e-2)
@test isapprox(K[2, 2], 6.244E-02, atol=0.002e-2)
@test isapprox(K[3, 3], 2.982E-01, atol=0.002e-1)
@test isapprox(K[4, 4], 1.349E-03, atol=0.001e-3)
@test isapprox(K[5, 5], 1.349E-03, atol=0.001e-3)
@test isapprox(K[6, 6], 9.120E-04, atol=0.003e-4)
@test isapprox(K[3, 5], 1.805E-02, atol=0.001e-2)
@test isapprox(K[2, 6], -7.529E-03, atol=0.003e-3) 

@test isapprox(sc[1], -1.206E-01, atol=0.001e-1) 
@test isapprox(sc[2], 0.0, atol=1e-6) 
@test isapprox(tc[1], -6.051E-02, atol=0.002e-2) 
@test isapprox(tc[2], 0.0, atol=1e-6) 


# A critical assessment of computer tools for calculating composite wind turbine blade properties
# ------ circular tube -------
E = 73e9
nu = 0.33
G = E/(2*(1 + nu))
tratio = 1.0/3
R = 0.3
t = tratio*2*R
circmat = Material(E, E, G, nu, 1.0)

nr = 20
nt = 100
r = range(R - t, R, length=nr)
theta = range(0.0, 2*pi, length=nt)

nodes = Vector{Node}(undef, nr*(nt-1))
elements = Vector{Element}(undef, (nr-1)*(nt-1))
let 
m = 1
for i = 1:nt-1
    for j = 1:nr
        nodes[m] = Node(r[j]*cos(theta[i]), r[j]*sin(theta[i]))
        m += 1
    end
end

n = 1
for i = 1:nt-1
    for j = 1:nr-1
        if i == nt-1
            ip = 0
        else
            ip = i
        end

        elements[n] = Element([nr*ip+j, nr*(i-1)+j, nr*(i-1)+j+1, nr*ip+j+1], circmat, 0.0)
        n += 1
    end
end
end

# plotmesh(nodes, elements)

S, sc, tc = compliance(nodes, elements)
K = inv(S)
K2 = reorder(K)

@test isapprox(K2[1, 1], 1.835e10, rtol=0.001)
@test isapprox(K2[2, 2], 4.682e9, rtol=0.001)
@test isapprox(K2[3, 3], 4.682e9, rtol=0.001)
@test isapprox(K2[4, 4], 3.519e8, rtol=0.03)
@test isapprox(K2[5, 5], 4.587e8, rtol=0.002)
@test isapprox(K2[6, 6], 4.587e8, rtol=0.002)

#  --- Generalized Timoshenko Theory of the Variational Asymptotic Beam Sectional Analysis ----
# ----- multi-layer composite pipe --------

# E1 = 141.963e9
# E2 = 9.79056e9
# nu12 = 0.42
# G12 = 59.9844e9
E1 = 20.59e6
E2 = 1.42e6
nu12 = 0.42
G12 = 0.87e6
rho = 1.0
pipemat = Material(E1, E2, G12, nu12, rho)

nx = 50
nt = 24
nr = 20
nodes = Vector{Node}(undef, (2*(nx-1) + 2*(nt-1))*nr)
elements = Vector{Element}(undef, (2*(nx-1) + 2*(nt-1))*(nr-1))

let
    # x1 = -50.8e-3/2
    # y1 = 7.62e-3
    # x2 = 50.8e-3/2
    # y2 = 25.4e-3/2
    x1 = -1.0
    y1 = 0.3
    x2 = 1.0
    y2 = 0.5
    
    x = range(x1, x2, length=nx)
    y = range(y1, y2, length=nr)

    m = 1
    for i = 1:nx
        for j = 1:nr
            nodes[m] = Node(x[i], y[j])
            m += 1
        end
    end

    n = 1
    for i = 1:nx-1
        for j = 1:nr-1
            if j <= (nr ÷ 2)
                theta = 90*pi/180
            else
                theta = 0.0
            end
            elements[n] = Element([nr*(i-1)+j, nr*(i)+j, nr*(i)+j+1, nr*(i-1)+j+1], pipemat, theta)
            n += 1
        end
    end

    I = nx-1

    r = y
    theta = reverse(range(-pi/2, pi/2, length=nt))
    for i = 2:nt
        for j = 1:nr
            nodes[m] = Node(x2 + r[j]*cos(theta[i]), r[j]*sin(theta[i]))
            m += 1
        end
    end

    for i = 1:nt-1
        for j = 1:nr-1
            if j <= (nr ÷ 2)
                theta = 45*pi/180
            else
                theta = -45*pi/180
            end
            elements[n] = Element([nr*(I+i-1)+j, nr*(I+i)+j, nr*(I+i)+j+1, nr*(I+i-1)+j+1], pipemat, theta)
            n += 1
        end
    end

    I += nt-1

    # x1 = 50.8e-3/2
    # y1 = -7.62e-3
    # x2 = -50.8e-3/2
    # y2 = -25.4e-3/2
    x1 = 1.0
    y1 = -0.3
    x2 = -1.0
    y2 = -0.5

    x = reverse(range(x2, x1, length=nx))
    y = reverse(range(y2, y1, length=nr))

    for i = 2:nx
        for j = 1:nr
            nodes[m] = Node(x[i], y[j])
            m += 1
        end
    end

    for i = 1:nx-1
        for j = 1:nr-1
            if j <= (nr ÷ 2)
                theta = 90*pi/180
            else
                theta = 0.0
            end
            elements[n] = Element([nr*(I+i-1)+j, nr*(I+i)+j, nr*(I+i)+j+1, nr*(I+i-1)+j+1], pipemat, theta)
            n += 1
        end
    end

    I += nx-1

    theta = reverse(range(pi/2, 3*pi/2, length=nt))
    for i = 2:nt-1
        for j = 1:nr
            nodes[m] = Node(x2 + r[j]*cos(theta[i]), r[j]*sin(theta[i]))
            m += 1
        end
    end

    for i = 1:nt-1
        for j = 1:nr-1
            if i == nt-1
                Ip = -i
            else
                Ip = I
            end
            if j <= (nr ÷ 2)
                theta = 45*pi/180
            else
                theta = -45*pi/180
            end
            elements[n] = Element([nr*(I+i-1)+j, nr*(Ip+i)+j, nr*(Ip+i)+j+1, nr*(I+i-1)+j+1], pipemat, theta)
            n += 1
        end
    end

end

# plotmesh(nodes, elements)

S, sc, tc = compliance(nodes, elements)
K = inv(S)
K2 = reorder(K)


# using BenchmarkTools
# cache = initializecache(nodes, elements)
# @btime compliance($nodes, $elements, $cache)
# @profview compliance(nodes, elements, cache)



@test isapprox(K2[1, 1], 1.03892e7, rtol=0.04)
@test isapprox(K2[2, 2], 7.85310e5, rtol=0.01)
@test isapprox(K2[3, 3], 3.29279e5, rtol=0.02)
@test isapprox(K2[1, 4], 9.84575e4, rtol=0.12)
@test isapprox(K2[2, 5], -8.21805e3, rtol=0.11)
@test isapprox(K2[3, 6], -5.20981e4, rtol=0.21)   
@test isapprox(K2[4, 4], 6.87275e5, rtol=0.01)
@test isapprox(K2[5, 5], 1.88238e6, rtol=0.04)
@test isapprox(K2[6, 6], 5.38987e6, rtol=0.03)

# println("K11 = ", round((K2[1, 1]/1.03892e7 - 1)*100, digits=2), "%")
# println("K22 = ", round((K2[2, 2]/7.85310e5 - 1)*100, digits=2), "%")
# println("K33 = ", round((K2[3, 3]/3.29279e5 - 1)*100, digits=2), "%")
# println("K14 = ", round((K2[1, 4]/9.84575e4 - 1)*100, digits=2), "%")
# println("K25 = ", round((K2[2, 5]/-8.21805e3 - 1)*100, digits=2), "%")
# println("K36 = ", round((K2[3, 6]/-5.20981e4 - 1)*100, digits=2), "%")
# println("K44 = ", round((K2[4, 4]/6.87275e5 - 1)*100, digits=2), "%")
# println("K55 = ", round((K2[5, 5]/1.88238e6 - 1)*100, digits=2), "%")
# println("K66 = ", round((K2[6, 6]/5.38987e6 - 1)*100, digits=2), "%")

