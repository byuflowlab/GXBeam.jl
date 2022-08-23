"""
    Material(E1, E2, E3, G12, G13, G23, nu12, nu13, nu23, rho)

General orthotropic material properties. 
1 is along main ply axis. 2 is transverse. 3 is normal to ply.
for a fiber orientation of zero, 1 is along the beam axis.

**Arguments**
- `E::float`: Young's modulus along 1st, 2nd and 3rd axes.
- `G::float`: shear moduli
- `nu::float`: Poisson's ratio.  ``nu_ij E_j = nu_ji E_i``
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

Base.eltype(::Material{TF}) where TF = TF
Base.eltype(::Type{Material{TF}}) where TF = TF

Material{TF}(m::Material) where {TF} = Material{TF}(m.E1, m.E2, m.E3, m.G12, m.G13, m.G23, 
    m.nu12, m.nu13, m.nu23, m.rho)
Base.convert(::Type{Material{TF}}, m::Material) where {TF} = Material{TF}(m)

"""
    Node(x, y)

A node in the finite element mesh at location x, y.  If assembled in a vector, the vector index corresponds to the node number.

**Arguments**
- `x::float`: x location of node in global coordinate system
- `y::float`: y location of node in global coordinate system
"""
struct Node{TF}
    x::TF
    y::TF
end

Base.eltype(::Node{TF}) where TF = TF
Base.eltype(::Type{Node{TF}}) where TF = TF

Node{TF}(n::Node) where {TF} = Node{TF}(n.x, n.y)
Base.convert(::Type{Node{TF}}, n::Node) where {TF} = Node{TF}(n)

"""
    MeshElement(nodenum, material, theta)

An element in the mesh, consisting of four ordered nodes, a material, and a fiber orientation.

**Arguments**
- `nodenum::Vector{integer}`: a vector of four node numbers corresponding the the four nodes defining this element (vector indices of the nodes). 
    Node order should be counterclockwise starting from the bottom left node using the local coordinate sytem (see figure).
- `material::Material`: material properties of this element
- `theta::float`: fiber orientation
"""
struct MeshElement{VI, TF}
    nodenum::VI
    material::Material{TF}
    theta::TF
end

Base.eltype(::MeshElement{VI, TF}) where {VI, TF} = TF
Base.eltype(::Type{MeshElement{VI, TF}}) where {VI, TF} = TF

MeshElement{VI,TF}(e::MeshElement) where {VI,TF} = MeshElement{VI,TF}(e.nodenum, e.material, e.theta)
Base.convert(::Type{MeshElement{VI,TF}}, e::MeshElement) where {VI,TF} = MeshElement{VI,TF}(e)

"""
internal cache so allocations happen only once upfront
"""
struct SectionCache{TM, TSM, TFM, TV}  # matrix, sparse matrix, float matrix, vector
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
    initialize_cache(nodes, elements, etype=Float64)

create cache.  set sizes of static matrices, and set sparsity patterns for those that are fixed.
"""
function initialize_cache(nodes, elements, etype=Float64)

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

    cache = SectionCache(Q, Ttheta, Tbeta, Z, S, N, SZ, SN, Bksi, Beta, dNM_dksi, dNM_deta, BN, Ae, Re, Ee, Ce, Le, Me, idx, A, R, E, C, L, M)

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
Reorder stiffness or compliance matrix from internal order to GXBeam order
"""
function reorder(K)  # reorder to GXBeam format
    idx = [3, 1, 2, 6, 4, 5]
    return K[idx, idx]
end

function linearsolve1(A, B1, AF)

    X1 = zeros(size(B1))
    _, n = size(B1)
    for j = 1:n
        X1[:, j] = AF\B1[:, j]
    end

    return X1
end

function linearsolve2(AM, B2)
    AM = factorize(Symmetric(sparse(AM)))

    X2 = zeros(size(B2))
    _, n = size(B2)
    for j = 1:n
        X2[:, j] = AM\B2[:, j]
    end

    return AM, X2
end

function linearsolve1(A::SparseMatrixCSC{<:ForwardDiff.Dual{T}}, B1, AF) where {T}

     # extract primal values
    #  Av = ForwardDiff.value.(A)
     B1v = ForwardDiff.value.(B1)
 
     # linear solve
    #  Av = factorize(Av)
     x1v = zeros(size(B1v))
     _, n = size(B1v)
     for j = 1:n
         x1v[:, j] = AF\B1v[:, j]
     end
 
     # extract dual values
     ap = ForwardDiff.partials.(A)
     m, n = size(A)
     d = length(ap[1, 1])
     Adot = zeros(m, n, d)  # TODO: remove these allocations
     for i = 1:m
         for j = 1:n
             Adot[i, j, :] .= ap[i, j].values
         end
     end    
     
     bp = ForwardDiff.partials.(B1)
     m, n = size(B1v)
     B1dot = zeros(m, n, d)  # TODO: and these
     for i = 1:m
         for j = 1:n
             B1dot[i, j, :] .= bp[i, j].values
         end
     end
 
     # analytic derivative of linear solve
     x1dot = zeros(m, n, d)  # TODO: and these
     for i = 1:d
         for j = 1:n
             x1dot[:, j, i] = AF \ (B1dot[:, j, i] - Adot[:, :, i] * x1v[:, j])
         end
     end
 
     # repack in dual
    x1 = ForwardDiff.Dual{T}.(x1v, ForwardDiff.Partials.(Tuple.(view(x1dot, i, j, :) for i = 1:m, j = 1:n)))
    
    return x1

end

function linearsolve2(A::SparseMatrixCSC{<:ForwardDiff.Dual{T}}, B2) where {T}

    # extract primal values
    Av = ForwardDiff.value.(A)
    B2v = B2  # no partials

    # linear solve
    AF = factorize(Symmetric(sparse(Av)))
    x2v = zeros(size(B2v))
    _, n = size(B2v)
    for j = 1:n
        x2v[:, j] = AF\B2v[:, j]
    end

    # extract dual values
    ap = ForwardDiff.partials.(A)
    m, n = size(Av)
    d = length(ap[1, 1])
    Adot = zeros(m, n, d)  # TODO: eliminate these allocations
    for i = 1:m
        for j = 1:n
            Adot[i, j, :] .= ap[i, j].values
        end
    end    

    m, n = size(B2v)
    
    # analytic derivative of linear solve
    x2dot = zeros(m, n, d)
    for i = 1:d
        for j = 1:n
            x2dot[:, j, i] = AF \ (-Adot[:, :, i] * x2v[:, j])
        end
    end

    # repack in dual
   x2 = ForwardDiff.Dual{T}.(x2v, ForwardDiff.Partials.(Tuple.(view(x2dot, i, j, :) for i = 1:m, j = 1:n)))

   return AF, x2

end


"""
    compliance_matrix(nodes, elements; cache=initialize_cache(nodes, elements), gxbeam_order=true)

Compute compliance matrix given the finite element mesh described by nodes and elements.

**Arguments**
- `nodes::Vector{Node}`: all the nodes in the mesh
- `elements::Vector{MeshElement}`: all the elements in the mesh
- `gxbeam_order::Bool`: true if output compliance matrix should be in GXBeam order or internal ordering
- `cache::SectionCache`: if number of nodes, number of elements, and connectivity of mesh stays the same (and you will be repeating calls)
    then you can should initialize cache yourself and pass in so you don't have to keep reconstructing it.

**Returns**
- `S::Matrix`: compliance matrix (about the shear center as long as gxbeam_order = true)
- `sc::Vector{float}`: x, y location of shear center
    (location where a transverse/shear force will not produce any torsion, i.e., beam will not twist)
- `tc::Vector{float}`: x, y location of tension center, aka elastic center, aka centroid 
    (location where an axial force will not produce any bending, i.e., beam will remain straight)
"""
function compliance_matrix(nodes, elements; cache=initialize_cache(nodes, elements, promote_type(eltype(eltype(nodes)), eltype(eltype(elements)))), gxbeam_order=true)

    TF = promote_type(eltype(eltype(nodes)), eltype(eltype(elements)))

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
    DT = spzeros(TF, 6, ndof)
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
    # AM = factorize(Symmetric(sparse(AM)))
    
    B2 = sparse([zeros(ndof, 6); Tr'; zeros(6, 6)])
    AF, X2 = linearsolve2(AM, B2)
    dX = X2[1:ndof, :]
    dY = X2[ndof+1:ndof+6, :]

    # solve second linear system
    Bsub1 = [C'-C  L;  # NOTE: error in manual should be C' - C
            -L' zeros(6, 6);  # NOTE: sign fix, error in BECAS documentation
            zeros(6, ndof+6)]
    Bsub2 = [zeros(ndof, 6); I; zeros(6, 6)]
    B1 = sparse(Bsub1*X2[1:end-6, :] + Bsub2)
    X1 = linearsolve1(AM, B1, AF)
    X = X1[1:ndof, :]
    Y = X1[ndof+1:ndof+6, :]

    # compliance matrix
    XY = [X; dX; Y]
    S = XY'*[E C R; C' M L; R' L' A]*XY

    xs = -S[6, 2]/S[6, 6]
    ys = S[6, 1]/S[6, 6]
    xt = (S[4, 4]*S[5, 3] - S[4, 5]*S[4, 3])/(S[4, 4]*S[5, 5] - S[4, 5]^2)
    yt = (-S[4, 3]*S[5, 5] + S[4, 5]*S[5, 3])/(S[4, 4]*S[5, 5] - S[4, 5]^2) 
    sc = [xs, ys]
    tc = [xt, yt]


    if gxbeam_order
        # change ordering to match gxbeam
        S = reorder(S)

        # move properties to be about the shear center (note x->y, y->z b.c. x is axial in derivation)
        P = [0.0 -ys xs
             ys   0   0
            -xs   0   0]
        Hinv = [I transpose(P); zeros(3, 3) I]
        HinvT = [I zeros(3, 3); P I]

        S = Hinv * S * HinvT
    end


    return S, sc, tc
end

function area_and_centroid_of_element(node)

    # shoelace formula for area
    A = 0.5 * (
        node[1].x * node[2].y - node[2].x * node[1].y + 
        node[2].x * node[3].y - node[3].x * node[2].y + 
        node[3].x * node[4].y - node[4].x * node[3].y + 
        node[4].x * node[1].y - node[1].x * node[4].y)

    # centroid of element
    xc = (node[1].x + node[2].x + node[3].x + node[4].x)/4
    yc = (node[1].y + node[2].y + node[3].y + node[4].y)/4

    return A, xc, yc
end


"""
    mass_matrix(nodes, elements)

Compute mass matrix for the structure using GXBeam ordering.
    
**Returns**
- `M::Matrix`: mass matrix
- `mc::Vector{float}`: x, y location of mass center
"""
function mass_matrix(nodes, elements)

    # --- find total mass and center of mass -----
    m = 0.0
    xm = 0.0
    ym = 0.0

    for elem in elements
        # extract nodes and density, compute area and centroid
        node = nodes[elem.nodenum]
        rho = elem.material.rho
        A, xc, yc = area_and_centroid_of_element(node)

        # mass and (numerator of) center of mass
        dm = rho * A
        m += dm
        xm += xc * dm
        ym += yc * dm
    end

    # center of mass
    xm /= m
    ym /= m

    # ----------- compute moments of inertia ---------
    Ixx = 0.0
    Iyy = 0.0
    Ixy = 0.0

    for elem in elements
        # extract nodes and density, compute area and centroid
        node = nodes[elem.nodenum]
        rho = elem.material.rho
        A, xc, yc = area_and_centroid_of_element(node)

        Ixx += (yc - ym)^2 * rho * A
        Iyy += (xc - xm)^2 * rho * A
        Ixy += (xc - xm) * (yc - ym) * rho * A
    end

    M = Symmetric([
        m 0.0 0 0 m*ym -m*xm
        0 m 0 -m*ym 0 0
        0 0 m m*xm 0 0
        0 -m*ym m*xm Ixx+Iyy 0 0
        m*ym 0 0 0 Ixx -Ixy
        -m*xm 0 0 0 -Ixy Iyy
    ])

    return M, [xm, ym]
end

"""
    plotmesh(nodes, elements, pyplot; plotnumbers=false)

plot nodes and elements for a quick visualization.
Need to pass in a PyPlot object as PyPlot is not loaded by this package.
"""
function plotmesh(nodes, elements, pyplot; plotnumbers=false)
    ne = length(elements)
    
    for i = 1:ne
        node = nodes[elements[i].nodenum]
        for i = 1:4
            iplus = i+1
            if iplus == 5
                iplus = 1 
            end
            pyplot.plot([node[i].x, node[iplus].x], [node[i].y, node[iplus].y], "k")
        end
        if plotnumbers
            barx = sum([n.x/4 for n in node])
            bary = sum([n.y/4 for n in node])
            pyplot.text(barx, bary, string(i), color="r")
        end
    end
    if plotnumbers
        nn = length(nodes)
        for i = 1:nn
            pyplot.text(nodes[i].x, nodes[i].y, string(i))
        end
    end
end



