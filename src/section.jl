"""
    Material(E1, E2, E3, G12, G13, G23, nu12, nu13, nu23, rho,
        S1t=1.0, S1c=1.0, S2t=1.0, S2c=1.0, S3t=1.0, S3c=1.0, S12=1.0, S13=1.0, S23=1.0)

Orthotropic material properties. Axis 1 is the main ply axis, axis 2 is the transverse
ply axis, and axis 3 is normal to the ply.  For a fiber orientation of zero, axis 1 is
along the beam axis.

# Arguments
- `Ei::float`: Young's modulus along ith axis.
- `Gij::float`: Shear moduli
- `nuij::float`: Poisson's ratio.  ``nu_ij E_j = nu_ji E_i``
- `rho::float`: Density
- `Sit::float`: Tensile strength in ith direction
- `Sic::float`: Compressive strength in ith direction
- `Sij::float`: Strength in ij direction
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
    S1t::TF
    S1c::TF
    S2t::TF
    S2c::TF
    S3t::TF
    S3c::TF
    S12::TF
    S13::TF
    S23::TF
end

# floating point type
Base.eltype(::Material{TF}) where TF = TF
Base.eltype(::Type{Material{TF}}) where TF = TF

# floating point type conversions
Material{TF}(m::Material) where {TF} = Material{TF}(
    m.E1, m.E2, m.E3, m.G12, m.G13, m.G23, m.nu12, m.nu13, m.nu23, m.rho,
    m.S1t, m.S1c, m.S2t, m.S2c, m.S3t, m.S3c, m.S12, m.S13, m.S23)
Base.convert(::Type{Material{TF}}, m::Material) where {TF} = Material{TF}(m)

# promote all inputs to the same type
function Material(E1, E2, E3, G12, G13, G23, nu12, nu13, nu23, rho,
    S1t, S1c, S2t, S2c, S3t, S3c, S12, S13, S23)

    return Material(promote(E1, E2, E3, G12, G13, G23, nu12, nu13, nu23, rho,
    S1t, S1c, S2t, S2c, S3t, S3c, S12, S13, S23)...)
end

# make strength properties optional
Material(E1, E2, E3, G12, G13, G23, nu12, nu13, nu23, rho) = Material(promote(E1, E2, E3, G12, G13, G23, nu12, nu13, nu23, rho, ones(9)...)...)
Material{TF}(E1, E2, E3, G12, G13, G23, nu12, nu13, nu23, rho) where TF = Material{TF}(E1, E2, E3, G12, G13, G23, nu12, nu13, nu23, rho, ones(TF, 9)...)


"""
    Node(x, y)

A node in the finite element mesh at location x, y.  If assembled in a vector, the vector
index is the node number.

# Arguments
- `x::float`: x-location of node in global coordinate system
- `y::float`: y-location of node in global coordinate system
"""
struct Node{TF}
    x::TF
    y::TF
end

# floating point type
Base.eltype(::Node{TF}) where TF = TF
Base.eltype(::Type{Node{TF}}) where TF = TF

# floating point type conversions
Node{TF}(n::Node) where {TF} = Node{TF}(n.x, n.y)
Base.convert(::Type{Node{TF}}, n::Node) where {TF} = Node{TF}(n)

# promote all inputs to the same type
Node(x, y) = Node(promote(x, y)...)


"""
    MeshElement(nodenum, material, theta)

An element in the mesh, consisting of four ordered nodes, a material, and a fiber orientation.

# Arguments
- `nodenum::Vector{integer}`: Node indices, defined counterclockwise starting from the
    bottom left node (as defined in the local coordinate system, see the figure).
- `material::Material`: Material properties
- `theta::float`: Fiber orientation
"""
struct MeshElement{TF}
    nodenum::Vector{Int}
    material::Material{TF}
    theta::TF
end

# floating point type
Base.eltype(::MeshElement{TF}) where {TF} = TF
Base.eltype(::Type{MeshElement{TF}}) where {TF} = TF

# floating point type conversions
MeshElement{TF}(e::MeshElement) where {TF} = MeshElement{TF}(e.nodenum, e.material, e.theta)
Base.convert(::Type{MeshElement{TF}}, e::MeshElement) where {TF} = MeshElement{TF}(e)

# promote all inputs to the same type
function MeshElement(nodenum, material, theta)
    TF = promote_type(eltype(Material), typeof(theta))
    material = convert(Material{TF}, material)
    theta = convert(TF, theta)
    return MeshElement{TF}(nodenum, material, theta)
end

"""
    SectionCache{TM, TSM, TMF, TSMF, TAF, TSAF}

Internal cache for section properties
"""
struct SectionCache{TM, TSM, TMF, TSMF, TAF, TSAF}  # matrix, sparse matrix, matrix of floats, array of floats
    # node mapping to global matrices
    idx::Vector{Int}
    # global matrices
    A::TM
    R::TM
    E::TSM
    C::TSM
    L::TM
    M::TSM
    # displacements
    X::TM
    Y::TM
    dX::TM
    dY::TM
    # system matrices
    Asys::TSM
    Bsys::TM
    Xsys::TM
    # values
    Av::TSMF
    Bv::TMF
    Xv::TMF
    # derivatives
    Adot::TSAF
    Bdot::TAF
    Xdot::TAF
end

function initialize_cache(nodes, elements)

    TF = promote_type(eltype(nodes[1]), eltype(elements[1]))

    return initialize_cache(TF, nodes, elements)
end

"""
    initialize_cache([TF,] nodes, elements)

Construct the cache for the section property calculations.
"""
function initialize_cache(TF::Type, nodes, elements)

    # 3 displacement dof per node
    ndof = 3 * length(nodes)

    # storage for node mapping to global matrices
    idx = zeros(Int64, 12)

    # global matrices
    A = zeros(TF, 6, 6)  # 6 x 6
    R = zeros(TF, ndof, 6)  #  nn*3 x 6
    E = spzeros(TF, ndof, ndof)  # nn*3 x nn*3
    C = spzeros(TF, ndof, ndof)  # nn*3 x nn*3
    L = zeros(TF, ndof, 6)  # nn*3 x 6
    M = spzeros(TF, ndof, ndof)  # nn*3 x nn*3

    # displacements
    X = zeros(TF, ndof, 6)
    Y = zeros(TF, 6, 6)
    dX = zeros(TF, ndof, 6)
    dY = zeros(TF, 6, 6)

    # system matrices
    Asys = spzeros(TF, ndof+12, ndof+12)
    Bsys = zeros(TF, ndof+12, 6)
    Xsys = zeros(TF, ndof+12, 6)

    # values
    Av = spzeros(ndof+12, ndof+12)
    Bv = zeros(ndof+12, 6)
    Xv = zeros(ndof+12, 6)

    # derivatives
    nderiv = TF <: ForwardDiff.Dual ? ForwardDiff.npartials(TF) : 0
    Adot = zeros(ndof+12, ndof+12, nderiv)
    Bdot = zeros(ndof+12, 6, nderiv)
    Xdot = zeros(ndof+12, 6, nderiv)

    return SectionCache(idx, A, R, E, C, L, M, X, Y, dX, dY, Asys, Bsys, Xsys, Av, Bv, Xv, Adot, Bdot, Xdot)
end


"""
    stiffness(material)

Construct constitutive matrix `Q` for the specified material (with the internal ordering).
"""
function stiffness(material)

    @unpack E1, E2, E3, nu12, nu13, nu23, G12, G13, G23 = material

    nu21 = nu12*E2/E1
    nu31 = nu13*E3/E1
    nu32 = nu23*E3/E2
    delta = 1.0 / (1 - nu12*nu21 - nu23*nu32 - nu13*nu31 - 2*nu21*nu32*nu13)

    Q66 = E1*(1 - nu23*nu32)*delta
    Q11 = E2*(1 - nu13*nu31)*delta
    Q22 = E3*(1 - nu12*nu21)*delta
    Q16 = E1*(nu21 + nu31*nu23)*delta
    Q12 = E2*(nu32 + nu12*nu31)*delta
    Q26 = E1*(nu31 + nu21*nu32)*delta
    Q44 = G12
    Q55 = G13
    Q33 = G23

    Q = @SMatrix [
        Q11 Q12   0   0   0 Q16;
        Q12 Q22   0   0   0 Q26;
        0    0   Q33  0   0   0;
        0    0    0  Q44  0   0;
        0    0    0   0  Q55  0;
        Q16 Q26   0   0   0 Q66;
    ]

    return Q
end


"""
    rotate_ply_to_element(Q, theta)

Rotate constitutive matrix `Q` by ply angle `theta`
"""
function rotate_ply_to_element(Q, theta)

    Ttheta = get_Ttheta(theta)

    return Ttheta * Q * Ttheta'
end

function get_Ttheta(theta)

    s, c = sincos(theta)

    c2 = c^2
    s2 = s^2
    sc = s*c

    Ttheta = @SMatrix [
        c2  0  0  2*sc 0  s2;
        0   1  0   0   0   0;
        0   0  c   0   s   0;
        -sc 0  0 c2-s2 0  sc;
        0   0 -s   0   c   0;
        s2  0  0 -2*sc 0  c2;
           ]

    return Ttheta
end

"""
    rotate_element_to_beam(Q, c, s)

Rotate constitutive matrix `Q` by element orientation `beta` where `c = cos(beta)` and
`s = sin(beta)`
"""
function rotate_element_to_beam(Q, c, s)

    Tbeta = get_Tbeta(c, s)

    return Tbeta * Q * Tbeta'
end

function get_Tbeta(c, s)

    c2 = c^2
    s2 = s^2
    sc = s*c

    Tbeta = @SMatrix [
        c2  s2 -2*sc   0  0  0;
        s2  c2  2*sc   0  0  0;
        sc -sc c2-s2   0  0  0;
        0   0    0     c -s  0;
        0   0    0     s  c  0;
        0   0    0     0  0  1;
    ]

    return Tbeta
end

"""
    elementQ(material, theta, cbeta, sbeta)

Calculate element constitutive matrix `Q` accounting for fiber orientation `theta` and
element orientation `beta` (where `cbeta = cos(beta)` and `sbeta=sin(beta)`)
"""
function elementQ(material, theta, cbeta, sbeta)
    Q = stiffness(material)
    Q = rotate_ply_to_element(Q, theta)
    Q = rotate_element_to_beam(Q, cbeta, sbeta)
    return Q
end


"""
    element_orientation(nodes)

Calculate the orientation `beta` of a single element using its nodes.  Return `cos(beta)`
and `sin(beta)`
"""
function element_orientation(nodes)

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

    return cbeta, sbeta
end

function element_matrices(ksi, eta, element, nodes)

    # shape functions
    mksi = 1 - ksi
    pksi = 1 + ksi
    meta = 1 - eta
    peta = 1 + eta

    N = SVector(
        mksi*meta/4,
        pksi*meta/4,
        pksi*peta/4,
        mksi*peta/4,
    )

    # x, y position
    x = 0.0
    y = 0.0
    for i = 1:4
        x += N[i]*nodes[i].x
        y += N[i]*nodes[i].y
    end

    # element orientation
    cbeta, sbeta = element_orientation(nodes)

    # element constitutive matrix
    Q = elementQ(element.material, element.theta, cbeta, sbeta)

    # translation and cross product
    SZ = @SMatrix [
        0  0  0  0   0   0
        0  0  0  0   0   0
        0  0  0  0   0   0
        1  0  0  0   0  -y
        0  1  0  0   0   x
        0  0  1  y  -x   0
    ]

    SN = @SMatrix [
        0     0    0    0    0    0    0    0    0    0    0      0;
        0     0    0    0    0    0    0    0    0    0    0      0;
        0     0    0    0    0    0    0    0    0    0    0      0;
        N[1]  0    0   N[2]  0    0   N[3]  0    0   N[4]  0      0;
        0    N[1]  0    0   N[2]  0    0   N[3]  0    0   N[4]    0;
        0     0   N[1]  0    0   N[2]  0    0   N[3]  0    0   N[4];
    ]

    # --- Jacobian --- #

    # derivative of `N` wrt `ksi`
    N_ksi = SVector(-meta/4, meta/4, peta/4, -peta/4)

    # derivative of `N` wrt `eta`
    N_eta = SVector(-mksi/4, -pksi/4, pksi/4, mksi/4)

    # derivative of `x, y` wrt `ksi, eta`
    x_ksi = 0.0
    x_eta = 0.0
    y_ksi = 0.0
    y_eta = 0.0
    for i = 1:4
        x_ksi += N_ksi[i]*nodes[i].x
        x_eta += N_eta[i]*nodes[i].x
        y_ksi += N_ksi[i]*nodes[i].y
        y_eta += N_eta[i]*nodes[i].y
    end

    J = @SMatrix [
        x_ksi y_ksi;
        x_eta y_eta;
    ]

    Jinv = inv(J)

    Bksi = @SMatrix [
        Jinv[1, 1]      0              0;
        0          Jinv[2, 1]          0;
        Jinv[2, 1] Jinv[1, 1]          0;
        0               0     Jinv[1, 1];
        0               0     Jinv[2, 1];
        0               0              0;
    ]

    Beta = @SMatrix [
        Jinv[1, 2]      0              0;
        0          Jinv[2, 2]          0;
        Jinv[2, 2] Jinv[1, 2]          0;
        0               0     Jinv[1, 2];
        0               0     Jinv[2, 2];
        0               0              0;
    ]

    # dNM_dksi = [N_ksi[1]*Matrix(1.0I, 3, 3) N_ksi[2]*I N_ksi[3]*I N_ksi[4]*I]
    NM_ksi = @SMatrix [
        N_ksi[1] 0 0 N_ksi[2] 0 0 N_ksi[3] 0 0 N_ksi[4] 0 0;
        0 N_ksi[1] 0 0 N_ksi[2] 0 0 N_ksi[3] 0 0 N_ksi[4] 0;
        0 0 N_ksi[1] 0 0 N_ksi[2] 0 0 N_ksi[3] 0 0 N_ksi[4];
    ]

    # NM_eta = [N_eta[1]*Matrix(1.0I, 3, 3) N_eta[2]*I N_eta[3]*I N_eta[4]*I]
    NM_eta = @SMatrix [
        N_eta[1] 0 0 N_eta[2] 0 0 N_eta[3] 0 0 N_eta[4] 0 0;
        0 N_eta[1] 0 0 N_eta[2] 0 0 N_eta[3] 0 0 N_eta[4] 0;
        0 0 N_eta[1] 0 0 N_eta[2] 0 0 N_eta[3] 0 0 N_eta[4];
    ]

    BN = Bksi * NM_ksi + Beta * NM_eta

    return Q, J, SZ, SN, BN
end

"""
    element_integrand(ksi, eta, element, nodes)

Compute the integrand for a single element with a given `ksi`, `eta`.
"""
function element_integrand(ksi, eta, element, nodes)

    Q, J, SZ, SN, BN = element_matrices(ksi, eta, element, nodes)

    # integrands
    detJ = det(J)

    Ae = SZ' * Q * SZ * detJ
    Re = BN' * Q * SZ * detJ
    Ee = BN' * Q * BN * detJ
    Ce = BN' * Q * SN * detJ  # use Giavottoa definition
    Le = SN' * Q * SZ * detJ
    Me = SN' * Q * SN * detJ

    return Ae, Re, Ee, Ce, Le, Me
end

"""
    element_submatrix(element, nodes)

Performs 2-point Gauss quadrature of one element
"""
function element_submatrix(element, nodes)

    Ae1, Re1, Ee1, Ce1, Le1, Me1 = element_integrand(1.0/sqrt(3), 1.0/sqrt(3), element, nodes)
    Ae2, Re2, Ee2, Ce2, Le2, Me2 = element_integrand(1.0/sqrt(3), -1.0/sqrt(3), element, nodes)
    Ae3, Re3, Ee3, Ce3, Le3, Me3 = element_integrand(-1.0/sqrt(3), 1.0/sqrt(3), element, nodes)
    Ae4, Re4, Ee4, Ce4, Le4, Me4 = element_integrand(-1.0/sqrt(3), -1.0/sqrt(3), element, nodes)

    Ae = Ae1 + Ae2 + Ae3 + Ae4
    Re = Re1 + Re2 + Re3 + Re4
    Ee = Ee1 + Ee2 + Ee3 + Ee4
    Ce = Ce1 + Ce2 + Ce3 + Ce4
    Le = Le1 + Le2 + Le3 + Le4
    Me = Me1 + Me2 + Me3 + Me4

    return Ae, Re, Ee, Ce, Le, Me
end

"""
    node2idx(nodenums)

Convenience function to map node numbers to locations in global matrix.
"""
function node2idx(nodenums)

    idx = Vector{Int64}(undef, length(nodenums)*3)

    return node2idx!(idx, nodenums)
end

"""
    node2idx!(idx, nodenums)

In-place version of [`node2idx`](@ref).
"""
function node2idx!(idx, nodenums)

    nn = length(nodenums)

    for i = 1:nn
        i1 = (i-1)*3+1
        i2 = i*3
        j1 = (nodenums[i]-1)*3+1
        j2 = nodenums[i]*3
        idx[i1:i2] = j1:j2
    end

    return idx
end

"""
    reorder(K)

Reorder stiffness or compliance matrix from internal order to GXBeam order
"""
function reorder(K)  # reorder to GXBeam format

    idx = [3, 1, 2, 6, 4, 5]

    return K[idx, idx]
end

"""
    linearsolve(A, B; kwargs...)

Linear solve which is overloaded with analytic derivatives.  Returns the factorized matrix
`Afact` and the result of the linear solve `X`.
"""
function linearsolve(A, B; Afact=factorize(A), kwargs...)

    X = Afact \ B

    return Afact, X
end

# overloaded for use with ForwardDiff
function linearsolve(A::AbstractMatrix{TF}, B; X = similar(B, TF), Afact=nothing,
    Av = similar(A, ForwardDiff.value(TF)), Bv = similar(B, ForwardDiff.value(TF)),
    Xv = similar(X, ForwardDiff.value(TF)), Adot = similar(A, size(A, 1), size(A, 2), ForwardDiff.npartials(TF)),
    Bdot = similar(B, size(B, 1), size(B, 2), ForwardDiff.npartials(TF)),
    Xdot = similar(X, size(X, 1), size(X, 2), ForwardDiff.npartials(TF))) where {TF<:ForwardDiff.Dual}

    # extract primal values
    Av .= ForwardDiff.value.(A)
    Bv .= ForwardDiff.value.(B)

    # compute factorization (if not provided)
    if isnothing(Afact)
        Afact = factorize(Av)
    end

    # perform linear solve
    Xv = Afact \ Bv

    # extract partials
    m, n = size(A)
    for i = 1:m
        for j = 1:n
            Adot[i, j, :] .= ForwardDiff.partials(A[i, j]).values
        end
    end

    m, n = size(B)
    for i = 1:m
        for j = 1:n
            Bdot[i, j, :] .= ForwardDiff.partials(B[i, j]).values
        end
    end

    nderiv = ForwardDiff.npartials(TF)
    for k = 1:nderiv
        # set right hand side to Bdot
        rhs = view(Bdot, :, :, k)
        # subtract Adot*Xv
        mul!(rhs, view(Adot, :, :, k), Xv, -1, 1)
        # calculate Xdot
        Xdot[:,:,k] .= Afact \ rhs
    end

    # repack in dual
    m, n = size(Xv)
    for i = 1:m
        for j = 1:n
            X[i,j] = TF(Xv[i,j], ForwardDiff.Partials(ntuple(k->Xdot[i,j,k], nderiv)))
        end
    end

    return Afact, X
end


"""
    compliance_matrix(nodes, elements; cache=initialize_cache(nodes, elements),
        gxbeam_order=true, shear_center=true)

Compute compliance matrix given a finite element mesh described by nodes and elements.

# Arguments
- `nodes`: Vector containing all the nodes in the mesh
- `elements`: Vector containing all the elements in the mesh
- `cache::SectionCache`: A pre-allocated cache which may be passed in to reduce allocations
    across multiple calls to this function when the number of nodes, number of elements, and
    connectivity remain the same.
- `gxbeam_order::Bool`: Indicates whether the compliance matrix should be provided in the
    order expected by GXBeam (rather than the internal ordering used by the section analysis)
- `shear_center::Bool`: Indicates whether the compliance matrix should be provided about the
    shear center

# Returns
- `S::Matrix`: compliance matrix
- `sc::Vector{float}`: x, y location of shear center (location where a transverse/shear
    force will not produce any torsion, i.e., beam will not twist)
- `tc::Vector{float}`: x, y location of tension center, aka elastic center, aka centroid
    (location where an axial force will not produce any bending, i.e., beam will remain
    straight)
"""
function compliance_matrix(nodes, elements; cache=initialize_cache(nodes, elements),
    gxbeam_order=true, shear_center=gxbeam_order)

    # problem dimensions
    ne = length(elements) # number of elements
    nn = length(nodes)    # number of nodes
    ndof = 3 * nn         # 3 displacement dof per node

    # reset global matrices
    A = cache.A .= 0
    R = cache.R .= 0
    E = cache.E .= 0
    C = cache.C .= 0
    L = cache.L .= 0
    M = cache.M .= 0

    # storage for element node indices
    idx = cache.idx

    # place element matrices in global matrices (scatter)
    for i = 1:ne
        # element nodes
        nodenum = elements[i].nodenum
        # location in global matrix
        node2idx!(idx, nodenum)
        # element submatrix
        Ae, Re, Ee, Ce, Le, Me = element_submatrix(elements[i], nodes[nodenum])
        # add to global matrix
        cache.A .+= Ae
        @views cache.R[idx, :] .+= Re
        @views cache.E[idx, idx] .+= Ee
        @views cache.C[idx, idx] .+= Ce
        @views cache.L[idx, :] .+= Le
        @views cache.M[idx, idx] .+= Me
    end

    # --- Construct and Solve First Linear System --- #

    # construct left hand side: Asys = [E R D; R' A 0 0; D' 0 0]
    Asys = cache.Asys .= 0
    Asys[1:ndof, 1:ndof] = E
    Asys[1:ndof, ndof+1:ndof+6] = R
    Asys[ndof+1:ndof+6, 1:ndof] = R'
    Asys[ndof+1:ndof+6, ndof+1:ndof+6] = A
    for i = 1:nn
        s = 3*(i-1)
        Asys[ndof+7, s+1] = Asys[s+1, ndof+7] = 1.0
        Asys[ndof+8, s+2] = Asys[s+2, ndof+8] = 1.0
        Asys[ndof+9, s+3] = Asys[s+3, ndof+9] = 1.0
        Asys[ndof+10, s+3] = Asys[s+3, ndof+10] = nodes[i].y
        Asys[ndof+11, s+3] = Asys[s+3, ndof+11] = -nodes[i].x
        Asys[ndof+12, s+1] = Asys[s+1, ndof+12] = -nodes[i].y
        Asys[ndof+12, s+2] = Asys[s+2, ndof+12] = nodes[i].x
    end

    # construct right hand side: Tr=zeros(6,6); Tr[1,5]=-1; Tr[2,4]=1; B2 = [0, Tr', 0];
    B2 = cache.Bsys .= 0
    B2[ndof+5, 1] = -1.0
    B2[ndof+4, 2] =  1.0

    # # solve linear system
    Afact, X2 = linearsolve(Asys, B2; X = cache.Xsys, Av = cache.Av, Bv = cache.Bv,
        Xv = cache.Xv, Adot = cache.Adot, Bdot = cache.Bdot, Xdot = cache.Xdot)

    # X2 = hcat([ImplicitAD.implicit_linear(Asys, B2[:,k]) for k = 1:6]...)

    # extract and save results
    dX = cache.dX .= view(X2, 1:ndof, :)
    dY = cache.dY .= view(X2, ndof+1:ndof+6, :)

    # --- Construct and Solve Second Linear System --- #

    # construct right hand side: B1 = [C'-C  L; -L' 0; 0 0]*[dX, dY] + [0, I, 0]
    # (note that there are a couple errors in the BECAS theory guide for this expression)

    B1 = cache.Bsys .= 0

    mul!(view(B1, 1:ndof, :), C', dX, 1, 1)
    mul!(view(B1, 1:ndof, :), C, dX, -1, 1)
    mul!(view(B1, 1:ndof, :), L, dY, 1, 1)

    B1[ndof+1, 1] = 1
    B1[ndof+2, 2] = 1
    B1[ndof+3, 3] = 1
    B1[ndof+4, 4] = 1
    B1[ndof+5, 5] = 1
    B1[ndof+6, 6] = 1
    mul!(view(B1, ndof+1:ndof+6, :), L', dX, -1, 1)

    # solve linear system
    _, X1 = linearsolve(Asys, B1; Afact = Afact, X = cache.Xsys, Av = cache.Av, Bv = cache.Bv,
        Xv = cache.Xv, Adot = cache.Adot, Bdot = cache.Bdot, Xdot = cache.Xdot)

    # X1 = hcat([ImplicitAD.implicit_linear(Asys, B1[:,k]) for k = 1:6]...)

    # extract and save results
    X = cache.X .= view(X1, 1:ndof, :)
    Y = cache.Y .= view(X1, ndof+1:ndof+6, :)

    # --- Compute Compliance Matrix --- #

    # The following is equivalent to the following matrix operation
    # S = hcat(X', dX', Y')*[E C R; C' M L; R' L' A]*vcat(X, dX, Y)

    # use state vector as temporary storage
    tmp1 = view(X1, 1:ndof, :)
    tmp2 = view(X1, ndof+1:ndof+6, :)

    S = zero(A)
    mul!(S, X', mul!(tmp1, E, X), 1, 1)
    mul!(S, X', mul!(tmp1, C, dX), 1, 1)
    mul!(S, mul!(tmp2, X', R), Y, 1, 1)
    mul!(S, dX', mul!(tmp1, C', X), 1, 1)
    mul!(S, dX', mul!(tmp1, M, dX), 1, 1)
    mul!(S, mul!(tmp2, dX', L), Y, 1, 1)
    mul!(S, Y', mul!(tmp2, R', X), 1, 1)
    mul!(S, Y', mul!(tmp2, L', dX), 1, 1)
    mul!(S, Y', mul!(tmp2, A, Y), 1, 1)

    # --- Find Shear and Tension Center --- #

    xs = -S[6, 2]/S[6, 6]
    ys = S[6, 1]/S[6, 6]
    xt = (S[4, 4]*S[5, 3] - S[4, 5]*S[4, 3])/(S[4, 4]*S[5, 5] - S[4, 5]^2)
    yt = (-S[4, 3]*S[5, 5] + S[4, 5]*S[5, 3])/(S[4, 4]*S[5, 5] - S[4, 5]^2)
    sc = [xs, ys]
    tc = [xt, yt]

    # --- Change Ordering to Match GXBeam --- #

    if shear_center
        P = [0 0 ys; 0 0 -xs; -ys xs 0]
        Hinv = [I transpose(P); zeros(3, 3) I]
        HinvT = [I zeros(3, 3); P I]
        S = Hinv * S * HinvT
    end

    if gxbeam_order
        S = reorder(S)
    end

    return S, sc, tc
end

"""
    area_and_centroid_of_element(node)

Compute area and centroid of an element specified by its four nodes
"""
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

# Returns
- `M::Matrix`: mass matrix
- `mc::Vector{float}`: x, y location of mass center
"""
function mass_matrix(nodes, elements)

    # --- find total mass and center of mass -----
    m = 0.0
    xm = 0.0
    ym = 0.0

    for element in elements
        # extract nodes and density, compute area and centroid
        node = nodes[element.nodenum]
        rho = element.material.rho
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

    for element in elements
        # extract nodes and density, compute area and centroid
        node = nodes[element.nodenum]
        rho = element.material.rho
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


"""
    strain_recovery(F, M, nodes, elements, cache; gxbeam_order=false)

Compute stresses and strains at each element in cross section.

# Arguments
- `F::Vector(3)`: force at this cross section in x, y, z directions
- `M::Vector(3)`: moment at this cross section in x, y, z directions
- `nodes::Vector{Node{TF}}`: all the nodes in the mesh
- `elements::Vector{MeshElement{TF}}`: all the elements in the mesh
- `cache::SectionCache`: needs to reuse data from the compliance solve
    (thus must initialize cache and pass it to both compliance and this function)

# Keyword Arguments
- `gxbeam_order=false::Bool`: if true, `F`` and `M` are assumed to be in the 
    local beam axis used by GXBeam (where the beam extends along the x-axis). This
    also returns beam stresses and strains in the axis order set by GXBeam 
    (e.g. axial stresses would correspond to the `xx` direction, or first index).

# Returns
- `strain_b::Vector(6, ne)`: strains in beam coordinate system for each element. order: xx, yy, zz, xy, xz, yz
    Note: this order (as well as those below) corresponds to the local beam reference frame if `gxbeam_order` 
    is set to `true`.
- `stress_b::Vector(6, ne)`: stresses in beam coordinate system for each element. order: xx, yy, zz, xy, xz, yz
- `strain_p::Vector(6, ne)`: strains in ply coordinate system for each element. order: 11, 22, 33, 12, 13, 23
- `stress_p::Vector(6, ne)`: stresses in ply coordinate system for each element. order: 11, 22, 33, 12, 13, 23
"""
function strain_recovery(F, M, nodes, elements, cache; gxbeam_order=false)

    # initialize outputs
    T = promote_type(eltype(F), eltype(M))
    ne = length(elements)
    epsilon_b = Matrix{T}(undef, 6, ne)
    sigma_b = Matrix{T}(undef, 6, ne)
    epsilon_p = Matrix{T}(undef, 6, ne)
    sigma_p = Matrix{T}(undef, 6, ne)

    # concatenate forces/moments
    if gxbeam_order
        new_idxs = [2,3,1]
        theta = vcat(SVector{3}(F[new_idxs]), SVector{3}(M[new_idxs])) # convert input loads from GXBeam local frame to stress recovery local frame
    else
    theta = vcat(SVector{3}(F), SVector{3}(M))
    end

    # indices into global matrices
    idx = cache.idx

    # save reordering index
    if gxbeam_order
        idx_b = [6, 1, 2, 5, 3, 4]   # zz, xx, yy, yz, xy, xz
    else
        idx_b = [1, 2, 6, 3, 4, 5]   # xx, yy, zz, xy, xz, yz
    end
    idx_p = [6, 1, 2, 4, 5, 3]   # 11, 22, 33, 12, 13, 23

    # iterate over elements
    @views for i = 1:ne

        # analyze this element
        element = elements[i]
        nodenum = element.nodenum
        material = element.material
        element_nodes = nodes[nodenum]

        # compute submatrices SZ, BN, SN (evaluated at center of element)
        Q, J, SZ, SN, BN = element_matrices(0.0, 0.0, element, element_nodes)

        # extract part of solution corresponding to this element
        node2idx!(idx, nodenum)
        Xe = SMatrix{12,6}(view(cache.X, idx, :))
        dXe = SMatrix{12,6}(view(cache.dX, idx, :))
        Y = SMatrix{6,6}(cache.Y)

        # element strains in beam reference frame
        epsilon_b[:, i] .= SZ*Y*theta + BN*Xe*theta + SN*dXe*theta

        # corresponding stress
        sigma_b[:, i] .= Q * epsilon_b[:, i]

        # element strains in ply reference frame
        cbeta, sbeta = element_orientation(element_nodes)
        Ttheta = get_Ttheta(element.theta)
        Tbeta = get_Tbeta(cbeta, sbeta)
        epsilon_p[:, i] .= Ttheta' * Tbeta' * epsilon_b[:, i]

        # corresponding stresses
        Q = stiffness(material)
        sigma_p[:, i] .= Q * epsilon_p[:, i]

        # reorder to a more conventional order
        epsilon_b[:, i] .= epsilon_b[idx_b, i]
        sigma_b[:, i] .= sigma_b[idx_b, i]
        epsilon_p[:, i] .= epsilon_p[idx_p, i]
        sigma_p[:, i] .= sigma_p[idx_p, i]
    end

    return epsilon_b, sigma_b, epsilon_p, sigma_p
end


"""
    plotsoln(nodes, elements, soln, pyplot)

plot stress/strain on mesh
soln could be any vector that is of length # of elements, e.g., sigma_b[3, :]
Need to pass in a PyPlot object as PyPlot is not loaded by this package.
"""
function plotsoln(nodes, elements, soln, pyplot)
    ne = length(elements)
    nn = length(nodes)

    # extract node points
    xpts = zeros(nn)
    ypts = zeros(nn)
    for i = 1:nn
        xpts[i] = nodes[i].x
        ypts[i] = nodes[i].y
    end

    # split quads into trianagles
    triangles = zeros(Int64, ne*2, 3)
    trisol = zeros(ne*2)

    for i = 1:ne
        nnum = elements[i].nodenum

        triangles[i*2-1, :] = nnum[1:3] .- 1
        triangles[i*2, :] = [nnum[1], nnum[3], nnum[4]] .- 1

        # same solution on both triangles (same quad element)
        trisol[2*i-1] = soln[i]
        trisol[2*i] = soln[i]
    end

    pyplot.tripcolor(xpts, ypts, trisol, triangles=triangles)
end


# function tsai_hill(sigma, strength)

#     (; S1t, S1c, S2t, S2c, S3t, S3c, S12, S13, S23) = strength

#     _, ne = size(sigma)
#     T = eltype(sigma)
#     failure = Vector{T}(undef, n)  # fails if > 1
#     s = Vector{T}(undef, 6)

#     for i = 1:ne
#         s .= sigma[:, i]

#         if s[1] >= 0.0
#             S1 = S1t
#         else
#             S1 = S1c
#         end
#         if s[2] >= 0.0
#             S2 = S2t
#         else
#             S2 = S2c
#         end
#         failure[i] = s[1]^2/S1^2 + s[2]^2/S2^2 + s[4]^2/S12^2 - s[1]*s[2]/S1^2
#     end

#     return failure
# end

"""
    tsai_wu(stress_p, elements)

Tsai Wu failure criteria

# Arguments
- `stress_p::vector(6, ne)`: stresses in ply coordinate system
- `elements::Vector{MeshElement{TF}}`: all the elements in the mesh

# Returns
- `failure::vector(ne)`: tsai-wu failure criteria for each element.  fails if >= 1
"""
function tsai_wu(stress_p, elements)

    ne = length(elements)
    T = eltype(stress_p)
    failure = Vector{T}(undef, ne)  # fails if > 1
    s = Vector{T}(undef, 6)

    @views for i = 1:ne
        m = elements[i].material
        s .= stress_p[:, i]
        failure[i] = s[1]^2/(m.S1t*m.S1c) +
                     s[2]^2/(m.S2t*m.S2c) +
                     s[3]^2/(m.S3t*m.S3c) +
                     s[4]^2/m.S12^2 +
                     s[5]^2/m.S13^2 +
                     s[6]^2/m.S23^2 +
                     s[1]*(1/m.S1t - 1/m.S1c) +
                     s[2]*(1/m.S2t - 1/m.S2c) +
                     s[3]*(1/m.S3t - 1/m.S3c) -
                     s[1]*s[2]/sqrt(m.S1t*m.S1c*m.S2t*m.S2c) -
                     s[1]*s[3]/sqrt(m.S1t*m.S1c*m.S3t*m.S3c) -
                     s[2]*s[3]/sqrt(m.S2t*m.S2c*m.S3t*m.S3c)
    end

    return failure
end
