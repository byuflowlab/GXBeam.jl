using LinearAlgebra: I, Symmetric, det, factorize
using SparseArrays: spzeros, sparse

"""
    Material(E1, E2, E2, nu12, nu13, nu23, G12, G13, G23, rho)

General orthotropic material properties. 
The "1" direction is along the beam axis for a fiber orientation of zero (theta=0).
"2" corresponds to the local x-direction and "3" the local y at theta=0. (see documentation for figures).  
If the two in-plane directions have the same stiffness properties then only one 
needs to be specified (e.g., compatible with the plane-stress assumption of CLT-based tools):
`Material(E1, E2, nu12, G12, rho)`

**Arguments**
- `E::float`: Young's modulus along 1st, 2nd and 3rd axes.
- `nu::float`: Poisson's ratio.  nu_ij E_j = nu_ji E_i
- `G::float`: shear moduli
- `rho::float`: density
"""
struct Material{TF}
    E1::TF
    E2::TF
    E3::TF
    nu12::TF
    nu13::TF
    nu23::TF
    G12::TF
    G13::TF
    G23::TF
    rho::TF
end

Material(E1, E2, nu12, G12, rho) = Material(E1, E2, E2, nu12, nu12, nu12, G12, G12, G12, rho)

"""
    Node(x, y, number)

A node in the finite element mesh at location x, y with a given index number.

**Arguments**
- `x::float`: x location of node in global coordinate system
- `y::float`: y location of node in global coordinate system
- `number::integer`: unique index of this node, elements use this indices
"""
struct Node{TF, TI}
    x::TF
    y::TF
    number::TI
end

"""
    Element(nodenum, material, theta)

An element in the mesh, consisting of four ordered nodes, a material, and a fiber orientation.

**Arguments**
- `nodenum::Vector{integer}`: a vector of four node numbers corresponding the the four nodes defining this element. 
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
Constituitive matrix of this material using the internal ordering.
"""
function stiffness(material) 
    E1 = material.E1; E2 = material.E2; E3 = material.E3
    nu12 = material.nu12; nu13 = material.nu13; nu23 = material.nu23
    G12 = material.G12; G13 = material.G13; G23 = material.G23

    nu21 = nu12*E2/E1
    nu31 = nu13*E3/E1
    nu32 = nu23*E3/E2
    delta = 1.0 / (1 - nu12*nu21 - nu23*nu32 - nu13*nu31 - 2*nu21*nu32*nu13)

    Q = zeros(eltype(E1), 6, 6)
    Q[6, 6] = E1*(1 - nu23*nu32)*delta
    Q[1, 1] = E2*(1 - nu13*nu31)*delta
    Q[2, 2] = E3*(1 - nu12*nu21)*delta
    Q[1, 6] = E1*(nu21 + nu31*nu23)*delta
    Q[2, 6] = E1*(nu31 + nu21*nu32)*delta
    Q[1, 2] = E2*(nu32 + nu12*nu31)*delta
    Q[4, 4] = G12
    Q[5, 5] = G13
    Q[3, 3] = G23

    return Symmetric(Q)
end

"""
Rotate constituitive matrix by ply angle
"""
function rotate_to_ply(K, theta)
    c = cos(theta)
    s = sin(theta)
    T = [c^2 0 0 2*s*c 0 s^2;
        0.0 1 0 0 0 0;
        0 0 c 0 s 0;
        -s*c 0 0 c^2-s^2 0 s*c;
        0 0 -s 0 c 0;
        s^2 0 0 -2*s*c 0 c^2]

    return T*K*T'
end

"""
Rotate constituitive matrix by element orientation where `c = cos(beta)` and `s = sin(beta)`
"""
function rotate_to_element(K, c, s)  # c = cos(beta), s = sin(beta)
    T = [c^2 s^2 -2*s*c 0 0 0;
          s^2 c^2 2*s*c 0 0 0;
          s*c -s*c c^2-s^2 0 0 0;
          0 0 0 c -s 0;
          0 0 0 s c 0;
          0.0 0 0 0 0 1]

    return T*K*T'
end


"""
Get element constituitive matrix accounting for fiber orientation and element orientation
"""
function elementQ(material, theta, cbeta, sbeta)
    Q1 = stiffness(material)  # material
    Q2 = rotate_to_ply(Q1, theta)
    Q3 = rotate_to_element(Q2, cbeta, sbeta)

    return Q3
end

"""
Compute the integrand for a single element with a given ksi, eta.
"""
function elementintegrand(ksi, eta, element, nodes)

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
    Z = [I [0.0 0 -y; 0 0 x; y -x 0]]  # translation and cross product
    S = [zeros(3, 3); I]
    Q = elementQ(element.material, element.theta, cbeta, sbeta)
    N = [N[1]*Matrix(1.0I, 3, 3) N[2]*I N[3]*I N[4]*I]
    SZ = S*Z
    SN = S*N

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
    Bksi = [Jinv[1, 1] 0 0;
            0 Jinv[2, 1] 0;
            Jinv[2, 1] Jinv[1, 1] 0;
            0 0 Jinv[1, 1];
            0 0 Jinv[2, 1];
            0 0 0]
    Beta = [Jinv[1, 2] 0 0;
            0 Jinv[2, 2] 0;
            Jinv[2, 2] Jinv[1, 2] 0;
            0 0 Jinv[1, 2];
            0 0 Jinv[2, 2];
            0 0 0]

    dNM_dksi = [dN_dksi[1]*Matrix(1.0I, 3, 3) dN_dksi[2]*I dN_dksi[3]*I dN_dksi[4]*I]
    dNM_deta = [dN_deta[1]*Matrix(1.0I, 3, 3) dN_deta[2]*I dN_deta[3]*I dN_deta[4]*I]

    BN = Bksi*dNM_dksi + Beta*dNM_deta
    
    # integrands
    A = SZ'*Q*SZ*detJ
    R = BN'*Q*SZ*detJ
    E = BN'*Q*BN*detJ
    # C = SN'*Q*BN*detJ  # NOTE: error in implementation sectiono of manual, theory is correct
    C = BN'*Q*SN*detJ  # use Giavottoa definition
    L = SN'*Q*SZ*detJ
    M = SN'*Q*SN*detJ

    return A, R, E, C, L, M
end

"""
2-point Gauss quadrature of one element
"""
function submatrix(element, nodes)

    # 2-point Gauss quadrature in both dimensions
    ksi = [1.0/sqrt(3), 1.0/sqrt(3), -1.0/sqrt(3), -1.0/sqrt(3)]
    eta = [1.0/sqrt(3), -1.0/sqrt(3), 1.0/sqrt(3), -1.0/sqrt(3)]

    A, R, E, C, L, M = elementintegrand(ksi[1], eta[1], element, nodes)
    for i = 2:4
        Ai, Ri, Ei, Ci, Li, Mi = elementintegrand(ksi[i], eta[i], element, nodes)
        A += Ai; R += Ri; E += Ei
        C += Ci; L += Li; M += Mi
    end
    
    return A, R, E, C, L, M
end

"""
Convenience function to map node numbers to locations in global matrix.
"""
function node2idx(nodes)
    nn = 4
    idx = Vector{Int64}(undef, nn*3)
    for i = 1:nn
        idx[((i-1)*3+1):i*3] = ((nodes[i]-1)*3+1):nodes[i]*3
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
function compliance(nodes, elements)

    # initialize
    ne = length(elements) # number of elements
    nn = length(nodes)  # number of nodes
    ndof = 3 * nn  # 3 displacement dof per node
    
    # TODO: benchmark sparsity

    # initialize
    etype = eltype(elements[1].theta)  # TODO add others
    A = zeros(etype, 6, 6)  # 6 x 6
    R = zeros(etype, ndof, 6)  #  nn*3 x 6
    E = spzeros(etype, ndof, ndof)  # nn*3 x nn*3
    C = spzeros(etype, ndof, ndof)  # nn*3 x nn*3
    L = zeros(etype, ndof, 6)  # nn*3 x 6
    M = spzeros(etype, ndof, ndof)  # nn*3 x nn*3

    # place element matrices in global matrices (scatter)
    for i = 1:ne
        nodenum = elements[i].nodenum
        Ae, Re, Ee, Ce, Le, Me = submatrix(elements[i], nodes[nodenum])
        idx = node2idx(nodenum)
        A += Ae
        R[idx, :] += Re
        E[idx, idx] += Ee
        C[idx, idx] += Ce
        L[idx, :] += Le
        M[idx, idx] += Me
    end

    # assemble displacement constraint matrix
    DT = spzeros(6, ndof)
    for i = 1:nn
        k = nodes[i].number    
        s = 3*(k-1)
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


"""
    Layer(material, t, theta)

A layer (could be one ply or many plys of same material).
A layup is a vector of layers.

**Arguments**
- `material::Material`: corresponding material
- `t::float`: thickness of ply 
- `theta::float`: fiber orientation (rad)
"""
struct Layer{TF}
    material::Material
    t::TF
    theta::TF
end



# --------- Unit Tests --------
using Test

E1 = rand()
E2 = rand()
E3 = rand()
nu12 = rand()
nu23 = rand()
nu13 = rand()
G12 = rand()
G13 = rand()
G23 = rand()
rho = 1.0
mat = Material(E1, E2, E3, nu12, nu13, nu23, G12, G13, G23, rho)
Q1 = stiffness(mat)

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
iso1 = Material(100.0, 100.0, 100.0, 0.2, 0.2, 0.2, 41.667, 41.667, 41.667, 1.0)
x = range(-0.05, 0.05, length=11)
y = range(-0.05, 0.05, length=11)

nodes = Vector{Node}(undef, 11*11)
elements = Vector{Element}(undef, 10*10)

let
m = 1
for i = 1:11
    for j = 1:11
        nodes[m] = Node(x[i], y[j], m)
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

# ne, _ = size(connectivity)
ne = length(elements)
figure()
for i = 1:ne
    # node = nodes[connectivity[i, :]]
    node = nodes[elements[i].nodenum]
    for i = 1:4
        text(node[i].x, node[i].y, string(node[i].number))
        iplus = i+1
        if iplus == 5
            iplus = 1 
        end
        plot([node[i].x, node[iplus].x], [node[i].y, node[iplus].y])
    end
    barx = sum([n.x/4 for n in node])
    bary = sum([n.y/4 for n in node])
    text(barx, bary, string(i))
end


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
iso2 = Material(100.0/alpha, 100.0/alpha, 0.2, 41.667/alpha, 1.0)  # note error in user guide for nu

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
ortho = Material(480.0, 120.0, 120.0, 0.19, 0.26, 0.19, 60.0, 50.0, 60.0, 1.0)

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
        nodes[m] = Node(r[j]*cos(theta[i]), r[j]*sin(theta[i]), m)
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

ne = length(elements)
nn = length(nodes)
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
    barx = sum([n.x/4 for n in node])
    bary = sum([n.y/4 for n in node])
    # text(barx, bary, string(i), color="r")
end
# for i = 1:nn
#     text(nodes[i].x, nodes[i].y, string(nodes[i].number))
# end

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
        nodes[m] = Node(r[j]*cos(theta[i]), r[j]*sin(theta[i]), m)
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

ne = length(elements)
nn = length(nodes)
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
    barx = sum([n.x/4 for n in node])
    bary = sum([n.y/4 for n in node])
    # text(barx, bary, string(i), color="r")
end
# for i = 1:nn
#     text(nodes[i].x, nodes[i].y, string(nodes[i].number))
# end

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
circmat = Material(E, E, nu, G, 1.0)

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
        nodes[m] = Node(r[j]*cos(theta[i]), r[j]*sin(theta[i]), m)
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

ne = length(elements)
nn = length(nodes)
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
    barx = sum([n.x/4 for n in node])
    bary = sum([n.y/4 for n in node])
    # text(barx, bary, string(i), color="r")
end
# for i = 1:nn
#     text(nodes[i].x, nodes[i].y, string(nodes[i].number))
# end

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
pipemat = Material(E1, E2, nu12, G12, rho)

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
            nodes[m] = Node(x[i], y[j], m)
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
            nodes[m] = Node(x2 + r[j]*cos(theta[i]), r[j]*sin(theta[i]), m)
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
            nodes[m] = Node(x[i], y[j], m)
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
            nodes[m] = Node(x2 + r[j]*cos(theta[i]), r[j]*sin(theta[i]), m)
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

ne = length(elements)
nn = length(nodes)
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
    barx = sum([n.x/4 for n in node])
    bary = sum([n.y/4 for n in node])
    # text(barx, bary, string(i), color="r")
end
# for i = 1:nn
#     text(nodes[i].x, nodes[i].y, string(nodes[i].number))
# end

S, sc, tc = compliance(nodes, elements)
K = inv(S)
K2 = reorder(K)

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
