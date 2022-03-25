using LinearAlgebra: I, Symmetric, det

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

struct Node{TF, TI}
    x::TF
    y::TF
    number::TI
end

struct Element{VI, TF}
    nodenum::VI
    material::Material
    theta::TF
    beta::TF
end

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

function rotate_to_element(K, beta)
    c = cos(beta)
    s = sin(beta)
    T = [c^2 s^2 -2*s*c 0 0 0;
          s^2 c^2 2*s*c 0 0 0;
          s*c -s*c c^2-s^2 0 0 0;
          0 0 0 c -s 0;
          0 0 0 s c 0;
          0.0 0 0 0 0 1]

    return T*K*T'
end

# function rotate_to_element(K, x, y)
#     dx = x[2] - x[1]
#     dy = y[2] - y[1]
#     ds = sqrt(dx^2 + dy^2)
#     c = dx/ds
#     s = dy/ds
#     T = [c^2 s^2 0 -2*s*c 0 0;
#           s^2 c^2 0 2*s*c 0 0;
#           0.0 0 1 0 0 0;
#           s*c -s*c 0 c^2-s^2 0 0;
#           0 0 0 0 c -s;
#           0 0 0 0 s c]

#     return T*K*T'
# end


# function transformation(Q, R)
#     Qnew = zeros(6, 6)
#     for i = 1:3
#         for j = i:3
#             for k = 1:3
#                 for l = 1:3
#                     Qnew[i, j] = R[i, k]*R[j, l]*Q[k, l]
#                 end
#             end
#         end
#     end
#     return Symmetric(Qnew)
# end


# Rmaterial(theta) = [cos(theta) 0.0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)]
# Rgeometry(beta) = [cos(beta) -sin(beta) 0.0; sin(beta) cos(beta) 0; 0 0 1]
# function Rgeometry(y1, y2, z2, z2)
#     dy = y2 - y1
#     dz = z2 - z1
#     ds = sqrt(dy^2 + dz^2)
#     cb = dx/ds
#     sb = dy/ds
#     return [cb -sb 0.0; sb cb 0; 0 0 1]
# end

# function reorder(Q)
#     idx = [2, 3, 6, 4, 5, 1]
#     return Q[idx, idx]
# end

function elementQ(material, theta, beta)
    Q1 = stiffness(material)  # material
    Q2 = rotate_to_ply(Q1, theta)
    Q3 = rotate_to_element(Q2, beta)

    return Q3
end


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
    
    # basic matrices
    Z = [I [0.0 0 -y; 0 0 x; y -x 0]]  # translation and cross product
    S = [zeros(3, 3); I]
    Q = elementQ(element.material, element.theta, element.beta)
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


function node2idx(nodes)
    nn = 4
    idx = Vector{Int64}(undef, nn*3)
    for i = 1:nn
        idx[((i-1)*3+1):i*3] = ((nodes[i]-1)*3+1):nodes[i]*3
    end
    return idx
end

function sectionprops(nodes, elements)

    # initialize
    ne = length(elements) # number of elements
    nn = length(nodes)  # number of nodes
    ndof = 3 * nn  # 3 displacement dof per node
    
    # initialize
    etype = eltype(elements[1].theta)  # TODO add others
    A = zeros(etype, 6, 6)  # 6 x 6
    R = zeros(etype, ndof, 6)  #  nn*3 x 6
    E = zeros(etype, ndof, ndof)  # nn*3 x nn*3
    C = zeros(etype, ndof, ndof)  # nn*3 x nn*3
    L = zeros(etype, ndof, 6)  # nn*3 x 6
    M = zeros(etype, ndof, ndof)  # nn*3 x nn*3

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
    DT = zeros(6, ndof)
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
    D = DT'

    # Tr matrix
    Tr = zeros(6, 6)
    Tr[1, 5] = -1.0
    Tr[2, 4] = 1.0

    # solve first linear system
    AM = [E R D;
        R' A zeros(6, 6);
        DT zeros(6, 12)]
    
    B2 = [zeros(ndof, 6); Tr'; zeros(6, 6)]
    X2 = zeros(ndof+12, 6)
    for i = 1:6
        X2[:, i] = AM\B2[:, i]  # TODO: save factorization
    end
    dX = X2[1:ndof, :]
    dY = X2[ndof+1:ndof+6, :]

    # solve second linear system
    Bsub1 = [C'-C  L;  # NOTE: error in manual should be C' - C
            -L' zeros(6, 6);  # NOTE: sign fix, error in BECAS documentation
            zeros(6, ndof+6)]
    Bsub2 = [zeros(ndof, 6); I; zeros(6, 6)]
    B1 = Bsub1*X2[1:end-6, :] + Bsub2
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

    # stiffness matrix
    K = inv(S)

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

    return S, K
end

function centers(S, z, L)
    xs = -(S[6, 2] + S[6, 4]*(L - z))/S[6, 6]
    ys = (S[6, 1] + S[6, 5]*(L - z))/S[6, 6]
    xt = (-S[4, 4]*S[5, 3] + S[4, 5]*S[4, 3])/(S[4, 4]*S[5, 5] - S[4, 5]^2)
    yt = (-S[4, 3]*S[5, 5] + S[4, 5]*S[5, 3])/(S[4, 4]*S[5, 5] - S[4, 5]^2) 

    return xs, ys, xt, yt
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
        elements[m] = Element([11*(i-1)+j, 11*(i)+j, 11*(i)+j+1, 11*(i-1)+j+1], iso1, 0.0, 0.0)
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


S, K = sectionprops(nodes, elements)


@test isapprox(K[1, 1], 3.4899e-1, atol=0.0001e-1)
@test isapprox(K[2, 2], 3.4899e-1, atol=0.0001e-1)
@test isapprox(K[3, 3], 1.0, atol=1e-4)
@test isapprox(K[4, 4], 8.3384e-4, atol=0.0001e-4)
@test isapprox(K[5, 5], 8.3384e-4, atol=0.0001e-4)
@test isapprox(K[6, 6], 5.9084e-4, atol=0.0001e-4)

xs, ys, xt, yt = centers(S, 0.0, 1.0)
@test isapprox(xs, 0.0, atol=1e-8)
@test isapprox(ys, 0.0, atol=1e-8)
@test isapprox(xt, 0.0, atol=1e-8)
@test isapprox(yt, 0.0, atol=1e-8)

# --------- Case 2 --------
alpha = 1e1
iso2 = Material(100.0/alpha, 100.0/alpha, 0.2, 41.667/alpha, 1.0)  # note error in user guide for nu

let
m = 1
for i = 1:10
    for j = 1:10
        if i <= 5
            elements[m] = Element([11*(i-1)+j, 11*(i)+j, 11*(i)+j+1, 11*(i-1)+j+1], iso1, 0.0, 0.0)
        else
            elements[m] = Element([11*(i-1)+j, 11*(i)+j, 11*(i)+j+1, 11*(i-1)+j+1], iso2, 0.0, 0.0)
        end
        m += 1
    end
end
end


S, K = sectionprops(nodes, elements)

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
            elements[m] = Element([11*(i-1)+j, 11*(i)+j, 11*(i)+j+1, 11*(i-1)+j+1], ortho, theta, 0.0)
            m += 1
        end
    end
end

S, K = sectionprops(nodes, elements)
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
            elements[m] = Element([11*(i-1)+j, 11*(i)+j, 11*(i)+j+1, 11*(i-1)+j+1], ortho, theta, 0.0)
            m += 1
        end
    end
end

S, K = sectionprops(nodes, elements)
@test isapprox(K[1, 1], 7.598E-01, atol=0.001e-1)
@test isapprox(K[2, 2], 4.129E-01, atol=0.001e-1)
@test isapprox(K[3, 3], 3.435E+00, atol=0.001e0)
@test isapprox(K[4, 4], 2.489E-03, atol=0.001e-3)
@test isapprox(K[5, 5], 2.274E-03, atol=0.001e-3)
@test isapprox(K[6, 6], 9.499E-04, atol=0.001e-4)
@test isapprox(K[1, 3], 7.387E-01, atol=0.001e-1)
@test isapprox(K[4, 6], -4.613E-04, atol=0.001e-4)