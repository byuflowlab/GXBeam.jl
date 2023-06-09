# Unit Tests for section.jl and afmesh.jl

using GXBeam, LinearAlgebra, Random, Test
using ForwardDiff, FiniteDiff

@testset "section properties: material stiffness matrix" begin

    RNG = MersenneTwister(1234)

    E1 = rand(RNG)
    E2 = rand(RNG)
    E3 = rand(RNG)
    G12 = rand(RNG)
    G13 = rand(RNG)
    G23 = rand(RNG)
    nu12 = rand(RNG)
    nu23 = rand(RNG)
    nu13 = rand(RNG)
    rho = 1.0
    mat = Material(E1, E2, E3, G12, G13, G23, nu12, nu13, nu23, rho)
    Q1 = GXBeam.stiffness(mat)

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

end


@testset "section properties: square cross sections" begin

    # BECAS User Guide
    #  Case 1: Square cross section of isotropic material - S1

    iso1 = Material(100.0, 100.0, 100.0, 41.667, 41.667, 41.667, 0.2, 0.2, 0.2, 1.0)
    x = range(-0.05, 0.05, length=11)
    y = range(-0.05, 0.05, length=11)

    nodes = Vector{Node{Float64}}(undef, 11*11)
    elements = Vector{MeshElement{Float64}}(undef, 10*10)

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
            elements[m] = MeshElement([11*(i-1)+j, 11*(i)+j, 11*(i)+j+1, 11*(i-1)+j+1], iso1, 0.0)
            m += 1
        end
    end

    end

    S, sc, tc = compliance_matrix(nodes, elements, gxbeam_order=false, shear_center=false)
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

    # --------- Case 2 (S2, two isotropic materials) --------
    alpha = 1e1
    iso2 = Material(100.0/alpha, 100.0/alpha, 100.0/alpha, 41.667/alpha, 41.667/alpha, 41.667/alpha, 0.2, 0.2, 0.2, 1.0)  # note error in user guide for nu

    let
    m = 1
    for i = 1:10
        for j = 1:10
            if i <= 5
                elements[m] = MeshElement([11*(i-1)+j, 11*(i)+j, 11*(i)+j+1, 11*(i-1)+j+1], iso1, 0.0)
            else
                elements[m] = MeshElement([11*(i-1)+j, 11*(i)+j, 11*(i)+j+1, 11*(i-1)+j+1], iso2, 0.0)
            end
            m += 1
        end
    end
    end


    S, sc, tc = compliance_matrix(nodes, elements, gxbeam_order=false, shear_center=false)
    K = inv(S)

    @test isapprox(K[1, 1], 1.28e-1, atol=0.01e-1)
    @test isapprox(K[2, 2], 1.92e-1, atol=0.01e-1)
    @test isapprox(K[3, 3], 5.50e-1, atol=0.01e-1)
    @test isapprox(K[4, 4], 4.59e-4, atol=0.01e-4)
    @test isapprox(K[5, 5], 4.59e-4, atol=0.01e-4)
    @test isapprox(K[6, 6], 2.77e-4, atol=0.01e-4)
    @test isapprox(K[2, 6], -3.93e-3, atol=0.01e-3)
    @test isapprox(K[3, 5], 1.13e-2, atol=0.01e-2)

    # ------ case 3: S3, orthotropic material-------
    ortho = Material(480.0, 120.0, 120.0, 60.0, 50.0, 60.0, 0.19, 0.26, 0.19, 1.0)

    theta = 0.0
    let
        m = 1
        for i = 1:10
            for j = 1:10
                elements[m] = MeshElement([11*(i-1)+j, 11*(i)+j, 11*(i)+j+1, 11*(i-1)+j+1], ortho, theta)
                m += 1
            end
        end
    end

    S, sc, tc = compliance_matrix(nodes, elements, gxbeam_order=false, shear_center=false)
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
                elements[m] = MeshElement([11*(i-1)+j, 11*(i)+j, 11*(i)+j+1, 11*(i-1)+j+1], ortho, theta)
                m += 1
            end
        end
    end

    S, sc, tc = compliance_matrix(nodes, elements, gxbeam_order=false, shear_center=false)
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
                elements[m] = MeshElement([11*(i-1)+j, 11*(i)+j, 11*(i)+j+1, 11*(i-1)+j+1], ortho, theta)
                m += 1
            end
        end
    end

    S, sc, tc = compliance_matrix(nodes, elements, gxbeam_order=false, shear_center=false)
    K = inv(S)
    @test isapprox(K[1, 1], 5.0202E-01, atol=0.0001e-1)
    @test isapprox(K[2, 2], 5.0406E-01, atol=0.0001e-1)
    @test isapprox(K[3, 3], 1.2, atol=0.001e0)
    @test isapprox(K[4, 4], 1.0004E-03, atol=0.0001e-3)
    @test isapprox(K[5, 5], 1.0002E-03, atol=0.0001e-3)
    @test isapprox(K[6, 6], 8.5081E-04, atol=0.0001e-4)

end

@testset "section properties: cylinder" begin
    # BECAS User Guide
    # --------- cylinder C1 isotropic -------
    R = 0.1
    t = 0.01
    nr = 4
    nt = 120
    r = range(R - t, R, length=nr)
    theta = range(0, 2*pi, length=nt)

    iso1 = Material(100.0, 100.0, 100.0, 41.667, 41.667, 41.667, 0.2, 0.2, 0.2, 1.0)

    nodes = Vector{Node{Float64}}(undef, nr*(nt-1))
    elements = Vector{MeshElement{Float64}}(undef, (nr-1)*(nt-1))
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

            elements[n] = MeshElement([nr*ip+j, nr*(i-1)+j, nr*(i-1)+j+1, nr*ip+j+1], iso1, 0.0)
            n += 1
        end
    end
    end

    # plotmesh(nodes, elements)

    S, sc, tc = compliance_matrix(nodes, elements, gxbeam_order=false, shear_center=false)
    K = inv(S)

    @test isapprox(K[1, 1], 1.249E-01, atol=0.001e-1/2)
    @test isapprox(K[2, 2], 1.249E-01, atol=0.001e-1/2)
    @test isapprox(K[3, 3], 5.965E-01, atol=0.0015e-1)
    @test isapprox(K[4, 4], 2.697E-03, atol=0.002e-3)  # mesh discretization is not the same
    @test isapprox(K[5, 5], 2.697E-03, atol=0.002e-3)
    @test isapprox(K[6, 6], 2.248E-03, atol=0.001e-3)


    # ---- C2, half cylinder isotropic ----

    R = 0.1
    t = 0.01
    nr = 4
    nt = 60
    r = range(R - t, R, length=nr)
    theta = range(pi/2, 3*pi/2, length=nt)

    nodes = Vector{Node{Float64}}(undef, nr*nt)
    elements = Vector{MeshElement{Float64}}(undef, (nr-1)*(nt-1))
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
            elements[n] = MeshElement([nr*i+j, nr*(i-1)+j, nr*(i-1)+j+1, nr*i+j+1], iso1, 0.0)
            n += 1
        end
    end
    end

    # plotmesh(nodes, elements)

    S, sc, tc = compliance_matrix(nodes, elements, gxbeam_order=false, shear_center=false)
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
end

@testset "section properties: circular tube" begin
    # A critical assessment of computer tools for calculating composite wind turbine blade properties, Chen, Yu, Capellaro
    # ------ circular tube -------
    E = 73e9
    nu = 0.33
    G = E/(2*(1 + nu))
    tratio = 1.0/3
    R = 0.3
    t = tratio*2*R
    rho = 2800.0
    circmat = Material(E, E, E, G, G, G, nu, nu, nu, rho)

    nr = 20
    nt = 100
    r = range(R - t, R, length=nr)
    theta = range(0.0, 2*pi, length=nt)

    nodes = Vector{Node{Float64}}(undef, nr*(nt-1))
    elements = Vector{MeshElement{Float64}}(undef, (nr-1)*(nt-1))
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

            elements[n] = MeshElement([nr*ip+j, nr*(i-1)+j, nr*(i-1)+j+1, nr*ip+j+1], circmat, 0.0)
            n += 1
        end
    end
    end

    # plotmesh(nodes, elements)

    S, sc, tc = compliance_matrix(nodes, elements)
    K = inv(S)

    @test isapprox(K[1, 1], 1.835e10, rtol=0.001)
    @test isapprox(K[2, 2], 4.682e9, rtol=0.001)
    @test isapprox(K[3, 3], 4.682e9, rtol=0.001)
    @test isapprox(K[4, 4], 3.519e8, rtol=0.03)
    @test isapprox(K[5, 5], 4.587e8, rtol=0.002)
    @test isapprox(K[6, 6], 4.587e8, rtol=0.002)

    M, mc = mass_matrix(nodes, elements)

    @test isapprox(M[1, 1], 7.037e2, rtol=0.001)
    @test isapprox(M[5, 5], 1.759e1, rtol=0.003)
    @test isapprox(M[6, 6], 1.759e1, rtol=0.003)
end


function composite_pipe()
      # E1 = 141.963e9
    # E2 = 9.79056e9
    # nu12 = 0.42
    # G12 = 59.9844e9
    E1 = 20.59e6
    E2 = 1.42e6
    nu12 = 0.42
    G12 = 0.87e6
    rho = 1.0
    pipemat = Material(E1, E2, E2, G12, G12, G12, nu12, nu12, nu12, rho)

    nx = 50
    nt = 24
    nr = 21
    nodes = Vector{Node{Float64}}(undef, (2*(nx-1) + 2*(nt-1))*nr)
    elements = Vector{MeshElement{Float64}}(undef, (2*(nx-1) + 2*(nt-1))*(nr-1))


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
            if j <= ((nr-1) ÷ 2)
                theta = 90*pi/180
            else
                theta = 0.0
            end
            elements[n] = MeshElement([nr*(i-1)+j, nr*(i)+j, nr*(i)+j+1, nr*(i-1)+j+1], pipemat, theta)
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
            if j <= ((nr-1) ÷ 2)
                theta = 45*pi/180
            else
                theta = -45*pi/180
            end
            elements[n] = MeshElement([nr*(I+i-1)+j, nr*(I+i)+j, nr*(I+i)+j+1, nr*(I+i-1)+j+1], pipemat, theta)
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
            if j <= ((nr-1) ÷ 2)
                theta = 90*pi/180
            else
                theta = 0.0
            end
            elements[n] = MeshElement([nr*(I+i-1)+j, nr*(I+i)+j, nr*(I+i)+j+1, nr*(I+i-1)+j+1], pipemat, theta)
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
            if j <= ((nr-1) ÷ 2)
                theta = 45*pi/180
            else
                theta = -45*pi/180
            end
            elements[n] = MeshElement([nr*(I+i-1)+j, nr*(Ip+i)+j, nr*(Ip+i)+j+1, nr*(I+i-1)+j+1], pipemat, theta)
            n += 1
        end
    end

    return nodes, elements
end

@testset "section properties: multi-layer composite pipe" begin

    #  --- Generalized Timoshenko Theory of the Variational Asymptotic Beam Sectional Analysis ----
    # multi-layer composite pipe
    # note that the previous paper has a composite pipe also, but the numbers are inconsistent.
    # see also preVABS documentation examples

    nodes, elements = composite_pipe()
    # plotmesh(nodes, elements)

    S, sc, tc = compliance_matrix(nodes, elements)
    K = inv(S)

    @test isapprox(K[1, 1], 1.03892e7, rtol=0.005)
    @test isapprox(K[2, 2], 7.85310e5, rtol=0.005)
    @test isapprox(K[3, 3], 3.29279e5, rtol=0.005)
    @test isapprox(K[1, 4], 9.84575e4, rtol=0.01)
    @test isapprox(K[2, 5], -8.21805e3, rtol=0.011)
    @test isapprox(K[3, 6], -5.20981e4, rtol=0.01)
    @test isapprox(K[4, 4], 6.87275e5, rtol=0.01)
    @test isapprox(K[5, 5], 1.88238e6, rtol=0.005)
    @test isapprox(K[6, 6], 5.38987e6, rtol=0.005)

    # println("K11 = ", round((K[1, 1]/1.03892e7 - 1)*100, digits=2), "%")
    # println("K22 = ", round((K[2, 2]/7.85310e5 - 1)*100, digits=2), "%")
    # println("K33 = ", round((K[3, 3]/3.29279e5 - 1)*100, digits=2), "%")
    # println("K14 = ", round((K[1, 4]/9.84575e4 - 1)*100, digits=2), "%")
    # println("K25 = ", round((K[2, 5]/-8.21805e3 - 1)*100, digits=2), "%")
    # println("K36 = ", round((K[3, 6]/-5.20981e4 - 1)*100, digits=2), "%")
    # println("K44 = ", round((K[4, 4]/6.87275e5 - 1)*100, digits=2), "%")
    # println("K55 = ", round((K[5, 5]/1.88238e6 - 1)*100, digits=2), "%")
    # println("K66 = ", round((K[6, 6]/5.38987e6 - 1)*100, digits=2), "%")
end

# borrowing from FLOWMath (just to avoid another dependency for this one off)
function findindex(xvec, x)

    n = length(xvec)
    i = searchsortedlast(real(xvec), real(x))

    # this version allows extrapolation
    if i == 0
        i = 1
    elseif i == n
        i = n - 1
    end

    return i
end

function linearinterp(xdata, ydata, x::Number)

    i = findindex(xdata, x)

    eta = (x - xdata[i]) / (xdata[i+1] - xdata[i])
    y = ydata[i] + eta*(ydata[i+1] - ydata[i])

    return y
end

linearinterp(xdata, ydata, x::AbstractVector) = linearinterp.(Ref(xdata), Ref(ydata), x)

@testset "strain recovery: multi-layer composite pipe" begin
    # Loss of Accuracy Using Smeared Properties in Composite Beam Modeling
    # Ning Liu, Purdue University
    # Uses same multi-layer composite pipe example from above (also with FEA comparisons)
    # also shows smeared properties, which are much less accurate

    nodes, elements = composite_pipe()
    cache = initialize_cache(nodes, elements)
    S, sc, tc = compliance_matrix(nodes, elements; cache)
    K = inv(S)

    F = [0.0; 0; 0]
    M = [-1000.0; 0; 0]
    epsilon_b, sigma_b, epsilon_p, sigma_p = strain_recovery(F, M, nodes, elements, cache; gxbeam_order=false)

    # using PyPlot
    # pygui(true); close("all")
    # figure()
    # plotsoln(nodes, elements, sigma_b[1, :], PyPlot)
    # colorbar()

    # figure()
    # plotsoln(nodes, elements, sigma_b[2, :], PyPlot)
    # colorbar()

    # figure()
    # plotsoln(nodes, elements, sigma_b[3, :], PyPlot)
    # colorbar()

    # figure()
    # ne = length(elements)
    # for i = 1:ne
    #     _, xc, yc = GXBeam.area_and_centroid_of_element(nodes[elements[i].nodenum])
    #     if abs(xc - 0.0) < 0.01 && yc > 0.0
    #         nn = nodes[elements[i].nodenum]
    #         for j = 1:4
    #             plot(nn[j].x, nn[j].y, "x")
    #             text(xc, yc, string(i))
    #         end
    #     end
    # end


    # grab elements at x = 0 from y = 0.3 -> 0.5
    idx = 481:500  # elements at x = 0 from y = 0.3 -> 0.5
    n = length(idx)
    yvec = zeros(n)
    s11 = zeros(n)
    s22 = zeros(n)
    for i = 1:n
        _, _, yvec[i] = GXBeam.area_and_centroid_of_element(nodes[elements[idx[i]].nodenum])
        s11[i] = sigma_b[3, idx[i]]
        s22[i] = sigma_b[1, idx[i]]
    end

    # data from paper
    data1 = [
    -1.0408340855860843e-17  -0.2233363719234286
    0.02504493708807669  -0.24156791248860665
    0.049970041941282226  -0.25979945305378427
    0.07501497902935894  -0.2734731084776677
    0.10005991611743567  -0.29170464904284577
    0.10005991611743564  -4.416590701914313
    0.12510485320551235  -4.690063810391979
    0.15002995805871785  -4.96809480401094
    0.17507489514679453  -5.241567912488606
    0.2001198322348712  -5.519598906107569
    ]

    # interpolate data onto my pts
    ydata = yvec .- 0.3
    s11interp = linearinterp(data1[:, 1], data1[:, 2], ydata)

    # figure()
    # plot(ydata, s11/1e3, ".")
    # plot(ydata, s11interp, "kx")
    # plot(data1[:, 1], data1[:, 2], "k--")
    # xlabel("x3")
    # ylabel("s11 (ksi)")


    data2 = [0.00011487650775416497  0.004832104832104944
    0.02504307869040781  0.06871416871416883
    0.04997128087306147  0.13226863226863234
    0.07501435956346927  0.19582309582309587
    0.09994256174612295  0.2597051597051597
    0.10005743825387711  -0.10196560196560189
    0.12498564043653071  -0.10491400491400493
    0.1500287191269386  -0.10786240786240764
    0.17484204480183801  -0.11146601146601148
    0.19999999999999993  -0.11408681408681415
    ]

    s22interp = linearinterp(data2[:, 1], data2[:, 2], ydata)

    # figure()
    # plot(ydata, s22/1e3, ".")
    # plot(ydata, s22interp, "kx")
    # plot(data2[:, 1], data2[:, 2], "k--")
    # xlabel("x3")
    # ylabel("s22 (ksi)")

    s11mine = s11/1e3
    s11fea = s11interp

    n = length(s11mine)
    for i = 1:n
        @test isapprox(s11mine[i], s11fea[i], rtol=0.015)
    end

    s22mine = s22/1e3
    s22fea = s22interp

    @test isapprox(s22mine[1], s22fea[1], atol=0.003) # close to zero
    @test isapprox(s22mine[2], s22fea[2], atol=0.002) # close to zero
    @test isapprox(s22mine[3], s22fea[3], atol=0.002) # close to zero
    @test isapprox(s22mine[4], s22fea[4], atol=0.002) # close to zero

    n = length(s22mine)
    for i = 5:n
        @test isapprox(s22mine[i], s22fea[i], rtol=0.01)
    end
end


@testset "section properties: composite wind turbine airfoil" begin
    # A critical assessment of computer tools for calculating composite wind turbine blade properties, Chen, Yu, Capellaro
    # ST1
    # See also: https://wenbinyugroup.github.io/ivabs/prevabs/contents/examples/example_airfoil.html

    xaf = [1.00000000, 0.99619582, 0.98515158, 0.96764209, 0.94421447, 0.91510964, 0.88074158, 0.84177999, 0.79894110, 0.75297076, 0.70461763, 0.65461515, 0.60366461, 0.55242353, 0.50149950, 0.45144530, 0.40276150, 0.35589801, 0.31131449, 0.26917194, 0.22927064, 0.19167283, 0.15672257, 0.12469599, 0.09585870, 0.07046974, 0.04874337, 0.03081405, 0.01681379, 0.00687971, 0.00143518, 0.00053606, 0.00006572, 0.00001249, 0.00023032, 0.00079945, 0.00170287, 0.00354717, 0.00592084, 0.01810144, 0.03471169, 0.05589286, 0.08132751, 0.11073805, 0.14391397, 0.18067874, 0.22089879, 0.26433734, 0.31062190, 0.35933893, 0.40999990, 0.46204424, 0.51483073, 0.56767889, 0.61998250, 0.67114514, 0.72054815, 0.76758733, 0.81168064, 0.85227225, 0.88883823, 0.92088961, 0.94797259, 0.96977487, 0.98607009, 0.99640466, 1.00000000]
    yaf = [0.00000000, 0.00017047, 0.00100213, 0.00285474, 0.00556001, 0.00906779, 0.01357364, 0.01916802, 0.02580144, 0.03334313, 0.04158593, 0.05026338, 0.05906756, 0.06766426, 0.07571157, 0.08287416, 0.08882939, 0.09329359, 0.09592864, 0.09626763, 0.09424396, 0.09023579, 0.08451656, 0.07727756, 0.06875796, 0.05918984, 0.04880096, 0.03786904, 0.02676332, 0.01592385, 0.00647946, 0.00370956, 0.00112514, -0.00046881, -0.00191488, -0.00329201, -0.00470585, -0.00688469, -0.00912202, -0.01720842, -0.02488211, -0.03226730, -0.03908459, -0.04503763, -0.04986836, -0.05338180, -0.05551392, -0.05636585, -0.05605816, -0.05472399, -0.05254383, -0.04969990, -0.04637175, -0.04264894, -0.03859653, -0.03433153, -0.02996944, -0.02560890, -0.02134397, -0.01726049, -0.01343567, -0.00993849, -0.00679919, -0.00402321, -0.00180118, -0.00044469, 0.00000000]

    uni = Material(37.00e9, 9.00e9, 9.00e9, 4.00e9, 4.00e9, 4.00e9, 0.28, 0.28, 0.28, 1.86e3)
    double = Material(10.30e9, 10.30e9, 10.30e9, 8.00e9, 8.00e9, 8.00e9, 0.30, 0.30, 0.30, 1.83e3)
    gelcoat = Material(1e1, 1e1, 1e1, 1.0, 1.0, 1.0, 0.30, 0.30, 0.30, 1.83e3)
    nexus = Material(10.30e9, 10.30e9, 10.30e9, 8.00e9, 8.00e9, 8.00e9, 0.30, 0.30, 0.30, 1.664e3)
    balsa = Material(0.01e9, 0.01e9, 0.01e9, 2e5, 2e5, 2e5, 0.30, 0.30, 0.30, 0.128e3)
    mat = [uni, double, gelcoat, nexus, balsa]

    chord = 1.9
    twist = 0.0*pi/180
    paxis = 0.4750 / chord
    xbreak = [0.0, 0.0041, 0.1147, 0.5366, 1.0]
    webloc = [0.15 0.5]

    idx = [3, 4, 2]
    t = [0.000381, 0.00051, 18*0.00053]
    theta = [0, 0, 20]*pi/180
    layup1 = Layer.(mat[idx], t, theta)
    idx = [3, 4, 2]
    t = [0.000381, 0.00051, 33*0.00053]
    theta = [0, 0, 20]*pi/180
    layup2 = Layer.(mat[idx], t, theta)
    idx = [3, 4, 2, 1, 5, 1, 2]
    t = [0.000381, 0.00051, 17*0.00053, 38*0.00053, 1*0.003125, 37*0.00053, 16*0.00053]
    theta = [0, 0, 20, 30, 0, 30, 20]*pi/180
    layup3 = Layer.(mat[idx], t, theta)
    idx = [3, 4, 2, 5, 2]
    t = [0.000381, 0.00051, 17*0.00053, 0.003125, 16*0.00053]
    theta = [0, 0, 20, 0, 0]*pi/180
    layup4 = Layer.(mat[idx], t, theta)

    idx = [1, 5, 1]
    t = [38*0.00053, 0.003125, 38*0.00053]
    theta = [0, 0, 0]*pi/180
    web = Layer.(mat[idx], t, theta)

    segments = [layup1, layup2, layup3, layup4]
    webs = [web, web]

    nodes, elements = afmesh(xaf, yaf, chord, twist, paxis, xbreak, webloc, segments, webs, ds=0.005, dt=0.01, wns=20)
    # nodes, elements = afmesh(xaf, yaf, chord, twist, paxis, xbreak, webloc, segments, webs, nt=[[1, 1, 7], [1, 1, 7], [1, 1, 1, 2, 1, 2, 1], [1, 1, 1, 1, 5]], wnt=[[1, 1, 1], [1, 1, 1]])

    S, sc, tc = compliance_matrix(nodes, elements)
    K = inv(S)

    @test isapprox(log10(K[1, 1]), log10(abs(2.389e9)), rtol=0.006)
    @test isapprox(log10(K[1, 2]), log10(abs(1.524e6)), rtol=0.05)
    @test isapprox(log10(K[2, 2]), log10(abs(4.334e8)), rtol=0.03)
    @test isapprox(log10(K[1, 3]), log10(abs(6.734e6)), rtol=0.09)
    @test isapprox(log10(-K[2, 3]), log10(abs(-3.741e6)), rtol=0.07)
    @test isapprox(log10(K[3, 3]), log10(abs(2.743e7)), rtol=0.03)
    @test isapprox(log10(-K[1, 4]), log10(abs(-3.382e7)), rtol=0.01)
    @test isapprox(log10(-K[2, 4]), log10(abs(-2.935e5)), rtol=0.03)
    @test isapprox(log10(-K[3, 4]), log10(abs(-4.592e4)), rtol=0.18)
    @test isapprox(log10(K[4, 4]), log10(abs(2.167e7)), rtol=0.01)
    @test isapprox(log10(-K[1, 5]), log10(abs(-2.627e7)), rtol=0.02)
    @test isapprox(log10(K[2, 5]), log10(abs(1.527e7)), rtol=0.03)
    # @test isapprox(log10(-K[3, 5]), log10(abs(-6.869e2)), rtol=0.5)
    @test isapprox(log10(-K[4, 5]), log10(abs(-6.279e4)), rtol=0.14)
    @test isapprox(log10(K[5, 5]), log10(abs(1.970e7)), rtol=0.005)
    @test isapprox(log10(-K[1, 6]), log10(abs(-4.736e8)), rtol=0.03)
    # @test isapprox(log10(K[2, 6]), log10(abs(3.835e5)), rtol=0.08)
    @test isapprox(log10(-K[3, 6]), log10(abs(-4.742e6)), rtol=0.02)
    @test isapprox(log10(K[4, 6]), log10(abs(1.430e6)), rtol=0.06)
    @test isapprox(log10(K[5, 6]), log10(abs(1.209e7)), rtol=0.02)
    @test isapprox(log10(K[6, 6]), log10(abs(4.406e8)), rtol=0.02)

    # println("K11 = ", round((K[1, 1]/2.389e9 - 1)*100, digits=2), "%")
    # println("K12 = ", round((K[1, 2]/1.524e6 - 1)*100, digits=2), "%")
    # println("K22 = ", round((K[2, 2]/4.334e8 - 1)*100, digits=2), "%")
    # println("K13 = ", round((K[1, 3]/6.734e6 - 1)*100, digits=2), "%")
    # println("K23 = ", round((K[2, 3]/-3.741e6 - 1)*100, digits=2), "%")
    # println("K33 = ", round((K[3, 3]/2.743e7 - 1)*100, digits=2), "%")
    # println("K14 = ", round((K[1, 4]/-3.382e7 - 1)*100, digits=2), "%")
    # println("K24 = ", round((K[2, 4]/-2.935e5 - 1)*100, digits=2), "%")
    # println("K34 = ", round((K[3, 4]/-4.592e4 - 1)*100, digits=2), "%")
    # println("K44 = ", round((K[4, 4]/2.167e7 - 1)*100, digits=2), "%")
    # println("K15 = ", round((K[1, 5]/-2.627e7 - 1)*100, digits=2), "%")
    # println("K25 = ", round((K[2, 5]/1.527e7 - 1)*100, digits=2), "%")
    # println("K35 = ", round((K[3, 5]/-6.869e2 - 1)*100, digits=2), "%")
    # println("K45 = ", round((K[4, 5]/-6.279e4 - 1)*100, digits=2), "%")
    # println("K55 = ", round((K[5, 5]/1.970e7 - 1)*100, digits=2), "%")
    # println("K16 = ", round((K[1, 6]/-4.736e8 - 1)*100, digits=2), "%")
    # println("K26 = ", round((K[2, 6]/3.835e5 - 1)*100, digits=2), "%")
    # println("K36 = ", round((K[3, 6]/-4.742e6 - 1)*100, digits=2), "%")
    # println("K46 = ", round((K[4, 6]/1.430e6 - 1)*100, digits=2), "%")
    # println("K56 = ", round((K[5, 6]/1.209e7 - 1)*100, digits=2), "%")
    # println("K66 = ", round((K[6, 6]/4.406e8 - 1)*100, digits=2), "%")

    # println("K11 = ", round((log10(K[1, 1])/log10(2.389e9) - 1)*100, digits=2), "%")
    # println("K12 = ", round((log10(K[1, 2])/log10(1.524e6) - 1)*100, digits=2), "%")
    # println("K22 = ", round((log10(K[2, 2])/log10(4.334e8) - 1)*100, digits=2), "%")
    # println("K13 = ", round((log10(K[1, 3])/log10(6.734e6) - 1)*100, digits=2), "%")
    # println("K23 = ", round((log10(-K[2, 3])/log10(3.741e6) - 1)*100, digits=2), "%")
    # println("K33 = ", round((log10(K[3, 3])/log10(2.743e7) - 1)*100, digits=2), "%")
    # println("K14 = ", round((log10(-K[1, 4])/log10(3.382e7) - 1)*100, digits=2), "%")
    # println("K24 = ", round((log10(-K[2, 4])/log10(2.935e5) - 1)*100, digits=2), "%")
    # println("K34 = ", round((log10(-K[3, 4])/log10(4.592e4) - 1)*100, digits=2), "%")
    # println("K44 = ", round((log10(K[4, 4])/log10(2.167e7) - 1)*100, digits=2), "%")
    # println("K15 = ", round((log10(-K[1, 5])/log10(2.627e7) - 1)*100, digits=2), "%")
    # println("K25 = ", round((log10(K[2, 5])/log10(1.527e7) - 1)*100, digits=2), "%")
    # println("K35 = ", round((log10(-K[3, 5])/log10(6.869e2) - 1)*100, digits=2), "%")
    # println("K45 = ", round((log10(-K[4, 5])/log10(6.279e4) - 1)*100, digits=2), "%")
    # println("K55 = ", round((log10(K[5, 5])/log10(1.970e7) - 1)*100, digits=2), "%")
    # println("K16 = ", round((log10(-K[1, 6])/log10(4.736e8) - 1)*100, digits=2), "%")
    # # println("K26 = ", round((log10(K[2, 6])/log10(3.835e5) - 1)*100, digits=2), "%")
    # println("K36 = ", round((log10(-K[3, 6])/log10(4.742e6) - 1)*100, digits=2), "%")
    # println("K46 = ", round((log10(K[4, 6])/log10(1.430e6) - 1)*100, digits=2), "%")
    # println("K56 = ", round((log10(K[5, 6])/log10(1.209e7) - 1)*100, digits=2), "%")
    # println("K66 = ", round((log10(K[6, 6])/log10(4.406e8) - 1)*100, digits=2), "%")


    M, mc = mass_matrix(nodes, elements)
    @test isapprox(M[1, 1], 258.053, rtol=0.01)
    @test isapprox(M[5, 5], 2.172, rtol=0.02)
    @test isapprox(M[6, 6], 46.418, rtol=0.03)
    # paper has it relative to pitch axis
    @test isapprox(mc[1] - paxis*chord, 0.2778, rtol=0.01)
    @test isapprox(mc[2], 0.02743, rtol=0.02)
    @test isapprox(tc[1] - paxis*chord, 0.233, rtol=0.02)
    @test isapprox(tc[2], 0.029, rtol=0.02)
    @test isapprox(sc[1], 0.031 + paxis*chord, rtol=0.22) # relative to leading edge since its close to pitch axis
    @test isapprox(sc[2], 0.040, rtol=0.05)

    Ixx = M[5, 5]
    Iyy = M[6, 6]
    Ixy = -M[5, 6]
    theta = 0.5 * atan(2*Ixy / (Iyy - Ixx))
    @test isapprox(theta*180/pi, -1.244, rtol=0.1)

end


@testset "strain recovery: composite wind turbine airfoil" begin
    # Loss of Accuracy Using Smeared Properties in Composite Beam Modeling
    # Ning Liu, Purdue University
    # Uses airfoil from above, but with different property values.

    xaf = [1.00000000, 0.99619582, 0.98515158, 0.96764209, 0.94421447, 0.91510964, 0.88074158, 0.84177999, 0.79894110, 0.75297076, 0.70461763, 0.65461515, 0.60366461, 0.55242353, 0.50149950, 0.45144530, 0.40276150, 0.35589801, 0.31131449, 0.26917194, 0.22927064, 0.19167283, 0.15672257, 0.12469599, 0.09585870, 0.07046974, 0.04874337, 0.03081405, 0.01681379, 0.00687971, 0.00143518, 0.00053606, 0.00006572, 0.00001249, 0.00023032, 0.00079945, 0.00170287, 0.00354717, 0.00592084, 0.01810144, 0.03471169, 0.05589286, 0.08132751, 0.11073805, 0.14391397, 0.18067874, 0.22089879, 0.26433734, 0.31062190, 0.35933893, 0.40999990, 0.46204424, 0.51483073, 0.56767889, 0.61998250, 0.67114514, 0.72054815, 0.76758733, 0.81168064, 0.85227225, 0.88883823, 0.92088961, 0.94797259, 0.96977487, 0.98607009, 0.99640466, 1.00000000]
    yaf = [0.00000000, 0.00017047, 0.00100213, 0.00285474, 0.00556001, 0.00906779, 0.01357364, 0.01916802, 0.02580144, 0.03334313, 0.04158593, 0.05026338, 0.05906756, 0.06766426, 0.07571157, 0.08287416, 0.08882939, 0.09329359, 0.09592864, 0.09626763, 0.09424396, 0.09023579, 0.08451656, 0.07727756, 0.06875796, 0.05918984, 0.04880096, 0.03786904, 0.02676332, 0.01592385, 0.00647946, 0.00370956, 0.00112514, -0.00046881, -0.00191488, -0.00329201, -0.00470585, -0.00688469, -0.00912202, -0.01720842, -0.02488211, -0.03226730, -0.03908459, -0.04503763, -0.04986836, -0.05338180, -0.05551392, -0.05636585, -0.05605816, -0.05472399, -0.05254383, -0.04969990, -0.04637175, -0.04264894, -0.03859653, -0.03433153, -0.02996944, -0.02560890, -0.02134397, -0.01726049, -0.01343567, -0.00993849, -0.00679919, -0.00402321, -0.00180118, -0.00044469, 0.00000000]

    uni = Material(5.3664e6, 1.3053e6, 1.3053e6, 5.8015e5, 5.8015e5, 5.8015e5, 0.28, 0.28, 0.28, 1.740449e-4)
    double = Material(1.4939e6, 1.4939e6, 1.4939e6, 1.1603e6, 1.1603e6, 1.1603e6, 0.30, 0.30, 0.30, 1.712378e-4)
    gelcoat = Material(1.4504e-3, 1.4504e-3, 1.4504e-3, 1.4504e-4, 1.4504e-4, 1.4504e-4, 0.30, 0.30, 0.30, 1.712378e-4)
    nexus = Material(1.4939e6, 1.4939e6, 1.4939e6, 1.1603e6, 1.1603e6, 1.1603e6, 0.30, 0.30, 0.30, 1.557047e-4)
    balsa = Material(1.4504e3, 1.4504e3, 1.4504e3, 2.9008e1, 2.9008e1, 2.9008e1, 0.30, 0.30, 0.30, 1.197729e-5)
    mat = [uni, double, gelcoat, nexus, balsa]

    chord = 74.8031
    twist = 0.0*pi/180
    paxis = 19.72747356 / chord
    xbreak = [0.0, 0.0041, 0.1147, 0.5366, 1.0]
    webloc = [0.15 0.5]

    idx = [3, 4, 2]
    t = [0.015, 0.02007874, 18*0.020866142]
    theta = [0, 0, 20]*pi/180
    layup1 = Layer.(mat[idx], t, theta)
    idx = [3, 4, 2]
    t = [0.015, 0.02007874, 33*0.020866142]
    theta = [0, 0, 20]*pi/180
    layup2 = Layer.(mat[idx], t, theta)
    idx = [3, 4, 2, 1, 5, 1, 2]
    t = [0.015, 0.02007874, 17*0.020866142, 38*0.020866142, 1*0.123031496, 37*0.020866142, 16*0.020866142]
    theta = [0, 0, 20, 30, 0, 30, 20]*pi/180
    layup3 = Layer.(mat[idx], t, theta)
    idx = [3, 4, 2, 5, 2]
    t = [0.015, 0.02007874, 17*0.020866142, 0.123031496, 16*0.020866142]
    theta = [0, 0, 20, 0, 0]*pi/180
    layup4 = Layer.(mat[idx], t, theta)

    idx = [1, 5, 1]
    t = [38*0.020866142, 2*0.061515748, 38*0.020866142]
    theta = [0, 0, 0]*pi/180
    web = Layer.(mat[idx], t, theta)

    segments = [layup1, layup2, layup3, layup4]
    webs = [web, web]

    nodes, elements = afmesh(xaf, yaf, chord, twist, paxis, xbreak, webloc, segments, webs, ds=0.01, dt=0.2, wns=20)#, ds=0.05, dt=0.01, wns=20)

    cache = initialize_cache(nodes, elements)
    S, sc, tc = compliance_matrix(nodes, elements; cache)
    K = inv(S)

    # figure()
    # ne = length(elements)
    # for i = 1:ne
    #     _, xc, yc = GXBeam.area_and_centroid_of_element(nodes[elements[i].nodenum])
    #     if abs(xc - (19.72747356+2.3725)) < 0.5 && yc > 0.0
    #         n = nodes[elements[i].nodenum]
    #         for j = 1:4
    #             plot(n[j].x, n[j].y, "x")
    #             text(xc, yc, string(i))
    #         end
    #     end
    # end

    F = [0.0; 0; 0]
    M = [0.0; 0; 1e6]
    epsilon_b, sigma_b, epsilon_p, sigma_p = strain_recovery(F, M, nodes, elements, cache; gxbeam_order=false)

    idx = 585:-1:571
    n = length(idx)
    x3vec = zeros(n)
    s11 = zeros(n)
    s22 = zeros(n)
    s12 = zeros(n)
    for i = 1:n
        _, _, x3vec[i] = GXBeam.area_and_centroid_of_element(nodes[elements[idx[i]].nodenum])
        s11[i] = sigma_b[3, idx[i]]
        s22[i] = sigma_b[1, idx[i]]
        s12[i] = sigma_b[5, idx[i]]
    end

    data11 = [0.0180411653073842  406.5894924309885
    0.12482948694456675  415.13802315227076
    0.22710651266109672  422.61798753339275
    0.33239062250525886  431.166518254675
    0.3319566379808254  -422.6179875333927
    0.3951603224558594  -443.98931433659857
    0.47943100878363254  -471.7720391807662
    0.583264485171124  -505.9661620658952
    0.6886021733516359  -540.1602849510242
    0.7924369891975358  -575.4229741763141
    0.8962718050434362  -610.6856634016029
    1.0001052814309277  -644.8797862867322
    1.1039387578184199  -679.0739091718616
    1.1045910740634781  0.5342831700800161
    1.2264322292981351  0.5342831700800161
    1.2266532399355787  -175.77916295636703
    1.312426798597963  -202.4933214603742
    1.3966988243841445  -231.3446126447019
    1.5200950908313546  -271.9501335707927
    1.6449982479884022  -314.6927871772042
    1.7699000656870412  -356.36687444345546
    1.89480322284409  -399.10952804986675
    2.0197050405427284  -440.78361531611824
    2.021628502817685  424.75512021371304
    2.123904189075806  433.30365093499523
    2.2487417127708254  442.9207479964378
    2.3735805759242536  451.46927871772004
    2.3753741107335635  220.6589492430985
    2.3964317363774406  221.72751558325876
    2.393700580682011  0.5342831700797888
    2.413255333991277  0.5342831700797888
    ]

    data22 = [0.001481042654028153  127.32714942430584
    0.12885071090047373  118.00254045514711
    0.2325236966824643  109.58167883523336
    0.33323459715639814  101.76177362196324
    0.3347156398104264  -1.2414072840613244
    0.41913507109004733  -11.762058266797624
    0.5242890995260663  -24.987902594537616
    0.6072274881516587  -35.80867597455273
    0.6916469194312798  -46.62962725758939
    0.7938388625592419  -59.85511577928645
    0.8975118483412323  -73.38108250430525
    1.0011848341232228  -86.30644862872356
    1.1063388625592419  -100.73349415766478
    1.1063388625592419  -0.13289355706416472
    1.2277843601895737  0.45311899577302484
    1.2292654028436019  40.99348163329208
    1.3758886255924172  22.056950315244137
    1.5017772511848344  6.125912642500339
    1.6039691943127963  -7.399876179497028
    1.7091232227488151  -20.625720507237133
    1.8142772511848342  -34.15186513527743
    1.916469194312796  -47.97795425757511
    2.021623222748815  -61.80439918591571
    2.0245853080568725  26.18323299602929
    2.0838270142180093  21.371312070364183
    2.166765402843602  14.754742894553374
    2.2719194312796214  6.634003671918379
    2.375592417061612  -1.1862573473947577
    2.375592417061612  230.34527418413685
    2.397808056872038  230.04230533851398
    2.397808056872038  -0.8886255924169859
    2.4170616113744083  -0.2903377310960309
    ]

    data12 = [
    0.017772511848341277  -577.8947368421051
    0.12440758293838872  -603.157894736842
    0.22956161137440767  -628.4210526315788
    0.33175355450236965  -652.6315789473684
    0.3347156398104265  -733.6842105263157
    0.3969194312796209  -754.7368421052631
    0.47985781990521326  -781.0526315789473
    0.5850118483412324  -816.8421052631579
    0.6886848341232227  -850.5263157894738
    0.7923578199052131  -884.2105263157894
    0.8960308056872036  -918.9473684210525
    1.0011848341232228  -952.6315789473683
    1.1033767772511847  -987.3684210526317
    1.1078199052132702  1.1368683772161603e-13
    1.2263033175355447  1.0526315789474552
    1.2277843601895735  -470.52631578947353
    1.294431279620853  -491.5789473684209
    1.3951421800947865  -525.2631578947368
    1.4988151658767772  -559.9999999999999
    1.6024881516587675  -595.7894736842104
    1.7076421800947865  -630.5263157894736
    1.8127962085308051  -665.2631578947367
    1.914988151658768  -698.9473684210525
    2.018661137440758  -733.6842105263155
    2.021623222748815  -656.8421052631578
    2.1238151658767768  -681.0526315789473
    2.249703791469194  -711.5789473684209
    2.375592417061611  -742.1052631578947
    2.374111374407583  -941.0526315789473
    2.396327014218009  -947.3684210526314
    2.396327014218009  3.157894736842138
    2.4140995260663507  2.1052631578947967
    ]

    # interpolate onto data points
    x3 = x3vec .- 4.77
    s11interp = linearinterp(data11[:, 1], data11[:, 2], x3)
    s22interp = linearinterp(data22[:, 1], data22[:, 2], x3)
    s12interp = linearinterp(data12[:, 1], data12[:, 2], x3)


    # figure()
    # plot(x3, s11, "x")
    # plot(data11[:, 1], data11[:, 2])
    # plot(x3, s11interp, "o")
    # xlabel("x3")
    # ylabel("s11")

    # figure()
    # plot(x3, s22, "x")
    # plot(data22[:, 1], data22[:, 2])
    # plot(x3, s22interp, "o")
    # xlabel("x3")
    # ylabel("s22")

    # figure()
    # plot(x3, s12, "x")
    # plot(data12[:, 1], data12[:, 2])
    # plot(x3, s12interp, "o")
    # xlabel("x3")
    # ylabel("s12")

    @test isapprox(s11[1], s11interp[1], rtol=0.06)
    @test isapprox(s11[2], s11interp[2], rtol=0.06)
    @test isapprox(s11[3], s11interp[3], atol=100.0) # these next 4 are offset
    @test isapprox(s11[4], s11interp[4], atol=110.0)
    @test isapprox(s11[5], s11interp[5], atol=110.0)
    @test isapprox(s11[6], s11interp[6], atol=120.0)
    @test isapprox(s11[7], s11interp[7], atol=1.0) # close to zero
    @test isapprox(s11[8], s11interp[8], atol=122.0)  #next 4 offset
    @test isapprox(s11[9], s11interp[9], atol=120.0)
    @test isapprox(s11[10], s11interp[10], atol=110.0)
    @test isapprox(s11[11], s11interp[11], atol=100.0)
    @test isapprox(s11[12], s11interp[12], rtol=0.08)
    @test isapprox(s11[13], s11interp[13], rtol=0.08)
    @test isapprox(s11[14], s11interp[14], rtol=0.06)
    @test isapprox(s11[15], s11interp[15], atol=1.0) # close to zero

    @test isapprox(s22[1], s22interp[1], rtol=0.12)
    @test isapprox(s22[2], s22interp[2], rtol=0.13)
    @test isapprox(s22[3], s22interp[3], atol=12.0)
    @test isapprox(s22[4], s22interp[4], atol=15.0)
    @test isapprox(s22[5], s22interp[5], atol=20.0)
    @test isapprox(s22[6], s22interp[6], atol=20.0)
    @test isapprox(s22[7], s22interp[7], atol=1.0) # close to zero
    @test isapprox(s22[8], s22interp[8], atol=20.0)
    @test isapprox(s22[9], s22interp[9], atol=20.0)
    @test isapprox(s22[10], s22interp[10], atol=15.0)
    @test isapprox(s22[11], s22interp[11], atol=15.0)
    @test isapprox(s22[12], s22interp[12], atol=15.0)
    @test isapprox(s22[13], s22interp[13], atol=15.0)
    @test isapprox(s22[14], s22interp[14], atol=32.0)
    @test isapprox(s22[15], s22interp[15], atol=1.0) # close to zero

    @test isapprox(s12[1], s12interp[1], rtol=0.10)
    @test isapprox(s12[2], s12interp[2], rtol=0.12)
    @test isapprox(s12[3], s12interp[3], rtol=0.13)
    @test isapprox(s12[4], s12interp[4], rtol=0.13)
    @test isapprox(s12[5], s12interp[5], rtol=0.13)
    @test isapprox(s12[6], s12interp[6], rtol=0.14)
    @test isapprox(s12[7], s12interp[7], atol=1.0) # close to zero
    @test isapprox(s12[8], s12interp[8], rtol=0.3)
    @test isapprox(s12[9], s12interp[9], rtol=0.25)
    @test isapprox(s12[10], s12interp[10], rtol=0.2)
    @test isapprox(s12[11], s12interp[11], rtol=0.16)
    @test isapprox(s12[12], s12interp[12], rtol=0.13)
    @test isapprox(s12[13], s12interp[13], rtol=0.1)
    @test isapprox(s12[14], s12interp[14], rtol=0.1)
    @test isapprox(s12[15], s12interp[15], atol=130) # close to zero

end

@testset "strain recovery: cantilever beam coupled with beam analysis" begin

    # Simple cantilever beam with tip load
    beam_length = 1.0
    num_beam_elements = 10

    # define material - aluminum
    E1 = 70e9
    E2 = 70e9
    E3 = 70e9
    G12 = 27e9
    G13 = 27e9
    G23 = 27e9
    nu12 = 0.3
    nu13 = 0.3
    nu23 = 0.3
    rho = 2.7e3

    material = GXBeam.Material(E1, E2, E3, G12, G13, G23,
                                    nu12, nu13, nu23, rho)

    # model square cross section mesh
    nx = 20
    ny = 20

    xs = range(0.0, stop=beam_length/10, length=nx+1)
    ys = range(0.0, stop=beam_length/10, length=ny+1)

    nodes = [GXBeam.Node(xs[i], ys[j]) for j in 1:ny+1 for i in 1:nx+1]
    elements = [GXBeam.MeshElement([i+(j-1)*(nx+1),
                                    i+1+(j-1)*(nx+1),
                                    nx+2+i+(j-1)*(nx+1),
                                    nx+1+i+(j-1)*(nx+1)], material, 0.0)
                for i in 1:nx for j in 1:ny]

    # get compliance, mass matrices
    cache = initialize_cache(nodes, elements)
    compliance = [GXBeam.compliance_matrix(nodes, elements; cache, gxbeam_order=true, shear_center=false)[1]
                    for i in 1:num_beam_elements]
    mass = [GXBeam.mass_matrix(nodes, elements)[1] for i in 1:num_beam_elements]

    # model cantilever beam in GXBeam
    xb = range(0.0, stop=beam_length, length=num_beam_elements+1)
    yb = zero(xb)
    zb = zero(xb)
    points = [[xb[i],yb[i],zb[i]] for i = 1:lastindex(xb)]

    start = 1:num_beam_elements
    stop = 2:num_beam_elements+1

    assembly = GXBeam.Assembly(points, start, stop; compliance, mass)

    # apply a point load at the tip
    Fx = 0.0
    Fy = 0.0
    Fz = -3.0

    prescribed_conditions = Dict(
        1 => GXBeam.PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
        num_beam_elements+1 => GXBeam.PrescribedConditions(Fx=Fx, Fy=Fy, Fz=Fz)
        )

    #solve GXBeam
    system, state, converged = static_analysis(assembly;
                                                prescribed_conditions,
                                                linear=true)

    # internal reaction loads at beam root
    F_GXBeam = -state.points[1].F
    M_GXBeam = -state.points[1].M

    # run GXBeam strain recovery
    strain_beam, stress_beam,
        strain_ply, stress_ply = GXBeam.strain_recovery(F_GXBeam, M_GXBeam, nodes, elements, cache;
                                                        gxbeam_order=true)

    #strains
    @test isapprox(strain_beam[1,1], -2.44191e-7, rtol=1e-4) #axial strain is in the first row now to match GXBeam coordinate system
    @test isapprox(strain_beam[1,400], 2.44191e-7, rtol=1e-4)
    @test isapprox(strain_beam[2,6], 3.47086e-8, rtol=1e-4)
    @test isapprox(strain_beam[3,3], 5.77984e-8, rtol=1e-4)
    @test isapprox(strain_beam[4,2], -1.51972e-8, rtol=1e-4)
    @test isapprox(strain_beam[5,5], -1.03921e-11, rtol=1e-4)
    @test isapprox(strain_beam[6,1], 4.22068e-9, rtol=1e-4)

    #stresses
    @test isapprox(stress_beam[1,1], -17090.9, rtol=1e-4) #axial stress is in the first row now to match GXBeam coordinate system
    @test isapprox(stress_beam[1,400], 17090.9, rtol=1e-4)
    @test isapprox(stress_beam[2,6], -0.018244, rtol=1e-4)
    @test isapprox(stress_beam[3,3], -2.78232, rtol=1e-4)
    @test isapprox(stress_beam[4,2], -410.326, rtol=1e-4)
    @test isapprox(stress_beam[5,5], -0.280588, rtol=1e-4)
    @test isapprox(stress_beam[6,1], 113.958, rtol=1e-4)

    end

function sectionwrapper(x)

    TF = eltype(x)

    xaf = TF[1.00000000, 0.99619582, 0.98515158, 0.96764209, 0.94421447, 0.91510964, 0.88074158, 0.84177999, 0.79894110, 0.75297076, 0.70461763, 0.65461515, 0.60366461, 0.55242353, 0.50149950, 0.45144530, 0.40276150, 0.35589801, 0.31131449, 0.26917194, 0.22927064, 0.19167283, 0.15672257, 0.12469599, 0.09585870, 0.07046974, 0.04874337, 0.03081405, 0.01681379, 0.00687971, 0.00143518, 0.00053606, 0.00006572, 0.00001249, 0.00023032, 0.00079945, 0.00170287, 0.00354717, 0.00592084, 0.01810144, 0.03471169, 0.05589286, 0.08132751, 0.11073805, 0.14391397, 0.18067874, 0.22089879, 0.26433734, 0.31062190, 0.35933893, 0.40999990, 0.46204424, 0.51483073, 0.56767889, 0.61998250, 0.67114514, 0.72054815, 0.76758733, 0.81168064, 0.85227225, 0.88883823, 0.92088961, 0.94797259, 0.96977487, 0.98607009, 0.99640466, 1.00000000]
    yaf = TF[0.00000000, 0.00017047, 0.00100213, 0.00285474, 0.00556001, 0.00906779, 0.01357364, 0.01916802, 0.02580144, 0.03334313, 0.04158593, 0.05026338, 0.05906756, 0.06766426, 0.07571157, 0.08287416, 0.08882939, 0.09329359, 0.09592864, 0.09626763, 0.09424396, 0.09023579, 0.08451656, 0.07727756, 0.06875796, 0.05918984, 0.04880096, 0.03786904, 0.02676332, 0.01592385, 0.00647946, 0.00370956, 0.00112514, -0.00046881, -0.00191488, -0.00329201, -0.00470585, -0.00688469, -0.00912202, -0.01720842, -0.02488211, -0.03226730, -0.03908459, -0.04503763, -0.04986836, -0.05338180, -0.05551392, -0.05636585, -0.05605816, -0.05472399, -0.05254383, -0.04969990, -0.04637175, -0.04264894, -0.03859653, -0.03433153, -0.02996944, -0.02560890, -0.02134397, -0.01726049, -0.01343567, -0.00993849, -0.00679919, -0.00402321, -0.00180118, -0.00044469, 0.00000000]

    xaf[15] = x[1]
    yaf[15] = x[2]

    uni = Material{TF}(37.00e9, 9.00e9, 9.00e9, 4.00e9, 4.00e9, 4.00e9, 0.28, 0.28, 0.28, 1.86e3)
    double = Material{TF}(x[3], 10.30e9, 10.30e9, 8.00e9, 8.00e9, 8.00e9, 0.30, 0.30, 0.30, 1.83e3)
    gelcoat = Material{TF}(1e1, 1e1, 1e1, 1.0, 1.0, 1.0, 0.30, 0.30, 0.30, 1.83e3)
    nexus = Material{TF}(10.30e9, 10.30e9, 10.30e9, 8.00e9, 8.00e9, 8.00e9, 0.30, 0.30, 0.30, 1.664e3)
    balsa = Material{TF}(0.01e9, 0.01e9, 0.01e9, 2e5, 2e5, 2e5, 0.30, 0.30, 0.30, 0.128e3)
    mat = [uni, double, gelcoat, nexus, balsa]

    chord = x[4]
    twist = x[5]
    paxis = 0.4750 / chord
    xbreak = [0.0, 0.0041, 0.1147, 0.5366, 1.0]
    webloc = [0.15 0.5]

    idx = [3, 4, 2]
    t = TF[0.000381, 0.00051, 18*0.00053]
    theta = TF[0, 0, 20]*pi/180
    layup1 = Layer.(mat[idx], t, theta)
    idx = [3, 4, 2]
    t = TF[0.000381, 0.00051, 33*0.00053]
    theta = TF[0, 0, 20]*pi/180
    layup2 = Layer.(mat[idx], t, theta)
    idx = [3, 4, 2, 1, 5, 1, 2]
    t = TF[0.000381, 0.00051, x[6], 38*0.00053, 1*0.003125, 37*0.00053, 16*0.00053]
    theta = TF[0, 0, x[7], 30, 0, 30, 20]*pi/180
    layup3 = Layer.(mat[idx], t, theta)
    idx = [3, 4, 2, 5, 2]
    t = TF[0.000381, 0.00051, 17*0.00053, 0.003125, 16*0.00053]
    theta = TF[0, 0, 20, 0, 0]*pi/180
    layup4 = Layer.(mat[idx], t, theta)

    idx = [1, 5, 1]
    t = TF[x[8], 0.003125, 38*0.00053]
    theta = TF[x[9], 0, 0]*pi/180
    web = Layer.(mat[idx], t, theta)

    segments = [layup1, layup2, layup3, layup4]

    webs = [web, web]

    nodes, elements = afmesh(xaf, yaf, chord, twist, paxis, xbreak, webloc, segments, webs)

    cache = initialize_cache(TF, nodes, elements)
    S, sc, tc = compliance_matrix(nodes, elements, cache=cache)
    M, mc = mass_matrix(nodes, elements)

    return vcat([S; M]...)
end

function sectionwrapper_nowebs(x)

    TF = eltype(x)

    xaf = TF[1.00000000, 0.99619582, 0.98515158, 0.96764209, 0.94421447, 0.91510964, 0.88074158, 0.84177999, 0.79894110, 0.75297076, 0.70461763, 0.65461515, 0.60366461, 0.55242353, 0.50149950, 0.45144530, 0.40276150, 0.35589801, 0.31131449, 0.26917194, 0.22927064, 0.19167283, 0.15672257, 0.12469599, 0.09585870, 0.07046974, 0.04874337, 0.03081405, 0.01681379, 0.00687971, 0.00143518, 0.00053606, 0.00006572, 0.00001249, 0.00023032, 0.00079945, 0.00170287, 0.00354717, 0.00592084, 0.01810144, 0.03471169, 0.05589286, 0.08132751, 0.11073805, 0.14391397, 0.18067874, 0.22089879, 0.26433734, 0.31062190, 0.35933893, 0.40999990, 0.46204424, 0.51483073, 0.56767889, 0.61998250, 0.67114514, 0.72054815, 0.76758733, 0.81168064, 0.85227225, 0.88883823, 0.92088961, 0.94797259, 0.96977487, 0.98607009, 0.99640466, 1.00000000]
    yaf = TF[0.00000000, 0.00017047, 0.00100213, 0.00285474, 0.00556001, 0.00906779, 0.01357364, 0.01916802, 0.02580144, 0.03334313, 0.04158593, 0.05026338, 0.05906756, 0.06766426, 0.07571157, 0.08287416, 0.08882939, 0.09329359, 0.09592864, 0.09626763, 0.09424396, 0.09023579, 0.08451656, 0.07727756, 0.06875796, 0.05918984, 0.04880096, 0.03786904, 0.02676332, 0.01592385, 0.00647946, 0.00370956, 0.00112514, -0.00046881, -0.00191488, -0.00329201, -0.00470585, -0.00688469, -0.00912202, -0.01720842, -0.02488211, -0.03226730, -0.03908459, -0.04503763, -0.04986836, -0.05338180, -0.05551392, -0.05636585, -0.05605816, -0.05472399, -0.05254383, -0.04969990, -0.04637175, -0.04264894, -0.03859653, -0.03433153, -0.02996944, -0.02560890, -0.02134397, -0.01726049, -0.01343567, -0.00993849, -0.00679919, -0.00402321, -0.00180118, -0.00044469, 0.00000000]

    xaf[15] = x[1]
    yaf[15] = x[2]

    uni = Material{TF}(37.00e9, 9.00e9, 9.00e9, 4.00e9, 4.00e9, 4.00e9, 0.28, 0.28, 0.28, 1.86e3)
    double = Material{TF}(x[3], 10.30e9, 10.30e9, 8.00e9, 8.00e9, 8.00e9, 0.30, 0.30, 0.30, 1.83e3)
    gelcoat = Material{TF}(1e1, 1e1, 1e1, 1.0, 1.0, 1.0, 0.30, 0.30, 0.30, 1.83e3)
    nexus = Material{TF}(10.30e9, 10.30e9, 10.30e9, 8.00e9, 8.00e9, 8.00e9, 0.30, 0.30, 0.30, 1.664e3)
    balsa = Material{TF}(0.01e9, 0.01e9, 0.01e9, 2e5, 2e5, 2e5, 0.30, 0.30, 0.30, 0.128e3)
    mat = [uni, double, gelcoat, nexus, balsa]

    chord = x[4]
    twist = x[5]
    paxis = 0.4750 / chord
    xbreak = [0.0, 0.0041, 0.1147, 0.5366, 1.0]
    webloc = []

    idx = [3, 4, 2]
    t = TF[0.000381, 0.00051, 18*0.00053]
    theta = TF[0, 0, 20]*pi/180
    layup1 = Layer.(mat[idx], t, theta)
    idx = [3, 4, 2]
    t = TF[0.000381, 0.00051, 33*0.00053]
    theta = TF[0, 0, 20]*pi/180
    layup2 = Layer.(mat[idx], t, theta)
    idx = [3, 4, 2, 1, 5, 1, 2]
    t = TF[0.000381, 0.00051, x[6], 38*0.00053, 1*0.003125, 37*0.00053, 16*0.00053]
    theta = TF[0, 0, x[7], 30, 0, 30, 20]*pi/180
    layup3 = Layer.(mat[idx], t, theta)
    idx = [3, 4, 2, 5, 2]
    t = TF[0.000381, 0.00051, 17*0.00053, 0.003125, 16*0.00053]
    theta = TF[0, 0, 20, 0, 0]*pi/180
    layup4 = Layer.(mat[idx], t, theta)

    segments = [layup1, layup2, layup3, layup4]

    webs = []

    nodes, elements = afmesh(xaf, yaf, chord, twist, paxis, xbreak, webloc, segments, webs)

    cache = initialize_cache(TF, nodes, elements)
    S, sc, tc = compliance_matrix(nodes, elements, cache=cache)
    M, mc = mass_matrix(nodes, elements)

    return vcat([S; M]...)
end

@testset "type stability" begin

    function checkstability()
        xvec = [0.5014995, 0.07571157, 10.30e9, 1.9, 0.0*pi/180, 17*0.00053, 20.0, 38*0.00053, 0.0]
        try
            @inferred Vector{Float64} sectionwrapper(xvec)
            return true
        catch err
            # println(err)
            return false
        end
    end

    @test checkstability()

end

@testset "Jacobian with Webs" begin

    # should use graph coloring b.c. plenty of sparsity, but dense is fine for purpose of this test.
    xvec = [0.5014995, 0.07571157, 10.30e9, 1.9, 0.0*pi/180, 17*0.00053, 20.0, 38*0.00053, 0.0]

    J1 = ForwardDiff.jacobian(sectionwrapper, xvec)
    J2 = FiniteDiff.finite_difference_jacobian(sectionwrapper, xvec, Val{:central})

    @test maximum(abs.(J1 .- J2)) < 1e-6

end

@testset "Jacobian without Webs" begin

    xvec = [0.5014995, 0.07571157, 10.30e9, 1.9, 0.0*pi/180, 17*0.00053, 20.0, 38*0.00053, 0.0]

    J1 = ForwardDiff.jacobian(sectionwrapper_nowebs, xvec)
    J2 = FiniteDiff.finite_difference_jacobian(sectionwrapper_nowebs, xvec, Val{:central})

    @test maximum(abs.(J1 .- J2)) < 1e-6

end