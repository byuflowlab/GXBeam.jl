# Unit Tests for section.jl and afmesh.jl

@testset "section properties: material stiffness matrix" begin
    
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
    cache = initialize_cache([Node(0.0, 0.0), Node(0.0, 0.0), Node(0.0, 0.0), Node(0.0, 0.0)], [MeshElement([1, 2, 3, 4], mat, 0.0)])
    GXBeam.stiffness!(mat, cache)
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

end


@testset "section properties: square cross sections" begin
    
    # BECAS User Guide
    #  Case 1: Square cross section of isotropic material - S1

    iso1 = Material(100.0, 100.0, 100.0, 41.667, 41.667, 41.667, 0.2, 0.2, 0.2, 1.0)
    x = range(-0.05, 0.05, length=11)
    y = range(-0.05, 0.05, length=11)

    nodes = Vector{Node{Float64}}(undef, 11*11)
    elements = Vector{MeshElement{Vector{Int64}, Float64}}(undef, 10*10)

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

    S, sc, tc = compliance_matrix(nodes, elements, gxbeam_order=false)
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


    S, sc, tc = compliance_matrix(nodes, elements, gxbeam_order=false)
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

    S, sc, tc = compliance_matrix(nodes, elements, gxbeam_order=false)
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

    S, sc, tc = compliance_matrix(nodes, elements, gxbeam_order=false)
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

    S, sc, tc = compliance_matrix(nodes, elements, gxbeam_order=false)
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
    elements = Vector{MeshElement{Vector{Int64},Float64}}(undef, (nr-1)*(nt-1))
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

    S, sc, tc = compliance_matrix(nodes, elements, gxbeam_order=false)
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
    elements = Vector{MeshElement{Vector{Int64},Float64}}(undef, (nr-1)*(nt-1))
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

    S, sc, tc = compliance_matrix(nodes, elements, gxbeam_order=false)
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
    elements = Vector{MeshElement{Vector{Int64},Float64}}(undef, (nr-1)*(nt-1))
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


@testset "section properties: multi-layer composite pipe" begin
    
    #  --- Generalized Timoshenko Theory of the Variational Asymptotic Beam Sectional Analysis ----
    # multi-layer composite pipe 
    # note that the previous paper has a composite pipe also, but the numbers are inconsistent.
    # see also preVABS documentation examples

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
    nr = 20
    nodes = Vector{Node{Float64}}(undef, (2*(nx-1) + 2*(nt-1))*nr)
    elements = Vector{MeshElement{Vector{Int64},Float64}}(undef, (2*(nx-1) + 2*(nt-1))*(nr-1))

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
                if j <= (nr ÷ 2)
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
                if j <= (nr ÷ 2)
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
                if j <= (nr ÷ 2)
                    theta = 45*pi/180
                else
                    theta = -45*pi/180
                end
                elements[n] = MeshElement([nr*(I+i-1)+j, nr*(Ip+i)+j, nr*(Ip+i)+j+1, nr*(I+i-1)+j+1], pipemat, theta)
                n += 1
            end
        end

    end

    # plotmesh(nodes, elements)

    S, sc, tc = compliance_matrix(nodes, elements)
    K = inv(S)

    @test isapprox(K[1, 1], 1.03892e7, rtol=0.04)
    @test isapprox(K[2, 2], 7.85310e5, rtol=0.01)
    @test isapprox(K[3, 3], 3.29279e5, rtol=0.02)
    @test isapprox(K[1, 4], 9.84575e4, rtol=0.12)
    @test isapprox(K[2, 5], -8.21805e3, rtol=0.11)
    @test isapprox(K[3, 6], -5.20981e4, rtol=0.21)   
    @test isapprox(K[4, 4], 6.87275e5, rtol=0.01)
    @test isapprox(K[5, 5], 1.88238e6, rtol=0.04)
    @test isapprox(K[6, 6], 5.38987e6, rtol=0.03)

    # println("K11 = ", round((K2[1, 1]/1.03892e7 - 1)*100, digits=2), "%")
    # println("K22 = ", round((K2[2, 2]/7.85310e5 - 1)*100, digits=2), "%")
    # println("K33 = ", round((K2[3, 3]/3.29279e5 - 1)*100, digits=2), "%")
    # println("K14 = ", round((K2[1, 4]/9.84575e4 - 1)*100, digits=2), "%")
    # println("K25 = ", round((K2[2, 5]/-8.21805e3 - 1)*100, digits=2), "%")
    # println("K36 = ", round((K2[3, 6]/-5.20981e4 - 1)*100, digits=2), "%")
    # println("K44 = ", round((K2[4, 4]/6.87275e5 - 1)*100, digits=2), "%")
    # println("K55 = ", round((K2[5, 5]/1.88238e6 - 1)*100, digits=2), "%")
    # println("K66 = ", round((K2[6, 6]/5.38987e6 - 1)*100, digits=2), "%")
end




@testset "section properties: composite wind turbine airfoil" begin
    # A critical assessment of computer tools for calculating composite wind turbine blade properties, Chen, Yu, Capellaro
    # ST1
    
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
    @test isapprox(log10(K[1, 2]), log10(abs(1.524e6)), rtol=0.03)
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
    @test isapprox(theta*180/pi, -1.244, rtol=0.09)
    
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
    
    S, sc, tc = compliance_matrix(nodes, elements)
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
            println(err)
            return false
        end
    end

    @test checkstability()

end

using ForwardDiff
using FiniteDiff

@testset "Jacobian" begin

    # should use graph coloring b.c. plenty of sparsity, but dense is fine for purpose of this test.
    xvec = [0.5014995, 0.07571157, 10.30e9, 1.9, 0.0*pi/180, 17*0.00053, 20.0, 38*0.00053, 0.0]
    
    J1 = ForwardDiff.jacobian(sectionwrapper, xvec)
    J2 = FiniteDiff.finite_difference_jacobian(sectionwrapper, xvec, Val{:central})

    @test maximum(abs.(J1 .- J2)) < 1e-6

end
