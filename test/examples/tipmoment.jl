using GXBeam, LinearAlgebra, Test

@testset "Nonlinear Analysis of a Cantilever Subjected to a Constant Moment" begin

    L = 12 # inches
    h = w = 1 # inches
    E = 30e6 # lb/in^4 Young's Modulus

    A = h*w
    Iyy = w*h^3/12
    Izz = w^3*h/12

    # bending moment (applied at end)
    # note that solutions for λ > 1.8 do not converge
    λ = [0.0, 0.4, 0.8, 1.2, 1.6, 1.8, 2.0]
    m = pi*E*Iyy/L
    M = λ*m

    # create points
    nelem = 16
    x = range(0, L, length=nelem+1)
    y = zero(x)
    z = zero(x)
    points = [[x[i],y[i],z[i]] for i = 1:length(x)]

    # index of endpoints for each beam element
    start = 1:nelem
    stop = 2:nelem+1

    # compliance matrix for each beam element
    compliance = fill(Diagonal([1/(E*A), 0, 0, 0, 1/(E*Iyy), 1/(E*Izz)]), nelem)

    # create assembly of interconnected nonlinear beams
    assembly = Assembly(points, start, stop, compliance=compliance)

    # pre-initialize system storage
    system = StaticSystem(assembly)

    # run an analysis for each prescribed bending moment
    states = Vector{AssemblyState{Float64}}(undef, length(M))
    for i = 1:length(M)

        # create dictionary of prescribed conditions
        prescribed_conditions = Dict(
            # fixed left side
            1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
            # moment on right side
            nelem+1 => PrescribedConditions(Mz = M[i])
        )

        # perform a static analysis
        static_analysis!(system, assembly, prescribed_conditions=prescribed_conditions)

        # post-process the results
        states[i] = AssemblyState(system, assembly, prescribed_conditions=prescribed_conditions)

    end

    # analytical solution (ρ = E*I/M)
    analytical(x, ρ) = ifelse(ρ == Inf, zeros(3), [ρ*sin(x/ρ)-x, ρ*(1-cos(x/ρ)), 0])

    # test element properties
    for i = 1:length(M)
        for ielem = 1:length(assembly.elements)
            xi = assembly.elements[ielem].x[1]
            u_a, v_a, w_a = analytical(xi, E*Iyy/M[i])
            @test isapprox(states[i].elements[ielem].u[1], u_a, atol=5e-2)
            @test isapprox(states[i].elements[ielem].u[2], v_a, atol=5e-2)
        end

        # test point properties
        for ipoint = 1:length(assembly.points)
            xi = assembly.points[ipoint][1]
            u_a, v_a, w_a = analytical(xi, E*Iyy/M[i])
            @test isapprox(states[i].points[ipoint].u[1], u_a, atol=5e-2)
            @test isapprox(states[i].points[ipoint].u[2], v_a, atol=5e-2)
        end
    end

    # pre-initialize system storage
    system = ExpandedSystem(assembly)

    # perform the same analysis for a constant mass matrix system
    states = Vector{AssemblyState{Float64}}(undef, length(M))
    for i = 1:length(M)

        prescribed_conditions = Dict(
            1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
            nelem+1 => PrescribedConditions(Mz = M[i])
        )

        steady_state_analysis!(system, assembly, 
            prescribed_conditions = prescribed_conditions,
            constant_mass_matrix = true)

        states[i] = AssemblyState(system, assembly, prescribed_conditions=prescribed_conditions)
    end

    for i = 1:length(M)
        # test element properties
        for ielem = 1:length(assembly.elements)
            xi = assembly.elements[ielem].x[1]
            u_a, v_a, w_a = analytical(xi, E*Iyy/M[i])
            @test isapprox(states[i].elements[ielem].u[1], u_a, atol=5e-2)
            @test isapprox(states[i].elements[ielem].u[2], v_a, atol=5e-2)
        end

        # test point properties
        for ipoint = 1:length(assembly.points)
            xi = assembly.points[ipoint][1]
            u_a, v_a, w_a = analytical(xi, E*Iyy/M[i])
            @test isapprox(states[i].points[ipoint].u[1], u_a, atol=5e-2)
            @test isapprox(states[i].points[ipoint].u[2], v_a, atol=5e-2)
        end
    end
end