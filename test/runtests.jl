using GXBeam
using LinearAlgebra
using DifferentialEquations
using Test
import Elliptic
using ForwardDiff

@testset "Math" begin
    
    c = rand(3)
    cdot = rand(3)

    # get_C_θ
    C_θ1, C_θ2, C_θ3 = GXBeam.get_C_θ(c)
    @test isapprox(C_θ1, ForwardDiff.derivative(c1 -> GXBeam.get_C([c1, c[2], c[3]]), c[1]))
    @test isapprox(C_θ2, ForwardDiff.derivative(c2 -> GXBeam.get_C([c[1], c2, c[3]]), c[2]))
    @test isapprox(C_θ3, ForwardDiff.derivative(c3 -> GXBeam.get_C([c[1], c[2], c3]), c[3]))

    # get_C_t_θ
    Cdot_θ1, Cdot_θ2, Cdot_θ3 = GXBeam.get_C_t_θ(c, cdot)
    @test isapprox(Cdot_θ1, ForwardDiff.derivative(c1 -> GXBeam.get_C_t([c1, c[2], c[3]], cdot), c[1]))
    @test isapprox(Cdot_θ2, ForwardDiff.derivative(c2 -> GXBeam.get_C_t([c[1], c2, c[3]], cdot), c[2]))
    @test isapprox(Cdot_θ3, ForwardDiff.derivative(c3 -> GXBeam.get_C_t([c[1], c[2], c3], cdot), c[3]))

    # get_C_t_θdot
    Cdot_θdot1, Cdot_θdot2, Cdot_θdot3 = GXBeam.get_C_t_θdot(c)
    @test isapprox(Cdot_θdot1, ForwardDiff.derivative(cdot1 -> GXBeam.get_C_t(c, [cdot1, cdot[2], cdot[3]]), cdot[1]))
    @test isapprox(Cdot_θdot2, ForwardDiff.derivative(cdot2 -> GXBeam.get_C_t(c, [cdot[1], cdot2, cdot[3]]), cdot[2]))
    @test isapprox(Cdot_θdot3, ForwardDiff.derivative(cdot3 -> GXBeam.get_C_t(c, [cdot[1], cdot[2], cdot3]), cdot[3]))
    
    # get_Q_θ
    Q_θ1, Q_θ2, Q_θ3 = GXBeam.get_Q_θ(c)
    @test isapprox(Q_θ1, ForwardDiff.derivative(c1 -> GXBeam.get_Q([c1, c[2], c[3]]), c[1]))
    @test isapprox(Q_θ2, ForwardDiff.derivative(c2 -> GXBeam.get_Q([c[1], c2, c[3]]), c[2]))
    @test isapprox(Q_θ3, ForwardDiff.derivative(c3 -> GXBeam.get_Q([c[1], c[2], c3]), c[3]))

    # get_Qinv_θ
    Qinv_θ1, Qinv_θ2, Qinv_θ3 = GXBeam.get_Qinv_θ(c)
    @test isapprox(Qinv_θ1, ForwardDiff.derivative(c1 -> GXBeam.get_Qinv([c1, c[2], c[3]]), c[1]))
    @test isapprox(Qinv_θ2, ForwardDiff.derivative(c2 -> GXBeam.get_Qinv([c[1], c2, c[3]]), c[2]))
    @test isapprox(Qinv_θ3, ForwardDiff.derivative(c3 -> GXBeam.get_Qinv([c[1], c[2], c3]), c[3]))

end

@testset "Jacobian and Mass Matrix Calculations" begin

    L = 60 # m

    # create points
    nelem = 1
    x = range(0, L, length=nelem+1)
    y = zero(x)
    z = zero(x)
    points = [[x[i],y[i],z[i]] for i = 1:length(x)]

    # index of endpoints of each beam element
    start = 1:nelem
    stop = 2:nelem+1

    # stiffness matrix for each beam element
    stiffness = fill(
        [2.389e9  1.524e6  6.734e6 -3.382e7 -2.627e7 -4.736e8
         1.524e6  4.334e8 -3.741e6 -2.935e5  1.527e7  3.835e5
         6.734e6 -3.741e6  2.743e7 -4.592e5 -6.869e5 -4.742e6
        -3.382e7 -2.935e5 -4.592e5  2.167e7 -6.279e5  1.430e6
        -2.627e7  1.527e7 -6.869e5 -6.279e5  1.970e7  1.209e7
        -4.736e8  3.835e5 -4.742e6  1.430e6  1.209e7  4.406e8],
        nelem)

    # mass matrix for each beam element
    mass = fill(
        [258.053      0.0        0.0      0.0      7.07839  -71.6871
           0.0      258.053      0.0     -7.07839  0.0        0.0
           0.0        0.0      258.053   71.6871   0.0        0.0
           0.0       -7.07839   71.6871  48.59     0.0        0.0
           7.07839    0.0        0.0      0.0      2.172      0.0
         -71.6871     0.0        0.0      0.0      0.0       46.418],
         nelem)

    # create assembly of interconnected nonlinear beams
    assembly = Assembly(points, start, stop; stiffness=stiffness, mass=mass)

    # prescribed conditions
    pcond = Dict(
        # fixed left side
        1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
        )

    # distributed loads
    dload = Dict()

    # point masses
    pmass = Dict(
        # point mass at the end of the beam
        nelem => PointMass(Symmetric(rand(6,6)))
    )

    # gravity vector
    gvec = rand(3)

    # --- Static System --- #
    static_system = System(assembly, true)

    force_scaling = static_system.force_scaling
    irow_point = static_system.irow_point
    irow_elem = static_system.irow_elem
    irow_elem1 = static_system.irow_elem1
    irow_elem2 = static_system.irow_elem2
    icol_point = static_system.icol_point
    icol_elem = static_system.icol_elem

    x = rand(length(static_system.x))
    J = similar(x, length(x), length(x))

    f = (x) -> GXBeam.system_residual!(similar(x), x, assembly, pcond, dload, pmass, gvec,
        force_scaling, irow_point, irow_elem, irow_elem1, irow_elem2,
        icol_point, icol_elem)

    GXBeam.system_jacobian!(J, x, assembly, pcond, dload, pmass, gvec, force_scaling,
        irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem)

    @test all(isapprox.(J, ForwardDiff.jacobian(f, x), atol=1e-10))

    # --- Dynamic System --- #
    system = System(assembly, false)

    force_scaling = system.force_scaling
    irow_point = system.irow_point
    irow_elem = system.irow_elem
    irow_elem1 = system.irow_elem1
    irow_elem2 = system.irow_elem2
    icol_point = system.icol_point
    icol_elem = system.icol_elem
    x0 = rand(3)
    v0 = rand(3)
    ω0 = rand(3)
    a0 = rand(3)
    α0 = rand(3)

    x = rand(length(system.x))
    J = similar(x, length(x), length(x))

    f = (x) -> GXBeam.system_residual!(similar(x), x, assembly, pcond, dload, pmass, gvec,
        force_scaling, irow_point, irow_elem, irow_elem1, irow_elem2,
        icol_point, icol_elem, x0, v0, ω0, a0, α0)

    GXBeam.system_jacobian!(J, x, assembly, pcond, dload, pmass, gvec, force_scaling,
        irow_point, irow_elem, irow_elem1, irow_elem2, icol_point,
        icol_elem, x0, v0, ω0, a0, α0)

    @test all(isapprox.(J, ForwardDiff.jacobian(f, x), atol=1e-10))

    # initial condition residual

    u0 = [rand(3) for ielem = 1:length(assembly.elements)]
    theta0 = [rand(3) for ielem = 1:length(assembly.elements)]
    udot0 = [rand(3) for ielem = 1:length(assembly.elements)]
    thetadot0 = [rand(3) for ielem = 1:length(assembly.elements)]

    x = rand(length(system.x))
    J = similar(x, length(x), length(x))

    f = (x) -> GXBeam.system_residual!(similar(x), x, assembly, pcond, dload, pmass, gvec,
        force_scaling, irow_point, irow_elem, irow_elem1, irow_elem2,
        icol_point, icol_elem, x0, v0, ω0, a0, α0, u0, theta0, udot0, thetadot0)

    GXBeam.system_jacobian!(J, x, assembly, pcond, dload, pmass, gvec, force_scaling,
        irow_point, irow_elem, irow_elem1, irow_elem2, icol_point,
        icol_elem, x0, v0, ω0, a0, α0, u0, theta0, udot0, thetadot0)

    @test all(isapprox.(J, ForwardDiff.jacobian(f, x), atol=1e-10))

    # time domain residual

    udot = [rand(3) for ielem = 1:length(assembly.elements)]
    θdot = [rand(3) for ielem = 1:length(assembly.elements)]
    Pdot = [rand(3) for ielem = 1:length(assembly.elements)]
    Hdot = [rand(3) for ielem = 1:length(assembly.elements)]
    dt = rand()

    x = rand(length(system.x))
    J = similar(x, length(x), length(x))

    f = (x) -> GXBeam.system_residual!(similar(x), x, assembly, pcond, dload, pmass, gvec,
        force_scaling, irow_point, irow_elem, irow_elem1, irow_elem2,
        icol_point, icol_elem, x0, v0, ω0, a0, α0, udot, θdot, Pdot, Hdot, dt)

    GXBeam.system_jacobian!(J, x, assembly, pcond, dload, pmass, gvec, force_scaling,
        irow_point, irow_elem, irow_elem1, irow_elem2, icol_point,
        icol_elem, x0, v0, ω0, a0, α0, udot, θdot, Pdot, Hdot, dt)

    @test all(isapprox.(J, ForwardDiff.jacobian(f, x), atol=1e-10))

    # system jacobian and mass matrix

    x = rand(length(system.x))
    dx = rand(length(x))
    J = similar(x, length(x), length(x))
    M = similar(x, length(x), length(x))

    fx = (x) -> GXBeam.dynamic_system_residual!(similar(x), x, dx, assembly, pcond, dload, 
        pmass, gvec, force_scaling, irow_point, irow_elem, irow_elem1, irow_elem2,
        icol_point, icol_elem, x0, v0, ω0, a0, α0)

    fdx = (dx) -> GXBeam.dynamic_system_residual!(similar(dx), x, dx, assembly, pcond, dload, 
        pmass, gvec, force_scaling, irow_point, irow_elem, irow_elem1, irow_elem2,
        icol_point, icol_elem, x0, v0, ω0, a0, α0)

    GXBeam.dynamic_system_jacobian!(J, x, dx, assembly, pcond, dload, pmass, gvec, force_scaling, 
        irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem,
        x0, v0, ω0, a0, α0)

    GXBeam.system_mass_matrix!(M, x, assembly, pmass, force_scaling, irow_point,
        irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem)

    @test all(isapprox.(J, ForwardDiff.jacobian(fx, x), atol=1e-10))

    @test all(isapprox.(M, ForwardDiff.jacobian(fdx, dx), atol=1e-10))

end

@testset "Linear Analysis of a Cantilever Partially Under a Uniform Distributed Load" begin

    nelem = 12

    # create points
    a = 0.3
    b = 0.7
    L = 1.0
    n1 = n3 = div(nelem, 3)
    n2 = nelem - n1 - n3
    x1 = range(0, a, length=n1+1)
    x2 = range(a, b, length=n2+1)
    x3 = range(b, L, length=n3+1)
    x = vcat(x1, x2[2:end], x3[2:end])
    y = zero(x)
    z = zero(x)
    points = [[x[i],y[i],z[i]] for i = 1:length(x)]

    # index of endpoints for each beam element
    start = 1:nelem
    stop = 2:nelem+1

    # create compliance matrix for each beam element
    EI = 1e9
    stiffness = fill(Diagonal([0, 0, 0, 0, EI, 0]), nelem)

    # create the assembly
    assembly = Assembly(points, start, stop, stiffness=stiffness)

    # set prescribed conditions (fixed right endpoint)
    prescribed_conditions = Dict(
        nelem+1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0)
    )

    # create distributed load
    q = 1000
    distributed_loads = Dict()
    for ielem in n1+1:n1+n2
        distributed_loads[ielem] = DistributedLoads(assembly, ielem; fz = (s) -> q)
    end

    system, converged = static_analysis(assembly, prescribed_conditions=prescribed_conditions,
        distributed_loads=distributed_loads, linear=true)

    state = AssemblyState(system, assembly, prescribed_conditions=prescribed_conditions)

    # analytical solution obtained using superposition
    initial_slope = -q/(6*EI)*((L-a)^3 - (L-b)^3)
    initial_deflection = q/(24*EI)*((L-a)^3*(3*L + a) - (L-b)^3*(3*L + b))
    analytical_M = function(x)
        if 0 < x <= a
            M = 0.0
        elseif a < x <= b
            M = q/2*(x-a)^2
        else
            M = q/2*((x-a)^2 - (x-b)^2)
        end
        return M
    end
    analytical_slope = function(x)
        slope = initial_slope
        if 0 < x <= a
            slope += 0.0
        elseif a < x <= b
            slope += q/(6*EI)*(x-a)^3
        else
            slope += q/(6*EI)*((x-a)^3 - (x-b)^3)
        end
        return slope
    end
    analytical_deflection = function(x)
        deflection = initial_deflection + initial_slope*x
        if 0 < x <= a
            deflection += 0.0
        elseif a < x <= b
            deflection += q/(24*EI)*(x-a)^4
        else
            deflection += q/(24*EI)*((x-a)^4 - (x-b)^4)
        end
        return deflection
    end

    # test element properties
    for i = 1:length(assembly.elements)
        xi = assembly.elements[i].x[1]
        @test isapprox(state.elements[i].u[3], analytical_deflection(xi), atol=1e-9)
        @test isapprox(state.elements[i].theta[2], -4*analytical_slope(xi)/4, atol=1e-9)
        @test isapprox(state.elements[i].M[2], -analytical_M(xi), atol=2)
    end

    # test point properties
    for i = 1:length(assembly.points)
        xi = assembly.points[i][1]
        @test isapprox(state.points[i].u[3], analytical_deflection(xi), atol=1e-8)
        @test isapprox(state.points[i].theta[2], -4*analytical_slope(xi)/4, atol=1e-7)
    end
end

@testset "Linear Analysis of a Beam Under a Linear Distributed Load" begin

    nelem = 16

    # create points
    L = 1
    x = range(0, L, length=nelem+1)
    y = zero(x)
    z = zero(x)
    points = [[x[i],y[i],z[i]] for i = 1:length(x)]

    # index of endpoints for each beam element
    start = 1:nelem
    stop = 2:nelem+1

    # create compliance matrix for each beam element
    EI = 1e7
    compliance = fill(Diagonal([0, 0, 0, 0, 1/EI, 0]), nelem)

    # create assembly
    assembly = Assembly(points, start, stop, compliance=compliance)

    # set prescribed conditions
    prescribed_conditions = Dict(
        # simply supported left endpoint
        1 => PrescribedConditions(uz=0),
        # clamped right endpoint
        nelem+1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0)
    )

    # create distributed load
    qmax = 1000
    distributed_loads = Dict()
    for i = 1:nelem
        distributed_loads[i] = DistributedLoads(assembly, i; s1=x[i],
            s2=x[i+1], fz = (s) -> qmax*s)
    end

    # solve system
    system, converged = static_analysis(assembly, prescribed_conditions=prescribed_conditions,
        distributed_loads=distributed_loads, linear=true)

    # post-process the results
    state = AssemblyState(system, assembly, prescribed_conditions=prescribed_conditions)

    # construct analytical solution
    analytical_deflection = (x) -> qmax*(1-x)^2/(120*EI)*(4 - 8*(1-x) + 5*(1-x)^2 - (1-x)^3)
    analytical_slope = (x) -> -qmax*(1-x)/(120*EI)*(8 - 24*(1-x) + 20*(1-x)^2 - 5*(1-x)^3)
    analytical_M = (x) -> qmax/120*(8 - 48*(1-x) + 60*(1-x)^2 - 20*(1-x)^3)

    # test element properties
    for i = 1:length(assembly.elements)
        xi = assembly.elements[i].x[1]
        @test isapprox(state.elements[i].u[3], analytical_deflection(xi), atol=1e-8)
        @test isapprox(state.elements[i].theta[2], -4*analytical_slope(xi)/4, atol=1e-7)
        @test isapprox(state.elements[i].M[2], -analytical_M(xi), atol=1)
    end

    # test point properties
    for i = 1:length(assembly.points)
        xi = assembly.points[i][1]
        @test isapprox(state.points[i].u[3], analytical_deflection(xi), atol=1e-8)
        @test isapprox(state.points[i].theta[2], -4*analytical_slope(xi)/4, atol=1e-8)
    end
end

@testset "Nonlinear Analysis of a Cantilever Subjected to a Constant Tip Load" begin

    L = 1
    EI = 1e6

    # shear force (applied at end)
    λ = 0:0.5:16
    p = EI/L^2
    P = λ*p

    # create points
    nelem = 16
    x = range(0, L, length=nelem+1)
    y = zero(x)
    z = zero(x)
    points = [[x[i],y[i],z[i]] for i = 1:length(x)]

    # index of endpoints of each beam element
    start = 1:nelem
    stop = 2:nelem+1

    # compliance matrix for each beam element
    compliance = fill(Diagonal([0, 0, 0, 0, 1/EI, 0]), nelem)

    # create assembly of interconnected nonlinear beams
    assembly = Assembly(points, start, stop, compliance=compliance)

    # pre-initialize system storage
    system = System(assembly, true)

    # run an analysis for each prescribed tip load
    states = Vector{AssemblyState{Float64}}(undef, length(P))
    for i = 1:length(P)

        # create dictionary of prescribed conditions
        prescribed_conditions = Dict(
            # fixed left side
            1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
            # shear force on right tip
            nelem+1 => PrescribedConditions(Fz = P[i])
        )

        # perform a static analysis
        static_analysis!(system, assembly, prescribed_conditions=prescribed_conditions)

        # post-process the results
        states[i] = AssemblyState(system, assembly, prescribed_conditions=prescribed_conditions)

    end

    # construct analytical solution
    δ = range(pi/4, pi/2, length=10^5)[2:end-1]

    k = @. cos(pi/4)/sin(δ)
    λ_a = @. (Elliptic.F(pi/2, k^2) - Elliptic.F(δ,  k^2))^2

    θ_a = @. 2*(pi/4 - acos(k))

    ξ_a = @. sqrt(2*sin(θ_a)/λ_a) .- 1

    η_a = @. 1-2/sqrt(λ_a)*(Elliptic.E(pi/2, k^2) - Elliptic.E(δ, k^2))

    # test tip displacements
    for i = 1:length(P)
        i_a = argmin(abs.(λ[i] .- λ_a))
        @test isapprox(states[i].points[end].u[1]/L, ξ_a[i_a], atol=1e-3)
        @test isapprox(states[i].points[end].u[3]/L, η_a[i_a], atol=1e-3)
        @test isapprox(states[i].points[end].theta[2], -4*tan(θ_a[i_a]/4), atol=1e-2)
    end
end

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
    system = System(assembly, true)

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
end

@testset "Nonlinear Analysis of the Bending of a Curved Beam in 3D Space" begin

    # problem constants
    R = 100
    L = R*pi/4 # inches
    h = w = 1 # inches
    E = 1e7 # psi Young's Modulus
    ν = 0.0
    G = E/(2*(1+ν))

    # beam starting point, frame, and curvature
    r = [0, 0, 0]
    frame = [0 -1 0; 1 0 0; 0 0 1]
    curvature = [0, 0, -1/R]

    # cross section properties
    A = h*w
    Ay = A
    Az = A
    Iyy = w*h^3/12
    Izz = w^3*h/12
    J = Iyy + Izz

    # discretize the beam
    nelem = 16
    ΔL, xp, xm, Cab = discretize_beam(L, r, nelem; frame=frame, curvature = curvature)

    # force
    P = 600 # lbs

    # index of left and right endpoints of each beam element
    start = 1:nelem
    stop = 2:nelem+1

    # compliance matrix for each beam element
    compliance = fill(Diagonal([1/(E*A), 1/(G*Ay), 1/(G*Az), 1/(G*J), 1/(E*Iyy), 1/(E*Izz)]), nelem)

    # create assembly of interconnected nonlinear beams
    assembly = Assembly(xp, start, stop, compliance=compliance, frames=Cab,
        lengths=ΔL, midpoints=xm)

    # create dictionary of prescribed conditions
    prescribed_conditions = Dict(
        # fixed left endpoint
        1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
        # force on right endpoint
        nelem+1 => PrescribedConditions(Fz=P)
    )

    # perform static analysis
    system, converged = static_analysis(assembly, prescribed_conditions=prescribed_conditions)

    # post-process results
    state = AssemblyState(system, assembly, prescribed_conditions=prescribed_conditions)

    # Results from "Large Displacement Analysis of Three-Dimensional Beam
    # Structures" by Bathe and Bolourch:
    # - Tip Displacement: [-13.4, -23.5, 53.4]

    # Note that these results are comparing computational solutions, rather than
    # the computational to the analytical solution, so some variation is expected.

    @test isapprox(state.points[end].u[1], -13.4, atol=0.2) # -13.577383726758564
    @test isapprox(state.points[end].u[2], -23.5, atol=0.1) # -23.545303336988038
    @test isapprox(state.points[end].u[3],  53.4, atol=0.1) #  53.45800757548929
end

@testset "Rotating Beam with a Swept Tip" begin
    sweep = 45 * pi/180
    rpm = 0:25:750

    # straight section of the beam
    L_b1 = 31.5 # inch
    r_b1 = [2.5, 0, 0]
    nelem_b1 = 13
    lengths_b1, xp_b1, xm_b1, Cab_b1 = discretize_beam(L_b1, r_b1, nelem_b1)

    # swept section of the beam
    L_b2 = 6 # inch
    r_b2 = [34, 0, 0]
    nelem_b2 = 3
    cs, ss = cos(sweep), sin(sweep)
    frame_b2 = [cs ss 0; -ss cs 0; 0 0 1]
    lengths_b2, xp_b2, xm_b2, Cab_b2 = discretize_beam(L_b2, r_b2, nelem_b2, frame=frame_b2)

    # combine elements and points into one array
    nelem = nelem_b1 + nelem_b2
    points = vcat(xp_b1, xp_b2[2:end])
    start = 1:nelem_b1 + nelem_b2
    stop = 2:nelem_b1 + nelem_b2 + 1
    lengths = vcat(lengths_b1, lengths_b2)
    midpoints = vcat(xm_b1, xm_b2)
    Cab = vcat(Cab_b1, Cab_b2)

    # cross section
    w = 1 # inch
    h = 0.063 # inch

    # material properties
    E = 1.06e7 # lb/in^2
    ν = 0.325
    ρ = 2.51e-4 # lb sec^2/in^4

    # shear and torsion correction factors
    ky = 1.2000001839588001
    kz = 14.625127919304001
    kt = 65.85255016982444

    A = h*w
    Iyy = w*h^3/12
    Izz = w^3*h/12
    J = Iyy + Izz

    # apply corrections
    Ay = A/ky
    Az = A/kz
    Jx = J/kt

    G = E/(2*(1+ν))

    compliance = fill(Diagonal([1/(E*A), 1/(G*Ay), 1/(G*Az), 1/(G*Jx), 1/(E*Iyy), 1/(E*Izz)]), nelem)

    mass = fill(Diagonal([ρ*A, ρ*A, ρ*A, ρ*J, ρ*Iyy, ρ*Izz]), nelem)

    # create assembly
    assembly = Assembly(points, start, stop, compliance=compliance, mass=mass, frames=Cab, lengths=lengths, midpoints=midpoints)

    # create dictionary of prescribed conditions
    prescribed_conditions = Dict(
        # root section is fixed
        1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0)
        )

    nonlinear_states = Vector{AssemblyState{Float64}}(undef, length(rpm))
    linear_states = Vector{AssemblyState{Float64}}(undef, length(rpm))
    for i = 1:length(rpm)
        # global frame rotation
        w0 = [0, 0, rpm[i]*(2*pi)/60]

        # perform nonlinear steady state analysis
        system, converged = steady_state_analysis(assembly,
            angular_velocity = w0,
            prescribed_conditions = prescribed_conditions)

        nonlinear_states[i] = AssemblyState(system, assembly, prescribed_conditions=prescribed_conditions)

        # perform linear steady state analysis
        system, converged = steady_state_analysis(assembly,
            angular_velocity = w0,
            prescribed_conditions = prescribed_conditions,
            linear = true)

        linear_states[i] = AssemblyState(system, assembly, prescribed_conditions=prescribed_conditions)
    end


    sweep = (0:2.5:45) * pi/180
    rpm = [0, 500, 750]
    nev = 30

    λ = Matrix{Vector{ComplexF64}}(undef, length(sweep), length(rpm))
    U = Matrix{Matrix{ComplexF64}}(undef, length(sweep), length(rpm))
    MV = Matrix{Matrix{ComplexF64}}(undef, length(sweep), length(rpm))
    state = Matrix{AssemblyState{Float64}}(undef, length(sweep), length(rpm))
    eigenstates = Matrix{Vector{AssemblyState{ComplexF64}}}(undef, length(sweep), length(rpm))
    for i = 1:length(sweep)
        local L_b1, r_b1, nelem_b1, lengths_b1 #hide
        local xp_b1, xm_b1, Cab_b1 #hide
        local cs, ss #hide
        local L_b2, r_b2, nelem_b2, frame_b2, lengths_b2 #hide
        local xp_b2, xm_b2, Cab_b2 #hide
        local nelem, points, start, stop #hide
        local lengths, midpoints, Cab, compliance, mass, assembly #hide

        # straight section of the beam
        L_b1 = 31.5 # inch
        r_b1 = [2.5, 0, 0]
        nelem_b1 = 20
        lengths_b1, xp_b1, xm_b1, Cab_b1 = discretize_beam(L_b1, r_b1, nelem_b1)

        # swept section of the beam
        L_b2 = 6 # inch
        r_b2 = [34, 0, 0]
        nelem_b2 = 20
        cs, ss = cos(sweep[i]), sin(sweep[i])
        frame_b2 = [cs ss 0; -ss cs 0; 0 0 1]
        lengths_b2, xp_b2, xm_b2, Cab_b2 = discretize_beam(L_b2, r_b2, nelem_b2, frame=frame_b2)

        # combine elements and points into one array
        nelem = nelem_b1 + nelem_b2
        points = vcat(xp_b1, xp_b2[2:end])
        start = 1:nelem_b1 + nelem_b2
        stop = 2:nelem_b1 + nelem_b2 + 1
        lengths = vcat(lengths_b1, lengths_b2)
        midpoints = vcat(xm_b1, xm_b2)
        Cab = vcat(Cab_b1, Cab_b2)

        compliance = fill(Diagonal([1/(E*A), 1/(G*Ay), 1/(G*Az), 1/(G*Jx), 1/(E*Iyy), 1/(E*Izz)]), nelem)

        mass = fill(Diagonal([ρ*A, ρ*A, ρ*A, ρ*J, ρ*Iyy, ρ*Izz]), nelem)

        # create assembly
        assembly = Assembly(points, start, stop, compliance=compliance, mass=mass, frames=Cab, lengths=lengths, midpoints=midpoints)

        # create system
        system = System(assembly, false)

        for j = 1:length(rpm)
            # global frame rotation
            w0 = [0, 0, rpm[j]*(2*pi)/60]

            # eigenvalues and (right) eigenvectors
            system, λ[i,j], V, converged = eigenvalue_analysis!(system, assembly,
                angular_velocity = w0,
                prescribed_conditions = prescribed_conditions,
                nev=nev)

            # corresponding left eigenvectors
            U[i,j] = left_eigenvectors(system, λ[i,j], V)

            # post-multiply mass matrix with right eigenvector matrix
            # (we use this later for correlating eigenvalues)
            MV[i,j] = system.M * V

            # process state and eigenstates
            state[i,j] = AssemblyState(system, assembly; prescribed_conditions=prescribed_conditions)
            eigenstates[i,j] = [AssemblyState(system, assembly, V[:,k];
                prescribed_conditions=prescribed_conditions) for k = 1:nev]
        end
    end


    # set previous left eigenvector matrix
    U_p = copy(U[1,1])

    for j = 1:length(rpm)
        for i = 1:length(sweep)
            # construct correlation matrix
            C = U_p*MV[i,j]

            # correlate eigenmodes
            perm, corruption = correlate_eigenmodes(C)

            # re-arrange eigenvalues and eigenvectors
            λ[i,j] = λ[i,j][perm]
            U[i,j] = U[i,j][perm,:]
            MV[i,j] = MV[i,j][:,perm]
            eigenstates[i,j] = eigenstates[i,j][perm]

            # update previous eigenvector matrix
            U_p .= U[i,j]
        end
        # update previous eigenvector matrix
        U_p .= U[1,j]
    end

    frequency = [[imag(λ[i,j][k])/(2*pi) for i = 1:length(sweep), j=1:length(rpm)] for k = 1:2:nev]

    indices = [1, 2, 4]
    experiment_rpm = [0, 500, 750]
    experiment_sweep = [0, 15, 30, 45]
    experiment_frequencies = [
        [1.4 1.8 1.7 1.6;
         10.2 10.1 10.2 10.2;
         14.8 14.4 14.9 14.7],
        [10.3 10.2 10.4 10.4;
         25.2 25.2 23.7 21.6;
         36.1 34.8 30.7 26.1],
        [27.7 27.2 26.6 24.8;
         47.0 44.4 39.3 35.1;
         62.9 55.9 48.6 44.8]
    ]

    for k = 1:length(experiment_frequencies)
        for j = 1:length(experiment_sweep)
            for i = 1:length(experiment_rpm)
                ii = argmin(abs.(rpm .- experiment_rpm[i]))
                jj = argmin(abs.(sweep*180/pi .- experiment_sweep[j]))
                kk = indices[k]
                @test isapprox(frequency[kk][jj,ii], experiment_frequencies[k][i,j], atol=1, rtol=0.1)
            end
        end
    end

    indices = [5, 7, 6]
    experiment_frequencies = [
        95.4 87.5 83.7 78.8;
        106.6 120.1 122.6 117.7;
        132.7 147.3 166.2 162.0
    ]

    for k = 1:size(experiment_frequencies, 1)
        for j = 1:length(experiment_sweep)
            ii = argmin(abs.(rpm .- 750))
            jj = argmin(abs.(sweep*180/pi .- experiment_sweep[j]))
            kk = indices[k]
            @test isapprox(frequency[kk][jj,ii], experiment_frequencies[k,j], rtol=0.1)
        end
    end
end

@testset "Nonlinear Dynamic Analysis of a Wind Turbine Blade" begin

    L = 60 # m

    # create points
    nelem = 10
    x = range(0, L, length=nelem+1)
    y = zero(x)
    z = zero(x)
    points = [[x[i],y[i],z[i]] for i = 1:length(x)]

    # index of endpoints of each beam element
    start = 1:nelem
    stop = 2:nelem+1

    # stiffness matrix for each beam element
    stiffness = fill(
        [2.389e9  1.524e6  6.734e6 -3.382e7 -2.627e7 -4.736e8
         1.524e6  4.334e8 -3.741e6 -2.935e5  1.527e7  3.835e5
         6.734e6 -3.741e6  2.743e7 -4.592e5 -6.869e5 -4.742e6
        -3.382e7 -2.935e5 -4.592e5  2.167e7 -6.279e5  1.430e6
        -2.627e7  1.527e7 -6.869e5 -6.279e5  1.970e7  1.209e7
        -4.736e8  3.835e5 -4.742e6  1.430e6  1.209e7  4.406e8],
        nelem)

    # mass matrix for each beam element
    mass = fill(
        [258.053      0.0        0.0      0.0      7.07839  -71.6871
           0.0      258.053      0.0     -7.07839  0.0        0.0
           0.0        0.0      258.053   71.6871   0.0        0.0
           0.0       -7.07839   71.6871  48.59     0.0        0.0
           7.07839    0.0        0.0      0.0      2.172      0.0
         -71.6871     0.0        0.0      0.0      0.0       46.418],
         nelem)

    # create assembly of interconnected nonlinear beams
    assembly = Assembly(points, start, stop; stiffness=stiffness, mass=mass)

    # simulation time
    tvec = 0:0.001:2.0

    # prescribed conditions
    prescribed_conditions = (t) -> begin
        Dict(
            # fixed left side
            1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
            # force on right side
            nelem+1 => PrescribedConditions(Fz=1e5*sin(20*t))
            )
    end

    system, history, converged = time_domain_analysis(assembly, tvec; prescribed_conditions=prescribed_conditions)

    @test converged
end

@testset "Nonlinear Static Analysis of a Joined-Wing" begin

    # Set endpoints of each beam
    p1 = [-7.1726, -12, -3.21539]
    p2 = [-5.37945, -9, -2.41154]
    p3 = [-3.5863, -6, -1.6077]
    p4 = [-1.79315, -3, -0.803848]
    p5 = [0, 0, 0]
    p6 = [7.1726, -12, 3.21539]

    # get transformation matrix for left beams

    # transformation from intermediate to global frame
    tmp1 = sqrt(p1[1]^2 + p1[2]^2)
    c1, s1 = -p1[1]/tmp1, -p1[2]/tmp1
    rot1 = [c1 -s1 0; s1 c1 0; 0 0 1]

    # transformation from local to intermediate frame
    tmp2 = sqrt(p1[1]^2 + p1[2]^2 + p1[3]^2)
    c2, s2 = tmp1/tmp2, -p1[3]/tmp2
    rot2 = [c2 0 -s2; 0 1 0; s2 0 c2]

    Cab_1 = rot1*rot2

    # get transformation matrix for right beam

    # transformation from intermediate frame to global frame
    tmp1 = sqrt(p6[1]^2 + p6[2]^2)
    c1, s1 = p6[1]/tmp1, p6[2]/tmp1
    rot1 = [c1 -s1 0; s1 c1 0; 0 0 1]

    # transformation from local beam frame to intermediate frame
    tmp2 = sqrt(p6[1]^2 + p6[2]^2 + p6[3]^2)
    c2, s2 = tmp1/tmp2, p6[3]/tmp2
    rot2 = [c2 0 -s2; 0 1 0; s2 0 c2]

    Cab_2 = rot1*rot2

    # beam 1
    L_b1 = norm(p2-p1)
    r_b1 = p1
    nelem_b1 = 5
    lengths_b1, xp_b1, xm_b1, Cab_b1 = discretize_beam(L_b1, r_b1, nelem_b1, frame=Cab_1)
    compliance_b1 = fill(Diagonal([1.05204e-9, 3.19659e-9, 2.13106e-8, 1.15475e-7, 1.52885e-7, 7.1672e-9]), nelem_b1)

    # beam 2
    L_b2 = norm(p3-p2)
    r_b2 = p2
    nelem_b2 = 5
    lengths_b2, xp_b2, xm_b2, Cab_b2 = discretize_beam(L_b2, r_b2, nelem_b2, frame=Cab_1)
    compliance_b2 = fill(Diagonal([1.24467e-9, 3.77682e-9, 2.51788e-8, 1.90461e-7, 2.55034e-7, 1.18646e-8]), nelem_b2)

    # beam 3
    L_b3 = norm(p4-p3)
    r_b3 = p3
    nelem_b3 = 5
    lengths_b3, xp_b3, xm_b3, Cab_b3 = discretize_beam(L_b3, r_b3, nelem_b3, frame=Cab_1)
    compliance_b3 = fill(Diagonal([1.60806e-9, 4.86724e-9, 3.24482e-8, 4.07637e-7, 5.57611e-7, 2.55684e-8]), nelem_b3)

    # beam 4
    L_b4 = norm(p5-p4)
    r_b4 = p4
    nelem_b4 = 5
    lengths_b4, xp_b4, xm_b4, Cab_b4 = discretize_beam(L_b4, r_b4, nelem_b4, frame=Cab_1)
    compliance_b4 = fill(Diagonal([2.56482e-9, 7.60456e-9, 5.67609e-8, 1.92171e-6, 2.8757e-6, 1.02718e-7]), nelem_b4)

    # beam 5
    L_b5 = norm(p6-p5)
    r_b5 = p5
    nelem_b5 = 20
    lengths_b5, xp_b5, xm_b5, Cab_b5 = discretize_beam(L_b5, r_b5, nelem_b5, frame=Cab_2)
    compliance_b5 = fill(Diagonal([2.77393e-9, 7.60456e-9, 1.52091e-7, 1.27757e-5, 2.7835e-5, 1.26026e-7]), nelem_b5)

    # combine elements and points into one array
    nelem = nelem_b1 + nelem_b2 + nelem_b3 + nelem_b4 + nelem_b5
    points = vcat(xp_b1, xp_b2[2:end], xp_b3[2:end], xp_b4[2:end], xp_b5[2:end])
    start = 1:nelem
    stop = 2:nelem + 1
    lengths = vcat(lengths_b1, lengths_b2, lengths_b3, lengths_b4, lengths_b5)
    midpoints = vcat(xm_b1, xm_b2, xm_b3, xm_b4, xm_b5)
    Cab = vcat(Cab_b1, Cab_b2, Cab_b3, Cab_b4, Cab_b5)
    compliance = vcat(compliance_b1, compliance_b2, compliance_b3, compliance_b4, compliance_b5)

    # create assembly
    assembly = Assembly(points, start, stop, compliance=compliance,
        frames=Cab, lengths=lengths, midpoints=midpoints)

    Fz = range(0, 70e3, length=141)

    # pre-allocate memory to reduce run-time
    system = System(assembly, true)

    linear_states = Vector{AssemblyState{Float64}}(undef, length(Fz))
    for i = 1:length(Fz)

        # create dictionary of prescribed conditions
        prescribed_conditions = Dict(
            # fixed endpoint on beam 1
            1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
            # force applied on point 4
            nelem_b1 + nelem_b2 + nelem_b3 + nelem_b4 + 1 => PrescribedConditions(Fz = Fz[i]),
            # fixed endpoint on last beam
            nelem+1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
        )

        _, converged = static_analysis!(system, assembly, prescribed_conditions=prescribed_conditions, linear=true)

        linear_states[i] = AssemblyState(system, assembly, prescribed_conditions=prescribed_conditions)

        @test converged
    end

    reset_state!(system)
    nonlinear_states = Vector{AssemblyState{Float64}}(undef, length(Fz))
    for i = 1:length(Fz)

        # create dictionary of prescribed conditions
        prescribed_conditions = Dict(
            # fixed endpoint on beam 1
            1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
            # force applied on point 4
            nelem_b1 + nelem_b2 + nelem_b3 + nelem_b4 + 1 => PrescribedConditions(Fz = Fz[i]),
            # fixed endpoint on last beam
            nelem+1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
        )

        _, converged = static_analysis!(system, assembly, prescribed_conditions=prescribed_conditions,
            reset_state = false)

        nonlinear_states[i] = AssemblyState(system, assembly;
            prescribed_conditions=prescribed_conditions)

        @test converged
    end

    reset_state!(system)
    nonlinear_follower_states = Vector{AssemblyState{Float64}}(undef, length(Fz))
    for i = 1:length(Fz)
        # create dictionary of prescribed conditions
        prescribed_conditions = Dict(
            # fixed endpoint on beam 1
            1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
            # force applied on point 4
            nelem_b1 + nelem_b2 + nelem_b3 + nelem_b4 + 1 => PrescribedConditions(Fz_follower = Fz[i]),
            # fixed endpoint on last beam
            nelem+1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
        )

        _, converged = static_analysis!(system, assembly, prescribed_conditions=prescribed_conditions,
            reset_state = false)

        nonlinear_follower_states[i] = AssemblyState(system, assembly;
            prescribed_conditions=prescribed_conditions)

        @test converged
    end
end

@testset "Nonlinear Dynamic Analysis of a Joined-Wing" begin

    # Set endpoints of each beam
    p1 = [0, 0, 0]
    p2 = [-7.1726, -12, -3.21539]
    p3 = [7.1726, -12,  3.21539]

    Cab_1 = [
    0.5         0.866025  0.0
    0.836516    -0.482963  0.258819
    0.224144     -0.12941   -0.965926
    ]

    Cab_2 = [
    0.5         0.866025  0.0
    -0.836516    0.482963 0.258819
    0.224144    -0.12941   0.965926
    ]

    # beam 1
    L_b1 = norm(p1-p2)
    r_b1 = p2
    nelem_b1 = 8
    lengths_b1, xp_b1, xm_b1, Cab_b1 = discretize_beam(L_b1, r_b1, nelem_b1, frame=Cab_1)

    # beam 2
    L_b2 = norm(p3-p1)
    r_b2 = p1
    nelem_b2 = 8
    lengths_b2, xp_b2, xm_b2, Cab_b2 = discretize_beam(L_b2, r_b2, nelem_b2, frame=Cab_2)

    # combine elements and points into one array
    nelem = nelem_b1 + nelem_b2
    points = vcat(xp_b1, xp_b2[2:end])
    start = 1:nelem
    stop = 2:nelem + 1
    lengths = vcat(lengths_b1, lengths_b2)
    midpoints = vcat(xm_b1, xm_b2)
    Cab = vcat(Cab_b1, Cab_b2)

    # assign all beams the same compliance and mass matrix
    compliance = fill(Diagonal([2.93944738387698e-10, 8.42991725049126e-10, 3.38313996669689e-08,
        4.69246721094557e-08, 6.79584100559513e-08, 1.37068861370898e-09]), nelem)
    mass = fill(Diagonal([4.86e-2, 4.86e-2, 4.86e-2,
        1.0632465e-2, 2.10195e-4, 1.042227e-2]), nelem)

    # create assembly
    assembly = Assembly(points, start, stop; compliance=compliance, mass=mass,
        frames=Cab, lengths=lengths, midpoints=midpoints)

    # time
    tvec = range(0, 0.04, length=1001)

    F_L = (t) -> begin
        if 0.0 <= t < 0.01
            1e6*t
        elseif 0.01 <= t < 0.02
            -1e6*(t-0.02)
        else
            zero(t)
        end
    end

    F_S = (t) -> begin
        if 0.0 <= t < 0.02
            5e3*(1-cos(pi*t/0.02))
        else
            1e4
        end
    end

    # assign boundary conditions and point load
    prescribed_conditions = (t) -> begin
        Dict(
        # fixed endpoint on beam 1
        1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
        # force applied on point 4
        nelem_b1 + 1 => PrescribedConditions(Fx=F_L(t), Fy=F_L(t), Fz=F_S(t)),
        # fixed endpoint on last beam
        nelem+1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
        )
    end

    system, history, converged = time_domain_analysis(assembly, tvec;
        prescribed_conditions=prescribed_conditions)

    @test converged
end

@testset "DifferentialEquations" begin

    L = 60 # m

    # create points
    nelem = 10
    x = range(0, L, length=nelem+1)
    y = zero(x)
    z = zero(x)
    points = [[x[i],y[i],z[i]] for i = 1:length(x)]

    # index of endpoints of each beam element
    start = 1:nelem
    stop = 2:nelem+1

    # stiffness matrix for each beam element
    stiffness = fill(
        [2.389e9  1.524e6  6.734e6 -3.382e7 -2.627e7 -4.736e8
         1.524e6  4.334e8 -3.741e6 -2.935e5  1.527e7  3.835e5
         6.734e6 -3.741e6  2.743e7 -4.592e5 -6.869e5 -4.742e6
        -3.382e7 -2.935e5 -4.592e5  2.167e7 -6.279e5  1.430e6
        -2.627e7  1.527e7 -6.869e5 -6.279e5  1.970e7  1.209e7
        -4.736e8  3.835e5 -4.742e6  1.430e6  1.209e7  4.406e8],
        nelem)

    # mass matrix for each beam element
    mass = fill(
        [258.053      0.0        0.0      0.0      7.07839  -71.6871
           0.0      258.053      0.0     -7.07839  0.0        0.0
           0.0        0.0      258.053   71.6871   0.0        0.0
           0.0       -7.07839   71.6871  48.59     0.0        0.0
           7.07839    0.0        0.0      0.0      2.172      0.0
         -71.6871     0.0        0.0      0.0      0.0       46.418],
         nelem)

    # create assembly of interconnected nonlinear beams
    assembly = Assembly(points, start, stop; stiffness=stiffness, mass=mass)

    # prescribed conditions
    prescribed_conditions = (t) -> begin
        Dict(
            # fixed left side
            1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
            # force on right side
            nelem+1 => PrescribedConditions(Fz=1e5*sin(20*t))
            )
    end

    # define simulation time
    tspan = (0.0, 2.0)

    # run initial condition analysis to get consistent set of initial conditions
    system, converged = initial_condition_analysis(assembly, tspan[1]; prescribed_conditions)

    # construct ODEProblem
    prob = ODEProblem(system, assembly, tspan; prescribed_conditions)

    # solve ODEProblem
    sol = solve(prob, Rodas4())

    # test that solution worked
    @test sol.t[end] == 2.0

    # construct DAEProblem
    prob = DAEProblem(system, assembly, tspan; prescribed_conditions)

    # solve DAEProblem
    sol = solve(prob, DABDF2())

    # test that solution worked
    @test sol.t[end] == 2.0
end

@testset "ForwardDiff" begin

    # Linear Analysis of a Beam Under a Linear Distributed Load

    function linear_analysis_test_with_AD(length) # this should affect just about everything

        nelem = 16

        # create points
        L = length[1]
        x = collect(range(0, L, length=nelem+1))
        y = zero(x)
        z = zero(x)

        points = [[x[i],y[i],z[i]] for i = 1:size(x,1)]

        # index of endpoints for each beam element
        start = 1:nelem
        stop = 2:nelem+1

        # create compliance matrix for each beam element
        EI = 1e7
        compliance = fill(Diagonal([0, 0, 0, 0, 1/EI, 0]), nelem)

        # create assembly
        assembly = Assembly(points, start, stop, compliance=compliance)

        # set prescribed conditions
        prescribed_conditions = Dict(
            # simply supported left endpoint
            1 => PrescribedConditions(uz=0),
            # clamped right endpoint
            nelem+1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0)
        )

        # create distributed load
        qmax = 1000
        distributed_loads = Dict()
        for i = 1:nelem
            distributed_loads[i] = DistributedLoads(assembly, i; s1=x[i],
                s2=x[i+1], fz = (s) -> qmax*s)
        end

        # solve system
        system, converged = static_analysis(assembly, prescribed_conditions=prescribed_conditions,
            distributed_loads=distributed_loads, linear=true)

        return system.x
    end

    # run FrowardDiff - no specific test, just make sure it runs fine
    J = ForwardDiff.jacobian(linear_analysis_test_with_AD, [1.0]) #length=1
end

@testset "Zero Mass Matrix" begin
    sweep = 45 * pi/180
    rpm = 750

    # straight section of the beam
    L_b1 = 31.5 # inch
    r_b1 = [2.5, 0, 0]
    nelem_b1 = 13
    lengths_b1, xp_b1, xm_b1, Cab_b1 = discretize_beam(L_b1, r_b1, nelem_b1)

    # swept section of the beam
    L_b2 = 6 # inch
    r_b2 = [34, 0, 0]
    nelem_b2 = 3
    cs, ss = cos(sweep), sin(sweep)
    frame_b2 = [cs ss 0; -ss cs 0; 0 0 1]
    lengths_b2, xp_b2, xm_b2, Cab_b2 = discretize_beam(L_b2, r_b2, nelem_b2, frame=frame_b2)

    # combine elements and points into one array
    nelem = nelem_b1 + nelem_b2
    points = vcat(xp_b1, xp_b2[2:end])
    start = 1:nelem_b1 + nelem_b2
    stop = 2:nelem_b1 + nelem_b2 + 1
    lengths = vcat(lengths_b1, lengths_b2)
    midpoints = vcat(xm_b1, xm_b2)
    Cab = vcat(Cab_b1, Cab_b2)

    # cross section
    w = 1 # inch
    h = 0.063 # inch

    # material properties
    E = 1.06e7 # lb/in^2
    ν = 0.325
    ρ = 2.51e-4 # lb sec^2/in^4

    # shear and torsion correction factors
    ky = 1.2000001839588001
    kz = 14.625127919304001
    kt = 65.85255016982444

    A = h*w
    Iyy = w*h^3/12
    Izz = w^3*h/12
    J = Iyy + Izz

    # apply corrections
    Ay = A/ky
    Az = A/kz
    Jx = J/kt

    G = E/(2*(1+ν))

    compliance = fill(Diagonal([1/(E*A), 1/(G*Ay), 1/(G*Az), 1/(G*Jx), 1/(E*Iyy), 1/(E*Izz)]), nelem)

    mass = fill(Diagonal(zeros(6)), nelem)

    # create assembly
    assembly = Assembly(points, start, stop, compliance=compliance, mass=mass, frames=Cab, lengths=lengths, midpoints=midpoints)

    # create dictionary of prescribed conditions
    prescribed_conditions = Dict(
        # root section is fixed
        1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0)
        )

    # set angular velocity vector
    w0 = [0, 0, rpm*(2*pi)/60]

    # perform nonlinear steady state analysis
    system, converged = steady_state_analysis(assembly,
        angular_velocity = w0,
        prescribed_conditions = prescribed_conditions)

    # test convergence
    @test converged
end

@testset "Zero Length Element" begin
    sweep = 45 * pi/180
    rpm = 750

    # straight section of the beam
    L_b1 = 31.5 # inch
    r_b1 = [2.5, 0, 0]
    nelem_b1 = 13
    lengths_b1, xp_b1, xm_b1, Cab_b1 = discretize_beam(L_b1, r_b1, nelem_b1)

    # zero length element between straight and swept sections
    L_b12 = 0
    r_b12 = [34, 0, 0]
    nelem_b12 = 1
    lengths_b12, xp_b12, xm_b12, Cab_b12 = discretize_beam(L_b12, r_b12, nelem_b12)

    # swept section of the beam
    L_b2 = 6 # inch
    r_b2 = [34, 0, 0]
    nelem_b2 = 3
    cs, ss = cos(sweep), sin(sweep)
    frame_b2 = [cs ss 0; -ss cs 0; 0 0 1]
    lengths_b2, xp_b2, xm_b2, Cab_b2 = discretize_beam(L_b2, r_b2, nelem_b2, frame=frame_b2)

    # combine elements and points into one array
    nelem = nelem_b1 + nelem_b12 + nelem_b2
    points = vcat(xp_b1, xp_b2[2:end]) # don't duplicate points
    lengths = vcat(lengths_b1, lengths_b12, lengths_b2)
    midpoints = vcat(xm_b1, xm_b12, xm_b2)
    Cab = vcat(Cab_b1, Cab_b12, Cab_b2)

    # specify connectivity
    start = vcat(1:nelem_b1+1, nelem_b1+1:nelem_b1+nelem_b2)
    stop = vcat(2:nelem_b1+1, nelem_b1+1:nelem_b1+nelem_b2+1)

    # cross section
    w = 1 # inch
    h = 0.063 # inch

    # material properties
    E = 1.06e7 # lb/in^2
    ν = 0.325
    ρ = 2.51e-4 # lb sec^2/in^4

    # shear and torsion correction factors
    ky = 1.2000001839588001
    kz = 14.625127919304001
    kt = 65.85255016982444

    A = h*w
    Iyy = w*h^3/12
    Izz = w^3*h/12
    J = Iyy + Izz

    # apply corrections
    Ay = A/ky
    Az = A/kz
    Jx = J/kt

    G = E/(2*(1+ν))

    compliance = fill(Diagonal([1/(E*A), 1/(G*Ay), 1/(G*Az), 1/(G*Jx), 1/(E*Iyy), 1/(E*Izz)]), nelem)

    mass = fill(Diagonal([ρ*A, ρ*A, ρ*A, ρ*J, ρ*Iyy, ρ*Izz]), nelem)

    # create assembly
    assembly = Assembly(points, start, stop, compliance=compliance, mass=mass, frames=Cab, lengths=lengths, midpoints=midpoints)

    # create dictionary of prescribed conditions
    prescribed_conditions = Dict(
        # root section is fixed
        1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0)
        )

    # set angular velocity vector
    w0 = [0, 0, rpm*(2*pi)/60]

    # perform nonlinear steady state analysis
    system, converged = steady_state_analysis(assembly,
        angular_velocity = w0,
        prescribed_conditions = prescribed_conditions)

    # test convergence
    @test converged
end

@testset "Element Gravitational Loads" begin

    # use arbitrary length
    ΔL = rand()

    # use random rotation matrix  
    CtCab = GXBeam.get_C(rand(3))

    # create random mass matrix
    μ = rand()
    xm2 = rand()
    xm3 = rand()
    i22 = rand()
    i33 = rand()
    i23 = rand()

    mass = [
        μ 0 0 0 μ*xm3 -μ*xm2; 
        0 μ 0 -μ*xm3 0 0; 
        0 0 μ μ*xm2 0 0; 
        0 -μ*xm3 μ*xm2 i22+i33 0 0; 
        μ*xm3 0 0 0 i22 0; 
        -μ*xm2 0 0 0 0 i33
        ]

    # use random gravity vector
    gvec = rand(3)

    # calculate integrated force and moment per unit length
    f = μ*gvec
    m = cross(CtCab*[0, xm2, xm3], f)
    f1 = f2 = ΔL*f/2
    m1 = m2 = ΔL*m/2

    # test against gravitational load function results
    mass11 = mass[1:3, 1:3]
    mass12 = mass[1:3, 4:6]
    f1t, f2t, m1t, m2t = GXBeam.element_gravitational_loads(ΔL, CtCab, mass11, mass12, gvec)

    @test isapprox(f1, f1t)
    @test isapprox(f2, f2t)
    @test isapprox(m1, m1t)
    @test isapprox(m2, m2t)

end