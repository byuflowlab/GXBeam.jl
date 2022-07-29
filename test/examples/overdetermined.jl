using GXBeam, LinearAlgebra, Test

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
        @test isapprox(state.elements[i].Mi[2], -analytical_M(xi), atol=1)
    end

    # test point properties
    for i = 1:length(assembly.points)
        xi = assembly.points[i][1]
        @test isapprox(state.points[i].u[3], analytical_deflection(xi), atol=1e-8)
        @test isapprox(state.points[i].theta[2], -4*analytical_slope(xi)/4, atol=1e-8)
    end

    # now check the state variables for a constant mass matrix system
    system, converged = steady_state_analysis(assembly, 
        prescribed_conditions=prescribed_conditions,
        distributed_loads=distributed_loads, 
        constant_mass_matrix=true,
        linear=true)

    state = AssemblyState(system, assembly, prescribed_conditions=prescribed_conditions)

    # test element properties
    for i = 1:length(assembly.elements)
        xi = assembly.elements[i].x[1]
        @test isapprox(state.elements[i].u[3], analytical_deflection(xi), atol=1e-8)
        @test isapprox(state.elements[i].theta[2], -4*analytical_slope(xi)/4, atol=1e-7)
        @test isapprox(state.elements[i].Mi[2], -analytical_M(xi), atol=1)
    end

    # test point properties
    for i = 1:length(assembly.points)
        xi = assembly.points[i][1]
        @test isapprox(state.points[i].u[3], analytical_deflection(xi), atol=1e-8)
        @test isapprox(state.points[i].theta[2], -4*analytical_slope(xi)/4, atol=1e-8)
    end
end