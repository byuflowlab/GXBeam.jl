using GXBeam, LinearAlgebra, Test

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

    system, state, converged = static_analysis(assembly, prescribed_conditions=prescribed_conditions,
        distributed_loads=distributed_loads, linear=true)

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
        @test isapprox(state.elements[i].Mi[2], -analytical_M(xi), atol=2)
    end

    # test point properties
    for i = 1:length(assembly.points)
        xi = assembly.points[i][1]
        @test isapprox(state.points[i].u[3], analytical_deflection(xi), atol=1e-8)
        @test isapprox(state.points[i].theta[2], -4*analytical_slope(xi)/4, atol=1e-7)
    end

    # # now check the state variables for a constant mass matrix system
    # system, state, converged = steady_state_analysis(assembly,
    #     prescribed_conditions = prescribed_conditions,
    #     distributed_loads = distributed_loads,
    #     constant_mass_matrix = true,
    #     linear = true)

    # for i = 1:length(assembly.elements)
    #     xi = assembly.elements[i].x[1]
    #     @test isapprox(state.elements[i].u[3], analytical_deflection(xi), atol=1e-9)
    #     @test isapprox(state.elements[i].theta[2], -4*analytical_slope(xi)/4, atol=1e-9)
    #     @test isapprox(state.elements[i].Mi[2], -analytical_M(xi), atol=2)
    # end

    # for i = 1:length(assembly.points)
    #     xi = assembly.points[i][1]
    #     @test isapprox(state.points[i].u[3], analytical_deflection(xi), atol=1e-8)
    #     @test isapprox(state.points[i].theta[2], -4*analytical_slope(xi)/4, atol=1e-7)
    # end

end
