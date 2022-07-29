using GXBeam, LinearAlgebra, Elliptic, Test

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
    system = StaticSystem(assembly)

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

    # perform the same analysis for a constant mass matrix system
    system = ExpandedSystem(assembly)

    states = Vector{AssemblyState{Float64}}(undef, length(P))
    for i = 1:length(P)

        prescribed_conditions = Dict(
            1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
            nelem+1 => PrescribedConditions(Fz = P[i])
        )

        steady_state_analysis!(system, assembly, 
            prescribed_conditions = prescribed_conditions,
            constant_mass_matrix = true)

        states[i] = AssemblyState(system, assembly, prescribed_conditions=prescribed_conditions)
    end

    # test tip displacements
    for i = 1:length(P)
        i_a = argmin(abs.(λ[i] .- λ_a))
        @test isapprox(states[i].points[end].u[1]/L, ξ_a[i_a], atol=1e-3)
        @test isapprox(states[i].points[end].u[3]/L, η_a[i_a], atol=1e-3)
        @test isapprox(states[i].points[end].theta[2], -4*tan(θ_a[i_a]/4), atol=1e-2)
    end

end