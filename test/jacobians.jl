using GXBeam, LinearAlgebra, Random, StaticArrays, ForwardDiff, Test

@testset "Rotation Parameter Jacobians" begin
    
    RNG = MersenneTwister(1234)

    θ = 1e3*rand(RNG, 3)
    Δθ = 1e3*rand(RNG, 3)

    # get_C_θ
    C_θ1, C_θ2, C_θ3 = GXBeam.get_C_θ(θ)
    @test isapprox(C_θ1, ForwardDiff.derivative(θ1 -> GXBeam.get_C([θ1, θ[2], θ[3]]), θ[1]))
    @test isapprox(C_θ2, ForwardDiff.derivative(θ2 -> GXBeam.get_C([θ[1], θ2, θ[3]]), θ[2]))
    @test isapprox(C_θ3, ForwardDiff.derivative(θ3 -> GXBeam.get_C([θ[1], θ[2], θ3]), θ[3]))

    # get_Q_θ
    Q_θ1, Q_θ2, Q_θ3 = GXBeam.get_Q_θ(θ)
    @test isapprox(Q_θ1, ForwardDiff.derivative(θ1 -> GXBeam.get_Q([θ1, θ[2], θ[3]]), θ[1]))
    @test isapprox(Q_θ2, ForwardDiff.derivative(θ2 -> GXBeam.get_Q([θ[1], θ2, θ[3]]), θ[2]))
    @test isapprox(Q_θ3, ForwardDiff.derivative(θ3 -> GXBeam.get_Q([θ[1], θ[2], θ3]), θ[3]))

    # get_Qinv_θ
    Qinv_θ1, Qinv_θ2, Qinv_θ3 = GXBeam.get_Qinv_θ(θ)
    @test isapprox(Qinv_θ1, ForwardDiff.derivative(θ1 -> GXBeam.get_Qinv([θ1, θ[2], θ[3]]), θ[1]))
    @test isapprox(Qinv_θ2, ForwardDiff.derivative(θ2 -> GXBeam.get_Qinv([θ[1], θ2, θ[3]]), θ[2]))
    @test isapprox(Qinv_θ3, ForwardDiff.derivative(θ3 -> GXBeam.get_Qinv([θ[1], θ[2], θ3]), θ[3]))

    # get_ΔQ
    ΔQ = GXBeam.get_ΔQ(θ, Δθ)
    @test isapprox(ΔQ, GXBeam.mul3(Q_θ1, Q_θ2, Q_θ3, Δθ))

    # get_ΔQ_θ
    ΔQ_θ1, ΔQ_θ2, ΔQ_θ3 = GXBeam.get_ΔQ_θ(θ, Δθ)
    @test isapprox(ΔQ_θ1, ForwardDiff.derivative(θ1 -> GXBeam.get_ΔQ([θ1, θ[2], θ[3]], Δθ), θ[1]))
    @test isapprox(ΔQ_θ2, ForwardDiff.derivative(θ2 -> GXBeam.get_ΔQ([θ[1], θ2, θ[3]], Δθ), θ[2]))
    @test isapprox(ΔQ_θ3, ForwardDiff.derivative(θ3 -> GXBeam.get_ΔQ([θ[1], θ[2], θ3], Δθ), θ[3]))

end

@testset "Prescribed Force Jacobians" begin

    RNG = MersenneTwister(1234)

    x = rand(RNG, 6)

    icol = 1

    force_scaling = 1.0

    prescribed_conditions = PrescribedConditions(;
        Fx = rand(RNG),
        Fy = rand(RNG),
        Fz = rand(RNG),
        Mx = rand(RNG),
        My = rand(RNG),
        Mz = rand(RNG),
        Fx_follower = rand(RNG),
        Fy_follower = rand(RNG),
        Fz_follower = rand(RNG),
        Mx_follower = rand(RNG),
        My_follower = rand(RNG),
        Mz_follower = rand(RNG))

    # point_displacement_jacobians
    u, θ = GXBeam.point_displacement(x, icol, prescribed_conditions)
    
    u_u, θ_θ = GXBeam.point_displacement_jacobians(prescribed_conditions)

    f = x -> vcat(GXBeam.point_displacement(x, icol, prescribed_conditions)...)

    dx = ForwardDiff.jacobian(f, x)

    @test isapprox(u_u, dx[1:3,1:3])
    @test isapprox(θ_θ, dx[4:6,4:6])

    # point_load_jacobians
    F, M = GXBeam.point_loads(x, icol, force_scaling, prescribed_conditions)

    F_θ, F_F, M_θ, M_M = GXBeam.point_load_jacobians(x, icol, force_scaling, prescribed_conditions)

    f = x -> vcat(GXBeam.point_loads(x, icol, force_scaling, prescribed_conditions)...)

    dx = ForwardDiff.jacobian(f, x)

    @test iszero(F_F)
    @test iszero(M_M)

    @test isapprox(F_θ, dx[1:3,4:6])
    @test isapprox(M_θ, dx[4:6,4:6])

end

@testset "Jacobian and Mass Matrix Calculations" begin

    RNG = MersenneTwister(1234)

    L = 60 # m

    # create points
    nelem = 2
    x = range(0, L, length=nelem+1)
    y = zero(x)
    z = zero(x)
    points = [[x[i],y[i],z[i]] for i in eachindex(x)]

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

    damping = fill(0.1*rand(RNG, 6), nelem)

    # create assembly of interconnected nonlinear beams
    assembly = Assembly(points, start, stop; stiffness=stiffness, mass=mass, damping=damping)

    # prescribed conditions
    pcond = Dict(
        # fixed left side
        1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, My=0, theta_z=0),
        )

    # distributed loads
    dload = Dict{Int,DistributedLoads{Float64}}()

    # point masses
    pmass = Dict(
        # point mass at the end of the beam
        nelem+1 => PointMass(Symmetric(1e3*rand(RNG,6,6)))
    )

    # gravity vector
    gvec = SVector{3}(1e3*rand(RNG, 3))

    # --- Static Analysis --- #

    system = StaticSystem(assembly)
    force_scaling = system.force_scaling
    indices = system.indices
    x = 1e2 .* rand(RNG, length(system.x))
    J = similar(x, length(x), length(x))
    two_dimensional = false

    f = (x) -> GXBeam.static_system_residual!(similar(x), x, indices, two_dimensional, force_scaling, 
        assembly, pcond, dload, pmass, gvec)

    GXBeam.static_system_jacobian!(J, x, indices, two_dimensional, force_scaling,
        assembly, pcond, dload, pmass, gvec)

    J_fd = ForwardDiff.jacobian(f, x)

    @test all(isapprox.(J, J_fd, atol=1e-10))

    # --- Steady State Analysis --- #

    system = DynamicSystem(assembly)
    GXBeam.update_body_acceleration_indices!(system, pcond)
    force_scaling = system.force_scaling
    indices = system.indices
    x = 1e2 .* rand(RNG, length(system.x))
    J = similar(x, length(x), length(x))

    structural_damping = true

    linear_velocity = SVector{3}(1e2*rand(RNG, 3))
    angular_velocity = SVector{3}(1e2*rand(RNG, 3))
    linear_acceleration = SVector{3}(1e2*rand(RNG, 3))
    angular_acceleration = SVector{3}(1e2*rand(RNG, 3))

    f = (x) -> GXBeam.steady_system_residual!(similar(x), x, indices, two_dimensional, 
        force_scaling, structural_damping, assembly, pcond, dload, pmass, gvec, 
        linear_velocity, angular_velocity, linear_acceleration, angular_acceleration)

    GXBeam.steady_system_jacobian!(J, x, indices, two_dimensional, force_scaling, 
        structural_damping, assembly, pcond, dload, pmass, gvec, 
        linear_velocity, angular_velocity, linear_acceleration, angular_acceleration)

    J_fd = ForwardDiff.jacobian(f, x)

    @test all(isapprox.(J, J_fd, atol=1e-10))

    # --- Initial Condition Analysis --- #

    u0 = [SVector{3}(rand(RNG, 3)) for ielem = 1:length(assembly.points)]
    theta0 = [SVector{3}(rand(RNG, 3)) for ielem = 1:length(assembly.points)]
    V0 = [SVector{3}(rand(RNG, 3)) for ielem = 1:length(assembly.points)]
    Omega0 = [SVector{3}(rand(RNG, 3)) for ielem = 1:length(assembly.points)]
    Vdot0 = [SVector{3}(rand(RNG, 3)) for ielem = 1:length(assembly.points)]
    Omegadot0 = [SVector{3}(rand(RNG, 3)) for ielem = 1:length(assembly.points)]

    x = rand(RNG, length(system.x))
    J = similar(x, length(x), length(x))
    M = similar(x, length(x), length(x))
    rate_vars1 = rand(RNG, Bool, length(x))
    rate_vars2 = rand(RNG, Bool, length(x))

    f = (x) -> GXBeam.initial_system_residual!(similar(x), x, indices, rate_vars1, 
        rate_vars2, two_dimensional, force_scaling, structural_damping, assembly, pcond, dload, pmass, gvec, 
        linear_velocity, angular_velocity, linear_acceleration, angular_acceleration,
        u0, theta0, V0, Omega0, Vdot0, Omegadot0)

    GXBeam.initial_system_jacobian!(J, x, indices, rate_vars1, rate_vars2, two_dimensional,
        force_scaling, structural_damping, assembly, pcond, dload, pmass, gvec, 
        linear_velocity, angular_velocity, linear_acceleration, angular_acceleration,
        u0, theta0, V0, Omega0, Vdot0, Omegadot0)

    J_fd = ForwardDiff.jacobian(f, x)

    @test all(isapprox.(J, J_fd, atol=1e-10))

    # --- Newmark Scheme Time-Marching Analysis --- #

    udot = [SVector{3}(rand(RNG, 3)) for ipoint = 1:length(assembly.points)]
    θdot = [SVector{3}(rand(RNG, 3)) for ipoint = 1:length(assembly.points)]
    Vdot = [SVector{3}(rand(RNG, 3)) for ipoint = 1:length(assembly.points)]
    Ωdot = [SVector{3}(rand(RNG, 3)) for ipoint = 1:length(assembly.points)]
    dt = rand(RNG)

    x = rand(RNG, length(system.x))
    J = similar(x, length(x), length(x))

    f = (x) -> GXBeam.newmark_system_residual!(similar(x), x, indices, two_dimensional, 
        force_scaling, structural_damping, assembly, pcond, dload, pmass, gvec, 
        linear_velocity, angular_velocity, udot, θdot, Vdot, Ωdot, dt)

    GXBeam.newmark_system_jacobian!(J, x, indices, two_dimensional, force_scaling, 
        structural_damping, assembly, pcond, dload, pmass, gvec, 
        linear_velocity, angular_velocity, udot, θdot, Vdot, Ωdot, dt)

    J_fd = ForwardDiff.jacobian(f, x)

    @test all(isapprox.(J, J_fd, atol=1e-10))

    # --- General Dynamic Analysis (3D) --- #

    dx = rand(RNG, length(system.x))
    x = rand(RNG, length(system.x))
    J = similar(x, length(x), length(x))
    M = similar(x, length(x), length(x))

    fx = (x) -> GXBeam.dynamic_system_residual!(similar(x), dx, x, indices, two_dimensional, 
        force_scaling, structural_damping, assembly, pcond, dload, pmass, gvec,
        linear_velocity, angular_velocity)

    fdx = (dx) -> GXBeam.dynamic_system_residual!(similar(dx), dx, x, indices, two_dimensional, 
        force_scaling, structural_damping, assembly, pcond, dload, pmass, gvec, 
        linear_velocity, angular_velocity)

    GXBeam.dynamic_system_jacobian!(J, dx, x, indices, two_dimensional, force_scaling, 
        structural_damping, assembly, pcond, dload, pmass, gvec, 
        linear_velocity, angular_velocity)

    GXBeam.system_mass_matrix!(M, x, indices, two_dimensional, force_scaling, assembly, pcond, pmass)

    J_fd = ForwardDiff.jacobian(fx, x)

    M_fd = ForwardDiff.jacobian(fdx, dx)

    @test all(isapprox.(J, J_fd, atol=1e-10))

    @test all(isapprox.(M, M_fd, atol=1e-10))

    # --- General Dynamic Analysis (2D) --- #

    two_dimensional = true

    GXBeam.dynamic_system_jacobian!(J, dx, x, indices, two_dimensional, force_scaling, 
        structural_damping, assembly, pcond, dload, pmass, gvec, 
        linear_velocity, angular_velocity)

    GXBeam.system_mass_matrix!(M, x, indices, two_dimensional, force_scaling, assembly, pcond, pmass)

    J_fd = ForwardDiff.jacobian(fx, x)

    M_fd = ForwardDiff.jacobian(fdx, dx)

    @test all(isapprox.(J, J_fd, atol=1e-10))

    @test all(isapprox.(M, M_fd, atol=1e-10))

    # --- Constant Mass Matrix (3D) --- #

    two_dimensional = false

    system = ExpandedSystem(assembly)
    force_scaling = system.force_scaling
    indices = system.indices
    dx = rand(RNG, length(system.x))
    x = rand(RNG, length(system.x))
    J = similar(x, length(x), length(x))
    M = similar(x, length(x), length(x))

    f = (x) -> GXBeam.expanded_steady_system_residual!(similar(x), x, indices, two_dimensional, force_scaling, 
        structural_damping, assembly, pcond, dload, pmass, gvec, 
        linear_velocity, angular_velocity, linear_acceleration, angular_acceleration)

    GXBeam.expanded_steady_system_jacobian!(J, x, indices, two_dimensional, force_scaling, 
        structural_damping, assembly, pcond, dload, pmass, gvec,         
        linear_velocity, angular_velocity, linear_acceleration, angular_acceleration)

    J_fd = ForwardDiff.jacobian(f, x)

    @test all(isapprox.(J, J_fd, atol=1e-10))

    fx = (x) -> GXBeam.expanded_dynamic_system_residual!(similar(x), dx, x, indices, two_dimensional, force_scaling, 
        structural_damping, assembly, pcond, dload, pmass, gvec, linear_velocity, angular_velocity)

    fdx = (dx) -> GXBeam.expanded_dynamic_system_residual!(similar(dx), dx, x, indices, two_dimensional, force_scaling, 
        structural_damping, assembly, pcond, dload, pmass, gvec, linear_velocity, angular_velocity)

    GXBeam.expanded_dynamic_system_jacobian!(J, dx, x, indices, two_dimensional, force_scaling, 
        structural_damping, assembly, pcond, dload, pmass, gvec, linear_velocity, angular_velocity)

    GXBeam.expanded_system_mass_matrix!(M, indices, two_dimensional, force_scaling, assembly, 
        pcond, pmass)

    J_fd = ForwardDiff.jacobian(fx, x)

    M_fd = ForwardDiff.jacobian(fdx, x)

    @test all(isapprox.(J, J_fd, atol=1e-10))

    @test all(isapprox.(M, M_fd, atol=1e-10))

    # -- Constant Mass Matrix (2D) --- #
    
    two_dimensional = true

    GXBeam.expanded_dynamic_system_jacobian!(J, dx, x, indices, two_dimensional, force_scaling, 
        structural_damping, assembly, pcond, dload, pmass, gvec, linear_velocity, angular_velocity)

    GXBeam.expanded_system_mass_matrix!(M, indices, two_dimensional, force_scaling, assembly, 
        pcond, pmass)

    J_fd = ForwardDiff.jacobian(fx, x)

    M_fd = ForwardDiff.jacobian(fdx, x)

    @test all(isapprox.(J, J_fd, atol=1e-10))

    @test all(isapprox.(M, M_fd, atol=1e-10))

end
