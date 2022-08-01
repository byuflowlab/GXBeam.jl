using GXBeam, DifferentialEquations, Elliptic, LinearAlgebra, ForwardDiff, Random, Test

@testset "Input and Output" begin

    RNG = MersenneTwister(1234)

    # --- Assembly --- #
    
    # beam length
    L = 1

    # number of elements
    nelem = 1

    # points
    x = range(0, L, length=nelem+1)
    y = zero(x)
    z = zero(x)
    points = [[x[i],y[i],z[i]] for i = 1:length(x)]

    # connectivity
    start = 1:nelem
    stop = 2:nelem+1

    # stiffness matrix
    stiffness = fill(
        [2.389e9  1.524e6  6.734e6 -3.382e7 -2.627e7 -4.736e8
         1.524e6  4.334e8 -3.741e6 -2.935e5  1.527e7  3.835e5
         6.734e6 -3.741e6  2.743e7 -4.592e5 -6.869e5 -4.742e6
        -3.382e7 -2.935e5 -4.592e5  2.167e7 -6.279e5  1.430e6
        -2.627e7  1.527e7 -6.869e5 -6.279e5  1.970e7  1.209e7
        -4.736e8  3.835e5 -4.742e6  1.430e6  1.209e7  4.406e8],
        nelem)

    # mass matrix
    mass = fill(
        [258.053      0.0        0.0      0.0      7.07839  -71.6871
           0.0      258.053      0.0     -7.07839  0.0        0.0
           0.0        0.0      258.053   71.6871   0.0        0.0
           0.0       -7.07839   71.6871  48.59     0.0        0.0
           7.07839    0.0        0.0      0.0      2.172      0.0
         -71.6871     0.0        0.0      0.0      0.0       46.418],
         nelem)

    # damping
    damping = fill(rand(RNG, 6), nelem)

    # assembly
    assembly = Assembly(points, start, stop; stiffness=stiffness, mass=mass, damping=damping)

    # --- Operating Conditions --- #

    # prescribed conditions
    prescribed_conditions = Dict(
        # fixed left side
        1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
        # loaded right side
        nelem+1 => PrescribedConditions(
            Fx=rand(RNG), Fy=rand(RNG), Fz=rand(RNG), 
            Mx=rand(RNG), My=rand(RNG), Mz=rand(RNG),
            Fx_follower=rand(RNG), Fy_follower=rand(RNG), Fz_follower=rand(RNG), 
            Mx_follower=rand(RNG), My_follower=rand(RNG), Mz_follower=rand(RNG),),
        )

    # distributed loads
    fx=rand(RNG); fy=rand(RNG); fz=rand(RNG)
    mx=rand(RNG); my=rand(RNG); mz=rand(RNG)
    fx_follower=rand(RNG); fy_follower=rand(RNG); fz_follower=rand(RNG)
    mx_follower=rand(RNG); my_follower=rand(RNG); mz_follower=rand(RNG)

    distributed_loads = Dict(
        ielem => DistributedLoads(assembly, ielem;
            fx=(s)->fx, fy=(s)->fy, fz=(s)->fz, 
            mx=(s)->mx, my=(s)->my, mz=(s)->mz,
            fx_follower=(s)->fx_follower, fy_follower=(s)->fy_follower, fz_follower=(s)->fz_follower, 
            mx_follower=(s)->mx_follower, my_follower=(s)->my_follower, mz_follower=(s)->mz_follower,
            ) for ielem = 1:nelem
    )

    # point masses
    point_masses = Dict(
        # point mass at the end of the beam
        nelem+1 => PointMass(Symmetric(rand(RNG, 6,6)))
    )
   
    # body frame linear velocity
    linear_velocity = rand(RNG, 3)
    
    # body frame angular velocity
    angular_velocity = rand(RNG, 3)
    
    # body frame linear acceleration
    linear_acceleration = rand(RNG, 3)
    
    # body frame angular acceleration
    angular_acceleration = rand(RNG, 3)

    # gravity vector
    gravity = rand(RNG, 3)

    # --- State Variables --- #

    # define linear and angular displacement
    u = [rand(RNG, 3) for i = 1:nelem+1]
    theta = [rand(RNG, 3) for i = 1:nelem+1]

    # fixed left side
    u[1] .= 0
    theta[1] .= 0 

    # define point velocities (in the deformed point frame)
    V_p = [rand(RNG, 3) for i = 1:nelem+1]
    Omega_p = [rand(RNG, 3) for i = 1:nelem+1]
    
    # define element velocities (in the deformed element frame)
    V_e = [rand(RNG, 3) for i = 1:nelem]
    Omega_e = [rand(RNG, 3) for i = 1:nelem]

    # define point velocities in the body frame
    V = [GXBeam.get_C(theta[i])'*V_p[i] for i = 1:nelem+1]
    Omega = [GXBeam.get_C(theta[i])'*Omega_p[i] for i = 1:nelem+1]

    # define external forces and moments
    F = [rand(RNG, 3) for i = 1:nelem+1]
    M = [rand(RNG, 3) for i = 1:nelem+1]
    
    # loaded right side
    F[nelem+1] .= prescribed_conditions[nelem+1].F + 
        GXBeam.get_C(theta[nelem+1])'*prescribed_conditions[nelem+1].Ff
    M[nelem+1] =  prescribed_conditions[nelem+1].M + 
        GXBeam.get_C(theta[nelem+1])'*prescribed_conditions[nelem+1].Mf
    
    # define resultant forces and moments
    F1 = [rand(RNG, 3) for i = 1:nelem]
    M1 = [rand(RNG, 3) for i = 1:nelem]
    F2 = [rand(RNG, 3) for i = 1:nelem]
    M2 = [rand(RNG, 3) for i = 1:nelem]
    
    # define element internal forces
    Fi = [F1[i]/2 + F2[i]/2 for i = 1:nelem]
    Mi = [M1[i]/2 + M2[i]/2 for i = 1:nelem]

    # --- Test Setting Static States --- #

    static_system = StaticSystem(assembly)

    set_state!(static_system, prescribed_conditions; 
        u, theta, F, M, Fi, Mi)

    static_states = AssemblyState(static_system, assembly; prescribed_conditions)

    for i = 1:nelem+1
        @test isapprox(static_states.points[i].F, F[i])
        @test isapprox(static_states.points[i].M, M[i])
        @test isapprox(static_states.points[i].u, u[i])
        @test isapprox(static_states.points[i].theta, theta[i])
    end

    for i = 1:nelem
        @test isapprox(static_states.elements[i].Fi, Fi[i])
        @test isapprox(static_states.elements[i].Mi, Mi[i])
    end

    # --- Test Setting Dynamic States --- #

    dynamic_system = DynamicSystem(assembly)

    set_state!(dynamic_system, prescribed_conditions; 
        u, theta, V, Omega, F, M, Fi, Mi)

    dynamic_states = AssemblyState(dynamic_system, assembly; prescribed_conditions)

    for i = 1:nelem+1
        @test isapprox(dynamic_states.points[i].F, F[i])
        @test isapprox(dynamic_states.points[i].M, M[i])
        @test isapprox(dynamic_states.points[i].u, u[i])
        @test isapprox(dynamic_states.points[i].theta, theta[i])
        @test isapprox(dynamic_states.points[i].V, V[i])
        @test isapprox(dynamic_states.points[i].Omega, Omega[i])
    end

    for i = 1:nelem
        @test isapprox(dynamic_states.elements[i].Fi, Fi[i])
        @test isapprox(dynamic_states.elements[i].Mi, Mi[i])
        @test isapprox(dynamic_states.points[i].V, V[i])
        @test isapprox(dynamic_states.points[i].Omega, Omega[i])
    end

    # --- Test Setting Expanded States --- #

    expanded_system = ExpandedSystem(assembly)

    set_state!(expanded_system, prescribed_conditions; 
        u, theta, F, M, F1, M1, F2, M2, V_p, Omega_p, V_e, Omega_e)    

    expanded_states = AssemblyState(expanded_system, assembly; prescribed_conditions)

    for i = 1:nelem+1
        @test isapprox(expanded_states.points[i].F, F[i])
        @test isapprox(expanded_states.points[i].M, M[i])
        @test isapprox(expanded_states.points[i].u, u[i])
        @test isapprox(expanded_states.points[i].theta, theta[i])
        @test isapprox(expanded_states.points[i].V, V[i])
        @test isapprox(expanded_states.points[i].Omega, Omega[i])
    end

    for i = 1:nelem
        @test isapprox(expanded_states.elements[i].Fi, Fi[i])
        @test isapprox(expanded_states.elements[i].Mi, Mi[i])
    end   

    # --- Test Copying Expanded States --- #

    # expanded -> static
    reset_state!(static_system)

    copy_state!(static_system, expanded_system, assembly;     
        prescribed_conditions=prescribed_conditions,
        distributed_loads=distributed_loads,
        point_masses=point_masses,
        linear_velocity=linear_velocity,
        angular_velocity=angular_velocity,
        linear_acceleration=linear_acceleration,
        angular_acceleration=angular_acceleration,
        gravity=gravity)

    static_states = AssemblyState(static_system, assembly; prescribed_conditions)

    for i = 1:nelem+1
        @test isapprox(static_states.points[i].F, F[i])
        @test isapprox(static_states.points[i].M, M[i])
        @test isapprox(static_states.points[i].u, u[i])
        @test isapprox(static_states.points[i].theta, theta[i])
    end

    for i = 1:nelem
        @test isapprox(static_states.elements[i].Fi, Fi[i])
        @test isapprox(static_states.elements[i].Mi, Mi[i])
    end

    # expanded -> dynamic
    reset_state!(dynamic_system)

    copy_state!(dynamic_system, expanded_system, assembly;     
        prescribed_conditions=prescribed_conditions,
        distributed_loads=distributed_loads,
        point_masses=point_masses,
        linear_velocity=linear_velocity,
        angular_velocity=angular_velocity,
        linear_acceleration=linear_acceleration,
        angular_acceleration=angular_acceleration,
        gravity=gravity)

    dynamic_states = AssemblyState(dynamic_system, assembly; prescribed_conditions)

    for i = 1:nelem+1
        @test isapprox(dynamic_states.points[i].F, F[i])
        @test isapprox(dynamic_states.points[i].M, M[i])
        @test isapprox(dynamic_states.points[i].u, u[i])
        @test isapprox(dynamic_states.points[i].theta, theta[i])
        @test isapprox(dynamic_states.points[i].V, V[i])
        @test isapprox(dynamic_states.points[i].Omega, Omega[i])
    end

    for i = 1:nelem
        @test isapprox(dynamic_states.elements[i].Fi, Fi[i])
        @test isapprox(dynamic_states.elements[i].Mi, Mi[i])
    end

    # expanded -> expanded
    copy_state!(expanded_system, deepcopy(expanded_system), assembly;     
        prescribed_conditions=prescribed_conditions,
        distributed_loads=distributed_loads,
        point_masses=point_masses,
        linear_velocity=linear_velocity,
        angular_velocity=angular_velocity,
        linear_acceleration=linear_acceleration,
        angular_acceleration=angular_acceleration,
        gravity=gravity)

    expanded_states = AssemblyState(expanded_system, assembly; prescribed_conditions)

    for i = 1:nelem+1
        @test isapprox(expanded_states.points[i].F, F[i])
        @test isapprox(expanded_states.points[i].M, M[i])
        @test isapprox(expanded_states.points[i].u, u[i])
        @test isapprox(expanded_states.points[i].theta, theta[i])
        @test isapprox(expanded_states.points[i].V, V[i])
        @test isapprox(expanded_states.points[i].Omega, Omega[i])
    end

    for i = 1:nelem
        @test isapprox(expanded_states.elements[i].Fi, Fi[i])
        @test isapprox(expanded_states.elements[i].Mi, Mi[i])
    end

    # --- Test Copying Dynamic States --- #

    copy_state!(static_system, dynamic_system, assembly;     
        prescribed_conditions=prescribed_conditions,
        distributed_loads=distributed_loads,
        point_masses=point_masses,
        linear_velocity=linear_velocity,
        angular_velocity=angular_velocity,
        linear_acceleration=linear_acceleration,
        angular_acceleration=angular_acceleration,
        gravity=gravity)

    static_states = AssemblyState(static_system, assembly; prescribed_conditions)

    for i = 1:nelem+1
        @test isapprox(static_states.points[i].F, F[i])
        @test isapprox(static_states.points[i].M, M[i])
        @test isapprox(static_states.points[i].u, u[i])
        @test isapprox(static_states.points[i].theta, theta[i])
    end

    for i = 1:nelem
        @test isapprox(static_states.elements[i].Fi, Fi[i])
        @test isapprox(static_states.elements[i].Mi, Mi[i])
    end

    # dynamic -> dynamic
    copy_state!(dynamic_system, deepcopy(dynamic_system), assembly;     
        prescribed_conditions=prescribed_conditions,
        distributed_loads=distributed_loads,
        point_masses=point_masses,
        linear_velocity=linear_velocity,
        angular_velocity=angular_velocity,
        linear_acceleration=linear_acceleration,
        angular_acceleration=angular_acceleration,
        gravity=gravity)

    dynamic_states = AssemblyState(dynamic_system, assembly; prescribed_conditions)

    for i = 1:nelem+1
        @test isapprox(dynamic_states.points[i].F, F[i])
        @test isapprox(dynamic_states.points[i].M, M[i])
        @test isapprox(dynamic_states.points[i].u, u[i])
        @test isapprox(dynamic_states.points[i].theta, theta[i])
        @test isapprox(dynamic_states.points[i].V, V[i])
        @test isapprox(dynamic_states.points[i].Omega, Omega[i])
    end

    for i = 1:nelem
        @test isapprox(dynamic_states.elements[i].Fi, Fi[i])
        @test isapprox(dynamic_states.elements[i].Mi, Mi[i])
    end

    # dynamic -> expanded
    copy_state!(expanded_system, dynamic_system, assembly;     
        prescribed_conditions=prescribed_conditions,
        distributed_loads=distributed_loads,
        point_masses=point_masses,
        linear_velocity=linear_velocity,
        angular_velocity=angular_velocity,
        linear_acceleration=linear_acceleration,
        angular_acceleration=angular_acceleration,
        gravity=gravity)

    expanded_states = AssemblyState(expanded_system, assembly; prescribed_conditions)

    for i = 1:nelem+1
        @test isapprox(expanded_states.points[i].F, F[i])
        @test isapprox(expanded_states.points[i].M, M[i])
        @test isapprox(expanded_states.points[i].u, u[i])
        @test isapprox(expanded_states.points[i].theta, theta[i])
        @test isapprox(expanded_states.points[i].V, V[i])
        @test isapprox(expanded_states.points[i].Omega, Omega[i])
    end

    for i = 1:nelem
        @test isapprox(expanded_states.elements[i].Fi, Fi[i])
        @test isapprox(expanded_states.elements[i].Mi, Mi[i])
    end

    # --- Test Copying Static States --- #

    # static -> static
    copy_state!(static_system, deepcopy(static_system), assembly;     
        prescribed_conditions=prescribed_conditions,
        distributed_loads=distributed_loads,
        point_masses=point_masses,
        linear_velocity=linear_velocity,
        angular_velocity=angular_velocity,
        linear_acceleration=linear_acceleration,
        angular_acceleration=angular_acceleration,
        gravity=gravity)

    static_states = AssemblyState(static_system, assembly; prescribed_conditions)

    for i = 1:nelem+1
        @test isapprox(static_states.points[i].F, F[i])
        @test isapprox(static_states.points[i].M, M[i])
        @test isapprox(static_states.points[i].u, u[i])
        @test isapprox(static_states.points[i].theta, theta[i])
    end

    for i = 1:nelem
        @test isapprox(static_states.elements[i].Fi, Fi[i])
        @test isapprox(static_states.elements[i].Mi, Mi[i])
    end

    # static -> dynamic
    copy_state!(dynamic_system, static_system, assembly;     
        prescribed_conditions=prescribed_conditions,
        distributed_loads=distributed_loads,
        point_masses=point_masses,
        linear_velocity=linear_velocity,
        angular_velocity=angular_velocity,
        linear_acceleration=linear_acceleration,
        angular_acceleration=angular_acceleration,
        gravity=gravity)

    dynamic_states = AssemblyState(dynamic_system, assembly; prescribed_conditions)

    for i = 1:nelem+1
        @test isapprox(dynamic_states.points[i].F, F[i])
        @test isapprox(dynamic_states.points[i].M, M[i])
        @test isapprox(dynamic_states.points[i].u, u[i])
        @test isapprox(dynamic_states.points[i].theta, theta[i])
    end

    for i = 1:nelem
        @test isapprox(dynamic_states.elements[i].Fi, Fi[i])
        @test isapprox(dynamic_states.elements[i].Mi, Mi[i])
    end

    # static -> expanded
    copy_state!(expanded_system, static_system, assembly;     
        prescribed_conditions=prescribed_conditions,
        distributed_loads=distributed_loads,
        point_masses=point_masses,
        linear_velocity=linear_velocity,
        angular_velocity=angular_velocity,
        linear_acceleration=linear_acceleration,
        angular_acceleration=angular_acceleration,
        gravity=gravity)

    expanded_states = AssemblyState(expanded_system, assembly; prescribed_conditions)

    for i = 1:nelem+1
        @test isapprox(expanded_states.points[i].F, F[i])
        @test isapprox(expanded_states.points[i].M, M[i])
        @test isapprox(expanded_states.points[i].u, u[i])
        @test isapprox(expanded_states.points[i].theta, theta[i])
    end

    for i = 1:nelem
        @test isapprox(expanded_states.elements[i].Fi, Fi[i])
        @test isapprox(expanded_states.elements[i].Mi, Mi[i])
    end

end

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

    damping = fill(rand(RNG, 6), nelem)

    # create assembly of interconnected nonlinear beams
    assembly = Assembly(points, start, stop; stiffness=stiffness, mass=mass, damping=damping)

    # prescribed conditions
    pcond = Dict(
        # fixed left side
        1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
        )

    # distributed loads
    dload = Dict{Int,DistributedLoads{Float64}}()

    # point masses
    pmass = Dict(
        # point mass at the end of the beam
        nelem+1 => PointMass(Symmetric(1e3*rand(RNG,6,6)))
    )

    # gravity vector
    gvec = 1e3*rand(RNG, 3)

    # --- Static Analysis --- #

    system = StaticSystem(assembly)
    force_scaling = system.force_scaling
    indices = system.indices
    x = 1e2 .* rand(RNG, length(system.x))
    J = similar(x, length(x), length(x))

    let

        f = (x) -> GXBeam.static_system_residual!(similar(x), x, indices, force_scaling, 
            assembly, pcond, dload, pmass, gvec)

        GXBeam.static_system_jacobian!(J, x, indices, force_scaling,
            assembly, pcond, dload, pmass, gvec)

        J_fd = ForwardDiff.jacobian(f, x)

        @test all(isapprox.(J, J_fd, atol=1e-10))

    end

    # --- Steady State Analysis --- #

    system = DynamicSystem(assembly)
    force_scaling = system.force_scaling
    indices = system.indices
    icol_accel = GXBeam.body_frame_acceleration_indices(system, pcond)
    x = 1e2 .* rand(RNG, length(system.x))
    J = similar(x, length(x), length(x))

    structural_damping = true

    ub_p = 1e2*rand(RNG, 3)
    θb_p = 1e2*rand(RNG, 3)
    vb_p = 1e2*rand(RNG, 3)
    ωb_p = 1e2*rand(RNG, 3)
    ab_p = 1e2*rand(RNG, 3)
    αb_p = 1e2*rand(RNG, 3)

    let

        f = (x) -> GXBeam.steady_state_system_residual!(similar(x), x, indices, icol_accel, 
            force_scaling, structural_damping, assembly, pcond, dload, pmass, gvec, 
            ub_p, θb_p, vb_p, ωb_p, ab_p, αb_p)

        GXBeam.steady_state_system_jacobian!(J, x, indices, icol_accel, force_scaling, 
            structural_damping, assembly, pcond, dload, pmass, gvec, ub_p, θb_p, 
            vb_p, ωb_p, ab_p, αb_p)

        J_fd = ForwardDiff.jacobian(f, x)

        @test all(isapprox.(J, J_fd, atol=1e-10))

    end

    # --- Initial Condition Analysis --- #

    u0 = [rand(RNG, 3) for ielem = 1:length(assembly.points)]
    theta0 = [rand(RNG, 3) for ielem = 1:length(assembly.points)]
    V0 = [rand(RNG, 3) for ielem = 1:length(assembly.points)]
    Omega0 = [rand(RNG, 3) for ielem = 1:length(assembly.points)]
    Vdot0 = [rand(RNG, 3) for ielem = 1:length(assembly.points)]
    Omegadot0 = [rand(RNG, 3) for ielem = 1:length(assembly.points)]

    x = rand(RNG, length(system.x))
    J = similar(x, length(x), length(x))
    M = similar(x, length(x), length(x))
    rate_vars1 = rand(RNG, Bool, length(x))
    rate_vars2 = rand(RNG, Bool, length(x))

    let

        f = (x) -> GXBeam.initial_condition_system_residual!(similar(x), x, indices, rate_vars1, rate_vars2, 
            icol_accel, force_scaling, structural_damping, assembly, pcond, dload, 
            pmass, gvec, ub_p, θb_p, vb_p, ωb_p, ab_p, αb_p, u0, theta0, V0, Omega0, Vdot0, Omegadot0)

        GXBeam.initial_condition_system_jacobian!(J, x, indices, rate_vars1, rate_vars2, icol_accel, 
            force_scaling, structural_damping, assembly, pcond, dload, pmass, gvec, 
            ub_p, θb_p, vb_p, ωb_p, ab_p, αb_p, u0, theta0, V0, Omega0, Vdot0, Omegadot0)

        J_fd = ForwardDiff.jacobian(f, x)

        @test all(isapprox.(J, J_fd, atol=1e-10))

    end

    # --- Newmark Scheme Time-Marching Analysis --- #

    ubdot = rand(RNG, 3)
    θbdot = rand(RNG, 3)
    vbdot = rand(RNG, 3)
    ωbdot = rand(RNG, 3)
    udot = [rand(RNG, 3) for ipoint = 1:length(assembly.points)]
    θdot = [rand(RNG, 3) for ipoint = 1:length(assembly.points)]
    Vdot = [rand(RNG, 3) for ipoint = 1:length(assembly.points)]
    Ωdot = [rand(RNG, 3) for ipoint = 1:length(assembly.points)]
    dt = rand(RNG)

    let

        f = (x) -> GXBeam.newmark_system_residual!(similar(x), x, indices, icol_accel, 
            force_scaling, structural_damping, assembly, pcond, dload, pmass, gvec, 
            ab_p, αb_p, ubdot, θbdot, vbdot, ωbdot, udot, θdot, Vdot, Ωdot, dt)

        GXBeam.newmark_system_jacobian!(J, x, indices, icol_accel, force_scaling, 
            structural_damping, assembly, pcond, dload, pmass, gvec, 
            ab_p, αb_p, ubdot, θbdot, vbdot, ωbdot, udot, θdot, Vdot, Ωdot, dt)

        J_fd = ForwardDiff.jacobian(f, x)

        @test all(isapprox.(J, J_fd, atol=1e-10))
    
    end

    # --- General Dynamic Analysis --- #

    dx = rand(RNG, length(system.x))
    x = rand(RNG, length(system.x))
    J = similar(x, length(x), length(x))
    M = similar(x, length(x), length(x))

    let

        fx = (x) -> GXBeam.dynamic_system_residual!(similar(x), dx, x, indices, icol_accel, 
            force_scaling, structural_damping, assembly, pcond, dload, pmass, gvec, ab_p, αb_p)
          
        GXBeam.dynamic_system_jacobian!(J, dx, x, indices, icol_accel, force_scaling, 
            structural_damping, assembly, pcond, dload, pmass, gvec, ab_p, αb_p)
    
        J_fd = ForwardDiff.jacobian(fx, x)

        @test all(isapprox.(J, J_fd, atol=1e-10))

    end

    let 
    
        fdx = (dx) -> GXBeam.dynamic_system_residual!(similar(dx), dx, x, indices, icol_accel, 
            force_scaling, structural_damping, assembly, pcond, dload, pmass, gvec, ab_p, αb_p)

        GXBeam.dynamic_system_mass_matrix!(M, x, indices, force_scaling, assembly, pcond, pmass)

        M_fd = ForwardDiff.jacobian(fdx, dx)

        @test all(isapprox.(M, M_fd, atol=1e-10))
    end

    # --- Constant Mass Matrix --- #

    system = ExpandedSystem(assembly)
    force_scaling = system.force_scaling
    indices = system.indices
    x = rand(RNG, length(system.x))
    J = similar(x, length(x), length(x))

    f = (x) -> GXBeam.expanded_steady_system_residual!(similar(x), x, indices, icol_accel, force_scaling, 
        structural_damping, assembly, pcond, dload, pmass, gvec, ub_p, θb_p, vb_p, ωb_p, ab_p, αb_p)

    GXBeam.expanded_steady_system_jacobian!(J, x, indices, icol_accel, force_scaling, 
        structural_damping, assembly, pcond, dload, pmass, gvec, ub_p, θb_p, vb_p, ωb_p, ab_p, αb_p)

    J_fd = ForwardDiff.jacobian(f, x)

    @test all(isapprox.(J, J_fd, atol=1e-10))

    f = (x) -> GXBeam.expanded_dynamic_system_residual!(similar(x), x, indices, icol_accel, force_scaling, 
        structural_damping, assembly, pcond, dload, pmass, gvec, ab_p, αb_p)

    GXBeam.expanded_dynamic_system_jacobian!(J, x, indices, icol_accel, force_scaling, 
        structural_damping, assembly, pcond, dload, pmass, gvec, ab_p, αb_p)

    J_fd = ForwardDiff.jacobian(f, x)

    @test all(isapprox.(J, J_fd, atol=1e-10))

end


# # jacobian tests
# # include("jacobians.jl")

# # cross-section tests
# include("section.jl")

# # examples
# include("examples/cantilever.jl")
# include("examples/overdetermined.jl")
# include("examples/tipforce.jl")
# include("examples/tipmoment.jl")
# include("examples/curved.jl")
# include("examples/rotating.jl")
# include("examples/wind-turbine-blade.jl")
# include("examples/static-joined-wing.jl")
# include("examples/dynamic-joined-wing.jl")

# # interfaces
# include("interfaces/diffeq.jl")
# include("interfaces/forwarddiff.jl")

# # issues
# include("issues/zeros.jl")
# include("issues/pointmass.jl")




