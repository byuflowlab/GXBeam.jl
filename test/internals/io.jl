using GXBeam, Random, Test

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