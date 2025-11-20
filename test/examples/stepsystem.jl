using GXBeam, LinearAlgebra, Test

@testset "Step system using wind turbine blade" begin

    L = 60 # m

    # create points
    nelem = 5
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
    tvec = 0:0.001:0.1

    # prescribed conditions
    prescribed_conditions = (t) -> begin
        Dict(
            # fixed left side
            1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
            # force on right side
            nelem+1 => PrescribedConditions(Fz=1e5*sin(20*t))
            )
    end

    system = GXBeam.DynamicSystem(assembly)

    nt = length(tvec)

    inittype = eltype(system.x)

    gxhistory = Array{GXBeam.AssemblyState{inittype, 
        Vector{GXBeam.PointState{inittype}},
        Vector{GXBeam.ElementState{inittype}}}}(undef, nt)

    converged = Array{Bool}(undef, nt)

    system, gxstate, constants, paug, xgx, converged[1] = GXBeam.initialize_system!(system, assembly, tvec; prescribed_conditions, structural_damping=true, reset_state=true) 
    gxhistory[1] = gxstate[1]
    
    for i in 2:nt
        system, gxhistory[i], constants, paug, xgx, converged[i] = GXBeam.step_system!(system, paug, xgx, constants, gxhistory[i-1], assembly, tvec, i; prescribed_conditions, structural_damping=true)
    end


    @test all(converged)

end