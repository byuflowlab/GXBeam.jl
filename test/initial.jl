using GXBeam, LinearAlgebra, UnPack, Test

@testset "Initial Condition Analysis" begin

    # straight section of the beam
    L_b1 = 31.5 # inch
    r_b1 = [0, 0, 0]
    nelem_b1 = 20
    lengths_b1, xp_b1, xm_b1, Cab_b1 = discretize_beam(L_b1, r_b1, nelem_b1)

    # swept section of the beam
    L_b2 = 6 # inch
    r_b2 = [31.5, 0, 0]
    nelem_b2 = 4
    frame_b2 = [
        0.707106781186548 0.707106781186548 0.0;
        -0.707106781186548 0.707106781186548 0.0;
        0.0               0.0               1.0]
    lengths_b2, xp_b2, xm_b2, Cab_b2 = discretize_beam(L_b2, r_b2, nelem_b2;
        frame = frame_b2)

    # combine elements and points into one array
    nelem = nelem_b1 + nelem_b2
    points = vcat(xp_b1, xp_b2[2:end])
    start = 1:nelem_b1 + nelem_b2
    stop = 2:nelem_b1 + nelem_b2 + 1
    lengths = vcat(lengths_b1, lengths_b2)
    midpoints = vcat(xm_b1, xm_b2)
    Cab = vcat(Cab_b1, Cab_b2)

    compliance = fill(Diagonal([
        1.4974543276E-06,
        4.7619054919E-06,
        5.8036221902E-05,
        3.1234387938E-03,
        4.5274507258E-03,
        1.7969451932E-05]), nelem)

    mass = fill(Diagonal([
        1.5813000000E-05,
        1.5813000000E-05,
        1.5813000000E-05,
        1.3229801498E-06,
        5.2301497500E-09,
        1.3177500000E-06]), nelem)

    # create assembly
    assembly = Assembly(points, start, stop;
        compliance = compliance,
        mass = mass,
        frames = Cab,
        lengths = lengths,
        midpoints = midpoints)

    # prescribed condition
    prescribed_conditions = Dict(
        1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
        nelem_b1 + 1 => PrescribedConditions(Fx=0, Fy=0, Fz=0, Mx=0, My=0, Mz=0)
    )

    # set linear and angular velocity
    linear_velocity = [0, 196.349540849362, 0]
    angular_velocity = [0, 0, 78.5398163397448]

    # perform time marching analysis
    system, state, converged = initial_condition_analysis(assembly, 0.0;
        structural_damping = false,
        linear_velocity = linear_velocity,
        angular_velocity = angular_velocity,
        prescribed_conditions = prescribed_conditions)

    # unpack system storage
    @unpack x, r, K, M, force_scaling, indices = system

    # set default inputs
    two_dimensional=false
    structural_damping=false
    distributed_loads=Dict{Int,DistributedLoads{Float64}}()
    point_masses=Dict{Int,PointMass{Float64}}()
    gravity=[0.0, 0.0, 0.0]

    # construct state vector
    x = set_state!(system.x, system, assembly, state; prescribed_conditions)

    # construct rate vector
    dx = similar(x) .= 0
    for ipoint = 1:length(assembly.points)
        icol = indices.icol_point[ipoint]
        udot_udot, θdot_θdot = GXBeam.point_displacement_jacobians(ipoint, prescribed_conditions)
        dx[icol:icol+2] = udot_udot*state.points[ipoint].udot
        dx[icol+3:icol+5] = θdot_θdot*state.points[ipoint].thetadot
        dx[icol+6:icol+8] = state.points[ipoint].Vdot
        dx[icol+9:icol+11] = state.points[ipoint].Omegadot
    end

    # compute residual
    GXBeam.dynamic_system_residual!(r, dx, x, indices, two_dimensional, force_scaling,
        structural_damping, assembly, prescribed_conditions, distributed_loads, point_masses,
        gravity, linear_velocity, angular_velocity)

    # check if residual is converged
    @test norm(r) <= 1e-9

end