using GXBeam, LinearAlgebra, UnPack, Test

@testset "Initial Condition Analysis" begin
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

    mass = fill(Diagonal([ρ*A, ρ*A, ρ*A, ρ*J, ρ*Iyy, ρ*Izz]), nelem)

    # create assembly
    assembly = Assembly(points, start, stop, compliance=compliance, mass=mass, frames=Cab, lengths=lengths, midpoints=midpoints)

    # create dictionary of prescribed conditions
    prescribed_conditions = Dict(
        # root section is fixed
        1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0)
        )

    # global frame rotation
    angular_velocity = [0, 0, rpm*(2*pi)/60]

    # current time
    t0 = 0.0

    # initialize time domain solution
    system, state, converged = initial_condition_analysis(assembly, t0;
        prescribed_conditions = prescribed_conditions,
        angular_velocity = angular_velocity)

    # unpack system storage
    @unpack x, r, K, M, force_scaling, indices = system

    # set default inputs
    two_dimensional=false
    structural_damping=true
    distributed_loads=Dict{Int,DistributedLoads{Float64}}()
    point_masses=Dict{Int,PointMass{Float64}}()
    gravity=[0.0, 0.0, 0.0]
    linear_velocity=[0.0, 0.0, 0.0]
    angular_velocity=angular_velocity

    # construct state vector
    x = set_state!(system.x, system, state, prescribed_conditions)

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