using GXBeam, LinearAlgebra, Test

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
        prescribed_conditions = prescribed_conditions,
        structural_damping = false)

    @test converged

    system, history, converged = time_domain_analysis(assembly, tvec;
        prescribed_conditions = prescribed_conditions,
        structural_damping = true)

    @test converged
end