using GXBeam, LinearAlgebra, Test

@testset "Nonlinear Static Analysis of a Joined-Wing" begin

    # Set endpoints of each beam
    p1 = [-7.1726, -12, -3.21539]
    p2 = [-5.37945, -9, -2.41154]
    p3 = [-3.5863, -6, -1.6077]
    p4 = [-1.79315, -3, -0.803848]
    p5 = [0, 0, 0]
    p6 = [7.1726, -12, 3.21539]

    # get transformation matrix for left beams

    # transformation from intermediate to global frame
    tmp1 = sqrt(p1[1]^2 + p1[2]^2)
    c1, s1 = -p1[1]/tmp1, -p1[2]/tmp1
    rot1 = [c1 -s1 0; s1 c1 0; 0 0 1]

    # transformation from local to intermediate frame
    tmp2 = sqrt(p1[1]^2 + p1[2]^2 + p1[3]^2)
    c2, s2 = tmp1/tmp2, -p1[3]/tmp2
    rot2 = [c2 0 -s2; 0 1 0; s2 0 c2]

    Cab_1 = rot1*rot2

    # get transformation matrix for right beam

    # transformation from intermediate frame to global frame
    tmp1 = sqrt(p6[1]^2 + p6[2]^2)
    c1, s1 = p6[1]/tmp1, p6[2]/tmp1
    rot1 = [c1 -s1 0; s1 c1 0; 0 0 1]

    # transformation from local beam frame to intermediate frame
    tmp2 = sqrt(p6[1]^2 + p6[2]^2 + p6[3]^2)
    c2, s2 = tmp1/tmp2, p6[3]/tmp2
    rot2 = [c2 0 -s2; 0 1 0; s2 0 c2]

    Cab_2 = rot1*rot2

    # beam 1
    L_b1 = norm(p2-p1)
    r_b1 = p1
    nelem_b1 = 5
    lengths_b1, xp_b1, xm_b1, Cab_b1 = discretize_beam(L_b1, r_b1, nelem_b1, frame=Cab_1)
    compliance_b1 = fill(Diagonal([1.05204e-9, 3.19659e-9, 2.13106e-8, 1.15475e-7, 1.52885e-7, 7.1672e-9]), nelem_b1)

    # beam 2
    L_b2 = norm(p3-p2)
    r_b2 = p2
    nelem_b2 = 5
    lengths_b2, xp_b2, xm_b2, Cab_b2 = discretize_beam(L_b2, r_b2, nelem_b2, frame=Cab_1)
    compliance_b2 = fill(Diagonal([1.24467e-9, 3.77682e-9, 2.51788e-8, 1.90461e-7, 2.55034e-7, 1.18646e-8]), nelem_b2)

    # beam 3
    L_b3 = norm(p4-p3)
    r_b3 = p3
    nelem_b3 = 5
    lengths_b3, xp_b3, xm_b3, Cab_b3 = discretize_beam(L_b3, r_b3, nelem_b3, frame=Cab_1)
    compliance_b3 = fill(Diagonal([1.60806e-9, 4.86724e-9, 3.24482e-8, 4.07637e-7, 5.57611e-7, 2.55684e-8]), nelem_b3)

    # beam 4
    L_b4 = norm(p5-p4)
    r_b4 = p4
    nelem_b4 = 5
    lengths_b4, xp_b4, xm_b4, Cab_b4 = discretize_beam(L_b4, r_b4, nelem_b4, frame=Cab_1)
    compliance_b4 = fill(Diagonal([2.56482e-9, 7.60456e-9, 5.67609e-8, 1.92171e-6, 2.8757e-6, 1.02718e-7]), nelem_b4)

    # beam 5
    L_b5 = norm(p6-p5)
    r_b5 = p5
    nelem_b5 = 20
    lengths_b5, xp_b5, xm_b5, Cab_b5 = discretize_beam(L_b5, r_b5, nelem_b5, frame=Cab_2)
    compliance_b5 = fill(Diagonal([2.77393e-9, 7.60456e-9, 1.52091e-7, 1.27757e-5, 2.7835e-5, 1.26026e-7]), nelem_b5)

    # combine elements and points into one array
    nelem = nelem_b1 + nelem_b2 + nelem_b3 + nelem_b4 + nelem_b5
    points = vcat(xp_b1, xp_b2[2:end], xp_b3[2:end], xp_b4[2:end], xp_b5[2:end])
    start = 1:nelem
    stop = 2:nelem + 1
    lengths = vcat(lengths_b1, lengths_b2, lengths_b3, lengths_b4, lengths_b5)
    midpoints = vcat(xm_b1, xm_b2, xm_b3, xm_b4, xm_b5)
    Cab = vcat(Cab_b1, Cab_b2, Cab_b3, Cab_b4, Cab_b5)
    compliance = vcat(compliance_b1, compliance_b2, compliance_b3, compliance_b4, compliance_b5)

    # create assembly
    assembly = Assembly(points, start, stop, compliance=compliance,
        frames=Cab, lengths=lengths, midpoints=midpoints)

    Fz = range(0, 70e3, length=141)

    # pre-allocate memory to reduce run-time
    system = StaticSystem(assembly)

    linear_states = Vector{AssemblyState{Float64}}(undef, length(Fz))
    for i = 1:length(Fz)

        # create dictionary of prescribed conditions
        prescribed_conditions = Dict(
            # fixed endpoint on beam 1
            1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
            # force applied on point 4
            nelem_b1 + nelem_b2 + nelem_b3 + nelem_b4 + 1 => PrescribedConditions(Fz = Fz[i]),
            # fixed endpoint on last beam
            nelem+1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
        )

        _, converged = static_analysis!(system, assembly, prescribed_conditions=prescribed_conditions, linear=true)

        linear_states[i] = AssemblyState(system, assembly; prescribed_conditions=prescribed_conditions)

        @test converged
    end

    reset_state!(system)
    nonlinear_states = Vector{AssemblyState{Float64}}(undef, length(Fz))
    for i = 1:length(Fz)

        # create dictionary of prescribed conditions
        prescribed_conditions = Dict(
            # fixed endpoint on beam 1
            1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
            # force applied on point 4
            nelem_b1 + nelem_b2 + nelem_b3 + nelem_b4 + 1 => PrescribedConditions(Fz = Fz[i]),
            # fixed endpoint on last beam
            nelem+1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
        )

        _, converged = static_analysis!(system, assembly, prescribed_conditions=prescribed_conditions,
            reset_state = false)

        nonlinear_states[i] = AssemblyState(system, assembly;
            prescribed_conditions=prescribed_conditions)

        @test converged
    end

    reset_state!(system)
    nonlinear_follower_states = Vector{AssemblyState{Float64}}(undef, length(Fz))
    for i = 1:length(Fz)
        # create dictionary of prescribed conditions
        prescribed_conditions = Dict(
            # fixed endpoint on beam 1
            1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
            # force applied on point 4
            nelem_b1 + nelem_b2 + nelem_b3 + nelem_b4 + 1 => PrescribedConditions(Fz_follower = Fz[i]),
            # fixed endpoint on last beam
            nelem+1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
        )

        _, converged = static_analysis!(system, assembly, prescribed_conditions=prescribed_conditions,
            reset_state = false)

        nonlinear_follower_states[i] = AssemblyState(system, assembly;
            prescribed_conditions=prescribed_conditions)

        @test converged
    end
end