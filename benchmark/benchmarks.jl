using GXBeam, LinearAlgebra, BenchmarkTools

function steady_benchmark()

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

    # create dictionary of prescribed conditions
    prescribed_conditions = Dict(
        # root section is fixed
        1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0)
        )

    # perform nonlinear steady state analysis
    system, state, converged = steady_state_analysis(assembly,
        linear_velocity = [0, 196.349540849362, 0],
        angular_velocity = [0, 0, 78.5398163397448],
        prescribed_conditions = prescribed_conditions)

    # test convergence
    @assert converged

    return state
end

function eigenvalue_benchmark()

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

    # create dictionary of prescribed conditions
    prescribed_conditions = Dict(
        # root section is fixed
        1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0)
        )

    # perform eigenvalue analysis
    system, λ, V, converged = eigenvalue_analysis(assembly; nev = 20,
        linear_velocity = [0, 196.349540849362, 0],
        angular_velocity = [0, 0, 78.5398163397448],
        prescribed_conditions = prescribed_conditions)

    # test convergence
    @assert converged

    return λ, V
end

function unsteady_benchmark()

    tvec = range(0, 0.5, length=2001)

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

    # time function
    tf1 = (t) -> t < 0.2 ? 0.1*(1 - sin(2*pi*(t/0.4 + 0.25))) : zero(t)

    # # prescribed condition
    # prescribed_conditions = (t) -> Dict(
    #     1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
    #     nelem_b1 + 1 => PrescribedConditions(Fx=0, Fy=100*tf1(t), Fz=0, Mx=0, My=0, Mz=0)
    # )

    # prescribed condition
    prescribed_conditions = Dict(
        1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
        nelem_b1 + 1 => PrescribedConditions(Fx=0, Fy=0, Fz=0, Mx=0, My=0, Mz=0)
    )

    # perform time marching analysis
    system, history, converged = time_domain_analysis(assembly, tvec; show_trace=true,
        linear_velocity = [0, 196.349540849362, 0],
        angular_velocity = [0, 0, 78.5398163397448],
        prescribed_conditions = prescribed_conditions)

    # test convergence
    @assert converged

    return history
end

println("Steady State Analysis: ")
@btime steady_benchmark();

# Results on a single i7-5500U CPU @ 2.40GHz core.
# GEBT: 13.722 ms
# GXBeam: 4.716 ms (2207 allocations: 2.79 MiB)

println("Eigenvalue Analysis: ")
@btime eigenvalue_benchmark();

# Results on a single i7-5500U CPU @ 2.40GHz core.
# GEBT: 33.712 ms
# GXBeam: 18.478 ms (2568 allocations: 5.16 MiB)

println("Time Marching Analysis: ")
@btime unsteady_benchmark();

# Results on a single i7-5500U CPU @ 2.40GHz core.
# GEBT: 26.870 s
# GXBeam: 9.019 s (851711 allocations: 4.30 GiB)
