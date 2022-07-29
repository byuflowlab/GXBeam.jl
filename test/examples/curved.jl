@testset "Nonlinear Analysis of the Bending of a Curved Beam in 3D Space" begin

    # problem constants
    R = 100
    L = R*pi/4 # inches
    h = w = 1 # inches
    E = 1e7 # psi Young's Modulus
    ν = 0.0
    G = E/(2*(1+ν))

    # beam starting point, frame, and curvature
    r = [0, 0, 0]
    frame = [0 -1 0; 1 0 0; 0 0 1]
    curvature = [0, 0, -1/R]

    # cross section properties
    A = h*w
    Ay = A
    Az = A
    Iyy = w*h^3/12
    Izz = w^3*h/12
    J = Iyy + Izz

    # discretize the beam
    nelem = 16
    ΔL, xp, xm, Cab = discretize_beam(L, r, nelem; frame=frame, curvature = curvature)

    # force
    P = 600 # lbs

    # index of left and right endpoints of each beam element
    start = 1:nelem
    stop = 2:nelem+1

    # compliance matrix for each beam element
    compliance = fill(Diagonal([1/(E*A), 1/(G*Ay), 1/(G*Az), 1/(G*J), 1/(E*Iyy), 1/(E*Izz)]), nelem)

    # create assembly of interconnected nonlinear beams
    assembly = Assembly(xp, start, stop, compliance=compliance, frames=Cab,
        lengths=ΔL, midpoints=xm)

    # create dictionary of prescribed conditions
    prescribed_conditions = Dict(
        # fixed left endpoint
        1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
        # force on right endpoint
        nelem+1 => PrescribedConditions(Fz=P)
    )

    # perform static analysis
    system, converged = static_analysis(assembly, prescribed_conditions=prescribed_conditions)

    # post-process results
    state = AssemblyState(system, assembly, prescribed_conditions=prescribed_conditions)

    # test tip deflections
    @test isapprox(state.points[end].u[1], -13.4, atol=0.2) # -13.577383726758564
    @test isapprox(state.points[end].u[2], -23.5, atol=0.1) # -23.545303336988038
    @test isapprox(state.points[end].u[3],  53.4, atol=0.1) #  53.45800757548929

    # Results from "Large Displacement Analysis of Three-Dimensional Beam
    # Structures" by Bathe and Bolourch:
    # - Tip Displacement: [-13.4, -23.5, 53.4]

    # Note that these results are comparing computational solutions, rather than
    # the computational to the analytical solution, so some variation is expected.

    # perform the same analysis for a constant mass matrix system
    system, converged = steady_state_analysis(assembly, 
        prescribed_conditions = prescribed_conditions,
        constant_mass_matrix = true)

    state = AssemblyState(system, assembly, prescribed_conditions=prescribed_conditions)
    
    # test tip deflections
    @test isapprox(state.points[end].u[1], -13.4, atol=0.2) # -13.577383726758564
    @test isapprox(state.points[end].u[2], -23.5, atol=0.1) # -23.545303336988038
    @test isapprox(state.points[end].u[3],  53.4, atol=0.1) #  53.45800757548929

end