using GXBeam, LinearAlgebra, Test

@testset "Zero Mass Matrix" begin
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

    mass = fill(Diagonal(zeros(6)), nelem)

    # create assembly
    assembly = Assembly(points, start, stop, compliance=compliance, mass=mass, frames=Cab, lengths=lengths, midpoints=midpoints)

    # create dictionary of prescribed conditions
    prescribed_conditions = Dict(
        # root section is fixed
        1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0)
        )

    # set angular velocity vector
    w0 = [0, 0, rpm*(2*pi)/60]

    # perform nonlinear steady state analysis
    system, converged = steady_state_analysis(assembly,
        angular_velocity = w0,
        prescribed_conditions = prescribed_conditions)

    # test convergence
    @test converged
end

@testset "Zero Length Element" begin
    sweep = 45 * pi/180
    rpm = 750

    # straight section of the beam
    L_b1 = 31.5 # inch
    r_b1 = [2.5, 0, 0]
    nelem_b1 = 13
    lengths_b1, xp_b1, xm_b1, Cab_b1 = discretize_beam(L_b1, r_b1, nelem_b1)

    # zero length element between straight and swept sections
    L_b12 = 0
    r_b12 = [34, 0, 0]
    nelem_b12 = 1
    lengths_b12, xp_b12, xm_b12, Cab_b12 = discretize_beam(L_b12, r_b12, nelem_b12)

    # swept section of the beam
    L_b2 = 6 # inch
    r_b2 = [34, 0, 0]
    nelem_b2 = 3
    cs, ss = cos(sweep), sin(sweep)
    frame_b2 = [cs ss 0; -ss cs 0; 0 0 1]
    lengths_b2, xp_b2, xm_b2, Cab_b2 = discretize_beam(L_b2, r_b2, nelem_b2, frame=frame_b2)

    # combine elements and points into one array
    nelem = nelem_b1 + nelem_b12 + nelem_b2
    points = vcat(xp_b1, xp_b2)
    lengths = vcat(lengths_b1, lengths_b12, lengths_b2)
    midpoints = vcat(xm_b1, xm_b12, xm_b2)
    Cab = vcat(Cab_b1, Cab_b12, Cab_b2)

    # specify connectivity
    start = 1:nelem
    stop = 2:nelem+1

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

    # set angular velocity vector
    w0 = [0, 0, rpm*(2*pi)/60]

    # perform nonlinear steady state analysis
    system, converged = steady_state_analysis(assembly,
        angular_velocity = w0,
        prescribed_conditions = prescribed_conditions)

    # test convergence
    @test converged
end

@testset "Massless Element Time Domain Initialization" begin

    comp = Symmetric([
        6.00001e-6  0.0         0.0          7.25923e-7  -8.1452e-7    0.0001;
        0.0         3.33333e-7  0.0          0.0          0.0          0.0;
        0.0         0.0         3.33333e-7   0.0          0.0          0.0;
        7.25923e-7  0.0         0.0          0.142898    -0.00285808   1.31466e-5;
       -8.1452e-7   0.0         0.0         -0.00285808   0.200057    -2.0263e-5;
        0.0001      0.0         0.0          1.31466e-5  -2.0263e-5    0.002;
        ])
    
    mass = Symmetric([
        0.02    0.0      0.0     0.0      -5.0e-7  -1.0e-7;
        0.0     0.02     0.0     5.0e-7    0.0      0.0001;
        0.0     0.0      0.02    1.0e-7   -0.0001   0.0;
        0.0     5.0e-7   1.0e-7  1.0e-5    1.0e-8   2.0e-10;
       -5.0e-7  0.0     -0.0001  1.0e-8    6.0e-7   9.0e-9;
       -1.0e-7  0.0001   0.0     2.0e-10   9.0e-9   1.0e-5;
        ])
    
    nodes = [[0,i,0] for i in 0:.1:1]
    
    nElements = length(nodes)-1
    start = 1:nElements
    stop =  2:nElements+1
    transformation = [[0 -1 0; 1 0 0; 0 0 1] for _ in 1:nElements];
    
    compliance = [comp for i in 1:nElements]
    
    pointmass = Dict(2 => PointMass(GXBeam.transform_properties(mass, transformation[2]')))
    for i in 4:2:nElements
        pointmass[i] = PointMass(GXBeam.transform_properties(mass, transformation[i]'))
    end
    pointmass[nElements] = PointMass(GXBeam.transform_properties(mass, transformation[nElements]')./2) # last lumped mass is half of the others, as it represents the last half of an element
    
    assembly = GXBeam.Assembly(nodes, start, stop, 
        compliance=compliance, 
        frames=transformation);
    
    prescribed_conditions = Dict(
        1 => GXBeam.PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
        length(nodes) => GXBeam.PrescribedConditions(Fz=30, My=-0.2))
    
    t0 = 0.0
    
    system, converged = initial_condition_analysis(assembly, t0;
                                                   prescribed_conditions = prescribed_conditions,
                                                   point_masses = pointmass);
    
    @test converged

end