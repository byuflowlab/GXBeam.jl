@testset "Point Masses" begin

    nodes = [[0,i,0] for i in 0:0.1:1]
    nelem = length(nodes)-1
    start = 1:nelem
    stop =  2:(nelem+1)

    stiff = [
        1.0e6  0.0    0.0    -0.5   -1.0  -50000.0
        0.0    3.0e6  0.0     0.0    0.0       0.0
        0.0    0.0    3.0e6   0.0    0.0       0.0
       -0.5    0.0    0.0     7.0    0.1      -0.02
       -1.0    0.0    0.0     0.1    5.0       0.1
   -50000.0    0.0    0.0    -0.02   0.1    3000.0
    ]

    mass = [
        0.02    0.0      0.0     0.0      -5.0e-7  -1.0e-7
        0.0     0.02     0.0     5.0e-7    0.0      0.0001
        0.0     0.0      0.02    1.0e-7   -0.0001   0.0
        0.0     5.0e-7   1.0e-7  1.0e-5    1.0e-8   2.0e-10
       -5.0e-7  0.0     -0.0001  1.0e-8    6.0e-7   9.0e-9
       -1.0e-7  0.0001   0.0     2.0e-10   9.0e-9   1.0e-5
    ]

    transformation = [0 -1 0; 1 0 0; 0 0 1]

    assembly = GXBeam.Assembly(nodes, start, stop; 
        frames = fill(transformation, nelem),
        stiffness = fill(stiff, nelem));

    pmass = GXBeam.transform_properties(mass, transformation)

    point_masses = Dict(1 => PointMass(pmass./2))
    for i = 2:nelem
        point_masses[i] = PointMass(pmass)
    end
    point_masses[nelem+1] = PointMass(pmass./2)

    prescribed_conditions = Dict(1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0))
    
    system, λ, V, converged = GXBeam.eigenvalue_analysis(assembly;
        prescribed_conditions = prescribed_conditions, 
        point_masses = point_masses,
        nev = 14);

    imagλ = imag(λ)
    isort = sortperm(abs.(imagλ))
    freq = imagλ[isort[1:2:10]]/(2*pi)

    # frequencies
    frequencies = [
        2.8182004347800804, 
        17.66611982731975, 
        27.978670985969078, 
        49.93431945680836, 
        66.07594270678581]

    @test isapprox(freq, frequencies)

end