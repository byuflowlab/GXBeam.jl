using GXBeam, LinearAlgebra, Test

@testset "Rotating Beam with a Swept Tip" begin
    sweep = 45 * pi/180
    rpm = 0:25:750

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

    nonlinear_states = Vector{AssemblyState{Float64}}(undef, length(rpm))
    linear_states = Vector{AssemblyState{Float64}}(undef, length(rpm))
    for i = 1:length(rpm)
        # global frame rotation
        w0 = [0, 0, rpm[i]*(2*pi)/60]

        # perform nonlinear steady state analysis
        system, converged = steady_state_analysis(assembly,
            angular_velocity = w0,
            prescribed_conditions = prescribed_conditions)

        nonlinear_states[i] = AssemblyState(system, assembly, prescribed_conditions=prescribed_conditions)

        # perform linear steady state analysis
        system, converged = steady_state_analysis(assembly,
            angular_velocity = w0,
            prescribed_conditions = prescribed_conditions,
            linear = true)

        linear_states[i] = AssemblyState(system, assembly, prescribed_conditions=prescribed_conditions)
    end

    sweep = (0:2.5:45) * pi/180
    rpm = [0, 500, 750]
    nev = 30

    λ = Matrix{Vector{ComplexF64}}(undef, length(sweep), length(rpm))
    U = Matrix{Matrix{ComplexF64}}(undef, length(sweep), length(rpm))
    MV = Matrix{Matrix{ComplexF64}}(undef, length(sweep), length(rpm))
    state = Matrix{AssemblyState{Float64}}(undef, length(sweep), length(rpm))
    eigenstates = Matrix{Vector{AssemblyState{ComplexF64}}}(undef, length(sweep), length(rpm))
    for i = 1:length(sweep)
        local L_b1, r_b1, nelem_b1, lengths_b1 #hide
        local xp_b1, xm_b1, Cab_b1 #hide
        local cs, ss #hide
        local L_b2, r_b2, nelem_b2, frame_b2, lengths_b2 #hide
        local xp_b2, xm_b2, Cab_b2 #hide
        local nelem, points, start, stop #hide
        local lengths, midpoints, Cab, compliance, mass, assembly #hide

        # straight section of the beam
        L_b1 = 31.5 # inch
        r_b1 = [2.5, 0, 0]
        nelem_b1 = 20
        lengths_b1, xp_b1, xm_b1, Cab_b1 = discretize_beam(L_b1, r_b1, nelem_b1)

        # swept section of the beam
        L_b2 = 6 # inch
        r_b2 = [34, 0, 0]
        nelem_b2 = 20
        cs, ss = cos(sweep[i]), sin(sweep[i])
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

        compliance = fill(Diagonal([1/(E*A), 1/(G*Ay), 1/(G*Az), 1/(G*Jx), 1/(E*Iyy), 1/(E*Izz)]), nelem)

        mass = fill(Diagonal([ρ*A, ρ*A, ρ*A, ρ*J, ρ*Iyy, ρ*Izz]), nelem)

        # create assembly
        assembly = Assembly(points, start, stop, compliance=compliance, mass=mass, frames=Cab, lengths=lengths, midpoints=midpoints)

        # create system
        system = DynamicSystem(assembly)

        for j = 1:length(rpm)
            # global frame rotation
            w0 = [0, 0, rpm[j]*(2*pi)/60]

            # eigenvalues and (right) eigenvectors
            system, λ[i,j], V, converged = eigenvalue_analysis!(system, assembly,
                angular_velocity = w0,
                prescribed_conditions = prescribed_conditions,
                nev=nev)

            # corresponding left eigenvectors
            U[i,j] = left_eigenvectors(system, λ[i,j], V)

            # post-multiply mass matrix with right eigenvector matrix
            # (we use this later for correlating eigenvalues)
            MV[i,j] = system.M * V

            # process state and eigenstates
            state[i,j] = AssemblyState(system, assembly; prescribed_conditions=prescribed_conditions)
            eigenstates[i,j] = [AssemblyState(system, assembly, V[:,k];
                prescribed_conditions=prescribed_conditions) for k = 1:nev]
        end
    end

    # set previous left eigenvector matrix
    U_p = copy(U[1,1])

    for j = 1:length(rpm)
        for i = 1:length(sweep)
            # construct correlation matrix
            C = U_p*MV[i,j]

            # correlate eigenmodes
            perm, corruption = correlate_eigenmodes(C)

            # re-arrange eigenvalues and eigenvectors
            λ[i,j] = λ[i,j][perm]
            U[i,j] = U[i,j][perm,:]
            MV[i,j] = MV[i,j][:,perm]
            eigenstates[i,j] = eigenstates[i,j][perm]

            # update previous eigenvector matrix
            U_p .= U[i,j]
        end
        # update previous eigenvector matrix
        U_p .= U[1,j]
    end

    frequency = [[imag(λ[i,j][k])/(2*pi) for i = 1:length(sweep), j=1:length(rpm)] for k = 1:2:nev]

    indices = [1, 2, 4]
    experiment_rpm = [0, 500, 750]
    experiment_sweep = [0, 15, 30, 45]
    experiment_frequencies = [
        [1.4 1.8 1.7 1.6;
         10.2 10.1 10.2 10.2;
         14.8 14.4 14.9 14.7],
        [10.3 10.2 10.4 10.4;
         25.2 25.2 23.7 21.6;
         36.1 34.8 30.7 26.1],
        [27.7 27.2 26.6 24.8;
         47.0 44.4 39.3 35.1;
         62.9 55.9 48.6 44.8]
    ]

    for k = 1:length(experiment_frequencies)
        for j = 1:length(experiment_sweep)
            for i = 1:length(experiment_rpm)
                ii = argmin(abs.(rpm .- experiment_rpm[i]))
                jj = argmin(abs.(sweep*180/pi .- experiment_sweep[j]))
                kk = indices[k]
                @test isapprox(frequency[kk][jj,ii], experiment_frequencies[k][i,j], atol=1, rtol=0.1)
            end
        end
    end

    indices = [5, 7, 6]
    experiment_frequencies = [
        95.4 87.5 83.7 78.8;
        106.6 120.1 122.6 117.7;
        132.7 147.3 166.2 162.0
    ]

    for k = 1:size(experiment_frequencies, 1)
        for j = 1:length(experiment_sweep)
            ii = argmin(abs.(rpm .- 750))
            jj = argmin(abs.(sweep*180/pi .- experiment_sweep[j]))
            kk = indices[k]
            @test isapprox(frequency[kk][jj,ii], experiment_frequencies[k,j], rtol=0.1)
        end
    end

    # # perform the same analysis for a constant mass matrix system
    # sweep = (0:2.5:45) * pi/180
    # rpm = [0, 500, 750]
    # nev = 30

    # λ = Matrix{Vector{ComplexF64}}(undef, length(sweep), length(rpm))
    # U = Matrix{Matrix{ComplexF64}}(undef, length(sweep), length(rpm))
    # MV = Matrix{Matrix{ComplexF64}}(undef, length(sweep), length(rpm))
    # state = Matrix{AssemblyState{Float64}}(undef, length(sweep), length(rpm))
    # eigenstates = Matrix{Vector{AssemblyState{ComplexF64}}}(undef, length(sweep), length(rpm))
    # for i = 1:length(sweep)
    #     local L_b1, r_b1, nelem_b1, lengths_b1 #hide
    #     local xp_b1, xm_b1, Cab_b1 #hide
    #     local cs, ss #hide
    #     local L_b2, r_b2, nelem_b2, frame_b2, lengths_b2 #hide
    #     local xp_b2, xm_b2, Cab_b2 #hide
    #     local nelem, points, start, stop #hide
    #     local lengths, midpoints, Cab, compliance, mass, assembly #hide

    #     # straight section of the beam
    #     L_b1 = 31.5 # inch
    #     r_b1 = [2.5, 0, 0]
    #     nelem_b1 = 20
    #     lengths_b1, xp_b1, xm_b1, Cab_b1 = discretize_beam(L_b1, r_b1, nelem_b1)

    #     # swept section of the beam
    #     L_b2 = 6 # inch
    #     r_b2 = [34, 0, 0]
    #     nelem_b2 = 20
    #     cs, ss = cos(sweep[i]), sin(sweep[i])
    #     frame_b2 = [cs ss 0; -ss cs 0; 0 0 1]
    #     lengths_b2, xp_b2, xm_b2, Cab_b2 = discretize_beam(L_b2, r_b2, nelem_b2, frame=frame_b2)

    #     # combine elements and points into one array
    #     nelem = nelem_b1 + nelem_b2
    #     points = vcat(xp_b1, xp_b2[2:end])
    #     start = 1:nelem_b1 + nelem_b2
    #     stop = 2:nelem_b1 + nelem_b2 + 1
    #     lengths = vcat(lengths_b1, lengths_b2)
    #     midpoints = vcat(xm_b1, xm_b2)
    #     Cab = vcat(Cab_b1, Cab_b2)

    #     compliance = fill(Diagonal([1/(E*A), 1/(G*Ay), 1/(G*Az), 1/(G*Jx), 1/(E*Iyy), 1/(E*Izz)]), nelem)

    #     mass = fill(Diagonal([ρ*A, ρ*A, ρ*A, ρ*J, ρ*Iyy, ρ*Izz]), nelem)

    #     # create assembly
    #     assembly = Assembly(points, start, stop, compliance=compliance, mass=mass, frames=Cab, lengths=lengths, midpoints=midpoints)

    #     # create system
    #     system = ExpandedSystem(assembly)

    #     for j = 1:length(rpm)
    #         # global frame rotation
    #         w0 = [0, 0, rpm[j]*(2*pi)/60]

    #         # eigenvalues and (right) eigenvectors
    #         system, λ[i,j], V, converged = eigenvalue_analysis!(system, assembly,
    #             angular_velocity = w0,
    #             prescribed_conditions = prescribed_conditions,
    #             constant_mass_matrix = true,
    #             nev=nev)
                
    #         # corresponding left eigenvectors
    #         U[i,j] = left_eigenvectors(system, λ[i,j], V)

    #         # post-multiply mass matrix with right eigenvector matrix
    #         # (we use this later for correlating eigenvalues)
    #         MV[i,j] = system.M * V

    #         # process state and eigenstates
    #         state[i,j] = AssemblyState(system, assembly; prescribed_conditions=prescribed_conditions)
    #         eigenstates[i,j] = [AssemblyState(system, assembly, V[:,k];
    #             prescribed_conditions=prescribed_conditions) for k = 1:nev]
    #     end
    # end

    # # set previous left eigenvector matrix
    # U_p = copy(U[1,1])

    # for j = 1:length(rpm)
    #     for i = 1:length(sweep)
    #         # construct correlation matrix
    #         C = U_p*MV[i,j]

    #         # correlate eigenmodes
    #         perm, corruption = correlate_eigenmodes(C)

    #         # re-arrange eigenvalues and eigenvectors
    #         λ[i,j] = λ[i,j][perm]
    #         U[i,j] = U[i,j][perm,:]
    #         MV[i,j] = MV[i,j][:,perm]
    #         eigenstates[i,j] = eigenstates[i,j][perm]

    #         # update previous eigenvector matrix
    #         U_p .= U[i,j]
    #     end
    #     # update previous eigenvector matrix
    #     U_p .= U[1,j]
    # end

    # frequency = [[imag(λ[i,j][k])/(2*pi) for i = 1:length(sweep), j=1:length(rpm)] for k = 1:2:nev]

    # indices = [1, 2, 4]
    # experiment_rpm = [0, 500, 750]
    # experiment_sweep = [0, 15, 30, 45]
    # experiment_frequencies = [
    #     [1.4 1.8 1.7 1.6;
    #      10.2 10.1 10.2 10.2;
    #      14.8 14.4 14.9 14.7],
    #     [10.3 10.2 10.4 10.4;
    #      25.2 25.2 23.7 21.6;
    #      36.1 34.8 30.7 26.1],
    #     [27.7 27.2 26.6 24.8;
    #      47.0 44.4 39.3 35.1;
    #      62.9 55.9 48.6 44.8]
    # ]

    # for k = 1:length(experiment_frequencies)
    #     for j = 1:length(experiment_sweep)
    #         for i = 1:length(experiment_rpm)
    #             ii = argmin(abs.(rpm .- experiment_rpm[i]))
    #             jj = argmin(abs.(sweep*180/pi .- experiment_sweep[j]))
    #             kk = indices[k]
    #             @test isapprox(frequency[kk][jj,ii], experiment_frequencies[k][i,j], atol=1, rtol=0.1)
    #         end
    #     end
    # end

    # indices = [5, 7, 6]
    # experiment_frequencies = [
    #     95.4 87.5 83.7 78.8;
    #     106.6 120.1 122.6 117.7;
    #     132.7 147.3 166.2 162.0
    # ]

    # for k = 1:size(experiment_frequencies, 1)
    #     for j = 1:length(experiment_sweep)
    #         ii = argmin(abs.(rpm .- 750))
    #         jj = argmin(abs.(sweep*180/pi .- experiment_sweep[j]))
    #         kk = indices[k]
    #         @test isapprox(frequency[kk][jj,ii], experiment_frequencies[k,j], rtol=0.1)
    #     end
    # end
end