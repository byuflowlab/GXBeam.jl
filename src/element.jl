# --- Element State Variable Helper Functions --- #

"""
    element_loads(x, ielem, icol_elem, force_scaling)

Extract the internal loads (`F`, `M`) for a beam element from the state variable vector.  
These loads are expressed in the deformed element frame.
"""
@inline function element_loads(x, ielem, icol_elem, force_scaling)

    icol = icol_elem[ielem]

    F = SVector(x[icol  ], x[icol+1], x[icol+2]) * force_scaling
    M = SVector(x[icol+3], x[icol+4], x[icol+5]) * force_scaling

    return F, M
end

"""
    expanded_element_loads(x, ielem, icol_elem, force_scaling)

Extract the internal loads of a beam element at its endpoints (`F1`, `M1`, `F2`, `M2`) 
from the state variable vector for a constant mass matrix system.
"""
@inline function expanded_element_loads(x, ielem, icol_elem, force_scaling)

    icol = icol_elem[ielem]

    F1 = SVector(x[icol  ], x[icol+1], x[icol+2]) * force_scaling
    M1 = SVector(x[icol+3], x[icol+4], x[icol+5]) * force_scaling
    F2 = SVector(x[icol+6], x[icol+7], x[icol+8]) * force_scaling
    M2 = SVector(x[icol+9], x[icol+10], x[icol+11]) * force_scaling

    return F1, M1, F2, M2
end

"""
    expanded_element_velocities(x, ielem, icol_elem)

Extract the velocities of a beam element from the state variable vector for a constant mass
matrix system.
"""
@inline function expanded_element_velocities(x, ielem, icol_elem)

    icol = icol_elem[ielem]

    V = SVector(x[icol+12], x[icol+13], x[icol+14])
    Ω = SVector(x[icol+15], x[icol+16], x[icol+17])

    return V, Ω
end

# --- Element Properties --- #

"""
    static_element_properties(x, indices, force_scaling, assembly, ielem, 
        prescribed_conditions, gravity)

Calculate/extract the element properties needed to construct the residual for a static 
analysis
"""
@inline function static_element_properties(x, indices, force_scaling, assembly, ielem, 
    prescribed_conditions, gravity)

    # unpack element parameters
    @unpack L, Cab, compliance, mass = assembly.elements[ielem]

    # scale compliance and mass matrices by the element length
    compliance *= L
    mass *= L

    # compliance submatrices (in the deformed element frame)
    S11 = compliance[SVector{3}(1:3), SVector{3}(1:3)]
    S12 = compliance[SVector{3}(1:3), SVector{3}(4:6)]
    S21 = compliance[SVector{3}(4:6), SVector{3}(1:3)]
    S22 = compliance[SVector{3}(4:6), SVector{3}(4:6)]

    # mass submatrices
    mass11 = mass[SVector{3}(1:3), SVector{3}(1:3)]
    mass12 = mass[SVector{3}(1:3), SVector{3}(4:6)]
    mass21 = mass[SVector{3}(4:6), SVector{3}(1:3)]
    mass22 = mass[SVector{3}(4:6), SVector{3}(4:6)]

    # linear and angular displacement
    u1, θ1 = point_displacement(x, assembly.start[ielem], indices.icol_point, prescribed_conditions)
    u2, θ2 = point_displacement(x, assembly.stop[ielem], indices.icol_point, prescribed_conditions)
    u = (u1 + u2)/2
    θ = (θ1 + θ2)/2

    # rotation parameter matrices
    C = get_C(θ)
    CtCab = C'*Cab
    Qinv = get_Qinv(θ)

    # forces and moments
    F, M = element_loads(x, ielem, indices.icol_elem, force_scaling)
    
    # strain and curvature
    γ = S11*F + S12*M
    κ = S21*F + S22*M

    # gravitational loads
    gvec = SVector{3}(gravity)

    return (; L, C, Cab, CtCab, Qinv, S11, S12, S21, S22, mass11, mass12, mass21, mass22, 
        u1, u2, θ1, θ2, u, θ, F, M, γ, κ, gvec)
end

"""
    steady_element_properties(x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, gravity, linear_velocity, angular_velocity, 
        linear_acceleration=(@SVector zeros(3)), angular_acceleration=(@SVector zeros(3)))

Calculate/extract the element properties needed to construct the residual for a steady 
state analysis
"""
@inline function steady_element_properties(x, indices, force_scaling, structural_damping, 
    assembly, ielem, prescribed_conditions, gravity, linear_velocity, angular_velocity, 
    linear_acceleration=(@SVector zeros(3)), angular_acceleration=(@SVector zeros(3)))

    properties = static_element_properties(x, indices, force_scaling, assembly, ielem, 
        prescribed_conditions, gravity)

    @unpack L, Cab, C, CtCab, Qinv, mass11, mass12, mass21, mass22, u1, u2, θ1, θ2, 
        u, θ, γ, κ = properties

    # rotation parameter matrices
    Q = get_Q(θ)

    # distance from the rotation center
    Δx = assembly.elements[ielem].x

    # body frame velocity (use prescribed values)
    vb, ωb = SVector{3}(linear_velocity), SVector{3}(angular_velocity)

    # body frame acceleration (use prescribed values)
    ab, αb = SVector{3}(linear_acceleration), SVector{3}(angular_acceleration)

    # linear and angular velocity
    V1, Ω1 = point_velocities(x, assembly.start[ielem], indices.icol_point)
    V2, Ω2 = point_velocities(x, assembly.stop[ielem], indices.icol_point)
    V = (V1 + V2)/2
    Ω = (Ω1 + Ω2)/2

    # linear and angular momentum
    P = CtCab*mass11*CtCab'*V + CtCab*mass12*CtCab'*Ω
    H = CtCab*mass21*CtCab'*V + CtCab*mass22*CtCab'*Ω

    # linear and angular displacement rates
    udot = @SVector zeros(3)
    θdot = @SVector zeros(3)

    # linear and angular acceleration
    Vdot = ab + cross(αb, Δx + u) + cross(ωb, udot)
    Ωdot = αb
   
    # linear and angular momentum rates
    Pdot = CtCab*mass11*CtCab'*Vdot + CtCab*mass12*CtCab'*Ωdot
    Hdot = CtCab*mass21*CtCab'*Vdot + CtCab*mass22*CtCab'*Ωdot

    # save properties
    properties = (; properties..., Q, Δx, vb, ωb, ab, αb, V1, V2, Ω1, Ω2, V, Ω, P, H, 
        udot, θdot, Vdot, Ωdot, Pdot, Hdot) 

    if structural_damping

        # damping coefficients
        μ = assembly.elements[ielem].mu

        # damping submatrices
        μ11 = @SMatrix [μ[1] 0 0; 0 μ[2] 0; 0 0 μ[3]]
        μ22 = @SMatrix [μ[4] 0 0; 0 μ[5] 0; 0 0 μ[6]]

        # rotation parameter matrices
        C1 = get_C(θ1)
        C2 = get_C(θ2)
        Qinv1 = get_Qinv(θ1)
        Qinv2 = get_Qinv(θ2)

        # distance from the reference location
        Δx1 = assembly.points[assembly.start[ielem]]
        Δx2 = assembly.points[assembly.stop[ielem]]

        # linear displacement rates 
        udot1 = V1 - vb - cross(ωb, Δx1 + u1)
        udot2 = V2 - vb - cross(ωb, Δx2 + u2)
        uedot = (udot1 + udot2)/2

        # angular displacement rates
        θdot1 = Qinv1*C1*(Ω1 - ωb)
        θdot2 = Qinv2*C2*(Ω2 - ωb)
        θedot = (θdot1 + θdot2)/2

        # change in linear and angular displacement
        Δu = u2 - u1
        Δθ = θ2 - θ1

        # change in linear and angular displacement rates
        Δudot = udot2 - udot1
        Δθdot = θdot2 - θdot1

        # ΔQ matrix (see structural damping theory)
        ΔQ = get_ΔQ(θ, Δθ, Q)

        # strain rates
        γdot = -CtCab'*tilde(Ω - ωb)*Δu + CtCab'*Δudot - L*CtCab'*tilde(Ω - ωb)*Cab*e1
        κdot = Cab'*Q*Δθdot + Cab'*ΔQ*θedot

        # adjust strains to account for strain rates
        γ -= μ11*γdot
        κ -= μ22*κdot

        # save new strains and structural damping properties
        properties = (; properties..., γ, κ, γdot, κdot,
            μ11, μ22, C1, C2, Qinv1, Qinv2, Δx1, Δx2, 
            udot1, udot2, uedot, θdot1, θdot2, θedot, Δu, Δθ, Δudot, Δθdot, ΔQ)
    end

    return properties
end

"""
    initial_element_properties(x, indices, rate_vars, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity, 
        linear_velocity, angular_velocity, linear_acceleration, angular_acceleration, u0, θ0, V0, Ω0, Vdot0, Ωdot0)

Calculate/extract the element properties needed to construct the residual for a time-domain
analysis initialization
"""
@inline function initial_element_properties(x, indices, rate_vars, 
    force_scaling, structural_damping, assembly, ielem, prescribed_conditions, gravity, 
    linear_velocity, angular_velocity, linear_acceleration, angular_acceleration, u0, θ0, V0, Ω0, Vdot0, Ωdot0)

    # unpack element parameters
    @unpack L, Cab, compliance, mass = assembly.elements[ielem]

    # scale compliance and mass matrices by the element length
    compliance *= L
    mass *= L

    # compliance submatrices
    S11 = compliance[SVector{3}(1:3), SVector{3}(1:3)]
    S12 = compliance[SVector{3}(1:3), SVector{3}(4:6)]
    S21 = compliance[SVector{3}(4:6), SVector{3}(1:3)]
    S22 = compliance[SVector{3}(4:6), SVector{3}(4:6)]

    # mass submatrices
    mass11 = mass[SVector{3}(1:3), SVector{3}(1:3)]
    mass12 = mass[SVector{3}(1:3), SVector{3}(4:6)]
    mass21 = mass[SVector{3}(4:6), SVector{3}(1:3)]
    mass22 = mass[SVector{3}(4:6), SVector{3}(4:6)]

    # linear and angular displacement
    u1, θ1 = initial_point_displacement(x, assembly.start[ielem], indices.icol_point, 
        prescribed_conditions, u0, θ0, rate_vars)
    u2, θ2 = initial_point_displacement(x, assembly.stop[ielem], indices.icol_point, 
        prescribed_conditions, u0, θ0, rate_vars)
    u = (u1 + u2)/2
    θ = (θ1 + θ2)/2

    # rotation parameter matrices
    C = get_C(θ)
    CtCab = C'*Cab
    Q = get_Q(θ)
    Qinv = get_Qinv(θ)

    # forces and moments
    F, M = element_loads(x, ielem, indices.icol_elem, force_scaling)

    # strain and curvature
    γ = S11*F + S12*M
    κ = S21*F + S22*M

    # distance from the rotation center
    Δx = assembly.elements[ielem].x
    
    # body frame velocity (use prescribed values)
    vb, ωb = SVector{3}(linear_velocity), SVector{3}(angular_velocity)

    # body frame acceleration (use prescribed values)
    ab, αb = SVector{3}(linear_acceleration), SVector{3}(angular_acceleration)

    # gravitational loads
    gvec = SVector{3}(gravity)

    # relative velocity
    V1 = SVector{3}(V0[assembly.start[ielem]])
    V2 = SVector{3}(V0[assembly.stop[ielem]])
    V = (V1 + V2)/2

    Ω1 = SVector{3}(Ω0[assembly.start[ielem]])
    Ω2 = SVector{3}(Ω0[assembly.stop[ielem]])
    Ω = (Ω1 + Ω2)/2

    # inertial velocity
    V += vb + cross(ωb, Δx + u)
    Ω += ωb

    # linear and angular momentum
    P = CtCab*mass11*CtCab'*V + CtCab*mass12*CtCab'*Ω
    H = CtCab*mass21*CtCab'*V + CtCab*mass22*CtCab'*Ω

    # linear and angular displacement rates
    u1dot, θ1dot = initial_point_displacement_rates(x, assembly.start[ielem], indices.icol_point)
    u2dot, θ2dot = initial_point_displacement_rates(x, assembly.stop[ielem], indices.icol_point)
    udot = (u1dot + u2dot)/2
    θdot = (θ1dot + θ2dot)/2

    # relative acceleration
    V1dot, Ω1dot = initial_point_velocity_rates(x, assembly.start[ielem], indices.icol_point, 
        prescribed_conditions, Vdot0, Ωdot0, rate_vars)
    V2dot, Ω2dot = initial_point_velocity_rates(x, assembly.stop[ielem], indices.icol_point, 
        prescribed_conditions, Vdot0, Ωdot0, rate_vars)
    Vdot = (V1dot + V2dot)/2
    Ωdot = (Ω1dot + Ω2dot)/2

    # inertial acceleration (excluding frame-rotation term)
    Vdot += ab + cross(αb, Δx) + cross(αb, u) + cross(ωb, udot)
    Ωdot += αb

    # linear and angular momentum rates
    CtCabdot = tilde(Ω - ωb)*CtCab

    Pdot = CtCab*mass11*CtCab'*Vdot + CtCab*mass12*CtCab'*Ωdot +
        CtCab*mass11*CtCabdot'*V + CtCab*mass12*CtCabdot'*Ω +
        CtCabdot*mass11*CtCab'*V + CtCabdot*mass12*CtCab'*Ω
    
    Hdot = CtCab*mass21*CtCab'*Vdot + CtCab*mass22*CtCab'*Ωdot +
        CtCab*mass21*CtCabdot'*V + CtCab*mass22*CtCabdot'*Ω +
        CtCabdot*mass21*CtCab'*V + CtCabdot*mass22*CtCab'*Ω

    # save properties
    properties = (; L, C, Cab, CtCab, Q, Qinv, S11, S12, S21, S22, mass11, mass12, mass21, mass22, 
        u1, u2, θ1, θ2, u, θ, F, M, γ, κ, gvec, Δx, vb, ωb, ab, αb, V1, V2, Ω1, Ω2, V, Ω, P, H, 
        u1dot, u2dot, θ1dot, θ2dot, udot, θdot, CtCabdot, Vdot, Ωdot, Pdot, Hdot) 

    if structural_damping 
        
        # damping coefficients
        μ = assembly.elements[ielem].mu

        # damping submatrices
        μ11 = @SMatrix [μ[1] 0 0; 0 μ[2] 0; 0 0 μ[3]]
        μ22 = @SMatrix [μ[4] 0 0; 0 μ[5] 0; 0 0 μ[6]]
    
        # change in linear and angular displacement
        Δu = u2 - u1
        Δθ = θ2 - θ1

        # change in linear and angular displacement rates
        Δudot = u2dot - u1dot
        Δθdot = θ2dot - θ1dot

        # ΔQ matrix (see structural damping theory)
        ΔQ = get_ΔQ(θ, Δθ, Q)

        # strain rates
        γdot = -CtCab'*tilde(Ω - ωb)*Δu + CtCab'*Δudot - L*CtCab'*tilde(Ω - ωb)*Cab*e1
        κdot = Cab'*Q*Δθdot + Cab'*ΔQ*θdot

        # adjust strains to account for strain rates
        γ -= μ11*γdot
        κ -= μ22*κdot

        # add structural damping properties
        properties = (; properties..., γ, κ, γdot, κdot,
            μ11, μ22, Δu, Δθ, Δudot, Δθdot, ΔQ)
    end

    return properties
end

"""
    newmark_element_properties(x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, gravity, linear_velocity, angular_velocity,
        Vdot_init, Ωdot_init, dt)

Calculate/extract the element properties needed to construct the residual for a newmark-
scheme time stepping analysis
"""
@inline function newmark_element_properties(x, indices, force_scaling, structural_damping, 
    assembly, ielem, prescribed_conditions, gravity, linear_velocity, angular_velocity, 
    Vdot_init, Ωdot_init, dt)

    properties = steady_element_properties(x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, gravity, linear_velocity, angular_velocity)

    @unpack ωb, C, Cab, CtCab, mass11, mass12, mass21, mass22, V1, V2, Ω1, Ω2, V, Ω = properties

    # linear and angular acceleration (including body frame motion)
    V1dot = 2/dt*V1 - SVector{3}(Vdot_init[assembly.start[ielem]])
    Ω1dot = 2/dt*Ω1 - SVector{3}(Ωdot_init[assembly.start[ielem]])

    V2dot = 2/dt*V2 - SVector{3}(Vdot_init[assembly.stop[ielem]])
    Ω2dot = 2/dt*Ω2 - SVector{3}(Ωdot_init[assembly.stop[ielem]])

    Vdot = (V1dot + V2dot)/2
    Ωdot = (Ω1dot + Ω2dot)/2

    # linear and angular momentum rates (including body frame motion)
    CtCabdot = tilde(Ω - ωb)*CtCab
    
    Pdot = CtCab*mass11*CtCab'*Vdot + CtCab*mass12*CtCab'*Ωdot +
        CtCab*mass11*CtCabdot'*V + CtCab*mass12*CtCabdot'*Ω +
        CtCabdot*mass11*CtCab'*V + CtCabdot*mass12*CtCab'*Ω
    
    Hdot = CtCab*mass21*CtCab'*Vdot + CtCab*mass22*CtCab'*Ωdot +
        CtCab*mass21*CtCabdot'*V + CtCab*mass22*CtCabdot'*Ω +
        CtCabdot*mass21*CtCab'*V + CtCabdot*mass22*CtCab'*Ω
   
    return (; properties..., CtCabdot, Vdot, Ωdot, Pdot, Hdot) 
end

"""
    dynamic_element_properties(dx, x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, gravity, linear_velocity, angular_velocity)

Calculate/extract the element properties needed to construct the residual for a dynamic 
analysis
"""
@inline function dynamic_element_properties(dx, x, indices, force_scaling, structural_damping, 
    assembly, ielem, prescribed_conditions, gravity, linear_velocity, angular_velocity)

    properties = steady_element_properties(x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, gravity, linear_velocity, angular_velocity)

    @unpack ωb, C, Cab, CtCab, mass11, mass12, mass21, mass22, V1, V2, Ω1, Ω2, V, Ω = properties

    # velocity rates
    V1dot, Ω1dot = point_velocities(dx, assembly.start[ielem], indices.icol_point)
    V2dot, Ω2dot = point_velocities(dx, assembly.stop[ielem], indices.icol_point)
    Vdot = (V1dot + V2dot)/2
    Ωdot = (Ω1dot + Ω2dot)/2

    # linear and angular momentum rates
    CtCabdot = tilde(Ω - ωb)*CtCab
    
    Pdot = CtCab*mass11*CtCab'*Vdot + CtCab*mass12*CtCab'*Ωdot +
        CtCab*mass11*CtCabdot'*V + CtCab*mass12*CtCabdot'*Ω +
        CtCabdot*mass11*CtCab'*V + CtCabdot*mass12*CtCab'*Ω
    
    Hdot = CtCab*mass21*CtCab'*Vdot + CtCab*mass22*CtCab'*Ωdot +
        CtCab*mass21*CtCabdot'*V + CtCab*mass22*CtCabdot'*Ω +
        CtCabdot*mass21*CtCab'*V + CtCabdot*mass22*CtCab'*Ω

    return (; properties..., CtCabdot, Vdot, Ωdot, Pdot, Hdot)
end

"""
    expanded_steady_element_properties(x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, gravity, linear_velocity, angular_velocity, 
        linear_acceleration=(@SVector zeros(3)), angular_acceleration=(@SVector zeros(3)))

Calculate/extract the element properties needed to construct the residual for a constant
mass matrix system
"""
@inline function expanded_steady_element_properties(x, indices, force_scaling, structural_damping, 
    assembly, ielem, prescribed_conditions, gravity, linear_velocity, angular_velocity, 
    linear_acceleration=(@SVector zeros(3)), angular_acceleration=(@SVector zeros(3)))

    # unpack element parameters
    @unpack L, Cab, compliance, mass = assembly.elements[ielem]

    # scale compliance and mass matrices by the element length
    compliance *= L
    mass *= L

    # compliance submatrices (in the deformed element frame)
    S11 = compliance[SVector{3}(1:3), SVector{3}(1:3)]
    S12 = compliance[SVector{3}(1:3), SVector{3}(4:6)]
    S21 = compliance[SVector{3}(4:6), SVector{3}(1:3)]
    S22 = compliance[SVector{3}(4:6), SVector{3}(4:6)]

    # mass submatrices
    mass11 = mass[SVector{3}(1:3), SVector{3}(1:3)]
    mass12 = mass[SVector{3}(1:3), SVector{3}(4:6)]
    mass21 = mass[SVector{3}(4:6), SVector{3}(1:3)]
    mass22 = mass[SVector{3}(4:6), SVector{3}(4:6)]

    # linear and angular displacement
    u1, θ1 = point_displacement(x, assembly.start[ielem], indices.icol_point, prescribed_conditions)
    u2, θ2 = point_displacement(x, assembly.stop[ielem], indices.icol_point, prescribed_conditions)
    u = (u1 + u2)/2
    θ = (θ1 + θ2)/2

    # rotation parameter matrices
    C = get_C(θ)
    C1 = get_C(θ1)
    C2 = get_C(θ2)
    CtCab = C'*Cab
    Q = get_Q(θ)
    Qinv = get_Qinv(θ)
    Qinv1 = get_Qinv(θ1)
    Qinv2 = get_Qinv(θ2)

    # forces and moments
    F1, M1, F2, M2 = expanded_element_loads(x, ielem, indices.icol_elem, force_scaling)
    F = (F1 + F2)/2
    M = (M1 + M2)/2

    # strain and curvature
    γ = S11*F + S12*M
    κ = S21*F + S22*M

    # gravitational loads
    gvec = SVector{3}(gravity)

    # distance from the rotation center
    Δx = assembly.elements[ielem].x
    Δx1 = assembly.points[assembly.start[ielem]]
    Δx2 = assembly.points[assembly.stop[ielem]]

    # body frame linear velocity
    vb = linear_velocity
    ωb = angular_velocity

    # body frame angular acceleration
    ab = linear_acceleration
    αb = angular_acceleration

    # linear and angular velocity
    V1, Ω1 = point_velocities(x, assembly.start[ielem], indices.icol_point)
    V2, Ω2 = point_velocities(x, assembly.stop[ielem], indices.icol_point)
    V, Ω = expanded_element_velocities(x, ielem, indices.icol_elem)

    # linear and angular momentum
    P = mass11*V + mass12*Ω
    H = mass21*V + mass22*Ω

    # linear and angular displacement rates
    udot1 = C1'*V1 - vb - cross(ωb, Δx1 + u1)
    udot2 = C2'*V2 - vb - cross(ωb, Δx2 + u2)
    udot = (udot1 + udot2)/2

    θdot1 = Qinv1*(Ω1 - C1*ωb)
    θdot2 = Qinv2*(Ω2 - C2*ωb)
    θdot = (θdot1 + θdot2)/2

    # linear and angular velocity rates
    Vdot = CtCab'*(ab + cross(αb, Δx) + cross(αb, u))
    Ωdot = CtCab'*αb
   
    # linear and angular momentum rates
    Pdot = mass11*Vdot + mass12*Ωdot
    Hdot = mass21*Vdot + mass22*Ωdot

    # save properties
    properties = (; L, C, C1, C2, Cab, CtCab, Q, Qinv, Qinv1, Qinv2, S11, S12, S21, S22, 
        mass11, mass12, mass21, mass22, u1, u2, θ1, θ2, u, θ, F1, F2, M1, M2, F, M, γ, κ, gvec, 
        Δx, Δx1, Δx2, vb, ωb, ab, αb, V1, V2, Ω1, Ω2, V, Ω, P, H, udot1, udot2, θdot1, θdot2, 
        udot, θdot, Vdot, Ωdot, Pdot, Hdot)

    if structural_damping

        # damping coefficients
        μ = assembly.elements[ielem].mu

        # damping submatrices
        μ11 = @SMatrix [μ[1] 0 0; 0 μ[2] 0; 0 0 μ[3]]
        μ22 = @SMatrix [μ[4] 0 0; 0 μ[5] 0; 0 0 μ[6]]

        # change in linear and angular displacement
        Δu = u2 - u1
        Δθ = θ2 - θ1

        # change in linear and angular displacement rates
        Δudot = udot2 - udot1
        Δθdot = θdot2 - θdot1

        # ΔQ matrix (see structural damping theory)
        ΔQ = get_ΔQ(θ, Δθ, Q)

        # strain rates
        γdot = -CtCab'*tilde(CtCab*Ω - ωb)*Δu + CtCab'*Δudot - L*CtCab'*tilde(CtCab*Ω - ωb)*Cab*e1
        κdot = Cab'*Q*Δθdot + Cab'*ΔQ*θdot

        # adjust strains to account for strain rates
        γ -= μ11*γdot
        κ -= μ22*κdot

        # add structural damping properties
        properties = (; properties..., γ, κ, γdot, κdot, μ11, μ22, Δu, Δθ, Δudot, Δθdot, ΔQ)
    end

    return properties
end

"""
    expanded_dynamic_element_properties(dx, x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, gravity, linear_velocity, angular_velocity)

Calculate/extract the element properties needed to construct the residual for a constant
mass matrix system
"""
function expanded_dynamic_element_properties(dx, x, indices, force_scaling, structural_damping, 
    assembly, ielem, prescribed_conditions, gravity, linear_velocity, angular_velocity)

    properties = expanded_steady_element_properties(x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity, 
        linear_velocity, angular_velocity)

    @unpack mass11, mass12, mass21, mass22 = properties

    # velocity rates
    Vdot, Ωdot = expanded_element_velocities(dx, ielem, indices.icol_elem)

    # linear and angular momentum rates
    Pdot = mass11*Vdot + mass12*Ωdot
    Hdot = mass21*Vdot + mass22*Ωdot

    return (; properties..., Vdot, Ωdot, Pdot, Hdot)
end

"""
    static_element_jacobian_properties(properties, x, indices, force_scaling, 
        assembly, ielem, prescribed_conditions, gravity)

Calculate/extract the element properties needed to calculate the jacobian entries 
corresponding to a beam element for a static analysis
"""
@inline function static_element_jacobian_properties(properties, x, indices, force_scaling, 
    assembly, ielem, prescribed_conditions, gravity)

    @unpack C, θ, S11, S12, S21, S22 = properties

    # linear and angular displacement
    u1_u1, θ1_θ1 = point_displacement_jacobians(assembly.start[ielem], prescribed_conditions)
    u2_u2, θ2_θ2 = point_displacement_jacobians(assembly.stop[ielem], prescribed_conditions)

    # rotation parameter matrices
    C_θ1, C_θ2, C_θ3 = get_C_θ(C, θ)
    Qinv_θ1, Qinv_θ2, Qinv_θ3 = get_Qinv_θ(θ)

    # strain and curvature
    γ_F, γ_M, κ_F, κ_M = S11, S12, S21, S22

    return (; properties..., u1_u1, u2_u2, θ1_θ1, θ2_θ2, 
        C_θ1, C_θ2, C_θ3, Qinv_θ1, Qinv_θ2, Qinv_θ3, γ_F, γ_M, κ_F, κ_M)
end

"""
    steady_element_jacobian_properties(properties, x, indices, 
        force_scaling, structural_damping, assembly, ielem, prescribed_conditions, 
        gravity)

Calculate/extract the element properties needed to calculate the jacobian entries 
corresponding to a element for a steady state analysis
"""
@inline function steady_element_jacobian_properties(properties, x, indices, 
    force_scaling, structural_damping, assembly, ielem, prescribed_conditions, 
    gravity)

    properties = static_element_jacobian_properties(properties, x, indices, force_scaling, 
        assembly, ielem, prescribed_conditions, gravity)

    @unpack Δx, αb, L, mass11, mass12, mass21, mass22, C, Cab, CtCab, Q, u, θ, V, Ω, Vdot, Ωdot, 
        C_θ1, C_θ2, C_θ3 = properties

    # rotation parameter matrices
    Q_θ1, Q_θ2, Q_θ3 = get_Q_θ(Q, θ)

    # linear and angular momentum
    P_θ = mul3(C_θ1', C_θ2', C_θ3', Cab*(mass11*CtCab'*V + mass12*CtCab'*Ω)) + 
        CtCab*mass11*Cab'*mul3(C_θ1, C_θ2, C_θ3, V) + 
        CtCab*mass12*Cab'*mul3(C_θ1, C_θ2, C_θ3, Ω)
    P_V = CtCab*mass11*CtCab'
    P_Ω = CtCab*mass12*CtCab'

    H_θ = mul3(C_θ1', C_θ2', C_θ3', Cab*(mass21*CtCab'*V + mass22*CtCab'*Ω)) + 
        CtCab*mass21*Cab'*mul3(C_θ1, C_θ2, C_θ3, V) + 
        CtCab*mass22*Cab'*mul3(C_θ1, C_θ2, C_θ3, Ω)
    H_V = CtCab*mass21*CtCab'
    H_Ω = CtCab*mass22*CtCab'

    # linear and angular acceleration
    Vdot_ab = I3
    Vdot_αb = -tilde(Δx + u)
    Vdot_u = tilde(αb)

    Ωdot_αb = I3

    # linear and angular momentum rates
    Pdot_Vdot = CtCab*mass11*CtCab'
    Pdot_Ωdot = CtCab*mass12*CtCab'
    Hdot_Vdot = CtCab*mass21*CtCab'
    Hdot_Ωdot = CtCab*mass22*CtCab'

    Pdot_ab = Pdot_Vdot*Vdot_ab
    Pdot_αb = Pdot_Vdot*Vdot_αb + Pdot_Ωdot*Ωdot_αb
    Hdot_ab = Hdot_Vdot*Vdot_ab
    Hdot_αb = Hdot_Vdot*Vdot_αb + Hdot_Ωdot*Ωdot_αb

    Pdot_u = Pdot_Vdot*Vdot_u
    Hdot_u = Hdot_Vdot*Vdot_u

    Pdot_θ = mul3(C_θ1', C_θ2', C_θ3', Cab*mass11*CtCab'*Vdot) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass12*CtCab'*Ωdot) +
        CtCab*mass11*Cab'*mul3(C_θ1, C_θ2, C_θ3, Vdot) + 
        CtCab*mass12*Cab'*mul3(C_θ1, C_θ2, C_θ3, Ωdot)
    Hdot_θ = mul3(C_θ1', C_θ2', C_θ3', Cab*mass21*CtCab'*Vdot) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass22*CtCab'*Ωdot) +
        CtCab*mass21*Cab'*mul3(C_θ1, C_θ2, C_θ3, Vdot) + 
        CtCab*mass22*Cab'*mul3(C_θ1, C_θ2, C_θ3, Ωdot)
    
    Pdot_V = @SMatrix zeros(3,3)
    Pdot_Ω = @SMatrix zeros(3,3)

    Hdot_V = @SMatrix zeros(3,3)
    Hdot_Ω = @SMatrix zeros(3,3)

    if structural_damping

        @unpack ωb, C1, C2, Qinv1, Qinv2, u1, u2, θ1, θ2, Ω1, Ω2, μ11, μ22, 
            uedot, θedot, Δx1, Δx2, Δu, Δθ, ΔQ, Δudot, Δθdot = properties

        # rotation parameter matrices
        C1_θ1, C1_θ2, C1_θ3 = get_C_θ(C1, θ1)
        C2_θ1, C2_θ2, C2_θ3 = get_C_θ(C2, θ2)
        Qinv1_θ1, Qinv1_θ2, Qinv1_θ3 = get_Qinv_θ(θ1)
        Qinv2_θ1, Qinv2_θ2, Qinv2_θ3 = get_Qinv_θ(θ2)

        # linear displacement rates
        udot1_u1 = -tilde(ωb)
        udot2_u2 = -tilde(ωb)
        udot1_V1 = I3
        udot2_V2 = I3

        # angular displacement rates
        θdot1_θ1 = mul3(Qinv1_θ1, Qinv1_θ2, Qinv1_θ3, C1*(Ω1 - ωb)) + Qinv1*mul3(C1_θ1, C1_θ2, C1_θ3, Ω1 - ωb)
        θdot2_θ2 = mul3(Qinv2_θ1, Qinv2_θ2, Qinv2_θ3, C2*(Ω2 - ωb)) + Qinv2*mul3(C2_θ1, C2_θ2, C2_θ3, Ω2 - ωb)
        
        θdot1_Ω1 = Qinv1*C1
        θdot2_Ω2 = Qinv2*C2

        θdot_θ1 = θdot1_θ1/2
        θdot_θ2 = θdot2_θ2/2

        θdot_Ω1 = θdot1_Ω1/2
        θdot_Ω2 = θdot2_Ω2/2

        # change in linear and angular displacement
        Δu_u1 = -I3
        Δu_u2 =  I3
        Δθ_θ1 = -I3
        Δθ_θ2 =  I3

        # change in linear and angular displacement rates
        Δudot_u1 = -udot1_u1
        Δudot_u2 =  udot2_u2

        Δudot_V1 = -udot1_V1
        Δudot_V2 =  udot2_V2

        Δθdot_θ1 = -θdot1_θ1
        Δθdot_θ2 =  θdot2_θ2

        Δθdot_Ω1 = -θdot1_Ω1
        Δθdot_Ω2 =  θdot2_Ω2   

        # ΔQ matrix (see structural damping theory)
        ΔQ_θ1, ΔQ_θ2, ΔQ_θ3 = get_ΔQ_θ(θ, Δθ, Q, Q_θ1, Q_θ2, Q_θ3)

        ΔQ_Δθ1 = mul3(Q_θ1, Q_θ2, Q_θ3, e1)
        ΔQ_Δθ2 = mul3(Q_θ1, Q_θ2, Q_θ3, e2)
        ΔQ_Δθ3 = mul3(Q_θ1, Q_θ2, Q_θ3, e3)

        # strain rates
        tmp = CtCab'*tilde(Ω - ωb)
        γdot_u1 = -tmp*Δu_u1 + CtCab'*Δudot_u1 
        γdot_u2 = -tmp*Δu_u2 + CtCab'*Δudot_u2 

        tmp = -Cab'*mul3(C_θ1, C_θ2, C_θ3, tilde(Ω - ωb)*Δu) + 
            Cab'*mul3(C_θ1, C_θ2, C_θ3, Δudot) - 
            L*Cab'*mul3(C_θ1, C_θ2, C_θ3, tilde(Ω - ωb)*Cab*e1)
        γdot_θ1 = 1/2*tmp
        γdot_θ2 = 1/2*tmp

        γdot_V1 = CtCab'*Δudot_V1
        γdot_V2 = CtCab'*Δudot_V2

        tmp = CtCab'*tilde(Δu) + L*CtCab'*tilde(Cab*e1)
        γdot_Ω1 = 1/2*tmp
        γdot_Ω2 = 1/2*tmp

        tmp1 = Cab'*mul3(Q_θ1, Q_θ2, Q_θ3, Δθdot)
        tmp2 = Cab'*mul3(ΔQ_θ1, ΔQ_θ2, ΔQ_θ3, θedot) 
        tmp3 = Cab'*mul3(ΔQ_Δθ1, ΔQ_Δθ2, ΔQ_Δθ3, θedot)
        κdot_θ1 = 1/2*tmp1 + 1/2*tmp2 + tmp3*Δθ_θ1 + Cab'*Q*Δθdot_θ1 + Cab'*ΔQ*θdot_θ1
        κdot_θ2 = 1/2*tmp1 + 1/2*tmp2 + tmp3*Δθ_θ2 + Cab'*Q*Δθdot_θ2 + Cab'*ΔQ*θdot_θ2  

        κdot_Ω1 = Cab'*Q*Δθdot_Ω1 + Cab'*ΔQ*θdot_Ω1
        κdot_Ω2 = Cab'*Q*Δθdot_Ω2 + Cab'*ΔQ*θdot_Ω2

        # adjust strains to account for strain rates
        
        γ_u1 = -μ11*γdot_u1
        γ_u2 = -μ11*γdot_u2

        γ_θ1 = -μ11*γdot_θ1
        γ_θ2 = -μ11*γdot_θ2
        
        γ_V1 = -μ11*γdot_V1
        γ_V2 = -μ11*γdot_V2
        
        γ_Ω1 = -μ11*γdot_Ω1
        γ_Ω2 = -μ11*γdot_Ω2
        
        κ_θ1 = -μ22*κdot_θ1
        κ_θ2 = -μ22*κdot_θ2
        
        κ_Ω1 = -μ22*κdot_Ω1
        κ_Ω2 = -μ22*κdot_Ω2

    else
        
        γ_u1 = @SMatrix zeros(3,3)
        γ_u2 = @SMatrix zeros(3,3)

        γ_θ1 = @SMatrix zeros(3,3)
        γ_θ2 = @SMatrix zeros(3,3)
        
        γ_V1 = @SMatrix zeros(3,3)
        γ_V2 = @SMatrix zeros(3,3)
        
        γ_Ω1 = @SMatrix zeros(3,3)
        γ_Ω2 = @SMatrix zeros(3,3)
        
        κ_θ1 = @SMatrix zeros(3,3)
        κ_θ2 = @SMatrix zeros(3,3)
        
        κ_Ω1 = @SMatrix zeros(3,3)
        κ_Ω2 = @SMatrix zeros(3,3)
         
    end

    return (; properties..., γ_u1, γ_u2, γ_θ1, γ_θ2, γ_V1, γ_V2, γ_Ω1, γ_Ω2, 
        κ_θ1, κ_θ2, κ_Ω1, κ_Ω2, P_θ, P_V, P_Ω, H_θ, H_V, H_Ω,
        Pdot_ab, Pdot_αb, Pdot_u, Pdot_θ, Pdot_V, Pdot_Ω, 
        Hdot_ab, Hdot_αb, Hdot_u, Hdot_θ, Hdot_V, Hdot_Ω) 
end

"""
    initial_element_jacobian_properties(properties, x, indices, rate_vars, 
        force_scaling, structural_damping, assembly, ielem, prescribed_conditions, gravity, 
        u0, θ0, V0, Ω0, Vdot0, Ωdot0)

Calculate/extract the element properties needed to calculate the jacobian entries 
corresponding to a element for a Newmark scheme time marching analysis
"""
@inline function initial_element_jacobian_properties(properties, x, indices, rate_vars, 
    force_scaling, structural_damping, assembly, ielem, prescribed_conditions, gravity, 
    u0, θ0, V0, Ω0, Vdot0, Ωdot0)

    @unpack Δx, ωb, αb, L, C, Cab, CtCab, CtCabdot, Q, S11, S12, S21, S22, mass11, mass12, mass21, mass22, 
        u, θ, V, Ω, Vdot, Ωdot = properties

    # linear and angular displacement
    u1_u1, θ1_θ1 = initial_point_displacement_jacobian(assembly.start[ielem], indices.icol_point, 
        prescribed_conditions, rate_vars)

    u2_u2, θ2_θ2 = initial_point_displacement_jacobian(assembly.stop[ielem], indices.icol_point, 
        prescribed_conditions, rate_vars)

    # rotation parameter matrices
    C_θ1, C_θ2, C_θ3 = get_C_θ(C, θ)
    Q_θ1, Q_θ2, Q_θ3 = get_Q_θ(Q, θ)
    Qinv_θ1, Qinv_θ2, Qinv_θ3 = get_Qinv_θ(θ)

    # strain and curvature
    γ_F, γ_M, κ_F, κ_M = S11, S12, S21, S22

    # linear and angular velocity (including body frame motion)
    V_u = tilde(ωb)

    # linear and angular momentum
    P_u = CtCab*mass11*CtCab'*V_u
    P_θ = mul3(C_θ1', C_θ2', C_θ3', Cab*(mass11*CtCab'*V + mass12*CtCab'*Ω)) + 
        CtCab*mass11*Cab'*mul3(C_θ1, C_θ2, C_θ3, V) + 
        CtCab*mass12*Cab'*mul3(C_θ1, C_θ2, C_θ3, Ω)

    H_u = CtCab*mass21*CtCab'*V_u
    H_θ = mul3(C_θ1', C_θ2', C_θ3', Cab*(mass21*CtCab'*V + mass22*CtCab'*Ω)) + 
        CtCab*mass21*Cab'*mul3(C_θ1, C_θ2, C_θ3, V) + 
        CtCab*mass22*Cab'*mul3(C_θ1, C_θ2, C_θ3, Ω)

    # linear and angular acceleration
    V1dot_V1dot, Ω1dot_Ω1dot = initial_point_velocity_rate_jacobian(assembly.start[ielem], 
        indices.icol_point, prescribed_conditions, rate_vars)

    V2dot_V2dot, Ω2dot_Ω2dot = initial_point_velocity_rate_jacobian(assembly.stop[ielem], 
        indices.icol_point, prescribed_conditions, rate_vars)
 
    # linear and angular acceleration (including body frame motion)
    Vdot_ab = I3
    Vdot_αb = -tilde(Δx + u)
    Vdot_u = tilde(αb)
    Vdot_udot = tilde(ωb)

    Ωdot_αb = I3

    # linear and angular momentum rates
    Pdot_V = CtCab*mass11*CtCabdot' + CtCabdot*mass11*CtCab'
    Pdot_Ω = CtCab*mass12*CtCabdot' + CtCabdot*mass12*CtCab'
    Pdot_Vdot = CtCab*mass11*CtCab'
    Pdot_Ωdot = CtCab*mass12*CtCab'
    Pdot_ab = Pdot_Vdot*Vdot_ab
    Pdot_αb = Pdot_Vdot*Vdot_αb + Pdot_Ωdot*Ωdot_αb
    Pdot_u = Pdot_V*V_u + Pdot_Vdot*Vdot_u
    Pdot_θ = mul3(C_θ1', C_θ2', C_θ3', Cab*mass11*CtCab'*Vdot) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass12*CtCab'*Ωdot) +
        CtCab*mass11*Cab'*mul3(C_θ1, C_θ2, C_θ3, Vdot) + 
        CtCab*mass12*Cab'*mul3(C_θ1, C_θ2, C_θ3, Ωdot) +
        tilde(Ω - ωb)*mul3(C_θ1', C_θ2', C_θ3', Cab*mass11*CtCab'*V) + 
        tilde(Ω - ωb)*mul3(C_θ1', C_θ2', C_θ3', Cab*mass12*CtCab'*Ω) +
        CtCabdot*mass11*Cab'*mul3(C_θ1, C_θ2, C_θ3, V) + 
        CtCabdot*mass12*Cab'*mul3(C_θ1, C_θ2, C_θ3, Ω) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass11*CtCabdot'*V) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass12*CtCabdot'*Ω) +
        -CtCab*mass11*Cab'*mul3(C_θ1, C_θ2, C_θ3, tilde(Ω - ωb)*V) + 
        -CtCab*mass12*Cab'*mul3(C_θ1, C_θ2, C_θ3, tilde(Ω - ωb)*Ω)
    Pdot_udot = Pdot_Vdot*Vdot_udot

    Hdot_V = CtCab*mass21*CtCabdot' + CtCabdot*mass21*CtCab'
    Hdot_Ω = CtCab*mass22*CtCabdot' + CtCabdot*mass22*CtCab'
    Hdot_Vdot = CtCab*mass21*CtCab'
    Hdot_Ωdot = CtCab*mass22*CtCab'
    Hdot_ab = Hdot_Vdot*Vdot_ab
    Hdot_αb = Hdot_Vdot*Vdot_αb + Hdot_Ωdot*Ωdot_αb
    Hdot_u = Hdot_V*V_u + Hdot_Vdot*Vdot_u
    Hdot_θ = mul3(C_θ1', C_θ2', C_θ3', Cab*mass21*CtCab'*Vdot) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass22*CtCab'*Ωdot) +
        CtCab*mass21*Cab'*mul3(C_θ1, C_θ2, C_θ3, Vdot) + 
        CtCab*mass22*Cab'*mul3(C_θ1, C_θ2, C_θ3, Ωdot) +
        tilde(Ω - ωb)*mul3(C_θ1', C_θ2', C_θ3', Cab*mass21*CtCab'*V) + 
        tilde(Ω - ωb)*mul3(C_θ1', C_θ2', C_θ3', Cab*mass22*CtCab'*Ω) +
        CtCabdot*mass21*Cab'*mul3(C_θ1, C_θ2, C_θ3, V) + 
        CtCabdot*mass22*Cab'*mul3(C_θ1, C_θ2, C_θ3, Ω) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass21*CtCabdot'*V) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass22*CtCabdot'*Ω) +
        -CtCab*mass21*Cab'*mul3(C_θ1, C_θ2, C_θ3, tilde(Ω - ωb)*V) + 
        -CtCab*mass22*Cab'*mul3(C_θ1, C_θ2, C_θ3, tilde(Ω - ωb)*Ω)
    Hdot_udot = Hdot_Vdot*Vdot_udot

    if structural_damping

        @unpack μ11, μ22, Δu, Δθ, ΔQ, udot, θdot, Δudot, Δθdot = properties

        # linear displacement rates
        u1dot_u1dot = I3
        u2dot_u2dot = I3

        # angular displacement rates
        θ1dot_θ1dot = I3
        θ2dot_θ2dot = I3
        
        θdot_θ1dot = 1/2*θ1dot_θ1dot
        θdot_θ2dot = 1/2*θ2dot_θ2dot

        # change in linear displacement
        Δu_u1 = -I3
        Δu_u2 =  I3

        # change in linear and angular displacement rates
        Δudot_u1dot = -u1dot_u1dot
        Δudot_u2dot =  u2dot_u2dot

        Δθdot_θ1dot = -θ1dot_θ1dot
        Δθdot_θ2dot =  θ2dot_θ2dot

        # ΔQ matrix (see structural damping theory)
        ΔQ_θ1, ΔQ_θ2, ΔQ_θ3 = get_ΔQ_θ(θ, Δθ, Q, Q_θ1, Q_θ2, Q_θ3)

        ΔQ_Δθ1 = mul3(Q_θ1, Q_θ2, Q_θ3, e1)
        ΔQ_Δθ2 = mul3(Q_θ1, Q_θ2, Q_θ3, e2)
        ΔQ_Δθ3 = mul3(Q_θ1, Q_θ2, Q_θ3, e3)

        # strain rates
        tmp = CtCab'*tilde(Ω - ωb)
        γdot_u1 = -tmp*Δu_u1
        γdot_u2 = -tmp*Δu_u2
        γdot_u1dot = CtCab'*Δudot_u1dot
        γdot_u2dot = CtCab'*Δudot_u2dot

        tmp = -Cab'*mul3(C_θ1, C_θ2, C_θ3, tilde(Ω - ωb)*Δu) + 
            Cab'*mul3(C_θ1, C_θ2, C_θ3, Δudot) - 
            L*Cab'*mul3(C_θ1, C_θ2, C_θ3, tilde(Ω - ωb)*Cab*e1)
        γdot_θ1 = 1/2*tmp
        γdot_θ2 = 1/2*tmp

        tmp1 = Cab'*mul3(Q_θ1, Q_θ2, Q_θ3, Δθdot)
        tmp2 = Cab'*mul3(ΔQ_θ1, ΔQ_θ2, ΔQ_θ3, θdot) 
        tmp3 = Cab'*mul3(ΔQ_Δθ1, ΔQ_Δθ2, ΔQ_Δθ3, θdot)
        κdot_θ1 = 1/2*tmp1 + 1/2*tmp2 - tmp3
        κdot_θ2 = 1/2*tmp1 + 1/2*tmp2 + tmp3  
        κdot_θ1dot = Cab'*Q*Δθdot_θ1dot + Cab'*ΔQ*θdot_θ1dot
        κdot_θ2dot = Cab'*Q*Δθdot_θ2dot + Cab'*ΔQ*θdot_θ2dot  

        # adjust strains to account for strain rates
        γ_u1 = -μ11*γdot_u1
        γ_u2 = -μ11*γdot_u2

        γ_θ1 = -μ11*γdot_θ1
        γ_θ2 = -μ11*γdot_θ2
        
        γ_u1dot = -μ11*γdot_u1dot
        γ_u2dot = -μ11*γdot_u2dot

        κ_θ1 = -μ22*κdot_θ1
        κ_θ2 = -μ22*κdot_θ2

        κ_θ1dot = -μ22*κdot_θ1dot
        κ_θ2dot = -μ22*κdot_θ2dot

    else
        
        γ_u1 = @SMatrix zeros(3,3)
        γ_u2 = @SMatrix zeros(3,3)

        γ_u1dot = @SMatrix zeros(3,3)
        γ_u2dot = @SMatrix zeros(3,3)

        γ_θ1 = @SMatrix zeros(3,3)
        γ_θ2 = @SMatrix zeros(3,3)
        
        κ_θ1 = @SMatrix zeros(3,3)
        κ_θ2 = @SMatrix zeros(3,3)
        
        κ_θ1dot = @SMatrix zeros(3,3)
        κ_θ2dot = @SMatrix zeros(3,3)

    end

    return (; properties..., u1_u1, u2_u2, θ1_θ1, θ2_θ2, V1dot_V1dot, V2dot_V2dot, 
        Ω1dot_Ω1dot, Ω2dot_Ω2dot, C_θ1, C_θ2, C_θ3, Qinv_θ1, Qinv_θ2, Qinv_θ3, 
        γ_u1, γ_u2, γ_θ1, γ_θ2, γ_u1dot, γ_u2dot, γ_F, γ_M, κ_θ1, κ_θ2, κ_θ1dot, κ_θ2dot, 
        κ_F, κ_M, V_u, P_u, P_θ, H_u, H_θ, 
        Pdot_ab, Pdot_αb, Pdot_u, Pdot_θ, Pdot_udot, Pdot_Vdot, Pdot_Ωdot,
        Hdot_ab, Hdot_αb, Hdot_u, Hdot_θ, Hdot_udot, Hdot_Vdot, Hdot_Ωdot)
end

"""
    newmark_element_jacobian_properties(properties, x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity, 
        Vdot_init, Ωdot_init, dt)

Calculate/extract the element properties needed to calculate the jacobian entries 
corresponding to a element for a Newmark scheme time marching analysis
"""
@inline function newmark_element_jacobian_properties(properties, x, indices, 
    force_scaling, structural_damping, assembly, ielem, prescribed_conditions, gravity, 
    Vdot_init, Ωdot_init, dt)

    properties = steady_element_jacobian_properties(properties, x, indices, 
        force_scaling, structural_damping, assembly, ielem, 
        prescribed_conditions, gravity)

    @unpack ωb, C, Cab, CtCab, CtCabdot, mass11, mass12, mass21, mass22, V, Ω, Vdot, Ωdot, 
        C_θ1, C_θ2, C_θ3 = properties

    # linear and angular momentum rates
    Pdot_u = @SMatrix zeros(3,3)

    Pdot_θ = mul3(C_θ1', C_θ2', C_θ3', Cab*mass11*CtCab'*Vdot) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass12*CtCab'*Ωdot) +
        CtCab*mass11*Cab'*mul3(C_θ1, C_θ2, C_θ3, Vdot) + 
        CtCab*mass12*Cab'*mul3(C_θ1, C_θ2, C_θ3, Ωdot) +
        tilde(Ω - ωb)*mul3(C_θ1', C_θ2', C_θ3', Cab*mass11*CtCab'*V) + 
        tilde(Ω - ωb)*mul3(C_θ1', C_θ2', C_θ3', Cab*mass12*CtCab'*Ω) +
        CtCabdot*mass11*Cab'*mul3(C_θ1, C_θ2, C_θ3, V) + 
        CtCabdot*mass12*Cab'*mul3(C_θ1, C_θ2, C_θ3, Ω) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass11*CtCabdot'*V) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass12*CtCabdot'*Ω) +
        -CtCab*mass11*Cab'*mul3(C_θ1, C_θ2, C_θ3, tilde(Ω - ωb)*V) + 
        -CtCab*mass12*Cab'*mul3(C_θ1, C_θ2, C_θ3, tilde(Ω - ωb)*Ω)

    Pdot_V = CtCab*mass11*CtCabdot' + CtCabdot*mass11*CtCab'

    Pdot_Ω = CtCab*mass12*CtCabdot' + CtCabdot*mass12*CtCab' +
        CtCab*mass11*CtCab'*tilde(V) + CtCab*mass12*CtCab'*tilde(Ω) +
        -tilde(CtCab*mass11*CtCab'*V) - tilde(CtCab*mass12*CtCab'*Ω)
    
    Pdot_V += 2/dt*CtCab*mass11*CtCab'
    
    Pdot_Ω += 2/dt*CtCab*mass12*CtCab'

    Hdot_u = @SMatrix zeros(3,3)

    Hdot_θ = mul3(C_θ1', C_θ2', C_θ3', Cab*mass21*CtCab'*Vdot) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass22*CtCab'*Ωdot) +
        CtCab*mass21*Cab'*mul3(C_θ1, C_θ2, C_θ3, Vdot) + 
        CtCab*mass22*Cab'*mul3(C_θ1, C_θ2, C_θ3, Ωdot) +
        tilde(Ω - ωb)*mul3(C_θ1', C_θ2', C_θ3', Cab*mass21*CtCab'*V) + 
        tilde(Ω - ωb)*mul3(C_θ1', C_θ2', C_θ3', Cab*mass22*CtCab'*Ω) +
        CtCabdot*mass21*Cab'*mul3(C_θ1, C_θ2, C_θ3, V) + 
        CtCabdot*mass22*Cab'*mul3(C_θ1, C_θ2, C_θ3, Ω) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass21*CtCabdot'*V) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass22*CtCabdot'*Ω) +
        -CtCab*mass21*Cab'*mul3(C_θ1, C_θ2, C_θ3, tilde(Ω - ωb)*V) + 
        -CtCab*mass22*Cab'*mul3(C_θ1, C_θ2, C_θ3, tilde(Ω - ωb)*Ω)

    Hdot_V = CtCab*mass21*CtCabdot' + CtCabdot*mass21*CtCab'

    Hdot_Ω = CtCab*mass22*CtCabdot' + CtCabdot*mass22*CtCab' +
        CtCab*mass21*CtCab'*tilde(V) + CtCab*mass22*CtCab'*tilde(Ω) +
        -tilde(CtCab*mass21*CtCab'*V) - tilde(CtCab*mass22*CtCab'*Ω)

    Hdot_V += 2/dt*CtCab*mass21*CtCab'

    Hdot_Ω += 2/dt*CtCab*mass22*CtCab'

    return (; properties..., Pdot_u, Pdot_θ, Pdot_V, Pdot_Ω, Hdot_u, Hdot_θ, Hdot_V, Hdot_Ω) 
end

"""
    dynamic_element_jacobian_properties(properties, dx, x, indices, 
        force_scaling, assembly, ielem, prescribed_conditions, gravity)

Calculate/extract the element properties needed to calculate the jacobian entries 
corresponding to a element for a dynamic analysis
"""
@inline function dynamic_element_jacobian_properties(properties, dx, x, indices, 
    force_scaling, structural_damping, assembly, ielem, prescribed_conditions, gravity)

    properties = steady_element_jacobian_properties(properties, x, indices, 
        force_scaling, structural_damping, assembly, ielem, 
        prescribed_conditions, gravity)

    @unpack ωb, C, Cab, CtCab, CtCabdot, mass11, mass12, mass21, mass22, V, Ω, Vdot, Ωdot, 
        C_θ1, C_θ2, C_θ3 = properties

    # linear and angular momentum rates

    Pdot_u = @SMatrix zeros(3,3)

    Pdot_θ = mul3(C_θ1', C_θ2', C_θ3', Cab*mass11*CtCab'*Vdot) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass12*CtCab'*Ωdot) +
        CtCab*mass11*Cab'*mul3(C_θ1, C_θ2, C_θ3, Vdot) + 
        CtCab*mass12*Cab'*mul3(C_θ1, C_θ2, C_θ3, Ωdot) +
        tilde(Ω - ωb)*mul3(C_θ1', C_θ2', C_θ3', Cab*mass11*CtCab'*V) + 
        tilde(Ω - ωb)*mul3(C_θ1', C_θ2', C_θ3', Cab*mass12*CtCab'*Ω) +
        CtCabdot*mass11*Cab'*mul3(C_θ1, C_θ2, C_θ3, V) + 
        CtCabdot*mass12*Cab'*mul3(C_θ1, C_θ2, C_θ3, Ω) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass11*CtCabdot'*V) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass12*CtCabdot'*Ω) +
        -CtCab*mass11*Cab'*mul3(C_θ1, C_θ2, C_θ3, tilde(Ω - ωb)*V) + 
        -CtCab*mass12*Cab'*mul3(C_θ1, C_θ2, C_θ3, tilde(Ω - ωb)*Ω)

    Pdot_V = CtCab*mass11*CtCabdot' + CtCabdot*mass11*CtCab'

    Pdot_Ω = CtCab*mass12*CtCabdot' + CtCabdot*mass12*CtCab' +
        CtCab*mass11*CtCab'*tilde(V) + CtCab*mass12*CtCab'*tilde(Ω) +
        -tilde(CtCab*mass11*CtCab'*V) - tilde(CtCab*mass12*CtCab'*Ω)
    
    Hdot_u = @SMatrix zeros(3,3)

    Hdot_θ = mul3(C_θ1', C_θ2', C_θ3', Cab*mass21*CtCab'*Vdot) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass22*CtCab'*Ωdot) +
        CtCab*mass21*Cab'*mul3(C_θ1, C_θ2, C_θ3, Vdot) + 
        CtCab*mass22*Cab'*mul3(C_θ1, C_θ2, C_θ3, Ωdot) +
        tilde(Ω - ωb)*mul3(C_θ1', C_θ2', C_θ3', Cab*mass21*CtCab'*V) + 
        tilde(Ω - ωb)*mul3(C_θ1', C_θ2', C_θ3', Cab*mass22*CtCab'*Ω) +
        CtCabdot*mass21*Cab'*mul3(C_θ1, C_θ2, C_θ3, V) + 
        CtCabdot*mass22*Cab'*mul3(C_θ1, C_θ2, C_θ3, Ω) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass21*CtCabdot'*V) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass22*CtCabdot'*Ω) +
        -CtCab*mass21*Cab'*mul3(C_θ1, C_θ2, C_θ3, tilde(Ω - ωb)*V) + 
        -CtCab*mass22*Cab'*mul3(C_θ1, C_θ2, C_θ3, tilde(Ω - ωb)*Ω)

    Hdot_V = CtCab*mass21*CtCabdot' + CtCabdot*mass21*CtCab'

    Hdot_Ω = CtCab*mass22*CtCabdot' + CtCabdot*mass22*CtCab' +
        CtCab*mass21*CtCab'*tilde(V) + CtCab*mass22*CtCab'*tilde(Ω) +
        -tilde(CtCab*mass21*CtCab'*V) - tilde(CtCab*mass22*CtCab'*Ω)

    return (; properties..., Pdot_u, Pdot_θ, Pdot_V, Pdot_Ω, Hdot_u, Hdot_θ, Hdot_V, Hdot_Ω) 
end
"""
    expanded_element_jacobian_properties(properties, x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity)
Calculate/extract the element properties needed to calculate the jacobian entries 
corresponding to a element for a steady state analysis
"""
@inline function expanded_steady_element_jacobian_properties(properties, x, indices, 
    force_scaling, structural_damping, assembly, ielem, prescribed_conditions, gravity)

    @unpack Δx, ωb, ab, αb, L, C1, C2, C, Cab, CtCab, Q, Qinv1, Qinv2, 
        mass11, mass12, mass21, mass22, S11, S12, S21, S22, u1, u2, θ1, θ2, V1, V2, Ω1, Ω2, 
        u, θ, V, Ω = properties

    # linear and angular displacement
    u1_u1, θ1_θ1 = point_displacement_jacobians(assembly.start[ielem], prescribed_conditions)
    u2_u2, θ2_θ2 = point_displacement_jacobians(assembly.stop[ielem], prescribed_conditions)

    # rotation parameter matrices
    C_θ1, C_θ2, C_θ3 = get_C_θ(C, θ)
    C1_θ1, C1_θ2, C1_θ3 = get_C_θ(C1, θ1)
    C2_θ1, C2_θ2, C2_θ3 = get_C_θ(C2, θ2)
    Q_θ1, Q_θ2, Q_θ3 = get_Q_θ(Q, θ)
    Qinv_θ1, Qinv_θ2, Qinv_θ3 = get_Qinv_θ(θ)
    Qinv1_θ1, Qinv1_θ2, Qinv1_θ3 = get_Qinv_θ(θ1)
    Qinv2_θ1, Qinv2_θ2, Qinv2_θ3 = get_Qinv_θ(θ2)

    # strain and curvature
    γ_F = S11
    γ_M = S12
    κ_F = S21
    κ_M = S22

    # linear displacement rates 
    udot1_u1 = -tilde(ωb)
    udot2_u2 = -tilde(ωb)

    udot1_θ1 = mul3(C1_θ1', C1_θ2', C1_θ3', V1)
    udot2_θ2 = mul3(C2_θ1', C2_θ2', C2_θ3', V2)

    udot1_V1 = C1'
    udot2_V2 = C2'

    udot_u1 = udot1_u1/2
    udot_u2 = udot2_u2/2

    udot_θ1 = udot1_θ1/2
    udot_θ2 = udot2_θ2/2

    udot_V1 = udot1_V1/2
    udot_V2 = udot2_V2/2

    # angular displacement rates
    θdot1_θ1 = mul3(Qinv1_θ1, Qinv1_θ2, Qinv1_θ3, Ω1 - C1*ωb) - Qinv1*mul3(C1_θ1, C1_θ2, C1_θ3, ωb)
    θdot2_θ2 = mul3(Qinv2_θ1, Qinv2_θ2, Qinv2_θ3, Ω2 - C2*ωb) - Qinv2*mul3(C2_θ1, C2_θ2, C2_θ3, ωb)
    
    θdot1_Ω1 = Qinv1
    θdot2_Ω2 = Qinv2

    θdot_θ1 = θdot1_θ1/2
    θdot_θ2 = θdot2_θ2/2

    θdot_Ω1 = θdot1_Ω1/2
    θdot_Ω2 = θdot2_Ω2/2

    # linear and angular momentum
    P_V = mass11
    P_Ω = mass12
    H_V = mass21
    H_Ω = mass22

    # linear and angular acceleration
    Vdot_ab = CtCab'
    Vdot_αb = -CtCab'*tilde(Δx + u)
    Vdot_u = CtCab'*tilde(αb)
    Vdot_θ = Cab'*mul3(C_θ1, C_θ2, C_θ3, ab + cross(αb, Δx) + cross(αb, u))

    Ωdot_αb = CtCab'
    Ωdot_θ = Cab'*mul3(C_θ1, C_θ2, C_θ3, αb)

    # linear and angular momentum rates
    Pdot_Vdot = mass11
    Pdot_Ωdot = mass12
    Hdot_Vdot = mass21
    Hdot_Ωdot = mass22

    Pdot_ab = Pdot_Vdot*Vdot_ab
    Pdot_αb = Pdot_Vdot*Vdot_αb + Pdot_Ωdot*Ωdot_αb
    Pdot_u = Pdot_Vdot*Vdot_u
    Pdot_θ = Pdot_Vdot*Vdot_θ + Pdot_Ωdot*Ωdot_θ

    Hdot_ab = Hdot_Vdot*Vdot_ab
    Hdot_αb = Hdot_Vdot*Vdot_αb + Hdot_Ωdot*Ωdot_αb
    Hdot_u = Hdot_Vdot*Vdot_u
    Hdot_θ = Hdot_Vdot*Vdot_θ + Hdot_Ωdot*Ωdot_θ

    if structural_damping

        @unpack C1, C2, Qinv1, Qinv2, u1, u2, θ1, θ2, Ω1, Ω2, μ11, μ22, udot, θdot, 
            Δx1, Δx2, Δu, Δθ, ΔQ, Δudot, Δθdot = properties

        # change in linear displacement
        Δu_u1 = -I3
        Δu_u2 =  I3

        # change in linear and angular displacement rates
        Δudot_u1 = -udot1_u1
        Δudot_u2 =  udot2_u2

        Δudot_θ1 = -udot1_θ1
        Δudot_θ2 =  udot2_θ2

        Δudot_V1 = -udot1_V1
        Δudot_V2 =  udot2_V2

        Δθdot_θ1 = -θdot1_θ1
        Δθdot_θ2 =  θdot2_θ2

        Δθdot_Ω1 = -θdot1_Ω1
        Δθdot_Ω2 =  θdot2_Ω2   

        # ΔQ matrix (see structural damping theory)
        ΔQ_θ1, ΔQ_θ2, ΔQ_θ3 = get_ΔQ_θ(θ, Δθ, Q, Q_θ1, Q_θ2, Q_θ3)

        ΔQ_Δθ1 = mul3(Q_θ1, Q_θ2, Q_θ3, e1)
        ΔQ_Δθ2 = mul3(Q_θ1, Q_θ2, Q_θ3, e2)
        ΔQ_Δθ3 = mul3(Q_θ1, Q_θ2, Q_θ3, e3)

        # strain rates
        tmp = CtCab'*tilde(C'*Ω - ωb)
        γdot_u1 = -tmp*Δu_u1 + CtCab'*Δudot_u1 
        γdot_u2 = -tmp*Δu_u2 + CtCab'*Δudot_u2 

        tmp = -Cab'*mul3(C_θ1, C_θ2, C_θ3, tilde(CtCab*Ω - ωb)*Δu) + 
            CtCab'*tilde(Δu)*mul3(C_θ1', C_θ2', C_θ3', Cab*Ω) + 
            Cab'*mul3(C_θ1, C_θ2, C_θ3, Δudot) - 
            L*Cab'*mul3(C_θ1, C_θ2, C_θ3, tilde(CtCab*Ω - ωb)*Cab*e1) + 
            L*CtCab'*tilde(Cab*e1)*mul3(C_θ1', C_θ2', C_θ3', Cab*Ω)
        γdot_θ1 = 1/2*tmp + CtCab'*Δudot_θ1
        γdot_θ2 = 1/2*tmp + CtCab'*Δudot_θ2

        γdot_V1 = CtCab'*Δudot_V1
        γdot_V2 = CtCab'*Δudot_V2

        γdot_Ω = CtCab'*tilde(Δu)*C' + L*CtCab'*tilde(Cab*e1)*C'

        tmp1 = Cab'*mul3(Q_θ1, Q_θ2, Q_θ3, Δθdot)
        tmp2 = Cab'*mul3(ΔQ_θ1, ΔQ_θ2, ΔQ_θ3, θdot) 
        tmp3 = Cab'*mul3(ΔQ_Δθ1, ΔQ_Δθ2, ΔQ_Δθ3, θdot)
        κdot_θ1 = 1/2*tmp1 + 1/2*tmp2 - tmp3 + Cab'*Q*Δθdot_θ1 + Cab'*ΔQ*θdot_θ1
        κdot_θ2 = 1/2*tmp1 + 1/2*tmp2 + tmp3 + Cab'*Q*Δθdot_θ2 + Cab'*ΔQ*θdot_θ2  

        κdot_Ω1 = Cab'*Q*Δθdot_Ω1 + Cab'*ΔQ*θdot_Ω1
        κdot_Ω2 = Cab'*Q*Δθdot_Ω2 + Cab'*ΔQ*θdot_Ω2

        # adjust strains to account for strain rates

        γ_u1 = -μ11*γdot_u1
        γ_u2 = -μ11*γdot_u2

        γ_θ1 = -μ11*γdot_θ1
        γ_θ2 = -μ11*γdot_θ2
        
        γ_V1 = -μ11*γdot_V1
        γ_V2 = -μ11*γdot_V2
        
        γ_Ω = -μ11*γdot_Ω

        κ_θ1 = -μ22*κdot_θ1
        κ_θ2 = -μ22*κdot_θ2
        
        κ_Ω1 = -μ22*κdot_Ω1
        κ_Ω2 = -μ22*κdot_Ω2

    else

        γ_u1 = @SMatrix zeros(3,3)
        γ_u2 = @SMatrix zeros(3,3)

        γ_θ1 = @SMatrix zeros(3,3)
        γ_θ2 = @SMatrix zeros(3,3)
        
        γ_V1 = @SMatrix zeros(3,3)
        γ_V2 = @SMatrix zeros(3,3)
        
        γ_Ω = @SMatrix zeros(3,3)

        κ_θ1 = @SMatrix zeros(3,3)
        κ_θ2 = @SMatrix zeros(3,3)
        
        κ_Ω1 = @SMatrix zeros(3,3)
        κ_Ω2 = @SMatrix zeros(3,3)
         
    end

    return (; properties..., u1_u1, u2_u2, θ1_θ1, θ2_θ2, 
        C1_θ1, C1_θ2, C1_θ3, C2_θ1, C2_θ2, C2_θ3, C_θ1, C_θ2, C_θ3, Qinv_θ1, Qinv_θ2, Qinv_θ3,
        γ_F, γ_M, γ_u1, γ_u2, γ_θ1, γ_θ2, γ_V1, γ_V2, γ_Ω, 
        κ_F, κ_M, κ_θ1, κ_θ2, κ_Ω1, κ_Ω2, P_V, P_Ω, H_V, H_Ω,
        udot_u1, udot_u2, udot_θ1, udot_θ2, udot_V1, udot_V2, θdot_θ1, θdot_θ2, θdot_Ω1, θdot_Ω2, 
        Pdot_ab, Pdot_αb, Pdot_u, Pdot_θ, Hdot_ab, Hdot_αb, Hdot_u, Hdot_θ,
        )
end

"""
    expanded_dynamic_element_jacobian_properties(properties, dx, x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity)

Calculate/extract the element properties needed to calculate the jacobian entries 
corresponding to a element for a dynamic analysis
"""
@inline function expanded_dynamic_element_jacobian_properties(properties, dx, x, indices, 
    force_scaling, structural_damping, assembly, ielem, prescribed_conditions, gravity)

    properties = expanded_steady_element_jacobian_properties(properties, x, indices, 
        force_scaling, structural_damping, assembly, ielem, prescribed_conditions, gravity)

    Pdot_u = @SMatrix zeros(3,3)
    Pdot_θ = @SMatrix zeros(3,3)
    Hdot_u = @SMatrix zeros(3,3)
    Hdot_θ = @SMatrix zeros(3,3)

    return (; properties..., Pdot_u, Pdot_θ, Hdot_u, Hdot_θ)
end

"""
    mass_matrix_element_jacobian_properties(x, indices, force_scaling, 
    assembly, ielem, prescribed_conditions)

Calculate/extract the element properties needed to calculate the mass matrix jacobian entries 
corresponding to a beam element
"""
@inline function mass_matrix_element_jacobian_properties(x, indices, force_scaling, 
    assembly, ielem, prescribed_conditions)

    # element properties
    @unpack L, Cab, compliance, mass = assembly.elements[ielem]
   
    # scale mass matrix by the element length
    mass *= L

    # mass submatrices
    mass11 = mass[SVector{3}(1:3), SVector{3}(1:3)]
    mass12 = mass[SVector{3}(1:3), SVector{3}(4:6)]
    mass21 = mass[SVector{3}(4:6), SVector{3}(1:3)]
    mass22 = mass[SVector{3}(4:6), SVector{3}(4:6)]

    # linear and angular displacement
    u1, θ1 = point_displacement(x, assembly.start[ielem], indices.icol_point, prescribed_conditions)
    u2, θ2 = point_displacement(x, assembly.stop[ielem], indices.icol_point, prescribed_conditions)
    u = (u1 + u2)/2
    θ = (θ1 + θ2)/2

    # transformation matrices
    C = get_C(θ)
    Ct = C'
    CtCab = Ct*Cab

    # linear and angular momentum rates
    Pdot_Vdot = CtCab*mass11*CtCab'
    Pdot_Ωdot = CtCab*mass12*CtCab'
    Hdot_Vdot = CtCab*mass21*CtCab'
    Hdot_Ωdot = CtCab*mass22*CtCab'

    return (; L, Pdot_Vdot, Pdot_Ωdot, Hdot_Vdot, Hdot_Ωdot)
end

"""
    expanded_mass_matrix_element_jacobian_properties(assembly, ielem, prescribed_conditions)
Calculate/extract the element properties needed to calculate the mass matrix jacobian entries 
corresponding to a beam element for a constant mass matrix system
"""
@inline function expanded_mass_matrix_element_jacobian_properties(assembly, ielem, prescribed_conditions)

    # element properties
    @unpack L, Cab, compliance, mass = assembly.elements[ielem]
   
    # scale mass matrix by the element length
    mass *= L

    # mass submatrices
    mass11 = mass[SVector{3}(1:3), SVector{3}(1:3)]
    mass12 = mass[SVector{3}(1:3), SVector{3}(4:6)]
    mass21 = mass[SVector{3}(4:6), SVector{3}(1:3)]
    mass22 = mass[SVector{3}(4:6), SVector{3}(4:6)]

    # linear and angular momentum rates
    Pdot_Vdot = mass11
    Pdot_Ωdot = mass12
    Hdot_Vdot = mass21
    Hdot_Ωdot = mass22

    return (; Pdot_Vdot, Pdot_Ωdot, Hdot_Vdot, Hdot_Ωdot)
end

# --- Compatability Residuals --- #

"""
    compatibility_residuals(properties)

Calculate the compatibility residuals for the beam element
"""
@inline function compatibility_residuals(properties)
   
    @unpack L, Cab, CtCab, Qinv, u1, u2, θ1, θ2, γ, κ = properties

    Δu = CtCab*(L*e1 + γ) - L*Cab*e1

    Δθ = Qinv*Cab*κ

    ru = u2 - u1 - Δu
    rθ = θ2 - θ1 - Δθ

    return (; ru, rθ)
end

@inline function static_compatibility_jacobians(properties)
   
    @unpack L, Cab, CtCab, Qinv, γ, κ, u1_u1, u2_u2, θ1_θ1, θ2_θ2, 
        γ_F, γ_M, κ_F, κ_M, C_θ1, C_θ2, C_θ3, Qinv_θ1, Qinv_θ2, Qinv_θ3 = properties

    ru_u1 = -u1_u1
    ru_u2 = u2_u2

    Δu_θ = mul3(C_θ1', C_θ2', C_θ3', Cab*(L*e1 + γ))
    ru_θ1 = -1/2*Δu_θ
    ru_θ2 = -1/2*Δu_θ
    
    Δu_F = CtCab*γ_F
    ru_F = -Δu_F
    
    Δu_M = CtCab*γ_M
    ru_M = -Δu_M

    Δθ_θ = mul3(Qinv_θ1, Qinv_θ2, Qinv_θ3, Cab*κ)
    rθ_θ1 = -I3 - 1/2*Δθ_θ
    rθ_θ2 = I3 - 1/2*Δθ_θ
    
    Δθ_F = Qinv*Cab*κ_F
    rθ_F = -Δθ_F
    
    Δθ_M = Qinv*Cab*κ_M
    rθ_M = -Δθ_M

    return (; ru_u1, ru_u2, ru_θ1, ru_θ2, ru_F, ru_M, rθ_θ1, rθ_θ2, rθ_F, rθ_M)
end

@inline function initial_compatibility_jacobians(properties)
   
    jacobians = static_compatibility_jacobians(properties)

    @unpack ru_u1, ru_u2, ru_θ1, ru_θ2, rθ_θ1, rθ_θ2 = jacobians

    @unpack Cab, CtCab, Qinv, γ_u1, γ_u2, γ_θ1, γ_θ2, γ_u1dot, γ_u2dot, κ_θ1, κ_θ2, κ_θ1dot, κ_θ2dot = properties

    ru_u1 -= CtCab*γ_u1
    ru_u2 -= CtCab*γ_u2

    ru_θ1 -= CtCab*γ_θ1
    ru_θ2 -= CtCab*γ_θ2
    
    ru_u1dot = -CtCab*γ_u1dot
    ru_u2dot = -CtCab*γ_u2dot

    rθ_θ1 -= Qinv*Cab*κ_θ1
    rθ_θ2 -= Qinv*Cab*κ_θ2
      
    rθ_θ1dot = -Qinv*Cab*κ_θ1dot
    rθ_θ2dot = -Qinv*Cab*κ_θ2dot

    return (; jacobians..., ru_u1, ru_u2, ru_θ1, ru_θ2, ru_u1dot, ru_u2dot, rθ_θ1, rθ_θ2, 
        rθ_θ1dot, rθ_θ2dot)
end

@inline function dynamic_compatibility_jacobians(properties)
   
    jacobians = static_compatibility_jacobians(properties)

    @unpack Cab, CtCab, Qinv, γ_u1, γ_u2, γ_θ1, γ_θ2, γ_V1, γ_V2, γ_Ω1, γ_Ω2, κ_θ1, κ_θ2, κ_Ω1, κ_Ω2 = properties

    @unpack ru_u1, ru_u2, ru_θ1, ru_θ2, rθ_θ1, rθ_θ2 = jacobians

    ru_u1 -= CtCab*γ_u1
    ru_u2 -= CtCab*γ_u2

    ru_θ1 -= CtCab*γ_θ1
    ru_θ2 -= CtCab*γ_θ2
    
    ru_V1 = -CtCab*γ_V1
    ru_V2 = -CtCab*γ_V2

    ru_Ω1 = -CtCab*γ_Ω1
    ru_Ω2 = -CtCab*γ_Ω2

    rθ_θ1 -= Qinv*Cab*κ_θ1
    rθ_θ2 -= Qinv*Cab*κ_θ2
    
    rθ_Ω1 = -Qinv*Cab*κ_Ω1
    rθ_Ω2 = -Qinv*Cab*κ_Ω2

    return (; jacobians..., ru_u1, ru_u2, ru_θ1, ru_θ2, ru_V1, ru_V2, ru_Ω1, ru_Ω2, 
        rθ_θ1, rθ_θ2, rθ_Ω1, rθ_Ω2)
end

@inline function expanded_compatibility_jacobians(properties)
   
    @unpack L, Cab, CtCab, Qinv, γ, κ, C_θ1, C_θ2, C_θ3, Qinv_θ1, Qinv_θ2, Qinv_θ3, 
        γ_u1, γ_u2, γ_θ1, γ_θ2, γ_V1, γ_V2, γ_Ω, γ_F, γ_M, 
        κ_θ1, κ_θ2, κ_Ω1, κ_Ω2, κ_F, κ_M = properties

    ru_u1 = -I3 - CtCab*γ_u1
    ru_u2 =  I3 - CtCab*γ_u2

    Δu_θ = mul3(C_θ1', C_θ2', C_θ3', Cab*(L*e1 + γ))
    ru_θ1 = -1/2*Δu_θ - CtCab*γ_θ1
    ru_θ2 = -1/2*Δu_θ - CtCab*γ_θ2
    
    ru_V1 = -CtCab*γ_V1
    ru_V2 = -CtCab*γ_V2

    ru_Ω = -CtCab*γ_Ω

    Δu_F = CtCab*γ_F
    ru_F1 = -1/2*Δu_F
    ru_F2 = -1/2*Δu_F

    Δu_M = CtCab*γ_M
    ru_M1 = -1/2*Δu_M
    ru_M2 = -1/2*Δu_M

    Δθ_θ = mul3(Qinv_θ1, Qinv_θ2, Qinv_θ3, Cab*κ)
    rθ_θ1 = -I3 - 1/2*Δθ_θ - Qinv*Cab*κ_θ1
    rθ_θ2 =  I3 - 1/2*Δθ_θ - Qinv*Cab*κ_θ2
    
    rθ_Ω1 = -Qinv*Cab*κ_Ω1
    rθ_Ω2 = -Qinv*Cab*κ_Ω2

    Δθ_F = Qinv*Cab*κ_F
    rθ_F1 = -1/2*Δθ_F
    rθ_F2 = -1/2*Δθ_F

    Δθ_M = Qinv*Cab*κ_M
    rθ_M1 = -1/2*Δθ_M
    rθ_M2 = -1/2*Δθ_M

    return (; ru_u1, ru_u2, ru_θ1, ru_θ2, ru_V1, ru_V2, ru_Ω, ru_F1, ru_F2, ru_M1, ru_M2, 
        rθ_θ1, rθ_θ2, rθ_Ω1, rθ_Ω2, rθ_F1, rθ_F2, rθ_M1, rθ_M2)
end

# --- Velocity Residuals --- #

"""
    expanded_element_velocity_residuals(properties)

Calculate the velocity residuals `rV` and `rΩ` of a beam element for a constant mass matrix 
system.
"""
@inline function expanded_element_velocity_residuals(properties)

    @unpack Δx, vb, ωb, C, Cab, CtCab, Qinv, u, V, Ω, udot, θdot = properties
    
    rV = CtCab*V - vb - cross(ωb, Δx) - cross(ωb, u) - udot
    rΩ = Qinv*(Cab*Ω - C*ωb) - θdot

    # @unpack CtCab, V, Ω, C1, V1, Ω1, C2, V2, Ω2 = properties
    # rV = CtCab*V - 1/2*(C1'*V1 + C2'*V2)
    # rΩ = CtCab*Ω - 1/2*(C1'*Ω1 + C2'*Ω2)

    return (; rV, rΩ)
end

"""
    expanded_element_velocity_jacobians(properties)
Calculate the jacobians of the element velocity residuals for a constant mass matrix system.
"""
@inline function expanded_element_velocity_jacobians(properties)

    @unpack ωb, C, Cab, CtCab, Qinv, V, Ω, C_θ1, C_θ2, C_θ3, Qinv_θ1, Qinv_θ2, Qinv_θ3, 
        udot_u1, udot_u2, udot_θ1, udot_θ2, udot_V1, udot_V2, 
        θdot_θ1, θdot_θ2, θdot_Ω1, θdot_Ω2 = properties
    
    tmp = -tilde(ωb)
    rV_u1 = 1/2*tmp - udot_u1
    rV_u2 = 1/2*tmp - udot_u2

    tmp = mul3(C_θ1', C_θ2', C_θ3', Cab*V)
    rV_θ1 = 1/2*tmp - udot_θ1
    rV_θ2 = 1/2*tmp - udot_θ2

    rV_V1 = -udot_V1
    rV_V2 = -udot_V2
    rV_V = CtCab

    tmp = mul3(Qinv_θ1, Qinv_θ2, Qinv_θ3, Cab*Ω - C*ωb) - Qinv*mul3(C_θ1, C_θ2, C_θ3, ωb)
    rΩ_θ1 = 1/2*tmp - θdot_θ1
    rΩ_θ2 = 1/2*tmp - θdot_θ2

    rΩ_Ω1 = -θdot_Ω1
    rΩ_Ω2 = -θdot_Ω2
    rΩ_Ω = Qinv*Cab

    # @unpack CtCab, V, Ω, C1, V1, Ω1, C2, V2, Ω2 = properties
    # rV = CtCab*V - 1/2*(C1'*V1 + C2'*V2)
    # rΩ = CtCab*Ω - 1/2*(C1'*Ω1 + C2'*Ω2)

    return (; rV_u1, rV_u2, rV_θ1, rV_θ2, rV_V1, rV_V2, rV_V, rΩ_θ1, rΩ_θ2, rΩ_Ω1, rΩ_Ω2, rΩ_Ω)
end

# --- Equilibrium Residuals --- #

"""
    expanded_element_equilibrium_residuals(properties, distributed_loads, ielem)

Calculate the equilibrium residuals of a beam element for a constant mass matrix system.
"""
@inline function expanded_element_equilibrium_residuals(properties, distributed_loads, ielem)
  
    @unpack L, Cab, CtCab, mass11, mass21, F1, F2, M1, M2, F, M, γ, κ, gvec,
        V, Ω, P, H, Pdot, Hdot = properties

    # initialize equilibrium residual
    rF = F2 - F1
    rM = M2 - M1

    # add loads due to internal forces/moments and stiffness
    rM += cross(L*e1 + γ, F)

    # add gravitational loads
    rF += mass11*CtCab'*gvec
    rM += mass21*CtCab'*gvec

    # add loads due to linear and angular momentum
    rF -= cross(Ω, P) + Pdot
    rM -= cross(Ω, H) + cross(V, P) + Hdot

    # add distributed loads
    if haskey(distributed_loads, ielem)
        dload = distributed_loads[ielem]
        f1 = CtCab'*dload.f1 + Cab'*dload.f1_follower
        f2 = CtCab'*dload.f2 + Cab'*dload.f2_follower
        rF += f1 + f2
        m1 = CtCab'*dload.m1 + Cab'*dload.m1_follower
        m2 = CtCab'*dload.m2 + Cab'*dload.m2_follower
        rM += m1 + m2
    end

    return (; rF, rM)
end

"""
    expanded_steady_element_equilibrium_jacobians(properties)

Calculate the jacobians of the element equilibrium residuals for a constant mass matrix system.
"""
@inline function expanded_steady_element_equilibrium_jacobians(properties, distributed_loads, ielem)

    jacobians = expanded_dynamic_element_equilibrium_jacobians(properties, distributed_loads, ielem)

    @unpack Pdot_ab, Pdot_αb, Hdot_ab, Hdot_αb = properties

    # add loads due to linear and angular momentum
    rF_ab = -Pdot_ab
    rF_αb = -Pdot_αb
    rM_ab = -Hdot_ab
    rM_αb = -Hdot_αb

    return (; jacobians..., rF_ab, rF_αb, rM_ab, rM_αb)
end


"""
    expanded_dynamic_element_equilibrium_jacobians(properties)
    
Calculate the jacobians of the element equilibrium residuals for a constant mass matrix system.
"""
@inline function expanded_dynamic_element_equilibrium_jacobians(properties, distributed_loads, ielem)

    @unpack L, C, Cab, CtCab, mass11, mass12, mass21, mass22, F1, F2, M1, M2, 
        V, Ω, P, H, F, M, γ, κ, gvec, γ_F, γ_M, γ_u1, γ_u2, γ_θ1, γ_θ2, γ_V1, γ_V2, γ_Ω, 
        P_V, P_Ω, H_V, H_Ω, C_θ1, C_θ2, C_θ3, Pdot_u, Pdot_θ, Hdot_u, Hdot_θ = properties

    # initialize equilibrium residual
    rF_F1 = -I3
    rF_F2 =  I3

    rM_M1 = -I3
    rM_M2 =  I3

    # add loads due to internal loads and stiffness
    
    tmp1 =  tilde(L*e1 + γ)
    tmp2 = -tilde(F)  

    rM_u1 = tmp2*γ_u1
    rM_u2 = tmp2*γ_u2

    rM_θ1 = tmp2*γ_θ1
    rM_θ2 = tmp2*γ_θ2

    rM_V1 = tmp2*γ_V1
    rM_V2 = tmp2*γ_V2

    rM_Ω = tmp2*γ_Ω

    rM_F1 = 1/2*(tmp1 + tmp2*γ_F)
    rM_F2 = 1/2*(tmp1 + tmp2*γ_F)

    rM_M1 += 1/2*tmp2*γ_M
    rM_M2 += 1/2*tmp2*γ_M

    # add gravitational loads
    
    tmp = mass11*Cab'*mul3(C_θ1, C_θ2, C_θ3, gvec)
    
    rF_θ1 = 1/2*tmp
    rF_θ2 = 1/2*tmp

    tmp = mass21*Cab'*mul3(C_θ1, C_θ2, C_θ3, gvec)
    rM_θ1 += 1/2*tmp
    rM_θ2 += 1/2*tmp

    # add loads due to linear and angular momentum
    rF_u1 = -1/2*Pdot_u
    rF_u2 = -1/2*Pdot_u
    
    rF_θ1 -= 1/2*Pdot_θ
    rF_θ2 -= 1/2*Pdot_θ

    rF_V = -tilde(Ω)*P_V

    rF_Ω = -tilde(Ω)*P_Ω + tilde(P)

    rM_u1 -= 1/2*Hdot_u
    rM_u2 -= 1/2*Hdot_u

    rM_θ1 -= 1/2*Hdot_θ
    rM_θ2 -= 1/2*Hdot_θ

    rM_V = -tilde(Ω)*H_V - tilde(V)*P_V + tilde(P)
    rM_Ω -= tilde(Ω)*H_Ω + tilde(V)*P_Ω - tilde(H)

    # add distributed loads
    if haskey(distributed_loads, ielem)
        dload = distributed_loads[ielem]

        tmp = Cab'*mul3(C_θ1, C_θ2, C_θ3, dload.f1 + dload.f2)
        rF_θ1 += 1/2*tmp
        rF_θ2 += 1/2*tmp

        tmp = Cab'*mul3(C_θ1, C_θ2, C_θ3, dload.m1 + dload.m2)
        rM_θ1 += 1/2*tmp
        rM_θ2 += 1/2*tmp

    end

    return (; rF_F1, rF_F2, rM_F1, rM_F2, rM_M1, rM_M2, rF_u1, rF_u2, rF_θ1, rF_θ2, rF_V, rF_Ω,  
        rM_u1, rM_u2, rM_θ1, rM_θ2, rM_V1, rM_V2, rM_V, rM_Ω)
end

"""
    expanded_mass_matrix_element_equilibrium_jacobians(properties)

Calculate the mass matrix jacobians for the resultant loads applied at each end of a 
beam element for a constant mass matrix system 
"""
@inline function expanded_mass_matrix_element_equilibrium_jacobians(properties)

    @unpack Pdot_Vdot, Pdot_Ωdot, Hdot_Vdot, Hdot_Ωdot = properties

    rF_Vdot = -Pdot_Vdot
    rF_Ωdot = -Pdot_Ωdot
    rM_Vdot = -Hdot_Vdot
    rM_Ωdot = -Hdot_Ωdot

    return (; rF_Vdot, rF_Ωdot, rM_Vdot, rM_Ωdot)
end

# --- Element Resultants --- #

"""
    static_element_resultants(properties, distributed_loads, ielem)

Calculate the resultant loads applied at each end of a beam element for a static analysis.
"""
@inline function static_element_resultants(properties, distributed_loads, ielem)
  
    @unpack L, C, Cab, CtCab, mass11, mass12, mass21, mass22, F, M, γ, κ, gvec = properties

    # add loads due to internal forces/moments and stiffness
    tmp = CtCab*F
    F1 = tmp
    F2 = tmp

    tmp1 = CtCab*M
    tmp2 = 1/2*CtCab*cross(L*e1 + γ, F)
    M1 = tmp1 + tmp2
    M2 = tmp1 - tmp2

    # add loads due to graviational acceleration
    tmp = -CtCab*mass11*CtCab'*gvec
    F1 -= 1/2*tmp
    F2 += 1/2*tmp

    tmp = -CtCab*mass21*CtCab'*gvec
    M1 -= 1/2*tmp
    M2 += 1/2*tmp

    # add distributed loads
    if haskey(distributed_loads, ielem)
        dload = distributed_loads[ielem]
        F1 += dload.f1 + C'*dload.f1_follower
        F2 -= dload.f2 + C'*dload.f2_follower
        M1 += dload.m1 + C'*dload.m1_follower
        M2 -= dload.m2 + C'*dload.m2_follower
    end  
   
    return (; F1, F2, M1, M2)
end

"""
    dynamic_element_resultants(properties, distributed_loads, ielem)

Calculate the resultant loads applied at each end of a beam element for a steady state 
analysis.
"""
@inline function dynamic_element_resultants(properties, distributed_loads, ielem)

    resultants = static_element_resultants(properties, distributed_loads, ielem)

    @unpack F1, F2, M1, M2 = resultants

    @unpack ωb, V, Ω, P, H, Pdot, Hdot = properties

    # add loads due to linear and angular momentum
    tmp = cross(ωb, P) + Pdot
    F1 -= 1/2*tmp
    F2 += 1/2*tmp

    tmp = cross(ωb, H) + cross(V, P) + Hdot
    M1 -= 1/2*tmp
    M2 += 1/2*tmp

    return (; resultants..., F1, F2, M1, M2)
end

"""
    expanded_element_resultants(properties)

Calculate the resultant loads applied at each end of a beam element for a constant mass
matrix system.
"""
@inline function expanded_element_resultants(properties)
   
    @unpack CtCab, C1, C2, F1, F2, M1, M2 = properties

    F1 = C1*CtCab*F1
    F2 = C2*CtCab*F2
    M1 = C1*CtCab*M1
    M2 = C2*CtCab*M2

    return (; F1, F2, M1, M2)
end

"""
    static_element_resultant_jacobians(properties, distributed_loads, ielem)

Calculate the jacobians for the resultant loads applied at each end of a beam element 
for a static analysis.
"""
@inline function static_element_resultant_jacobians(properties, distributed_loads, ielem)

    @unpack L, Cab, CtCab, mass11, mass12, mass21, mass22, F, M, γ, κ, gvec, 
        C_θ1, C_θ2, C_θ3, γ_F, γ_M = properties

    # add loads due to internal loads and stiffness
    
    tmp = mul3(C_θ1', C_θ2', C_θ3', Cab*F)

    F1_θ1 = 1/2*tmp
    F2_θ1 = 1/2*tmp
    
    F1_θ2 = 1/2*tmp
    F2_θ2 = 1/2*tmp

    F1_F = CtCab
    F2_F = CtCab

    tmp1 = mul3(C_θ1', C_θ2', C_θ3', Cab*M)
    tmp2 = 1/2*mul3(C_θ1', C_θ2', C_θ3', Cab*cross(L*e1 + γ, F))
    
    M1_θ1 = 1/2*(tmp1 + tmp2)
    M2_θ1 = 1/2*(tmp1 - tmp2)

    M1_θ2 = 1/2*(tmp1 + tmp2)
    M2_θ2 = 1/2*(tmp1 - tmp2)

    tmp = 1/2*CtCab*(tilde(L*e1 + γ) - tilde(F)*γ_F)
    M1_F = tmp
    M2_F = -tmp

    tmp = 1/2*CtCab*tilde(F)*γ_M
    M1_M = CtCab - tmp
    M2_M = CtCab + tmp
    
    # add gravitational loads
    
    tmp = -1/2*(
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass11*CtCab'*gvec) + 
        CtCab*mass11*Cab'*mul3(C_θ1, C_θ2, C_θ3, gvec)
    )
    
    F1_θ1 -= 1/2*tmp
    F2_θ1 += 1/2*tmp
    
    F1_θ2 -= 1/2*tmp
    F2_θ2 += 1/2*tmp

    tmp = -1/2*(
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass21*CtCab'*gvec) + 
        CtCab*mass21*Cab'*mul3(C_θ1, C_θ2, C_θ3, gvec)
    )
    
    M1_θ1 -= 1/2*tmp
    M2_θ1 += 1/2*tmp

    M1_θ2 -= 1/2*tmp
    M2_θ2 += 1/2*tmp

    # add distributed loads
    if haskey(distributed_loads, ielem)
        dload = distributed_loads[ielem]

        tmp = mul3(C_θ1', C_θ2', C_θ3', dload.f1_follower)
        F1_θ1 += 1/2*tmp
        F1_θ2 += 1/2*tmp

        tmp = mul3(C_θ1', C_θ2', C_θ3', dload.f2_follower)
        F2_θ1 -= 1/2*tmp
        F2_θ2 -= 1/2*tmp

        tmp = mul3(C_θ1', C_θ2', C_θ3', dload.m1_follower)
        M1_θ1 += 1/2*tmp
        M1_θ2 += 1/2*tmp

        tmp = mul3(C_θ1', C_θ2', C_θ3', dload.m2_follower)
        M2_θ1 -= 1/2*tmp
        M2_θ2 -= 1/2*tmp
    end

    return (; F1_θ1, F1_θ2, F1_F, F2_θ1, F2_θ2, F2_F,
        M1_θ1, M1_θ2, M1_F, M1_M, M2_θ1, M2_θ2, M2_F, M2_M)
end

"""
    steady_element_resultant_jacobians(properties, distributed_loads, ielem)

Calculate the jacobians for the resultant loads applied at each end of a beam element 
for a steady state analysis.
"""
@inline function steady_element_resultant_jacobians(properties, distributed_loads, ielem)

    jacobians = dynamic_element_resultant_jacobians(properties, distributed_loads, ielem)

    @unpack Pdot_ab, Pdot_αb, Hdot_ab, Hdot_αb = properties

    # add loads due to linear and angular momentum

    tmp = 1/2*Pdot_ab

    F1_ab = -tmp
    F2_ab = tmp

    tmp = 1/2*Pdot_αb

    F1_αb = -tmp
    F2_αb = tmp

    tmp = 1/2*Hdot_ab

    M1_ab = -tmp
    M2_ab = tmp

    tmp = 1/2*Hdot_αb

    M1_αb = -tmp
    M2_αb = tmp

    return (; jacobians..., F1_ab, F1_αb, F2_ab, F2_αb, M1_ab, M1_αb, M2_ab, M2_αb)
end
    
"""
    initial_element_resultant_jacobians(properties, distributed_loads, ielem)

Calculate the jacobians for the resultant loads applied at each end of V beam element 
for the initialization of a time domain analysis.
"""
@inline function initial_element_resultant_jacobians(properties, distributed_loads, ielem)

    jacobians = static_element_resultant_jacobians(properties, distributed_loads, ielem)

    @unpack ωb, L, CtCab, mass11, mass12, mass21, mass22, F, V, Ω, P, H,  
        γ_u1, γ_u2, γ_θ1, γ_θ2, γ_u1dot, γ_u2dot, V_u, P_u, P_θ, H_u, H_θ,  
        Pdot_ab, Pdot_αb, Pdot_u, Pdot_θ, Pdot_udot, Pdot_Vdot, Pdot_Ωdot,
        Hdot_ab, Hdot_αb, Hdot_u, Hdot_θ, Hdot_udot, Hdot_Vdot, Hdot_Ωdot = properties

    @unpack F1_θ1, F1_θ2, F2_θ1, F2_θ2, M1_θ1, M1_θ2, M2_θ1, M2_θ2 = jacobians

    # add loads due to internal forces/moments and stiffness
    tmp = 1/2*CtCab*tilde(F)
    
    M1_u1 = -tmp*γ_u1
    M2_u1 = tmp*γ_u1

    M1_u2 = -tmp*γ_u2
    M2_u2 =  tmp*γ_u2

    M1_θ1 -= tmp*γ_θ1
    M2_θ1 += tmp*γ_θ1

    M1_θ2 -= tmp*γ_θ2
    M2_θ2 += tmp*γ_θ2

    M1_u1dot = -tmp*γ_u1dot
    M2_u1dot =  tmp*γ_u1dot

    M1_u2dot = -tmp*γ_u2dot
    M2_u2dot =  tmp*γ_u2dot

    # add loads due to linear and angular momentum

    tmp = 1/2*Pdot_ab

    F1_ab = -tmp
    F2_ab = tmp

    tmp = 1/2*Pdot_αb

    F1_αb = -tmp
    F2_αb = tmp

    tmp = 1/2*(tilde(ωb)*P_u + Pdot_u)

    F1_u1 = -1/2*tmp
    F2_u1 = 1/2*tmp
    
    F1_u2 = -1/2*tmp
    F2_u2 = 1/2*tmp

    tmp = 1/2*(tilde(ωb)*P_θ + Pdot_θ)

    F1_θ1 -= 1/2*tmp
    F2_θ1 += 1/2*tmp
    
    F1_θ2 -= 1/2*tmp
    F2_θ2 += 1/2*tmp

    tmp = 1/2*Pdot_udot

    F1_u1dot =  -1/2*tmp
    F2_u1dot =  1/2*tmp

    F1_u2dot = -1/2*tmp
    F2_u2dot =  1/2*tmp

    tmp = 1/2*Pdot_Vdot

    F1_V1dot =  -1/2*tmp
    F2_V1dot =  1/2*tmp

    F1_V2dot = -1/2*tmp
    F2_V2dot =  1/2*tmp

    tmp = 1/2*Pdot_Ωdot

    F1_Ω1dot = -1/2*tmp
    F2_Ω1dot =  1/2*tmp

    F1_Ω2dot = -1/2*tmp
    F2_Ω2dot =  1/2*tmp

    tmp = 1/2*Hdot_ab

    M1_ab = -tmp
    M2_ab = tmp

    tmp = 1/2*Hdot_αb

    M1_αb = -tmp
    M2_αb = tmp

    tmp = 1/2*(tilde(ωb)*H_u + tilde(V)*P_u - tilde(P)*V_u + Hdot_u)

    M1_u1 -= 1/2*tmp
    M2_u1 += 1/2*tmp

    M1_u2 -= 1/2*tmp
    M2_u2 += 1/2*tmp

    tmp = 1/2*(tilde(ωb)*H_θ + tilde(V)*P_θ + Hdot_θ)
    
    M1_θ1 -= 1/2*tmp
    M2_θ1 += 1/2*tmp

    M1_θ2 -= 1/2*tmp
    M2_θ2 += 1/2*tmp

    tmp = 1/2*Hdot_udot

    M1_u1dot -= 1/2*tmp
    M2_u1dot += 1/2*tmp

    M1_u2dot -= 1/2*tmp
    M2_u2dot += 1/2*tmp

    tmp = 1/2*Hdot_Vdot

    M1_V1dot = -1/2*tmp
    M2_V1dot =  1/2*tmp

    M1_V2dot = -1/2*tmp
    M2_V2dot =  1/2*tmp

    tmp = 1/2*Hdot_Ωdot

    M1_Ω1dot = -1/2*tmp
    M2_Ω1dot =  1/2*tmp

    M1_Ω2dot = -1/2*tmp
    M2_Ω2dot =  1/2*tmp

    return (; jacobians..., 
        F1_ab, F1_αb, F1_u1, F1_u2, F1_θ1, F1_θ2, F1_u1dot, F1_u2dot, F1_V1dot, F1_V2dot, F1_Ω1dot, F1_Ω2dot, 
        F2_ab, F2_αb, F2_u1, F2_u2, F2_θ1, F2_θ2, F2_u1dot, F2_u2dot, F2_V1dot, F2_V2dot, F2_Ω1dot, F2_Ω2dot, 
        M1_ab, M1_αb, M1_u1, M1_u2, M1_θ1, M1_θ2, M1_u1dot, M1_u2dot, M1_V1dot, M1_V2dot, M1_Ω1dot, M1_Ω2dot, 
        M2_ab, M2_αb, M2_u1, M2_u2, M2_θ1, M2_θ2, M2_u1dot, M2_u2dot, M2_V1dot, M2_V2dot, M2_Ω1dot, M2_Ω2dot)
end

"""
    dynamic_element_resultant_jacobians(properties, distributed_loads, ielem)

Calculate the jacobians for the resultant loads applied at each end of a beam element 
for a dynamic analysis.
"""
@inline function dynamic_element_resultant_jacobians(properties, distributed_loads, ielem)

    jacobians = static_element_resultant_jacobians(properties, distributed_loads, ielem)

    @unpack ωb, L, CtCab, mass11, mass12, mass21, mass22, F, V, Ω, P, H,  
        γ_u1, γ_u2, γ_θ1, γ_θ2, γ_V1, γ_V2, γ_Ω1, γ_Ω2, P_θ, P_V, P_Ω, H_θ, H_V, H_Ω,
        Pdot_u, Pdot_θ, Pdot_V, Pdot_Ω, Hdot_u, Hdot_θ, Hdot_V, Hdot_Ω = properties

    @unpack F1_θ1, F1_θ2, F2_θ1, F2_θ2, M1_θ1, M1_θ2, M2_θ1, M2_θ2 = jacobians

    # add loads due to internal forces/moments and stiffness
    tmp = 1/2*CtCab*tilde(F)
    
    M1_u1 = -tmp*γ_u1
    M2_u1 = tmp*γ_u1

    M1_u2 = -tmp*γ_u2
    M2_u2 =  tmp*γ_u2

    M1_θ1 -= tmp*γ_θ1
    M2_θ1 += tmp*γ_θ1

    M1_θ2 -= tmp*γ_θ2
    M2_θ2 += tmp*γ_θ2

    M1_V1 = -tmp*γ_V1
    M2_V1 =  tmp*γ_V1

    M1_V2 = -tmp*γ_V2
    M2_V2 =  tmp*γ_V2

    M1_Ω1 = -tmp*γ_Ω1
    M2_Ω1 =  tmp*γ_Ω1

    M1_Ω2 = -tmp*γ_Ω2
    M2_Ω2 =  tmp*γ_Ω2
    
    # add loads due to linear and angular momentum

    tmp = 1/2*Pdot_u

    F1_u1 = -1/2*tmp
    F2_u1 = 1/2*tmp
    
    F1_u2 = -1/2*tmp
    F2_u2 = 1/2*tmp

    tmp = 1/2*(tilde(ωb)*P_θ + Pdot_θ)

    F1_θ1 -= 1/2*tmp
    F2_θ1 += 1/2*tmp
    
    F1_θ2 -= 1/2*tmp
    F2_θ2 += 1/2*tmp

    tmp = 1/2*(tilde(ωb)*P_V + Pdot_V)

    F1_V1 = -1/2*tmp
    F2_V1 = 1/2*tmp

    F1_V2 = -1/2*tmp
    F2_V2 = 1/2*tmp

    tmp = 1/2*(tilde(ωb)*P_Ω + Pdot_Ω)

    F1_Ω1 = -1/2*tmp
    F2_Ω1 = 1/2*tmp

    F1_Ω2 = -1/2*tmp
    F2_Ω2 = 1/2*tmp

    tmp = 1/2*Hdot_u

    M1_u1 -= 1/2*tmp
    M2_u1 += 1/2*tmp

    M1_u2 -= 1/2*tmp
    M2_u2 += 1/2*tmp

    tmp = 1/2*(tilde(ωb)*H_θ + tilde(V)*P_θ + Hdot_θ)
    
    M1_θ1 -= 1/2*tmp
    M2_θ1 += 1/2*tmp

    M1_θ2 -= 1/2*tmp
    M2_θ2 += 1/2*tmp

    tmp = 1/2*(tilde(ωb)*H_V + tilde(V)*P_V - tilde(P) + Hdot_V)
    
    M1_V1 -= 1/2*tmp
    M2_V1 += 1/2*tmp

    M1_V2 -= 1/2*tmp
    M2_V2 += 1/2*tmp

    tmp = 1/2*(tilde(ωb)*H_Ω + tilde(V)*P_Ω + Hdot_Ω)
    
    M1_Ω1 -= 1/2*tmp
    M2_Ω1 += 1/2*tmp

    M1_Ω2 -= 1/2*tmp
    M2_Ω2 += 1/2*tmp

    return (; jacobians..., 
        F1_u1, F1_u2, F1_θ1, F1_θ2, F1_V1, F1_V2, F1_Ω1, F1_Ω2, 
        F2_u1, F2_u2, F2_θ1, F2_θ2, F2_V1, F2_V2, F2_Ω1, F2_Ω2, 
        M1_u1, M1_u2, M1_θ1, M1_θ2, M1_V1, M1_V2, M1_Ω1, M1_Ω2, 
        M2_u1, M2_u2, M2_θ1, M2_θ2, M2_V1, M2_V2, M2_Ω1, M2_Ω2)

end

"""
    expanded_element_resultant_jacobians(properties)

Calculate the jacobians for the resultant loads applied at each end of a beam element 
for a constant mass matrix system.
"""
@inline function expanded_element_resultant_jacobians(properties)
   
    @unpack Cab, C1, C2, CtCab, F1, F2, M1, M2, C1_θ1, C1_θ2, C1_θ3, 
        C2_θ1, C2_θ2, C2_θ3, C_θ1, C_θ2, C_θ3 = properties

    F1_θ = C1*mul3(C_θ1', C_θ2', C_θ3', Cab*F1)
    F1_θ1 = 1/2*F1_θ + mul3(C1_θ1, C1_θ2, C1_θ3, CtCab*F1)
    F1_θ2 = 1/2*F1_θ
    F1_F1 = C1*CtCab

    F2_θ = C2*mul3(C_θ1', C_θ2', C_θ3', Cab*F2)
    F2_θ1 = 1/2*F2_θ
    F2_θ2 = 1/2*F2_θ + mul3(C2_θ1, C2_θ2, C2_θ3, CtCab*F2) 
    F2_F2 = C2*CtCab

    M1_θ = C1*mul3(C_θ1', C_θ2', C_θ3', Cab*M1)
    M1_θ1 = 1/2*M1_θ + mul3(C1_θ1, C1_θ2, C1_θ3, CtCab*M1) 
    M1_θ2 = 1/2*M1_θ
    M1_M1 = C1*CtCab

    M2_θ = C2*mul3(C_θ1', C_θ2', C_θ3', Cab*M2)
    M2_θ1 = 1/2*M2_θ
    M2_θ2 = 1/2*M2_θ + mul3(C2_θ1, C2_θ2, C2_θ3, CtCab*M2)
    M2_M2 = C2*CtCab

    return (; F1_θ1, F1_θ2, F1_F1, F2_θ1, F2_θ2, F2_F2, 
        M1_θ1, M1_θ2, M1_M1, M2_θ1, M2_θ2, M2_M2)
end

"""
    mass_matrix_element_resultant_jacobians(properties)

Calculate the mass matrix jacobians for the resultant loads applied at each end of a beam element 
"""
@inline function mass_matrix_element_resultant_jacobians(properties)

    @unpack Pdot_Vdot, Pdot_Ωdot, Hdot_Vdot, Hdot_Ωdot = properties

    tmp = 1/2*Pdot_Vdot

    F1_V1dot = -1/2*tmp
    F2_V1dot =  1/2*tmp

    F1_V2dot = -1/2*tmp
    F2_V2dot =  1/2*tmp

    tmp = 1/2*Pdot_Ωdot

    F1_Ω1dot = -1/2*tmp
    F2_Ω1dot =  1/2*tmp

    F1_Ω2dot = -1/2*tmp
    F2_Ω2dot =  1/2*tmp

    tmp = 1/2*Hdot_Vdot

    M1_V1dot = -1/2*tmp
    M2_V1dot =  1/2*tmp

    M1_V2dot = -1/2*tmp
    M2_V2dot =  1/2*tmp

    tmp = 1/2*Hdot_Ωdot

    M1_Ω1dot = -1/2*tmp
    M2_Ω1dot =  1/2*tmp

    M1_Ω2dot = -1/2*tmp
    M2_Ω2dot =  1/2*tmp

    return (; 
        F1_V1dot, F1_V2dot, F1_Ω1dot, F1_Ω2dot, 
        F2_V1dot, F2_V2dot, F2_Ω1dot, F2_Ω2dot, 
        M1_V1dot, M1_V2dot, M1_Ω1dot, M1_Ω2dot, 
        M2_V1dot, M2_V2dot, M2_Ω1dot, M2_Ω2dot)
end

# --- Insert Residual --- #

"""
    insert_element_residuals!(resid, indices, force_scaling, assembly, ielem, 
        compatibility, resultants)

Insert the residual entries corresponding to a beam element into the system residual vector.
"""
@inline function insert_element_residuals!(resid, indices, force_scaling, assembly, ielem, 
    compatibility, resultants)

    @unpack ru, rθ = compatibility
    @unpack F1, F2, M1, M2 = resultants

    # compatibility equations
    irow = indices.irow_elem[ielem]
    resid[irow:irow+2] .= ru
    resid[irow+3:irow+5] .= rθ

    # equilibrium equations for the start of the beam element
    irow = indices.irow_point[assembly.start[ielem]]
    @views resid[irow:irow+2] .-= F1 ./ force_scaling
    @views resid[irow+3:irow+5] .-= M1 ./ force_scaling

    # equilibrium equations for the end of the beam element
    irow = indices.irow_point[assembly.stop[ielem]]
    @views resid[irow:irow+2] .+= F2 ./ force_scaling
    @views resid[irow+3:irow+5] .+= M2 ./ force_scaling

    return resid
end

"""
    insert_expanded_element_residuals!(resid, indices, force_scaling, assembly, ielem, 
        compatibility, velocities, equilibrium, resultants)

Insert the residual entries corresponding to a beam element into the system residual vector
for a constant mass matrix system.
"""
@inline function insert_expanded_element_residuals!(resid, indices, force_scaling, assembly, 
    ielem, compatibility, velocities, equilibrium, resultants)

    @unpack ru, rθ = compatibility
    @unpack rV, rΩ = velocities
    @unpack rF, rM = equilibrium
    @unpack F1, F2, M1, M2 = resultants

    irow = indices.irow_elem[ielem]
    
    # compatibility equations for the start of the beam element
    resid[irow:irow+2] .= ru
    resid[irow+3:irow+5] .= rθ

    # velocity residuals
    resid[irow+6:irow+8] .= rV
    resid[irow+9:irow+11] .= rΩ

    # equilibrium residuals for the beam element
    resid[irow+12:irow+14] .= rF ./ force_scaling
    resid[irow+15:irow+17] .= rM ./ force_scaling

    # equilibrium equations for the start of the beam element
    irow = indices.irow_point[assembly.start[ielem]]
    @views resid[irow+6:irow+8] .-= F1 ./ force_scaling
    @views resid[irow+9:irow+11] .-= M1 ./ force_scaling

    # equilibrium equations for the end of the beam element
    irow = indices.irow_point[assembly.stop[ielem]]
    @views resid[irow+6:irow+8] .+= F2 ./ force_scaling
    @views resid[irow+9:irow+11] .+= M2 ./ force_scaling

    return resid
end

@inline function insert_static_element_jacobians!(jacob, indices, force_scaling, 
    assembly, ielem, properties, compatibility, resultants)

    @unpack u1_u1, u2_u2, θ1_θ1, θ2_θ2 = properties

    @unpack ru_u1, ru_u2, ru_θ1, ru_θ2, ru_F, ru_M, 
                          rθ_θ1, rθ_θ2, rθ_F, rθ_M = compatibility

    @unpack F1_θ1, F1_θ2, F1_F,
            F2_θ1, F2_θ2, F2_F,
            M1_θ1, M1_θ2, M1_F, M1_M,
            M2_θ1, M2_θ2, M2_F, M2_M = resultants

    icol = indices.icol_elem[ielem]
    icol1 = indices.icol_point[assembly.start[ielem]]
    icol2 = indices.icol_point[assembly.stop[ielem]]

    # compatibility equations
    irow = indices.irow_elem[ielem]

    jacob[irow:irow+2, icol1:icol1+2] .= ru_u1*u1_u1
    jacob[irow:irow+2, icol1+3:icol1+5] .= ru_θ1*θ1_θ1

    jacob[irow:irow+2, icol2:icol2+2] .= ru_u2*u2_u2
    jacob[irow:irow+2, icol2+3:icol2+5] .= ru_θ2*θ2_θ2

    jacob[irow:irow+2, icol:icol+2] .= ru_F .* force_scaling
    jacob[irow:irow+2, icol+3:icol+5] .= ru_M .* force_scaling

    jacob[irow+3:irow+5, icol1+3:icol1+5] .= rθ_θ1*θ1_θ1

    jacob[irow+3:irow+5, icol2+3:icol2+5] .= rθ_θ2*θ2_θ2

    jacob[irow+3:irow+5, icol:icol+2] .= rθ_F .* force_scaling
    jacob[irow+3:irow+5, icol+3:icol+5] .= rθ_M .* force_scaling

    # equilibrium equations for the start of the beam element
    irow1 = indices.irow_point[assembly.start[ielem]]

    @views jacob[irow1:irow1+2, icol1+3:icol1+5] .-= F1_θ1*θ1_θ1 ./ force_scaling
    @views jacob[irow1:irow1+2, icol2+3:icol2+5] .-= F1_θ2*θ2_θ2 ./ force_scaling

    jacob[irow1:irow1+2, icol:icol+2] .= -F1_F

    @views jacob[irow1+3:irow1+5, icol1+3:icol1+5] .-= M1_θ1*θ1_θ1 ./ force_scaling
    @views jacob[irow1+3:irow1+5, icol2+3:icol2+5] .-= M1_θ2*θ2_θ2 ./ force_scaling

    jacob[irow1+3:irow1+5, icol:icol+2] .= -M1_F
    jacob[irow1+3:irow1+5, icol+3:icol+5] .= -M1_M

    # equilibrium equations for the end of the beam element
    irow2 = indices.irow_point[assembly.stop[ielem]]
    
    @views jacob[irow2:irow2+2, icol1+3:icol1+5] .+= F2_θ1*θ1_θ1 ./ force_scaling
    @views jacob[irow2:irow2+2, icol2+3:icol2+5] .+= F2_θ2*θ2_θ2 ./ force_scaling

    jacob[irow2:irow2+2, icol:icol+2] .= F2_F

    @views jacob[irow2+3:irow2+5, icol1+3:icol1+5] .+= M2_θ1*θ1_θ1 ./ force_scaling
    @views jacob[irow2+3:irow2+5, icol2+3:icol2+5] .+= M2_θ2*θ2_θ2 ./ force_scaling

    jacob[irow2+3:irow2+5, icol:icol+2] .= M2_F
    jacob[irow2+3:irow2+5, icol+3:icol+5] .= M2_M

    return jacob
end

@inline function insert_steady_element_jacobians!(jacob, indices, force_scaling,
    assembly, ielem, properties, compatibility, resultants)

    insert_dynamic_element_jacobians!(jacob, indices, force_scaling,
        assembly, ielem, properties, compatibility, resultants)

    @unpack F1_ab, F1_αb, F2_ab, F2_αb, M1_ab, M1_αb, M2_ab, M2_αb = resultants

    # jacobians of prescribed linear and angular acceleration
    irow1 = indices.irow_point[assembly.start[ielem]]
    irow2 = indices.irow_point[assembly.stop[ielem]]
    icol = indices.icol_body

    for i = 1:3
        if !iszero(icol[i])
            jacob[irow1:irow1+2, icol[i]] .-= F1_ab[:,i] ./ force_scaling
            jacob[irow1+3:irow1+5, icol[i]] .-= M1_ab[:,i] ./ force_scaling
            jacob[irow2:irow2+2, icol[i]] .+= F2_ab[:,i] ./ force_scaling
            jacob[irow2+3:irow2+5, icol[i]] .+= M2_ab[:,i] ./ force_scaling
        end
    end

    for i = 4:6
        if !iszero(icol[i])
            jacob[irow1:irow1+2, icol[i]] .-= F1_αb[:,i-3] ./ force_scaling
            jacob[irow1+3:irow1+5, icol[i]] .-= M1_αb[:,i-3] ./ force_scaling
            jacob[irow2:irow2+2, icol[i]] .+= F2_αb[:,i-3] ./ force_scaling
            jacob[irow2+3:irow2+5, icol[i]] .+= M2_αb[:,i-3] ./ force_scaling
        end
    end

    return jacob
end


@inline function insert_initial_element_jacobians!(jacob, indices, force_scaling, 
    assembly, ielem, properties, compatibility, resultants)

    insert_static_element_jacobians!(jacob, indices, force_scaling, assembly, ielem, 
        properties, compatibility, resultants)

    @unpack u1_u1, u2_u2, θ1_θ1, θ2_θ2, V1dot_V1dot, V2dot_V2dot, Ω1dot_Ω1dot, Ω2dot_Ω2dot = properties

    @unpack ru_u1dot, ru_u2dot, rθ_θ1dot, rθ_θ2dot = compatibility

    @unpack F1_ab, F1_αb, F1_u1, F1_u2, F1_u1dot, F1_u2dot, F1_V1dot, F1_V2dot, F1_Ω1dot, F1_Ω2dot,
            F2_ab, F2_αb, F2_u1, F2_u2, F2_u1dot, F2_u2dot, F2_V1dot, F2_V2dot, F2_Ω1dot, F2_Ω2dot,
            M1_ab, M1_αb, M1_u1, M1_u2, M1_θ1, M1_θ2, M1_u1dot, M1_u2dot, M1_V1dot, M1_V2dot, M1_Ω1dot, M1_Ω2dot, 
            M2_ab, M2_αb, M2_u1, M2_u2, M2_θ1, M2_θ2, M2_u1dot, M2_u2dot, M2_V1dot, M2_V2dot, M2_Ω1dot, M2_Ω2dot = resultants

    icol1 = indices.icol_point[assembly.start[ielem]]
    icol2 = indices.icol_point[assembly.stop[ielem]]

    # compatibility equations
    irow = indices.irow_elem[ielem]

    jacob[irow:irow+2, icol1+6:icol1+8] .= ru_u1dot
    jacob[irow:irow+2, icol2+6:icol2+8] .= ru_u2dot

    jacob[irow+3:irow+5, icol1+9:icol1+11] .= rθ_θ1dot
    jacob[irow+3:irow+5, icol2+9:icol2+11] .= rθ_θ2dot

    # equilibrium equations for the start of the beam element
    irow1 = indices.irow_point[assembly.start[ielem]]

    @views jacob[irow1:irow1+2, icol1:icol1+2] .-= F1_u1*u1_u1 ./ force_scaling
    @views jacob[irow1:irow1+2, icol1:icol1+2] .-= F1_V1dot*V1dot_V1dot ./ force_scaling
    @views jacob[irow1:irow1+2, icol1+3:icol1+5] .-= F1_Ω1dot*Ω1dot_Ω1dot ./ force_scaling
    @views jacob[irow1:irow1+2, icol1+6:icol1+8] .-= F1_u1dot ./ force_scaling

    @views jacob[irow1:irow1+2, icol2:icol2+2] .-= F1_u2*u2_u2 ./ force_scaling
    @views jacob[irow1:irow1+2, icol2:icol2+2] .-= F1_V2dot*V2dot_V2dot ./ force_scaling
    @views jacob[irow1:irow1+2, icol2+3:icol2+5] .-= F1_Ω2dot*Ω2dot_Ω2dot ./ force_scaling
    @views jacob[irow1:irow1+2, icol2+6:icol2+8] .-= F1_u2dot ./ force_scaling

    @views jacob[irow1+3:irow1+5, icol1:icol1+2] .-= M1_u1*u1_u1 ./ force_scaling
    @views jacob[irow1+3:irow1+5, icol1:icol1+2] .-= M1_V1dot*V1dot_V1dot ./ force_scaling
    @views jacob[irow1+3:irow1+5, icol1+3:icol1+5] .-= M1_Ω1dot*Ω1dot_Ω1dot ./ force_scaling
    @views jacob[irow1+3:irow1+5, icol1+6:icol1+8] .-= M1_u1dot ./ force_scaling

    @views jacob[irow1+3:irow1+5, icol2:icol2+2] .-= M1_u2*u2_u2 ./ force_scaling
    @views jacob[irow1+3:irow1+5, icol2:icol2+2] .-= M1_V2dot*V2dot_V2dot ./ force_scaling
    @views jacob[irow1+3:irow1+5, icol2+3:icol2+5] .-= M1_Ω2dot*Ω2dot_Ω2dot ./ force_scaling
    @views jacob[irow1+3:irow1+5, icol2+6:icol2+8] .-= M1_u2dot ./ force_scaling

    # equilibrium equations for the end of the beam element
    irow2 = indices.irow_point[assembly.stop[ielem]]

    @views jacob[irow2:irow2+2, icol1:icol1+2] .+= F2_u1*u1_u1 ./ force_scaling
    @views jacob[irow2:irow2+2, icol1:icol1+2] .+= F2_V1dot*V1dot_V1dot ./ force_scaling
    @views jacob[irow2:irow2+2, icol1+3:icol1+5] .+= F2_Ω1dot*Ω1dot_Ω1dot ./ force_scaling
    @views jacob[irow2:irow2+2, icol1+6:icol1+8] .+= F2_u1dot ./ force_scaling

    @views jacob[irow2:irow2+2, icol2:icol2+2] .+= F2_u2*u2_u2 ./ force_scaling
    @views jacob[irow2:irow2+2, icol2:icol2+2] .+= F2_V2dot*V2dot_V2dot ./ force_scaling
    @views jacob[irow2:irow2+2, icol2+3:icol2+5] .+= F2_Ω2dot*Ω2dot_Ω2dot ./ force_scaling
    @views jacob[irow2:irow2+2, icol2+6:icol2+8] .+= F2_u2dot ./ force_scaling

    @views jacob[irow2+3:irow2+5, icol1:icol1+2] .+= M2_u1*u1_u1 ./ force_scaling
    @views jacob[irow2+3:irow2+5, icol1:icol1+2] .+= M2_V1dot*V1dot_V1dot ./ force_scaling
    @views jacob[irow2+3:irow2+5, icol1+3:icol1+5] .+= M2_Ω1dot*Ω1dot_Ω1dot ./ force_scaling
    @views jacob[irow2+3:irow2+5, icol1+6:icol1+8] .+= M2_u1dot ./ force_scaling

    @views jacob[irow2+3:irow2+5, icol2:icol2+2] .+= M2_u2*u2_u2 ./ force_scaling
    @views jacob[irow2+3:irow2+5, icol2:icol2+2] .+= M2_V2dot*V2dot_V2dot ./ force_scaling
    @views jacob[irow2+3:irow2+5, icol2+3:icol2+5] .+= M2_Ω2dot*Ω2dot_Ω2dot ./ force_scaling
    @views jacob[irow2+3:irow2+5, icol2+6:icol2+8] .+= M2_u2dot ./ force_scaling

    # jacobians of prescribed linear and angular acceleration
    icol = indices.icol_body
    for i = 1:3
        if !iszero(icol[i])
            jacob[irow1:irow1+2, icol[i]] .-= F1_ab[:,i] ./ force_scaling
            jacob[irow1+3:irow1+5, icol[i]] .-= M1_ab[:,i] ./ force_scaling
            jacob[irow2:irow2+2, icol[i]] .+= F2_ab[:,i] ./ force_scaling
            jacob[irow2+3:irow2+5, icol[i]] .+= M2_ab[:,i] ./ force_scaling
        end
    end

    for i = 4:6
        if !iszero(icol[i])
            jacob[irow1:irow1+2, icol[i]] .-= F1_αb[:,i-3] ./ force_scaling
            jacob[irow1+3:irow1+5, icol[i]] .-= M1_αb[:,i-3] ./ force_scaling
            jacob[irow2:irow2+2, icol[i]] .+= F2_αb[:,i-3] ./ force_scaling
            jacob[irow2+3:irow2+5, icol[i]] .+= M2_αb[:,i-3] ./ force_scaling
        end
    end

    return jacob
end

@inline function insert_dynamic_element_jacobians!(jacob, indices, force_scaling,
    assembly, ielem, properties, compatibility, resultants)

    insert_static_element_jacobians!(jacob, indices, force_scaling, assembly, ielem, 
        properties, compatibility, resultants)

    @unpack u1_u1, u2_u2 = properties

    @unpack ru_V1, ru_V2, ru_Ω1, ru_Ω2, rθ_Ω1, rθ_Ω2 = compatibility

    @unpack F1_u1, F1_u2, F1_V1, F1_V2, F1_Ω1, F1_Ω2,
            F2_u1, F2_u2, F2_V1, F2_V2, F2_Ω1, F2_Ω2,
            M1_u1, M1_u2, M1_V1, M1_V2, M1_Ω1, M1_Ω2,
            M2_u1, M2_u2, M2_V1, M2_V2, M2_Ω1, M2_Ω2 = resultants   

    icol1 = indices.icol_point[assembly.start[ielem]]
    icol2 = indices.icol_point[assembly.stop[ielem]]

    # compatibility equations
    irow = indices.irow_elem[ielem]

    jacob[irow:irow+2, icol1+6:icol1+8] .= ru_V1
    jacob[irow:irow+2, icol1+9:icol1+11] .= ru_Ω1

    jacob[irow:irow+2, icol2+6:icol2+8] .= ru_V2
    jacob[irow:irow+2, icol2+9:icol2+11] .= ru_Ω2

    jacob[irow+3:irow+5, icol1+9:icol1+11] .= rθ_Ω1

    jacob[irow+3:irow+5, icol2+9:icol2+11] .= rθ_Ω2

    # equilibrium equations for the start of the beam element
    irow1 = indices.irow_point[assembly.start[ielem]]

    @views jacob[irow1:irow1+2, icol1:icol1+2] .-= F1_u1*u1_u1 ./ force_scaling
    @views jacob[irow1:irow1+2, icol2:icol2+2] .-= F1_u2*u2_u2 ./ force_scaling

    @views jacob[irow1:irow1+2, icol1+6:icol1+8] .-= F1_V1 ./ force_scaling
    @views jacob[irow1:irow1+2, icol1+9:icol1+11] .-= F1_Ω1 ./ force_scaling

    @views jacob[irow1:irow1+2, icol2+6:icol2+8] .-= F1_V2 ./ force_scaling
    @views jacob[irow1:irow1+2, icol2+9:icol2+11] .-= F1_Ω2 ./ force_scaling

    @views jacob[irow1+3:irow1+5, icol1:icol1+2] .-= M1_u1*u1_u1 ./ force_scaling
    @views jacob[irow1+3:irow1+5, icol2:icol2+2] .-= M1_u2*u2_u2 ./ force_scaling

    @views jacob[irow1+3:irow1+5, icol1+6:icol1+8] .-= M1_V1 ./ force_scaling
    @views jacob[irow1+3:irow1+5, icol1+9:icol1+11] .-= M1_Ω1 ./ force_scaling

    @views jacob[irow1+3:irow1+5, icol2+6:icol2+8] .-= M1_V2 ./ force_scaling
    @views jacob[irow1+3:irow1+5, icol2+9:icol2+11] .-= M1_Ω2 ./ force_scaling

    # equilibrium equations for the end of the beam element
    irow2 = indices.irow_point[assembly.stop[ielem]]

    @views jacob[irow2:irow2+2, icol1:icol1+2] .+= F2_u1*u1_u1 ./ force_scaling
    @views jacob[irow2:irow2+2, icol2:icol2+2] .+= F2_u2*u2_u2 ./ force_scaling

    @views jacob[irow2:irow2+2, icol1+6:icol1+8] .+= F2_V1 ./ force_scaling
    @views jacob[irow2:irow2+2, icol1+9:icol1+11] .+= F2_Ω1 ./ force_scaling

    @views jacob[irow2:irow2+2, icol2+6:icol2+8] .+= F2_V2 ./ force_scaling
    @views jacob[irow2:irow2+2, icol2+9:icol2+11] .+= F2_Ω2 ./ force_scaling

    @views jacob[irow2+3:irow2+5, icol1:icol1+2] .+= M2_u1*u1_u1 ./ force_scaling
    @views jacob[irow2+3:irow2+5, icol2:icol2+2] .+= M2_u2*u2_u2 ./ force_scaling

    @views jacob[irow2+3:irow2+5, icol1+6:icol1+8] .+= M2_V1 ./ force_scaling
    @views jacob[irow2+3:irow2+5, icol1+9:icol1+11] .+= M2_Ω1 ./ force_scaling

    @views jacob[irow2+3:irow2+5, icol2+6:icol2+8] .+= M2_V2 ./ force_scaling
    @views jacob[irow2+3:irow2+5, icol2+9:icol2+11] .+= M2_Ω2 ./ force_scaling

    return jacob
end

@inline function insert_expanded_steady_element_jacobians!(jacob, indices, force_scaling,
    assembly, ielem, properties, compatibility, velocities, equilibrium, resultants)

    insert_expanded_dynamic_element_jacobians!(jacob, indices, force_scaling,
        assembly, ielem, properties, compatibility, velocities, equilibrium, resultants)

    @unpack rF_ab, rF_αb, rM_ab, rM_αb = equilibrium

    # NOTE: We have to switch the order of the equations here in order to match the indices
    # of the differential variables with their equations.  This is done for compatibility
    # with the DiffEqSensitivity package.

    # jacobians wrt prescribed linear and angular acceleration
    irow = indices.irow_elem[ielem]
    icol = indices.icol_body

    for i = 1:3
        if !iszero(icol[i])
            jacob[irow+12:irow+14, icol[i]] .= rF_ab[:,i] ./ force_scaling
            jacob[irow+15:irow+17, icol[i]] .= rM_ab[:,i] ./ force_scaling
        end
    end

    for i = 4:6
        if !iszero(icol[i])
            jacob[irow+12:irow+14, icol[i]] .= rF_αb[:,i-3] ./ force_scaling
            jacob[irow+15:irow+17, icol[i]] .= rM_αb[:,i-3] ./ force_scaling
        end
    end

    return jacob
end

@inline function insert_expanded_dynamic_element_jacobians!(jacob, indices, force_scaling,
    assembly, ielem, properties, compatibility, velocities, equilibrium, resultants)

    @unpack u1_u1, u2_u2, θ1_θ1, θ2_θ2 = properties

    @unpack ru_u1, ru_u2, ru_θ1, ru_θ2, ru_V1, ru_V2, ru_Ω, ru_F1, ru_F2, ru_M1, ru_M2, 
            rθ_θ1, rθ_θ2, rθ_Ω1, rθ_Ω2, rθ_F1, rθ_F2, rθ_M1, rθ_M2 = compatibility
    
    @unpack rV_u1, rV_u2, rV_θ1, rV_θ2, rV_V1, rV_V2, rV_V, rΩ_θ1, rΩ_θ2, rΩ_Ω1, rΩ_Ω2, rΩ_Ω = velocities

    @unpack rF_u1, rF_u2, rF_θ1, rF_θ2, rF_F1, rF_F2, rF_V, rF_Ω,  
            rM_u1, rM_u2, rM_θ1, rM_θ2, rM_F1, rM_F2, rM_M1, rM_M2, rM_V1, rM_V2, rM_V, rM_Ω = equilibrium
    
    @unpack F1_θ1, F1_θ2, F1_F1, F2_θ1, F2_θ2, F2_F2, 
            M1_θ1, M1_θ2, M1_M1, M2_θ1, M2_θ2, M2_M2 = resultants

    icol = indices.icol_elem[ielem]
    icol1 = indices.icol_point[assembly.start[ielem]]
    icol2 = indices.icol_point[assembly.stop[ielem]]

    irow = indices.irow_elem[ielem]

    # NOTE: We have to switch the order of the equations here in order to match the indices
    # of the differential variables with their equations.  This is done for compatibility
    # with the DiffEqSensitivity package.

    # element compatibility residuals

    jacob[irow:irow+2, icol1:icol1+2] .= ru_u1*u1_u1
    jacob[irow:irow+2, icol2:icol2+2] .= ru_u2*u2_u2
    jacob[irow:irow+2, icol1+3:icol1+5] .= ru_θ1*θ1_θ1
    jacob[irow:irow+2, icol2+3:icol2+5] .= ru_θ2*θ2_θ2
    jacob[irow:irow+2, icol1+6:icol1+8] .= ru_V1
    jacob[irow:irow+2, icol2+6:icol2+8] .= ru_V2
    jacob[irow:irow+2, icol:icol+2] .= ru_F1 .* force_scaling
    jacob[irow:irow+2, icol+3:icol+5] .= ru_M1 .* force_scaling
    jacob[irow:irow+2, icol+6:icol+8] .= ru_F2 .* force_scaling
    jacob[irow:irow+2, icol+9:icol+11] .= ru_M2 .* force_scaling
    jacob[irow:irow+2, icol+15:icol+17] .= ru_Ω

    jacob[irow+3:irow+5, icol1+3:icol1+5] .= rθ_θ1*θ1_θ1
    jacob[irow+3:irow+5, icol2+3:icol2+5] .= rθ_θ2*θ2_θ2
    jacob[irow+3:irow+5, icol1+9:icol1+11] .= rθ_Ω1
    jacob[irow+3:irow+5, icol2+9:icol2+11] .= rθ_Ω2
    jacob[irow+3:irow+5, icol:icol+2] .= rθ_F1 .* force_scaling
    jacob[irow+3:irow+5, icol+3:icol+5] .= rθ_M1 .* force_scaling
    jacob[irow+3:irow+5, icol+6:icol+8] .= rθ_F2 .* force_scaling
    jacob[irow+3:irow+5, icol+9:icol+11] .= rθ_M2 .* force_scaling

    # velocity residuals

    jacob[irow+6:irow+8, icol1:icol1+2] .= rV_u1*u1_u1
    jacob[irow+6:irow+8, icol2:icol2+2] .= rV_u2*u2_u2
    jacob[irow+6:irow+8, icol1+3:icol1+5] .= rV_θ1*θ1_θ1
    jacob[irow+6:irow+8, icol2+3:icol2+5] .= rV_θ2*θ2_θ2
    jacob[irow+6:irow+8, icol1+6:icol1+8] .= rV_V1
    jacob[irow+6:irow+8, icol2+6:icol2+8] .= rV_V2
    jacob[irow+6:irow+8, icol+12:icol+14] .= rV_V

    jacob[irow+9:irow+11, icol1+3:icol1+5] .= rΩ_θ1*θ1_θ1
    jacob[irow+9:irow+11, icol2+3:icol2+5] .= rΩ_θ2*θ2_θ2
    jacob[irow+9:irow+11, icol1+9:icol1+11] .= rΩ_Ω1
    jacob[irow+9:irow+11, icol2+9:icol2+11] .= rΩ_Ω2
    jacob[irow+9:irow+11, icol+15:icol+17] .= rΩ_Ω

    # element equilibrium residuals
    jacob[irow+12:irow+14, icol1:icol1+2] .= rF_u1*u1_u1 ./ force_scaling
    jacob[irow+12:irow+14, icol2:icol2+2] .= rF_u2*u2_u2 ./ force_scaling
    jacob[irow+12:irow+14, icol1+3:icol1+5] .= rF_θ1*θ1_θ1 ./ force_scaling
    jacob[irow+12:irow+14, icol2+3:icol2+5] .= rF_θ2*θ2_θ2 ./ force_scaling
    jacob[irow+12:irow+14, icol:icol+2] .= rF_F1
    jacob[irow+12:irow+14, icol+6:icol+8] .= rF_F2
    jacob[irow+12:irow+14, icol+12:icol+14] .= rF_V ./ force_scaling
    jacob[irow+12:irow+14, icol+15:icol+17] .= rF_Ω ./ force_scaling

    jacob[irow+15:irow+17, icol1:icol1+2] .= rM_u1*u1_u1 ./ force_scaling
    jacob[irow+15:irow+17, icol2:icol2+2] .= rM_u2*u2_u2 ./ force_scaling
    jacob[irow+15:irow+17, icol1+3:icol1+5] .= rM_θ1*θ1_θ1 ./ force_scaling
    jacob[irow+15:irow+17, icol2+3:icol2+5] .= rM_θ2*θ2_θ2 ./ force_scaling
    jacob[irow+15:irow+17, icol1+6:icol1+8] .= rM_V1 ./ force_scaling
    jacob[irow+15:irow+17, icol2+6:icol2+8] .= rM_V2 ./ force_scaling
    jacob[irow+15:irow+17, icol:icol+2] .= rM_F1
    jacob[irow+15:irow+17, icol+3:icol+5] .= rM_M1
    jacob[irow+15:irow+17, icol+6:icol+8] .= rM_F2
    jacob[irow+15:irow+17, icol+9:icol+11] .= rM_M2
    jacob[irow+15:irow+17, icol+12:icol+14] .= rM_V ./ force_scaling
    jacob[irow+15:irow+17, icol+15:icol+17] .= rM_Ω ./ force_scaling

    # equilibrium equations for the start of the beam element
    irow = indices.irow_point[assembly.start[ielem]]
    @views jacob[irow+6:irow+8, icol1+3:icol1+5] .-= F1_θ1*θ1_θ1 ./ force_scaling
    @views jacob[irow+6:irow+8, icol2+3:icol2+5] .-= F1_θ2*θ2_θ2 ./ force_scaling
    @views jacob[irow+6:irow+8, icol:icol+2] .-= F1_F1
    @views jacob[irow+9:irow+11, icol1+3:icol1+5] .-= M1_θ1*θ1_θ1 ./ force_scaling
    @views jacob[irow+9:irow+11, icol2+3:icol2+5] .-= M1_θ2*θ2_θ2 ./ force_scaling
    @views jacob[irow+9:irow+11, icol+3:icol+5] .-= M1_M1

    # equilibrium equations for the end of the beam element
    irow = indices.irow_point[assembly.stop[ielem]]
    @views jacob[irow+6:irow+8, icol1+3:icol1+5] .+= F2_θ1*θ1_θ1 ./ force_scaling
    @views jacob[irow+6:irow+8, icol2+3:icol2+5] .+= F2_θ2*θ2_θ2 ./ force_scaling
    @views jacob[irow+6:irow+8, icol+6:icol+8] .+= F2_F2
    @views jacob[irow+9:irow+11, icol1+3:icol1+5] .+= M2_θ1*θ1_θ1 ./ force_scaling
    @views jacob[irow+9:irow+11, icol2+3:icol2+5] .+= M2_θ2*θ2_θ2 ./ force_scaling
    @views jacob[irow+9:irow+11, icol+9:icol+11] .+= M2_M2

    return jacob
end

@inline function insert_mass_matrix_element_jacobians!(jacob, gamma, indices, two_dimensional, force_scaling,
    assembly, ielem, resultants)

    @unpack F1_V1dot, F1_V2dot, F1_Ω1dot, F1_Ω2dot, 
            F2_V1dot, F2_V2dot, F2_Ω1dot, F2_Ω2dot, 
            M1_V1dot, M1_V2dot, M1_Ω1dot, M1_Ω2dot, 
            M2_V1dot, M2_V2dot, M2_Ω1dot, M2_Ω2dot = resultants

    icol1 = indices.icol_point[assembly.start[ielem]]
    icol2 = indices.icol_point[assembly.stop[ielem]]
    
    if two_dimensional
        lmask = SVector(1, 1, 0)
        amask = SVector(0, 0, 1)
    else
        lmask = SVector(1, 1, 1)
        amask = SVector(1, 1, 1)
    end

    # equilibrium equations for the start of the beam element
    irow1 = indices.irow_point[assembly.start[ielem]]

    @views jacob[irow1:irow1+2, icol1+6:icol1+8] .-= lmask .* F1_V1dot .* gamma ./ force_scaling
    @views jacob[irow1:irow1+2, icol1+9:icol1+11] .-= lmask .* F1_Ω1dot .* gamma ./ force_scaling

    @views jacob[irow1:irow1+2, icol2+6:icol2+8] .-= lmask .* F1_V2dot .* gamma ./ force_scaling
    @views jacob[irow1:irow1+2, icol2+9:icol2+11] .-= lmask .* F1_Ω2dot .* gamma ./ force_scaling

    @views jacob[irow1+3:irow1+5, icol1+6:icol1+8] .-= amask .* M1_V1dot .* gamma ./ force_scaling
    @views jacob[irow1+3:irow1+5, icol1+9:icol1+11] .-= amask .* M1_Ω1dot .* gamma ./ force_scaling

    @views jacob[irow1+3:irow1+5, icol2+6:icol2+8] .-= amask .* M1_V2dot .* gamma ./ force_scaling
    @views jacob[irow1+3:irow1+5, icol2+9:icol2+11] .-= amask .* M1_Ω2dot .* gamma ./ force_scaling

    # equilibrium equations for the end of the beam element
    irow2 = indices.irow_point[assembly.stop[ielem]]

    @views jacob[irow2:irow2+2, icol1+6:icol1+8] .+= lmask .* F2_V1dot .* gamma ./ force_scaling
    @views jacob[irow2:irow2+2, icol1+9:icol1+11] .+= lmask .* F2_Ω1dot .* gamma ./ force_scaling

    @views jacob[irow2:irow2+2, icol2+6:icol2+8] .+= lmask .* F2_V2dot .* gamma ./ force_scaling
    @views jacob[irow2:irow2+2, icol2+9:icol2+11] .+= lmask .* F2_Ω2dot .* gamma ./ force_scaling

    @views jacob[irow2+3:irow2+5, icol1+6:icol1+8] .+= amask .* M2_V1dot .* gamma ./ force_scaling
    @views jacob[irow2+3:irow2+5, icol1+9:icol1+11] .+= amask .* M2_Ω1dot .* gamma ./ force_scaling

    @views jacob[irow2+3:irow2+5, icol2+6:icol2+8] .+= amask .* M2_V2dot .* gamma ./ force_scaling
    @views jacob[irow2+3:irow2+5, icol2+9:icol2+11] .+= amask .* M2_Ω2dot .* gamma ./ force_scaling

    return jacob
end

@inline function insert_expanded_mass_matrix_element_jacobians!(jacob, gamma, indices, 
    two_dimensional, force_scaling, assembly, ielem, equilibrium)

    @unpack rF_Vdot, rF_Ωdot, rM_Vdot, rM_Ωdot = equilibrium

    irow = indices.irow_elem[ielem] 
    icol = indices.icol_elem[ielem]

    if two_dimensional
        lmask = SVector(1, 1, 0)
        amask = SVector(0, 0, 1)
    else
        lmask = SVector(1, 1, 1)
        amask = SVector(1, 1, 1)
    end

    # NOTE: We have to switch the order of the equations here in order to match the indices
    # of the differential variables with their equations.  This is done for compatibility
    # with the DiffEqSensitivity package.

    # equilibrium residuals
    @views jacob[irow+12:irow+14, icol+12:icol+14] .+= lmask .* rF_Vdot .* gamma ./ force_scaling
    @views jacob[irow+12:irow+14, icol+15:icol+17] .+= lmask .* rF_Ωdot .* gamma ./ force_scaling

    @views jacob[irow+15:irow+17, icol+12:icol+14] .+= amask .* rM_Vdot .* gamma ./ force_scaling
    @views jacob[irow+15:irow+17, icol+15:icol+17] .+= amask .* rM_Ωdot .* gamma ./ force_scaling

    return jacob
end


# --- Element Residual --- #

"""
    static_element_residual!(resid, x, indices, force_scaling, assembly, ielem,  
        prescribed_conditions, distributed_loads, gravity)

Calculate and insert the residual entries corresponding to a beam element for a static 
analysis into the system residual vector.
"""
@inline function static_element_residual!(resid, x, indices, force_scaling, assembly, ielem,  
    prescribed_conditions, distributed_loads, gravity)

    properties = static_element_properties(x, indices, force_scaling, assembly, ielem, 
        prescribed_conditions, gravity)

    compatibility = compatibility_residuals(properties)

    resultants = static_element_resultants(properties, distributed_loads, ielem)

    insert_element_residuals!(resid, indices, force_scaling, assembly, ielem, 
        compatibility, resultants)

    return resid
end

"""
    steady_element_residual!(resid, x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, distributed_loads, gravity, 
        linear_velocity, angular_velocity, linear_acceleration, angular_acceleration)

Calculate and insert the residual entries corresponding to a beam element for a steady state 
analysis into the system residual vector.
"""
@inline function steady_element_residual!(resid, x, indices, force_scaling, 
    structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, gravity, 
    linear_velocity, angular_velocity, linear_acceleration, angular_acceleration)

    properties = steady_element_properties(x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity, 
        linear_velocity, angular_velocity, linear_acceleration, angular_acceleration)

    compatibility = compatibility_residuals(properties)

    resultants = dynamic_element_resultants(properties, distributed_loads, ielem)

    insert_element_residuals!(resid, indices, force_scaling, assembly, ielem, 
        compatibility, resultants)

    return resid
end

"""
    initial_element_residual!(resid, x, indices, rate_vars, 
        force_scaling, structural_damping, assembly, ielem, prescribed_conditions, 
        distributed_loads, gravity, linear_velocity, angular_velocity, 
        linear_acceleration, angular_acceleration, u0, θ0, V0, Ω0, Vdot0, Ωdot0)

Calculate and insert the residual entries corresponding to a beam element for the 
initialization of a time domain simulation into the system residual vector.
"""
@inline function initial_element_residual!(resid, x, indices, rate_vars, 
    force_scaling, structural_damping, assembly, ielem, prescribed_conditions, 
    distributed_loads, gravity, linear_velocity, angular_velocity, 
    linear_acceleration, angular_acceleration, u0, θ0, V0, Ω0, Vdot0, Ωdot0)

    properties = initial_element_properties(x, indices, rate_vars, 
        force_scaling, structural_damping, assembly, ielem, prescribed_conditions, gravity, 
        linear_velocity, angular_velocity, linear_acceleration, angular_acceleration, 
        u0, θ0, V0, Ω0, Vdot0, Ωdot0)

    compatibility = compatibility_residuals(properties)

    resultants = dynamic_element_resultants(properties, distributed_loads, ielem)

    insert_element_residuals!(resid, indices, force_scaling, assembly, ielem, 
        compatibility, resultants)

    return resid
end

"""
    newmark_element_residual!(resid, x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, distributed_loads, gravity, 
        linear_velocity, angular_velocity, Vdot_init, Ωdot_init, dt)

Calculate and insert the residual entries corresponding to a beam element for a 
newmark-scheme time marching analysis into the system residual vector.
"""
@inline function newmark_element_residual!(resid, x, indices, force_scaling, structural_damping, 
    assembly, ielem, prescribed_conditions, distributed_loads, gravity, 
    linear_velocity, angular_velocity, Vdot_init, Ωdot_init, dt)

    properties = newmark_element_properties(x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, gravity, linear_velocity, angular_velocity, 
        Vdot_init, Ωdot_init, dt)

    compatibility = compatibility_residuals(properties)

    resultants = dynamic_element_resultants(properties, distributed_loads, ielem)

    insert_element_residuals!(resid, indices, force_scaling, assembly, ielem, 
        compatibility, resultants)

    return resid
end

"""
    dynamic_element_residual!(resid, dx, x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, distributed_loads, gravity, 
        linear_velocity, angular_velocity)

Calculate and insert the residual entries corresponding to a beam element for a dynamic
analysis into the system residual vector.
"""
@inline function dynamic_element_residual!(resid, dx, x, indices, force_scaling, structural_damping, 
    assembly, ielem, prescribed_conditions, distributed_loads, gravity, linear_velocity, angular_velocity)

    properties = dynamic_element_properties(dx, x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, gravity, linear_velocity, angular_velocity)

    compatibility = compatibility_residuals(properties)

    resultants = dynamic_element_resultants(properties, distributed_loads, ielem)

    insert_element_residuals!(resid, indices, force_scaling, assembly, ielem, 
        compatibility, resultants)

    return resid
end

"""
    expanded_steady_element_residual!(resid, x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, distributed_loads, gravity,
        linear_velocity, angular_velocity, linear_acceleration, angular_acceleration)

Calculate and insert the residual entries corresponding to a beam element for a constant
mass matrix system into the system residual vector.
"""
@inline function expanded_steady_element_residual!(resid, x, indices, force_scaling, 
    structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, gravity, 
    linear_velocity, angular_velocity, linear_acceleration, angular_acceleration)

    properties = expanded_steady_element_properties(x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, gravity, linear_velocity, angular_velocity, 
        linear_acceleration, angular_acceleration)

    compatibility = compatibility_residuals(properties)

    resultants = expanded_element_resultants(properties)

    velocities = expanded_element_velocity_residuals(properties)

    equilibrium = expanded_element_equilibrium_residuals(properties, distributed_loads, ielem)

    insert_expanded_element_residuals!(resid, indices, force_scaling, assembly, ielem, 
        compatibility, velocities, equilibrium, resultants)

    return resid
end

"""
    expanded_dynamic_element_residual!(resid, dx, x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, 
        gravity, linear_velocity, angular_velocity)

Calculate and insert the residual entries corresponding to a beam element for a constant
mass matrix system into the system residual vector.
"""
@inline function expanded_dynamic_element_residual!(resid, dx, x, indices, force_scaling, 
    structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, 
    gravity, linear_velocity, angular_velocity)

    properties = expanded_dynamic_element_properties(dx, x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, gravity, linear_velocity, angular_velocity)

    compatibility = compatibility_residuals(properties)

    resultants = expanded_element_resultants(properties)

    velocities = expanded_element_velocity_residuals(properties)

    equilibrium = expanded_element_equilibrium_residuals(properties, distributed_loads, ielem)

    insert_expanded_element_residuals!(resid, indices, force_scaling, assembly, ielem, 
        compatibility, velocities, equilibrium, resultants)

    return resid
end

"""
    static_element_jacobian!(jacob, x, indices, force_scaling, 
        assembly, ielem, prescribed_conditions, distributed_loads, gravity)

Calculate and insert the jacobian entries corresponding to a beam element for a static 
analysis into the system jacobian matrix.
"""
@inline function static_element_jacobian!(jacob, x, indices, force_scaling, 
    assembly, ielem, prescribed_conditions, distributed_loads, gravity)

    properties = static_element_properties(x, indices, force_scaling, assembly, ielem, 
        prescribed_conditions, gravity)

    properties = static_element_jacobian_properties(properties, x, indices, force_scaling, 
        assembly, ielem, prescribed_conditions, gravity)

    compatibility = static_compatibility_jacobians(properties)

    resultants = static_element_resultant_jacobians(properties, distributed_loads, ielem)

    insert_static_element_jacobians!(jacob, indices, force_scaling, assembly, ielem, 
        properties, compatibility, resultants)

    return jacob
end

"""
    steady_element_jacobian!(jacob, x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, 
        gravity, linear_velocity, angular_velocity, linear_acceleration, angular_acceleration)

Calculate and insert the jacobian entries corresponding to a beam element for a steady state 
analysis into the system jacobian matrix.
"""
@inline function steady_element_jacobian!(jacob, x, indices, 
    force_scaling, structural_damping, assembly, ielem, prescribed_conditions, 
    distributed_loads, gravity, linear_velocity, angular_velocity, linear_acceleration, angular_acceleration)

    properties = steady_element_properties(x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity, 
        linear_velocity, angular_velocity, linear_acceleration, angular_acceleration)

    properties = steady_element_jacobian_properties(properties, x, indices, 
        force_scaling, structural_damping, assembly, ielem, 
        prescribed_conditions, gravity)

    compatibility = dynamic_compatibility_jacobians(properties)

    resultants = steady_element_resultant_jacobians(properties, distributed_loads, ielem)

    insert_steady_element_jacobians!(jacob, indices, force_scaling, assembly, ielem, 
        properties, compatibility, resultants)

    return jacob
end

"""
    initial_element_jacobian!(jacob, x, indices, rate_vars, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, 
        gravity, ub, θb, vb, ωb, ab, αb, u0, θ0, V0, Ω0, Vdot0, Ωdot0)

Calculate and insert the jacobian entries corresponding to a beam element for the 
initialization of a time domain analysis into the system jacobian matrix.
"""
@inline function initial_element_jacobian!(jacob, x, indices, rate_vars, force_scaling, 
    structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, gravity, 
    linear_velocity, angular_velocity, linear_acceleration, angular_acceleration,
    u0, θ0, V0, Ω0, Vdot0, Ωdot0)

    properties = initial_element_properties(x, indices, rate_vars, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity, 
        linear_velocity, angular_velocity, linear_acceleration, angular_acceleration, 
        u0, θ0, V0, Ω0, Vdot0, Ωdot0)

    properties = initial_element_jacobian_properties(properties, x, indices, 
        rate_vars, force_scaling, structural_damping, assembly, ielem, prescribed_conditions, 
        gravity, u0, θ0, V0, Ω0, Vdot0, Ωdot0)

    compatibility = initial_compatibility_jacobians(properties)

    resultants = initial_element_resultant_jacobians(properties, distributed_loads, ielem)

    insert_initial_element_jacobians!(jacob, indices, force_scaling, assembly, ielem, 
        properties, compatibility, resultants)

    return jacob
end

"""
    newmark_element_jacobian!(jacob, x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, 
        gravity, linear_velocity, angular_velocity, Vdot_init, Ωdot_init, dt)

Calculate and insert the jacobian entries corresponding to a beam element for a Newmark-scheme
time marching analysis into the system jacobian matrix.
"""
@inline function newmark_element_jacobian!(jacob, x, indices, force_scaling, 
    structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, gravity, 
    linear_velocity, angular_velocity, Vdot_init, Ωdot_init, dt)

    properties = newmark_element_properties(x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity, 
        linear_velocity, angular_velocity, Vdot_init, Ωdot_init, dt)

    properties = newmark_element_jacobian_properties(properties, x, indices, 
        force_scaling, structural_damping, assembly, ielem, prescribed_conditions, gravity, 
        Vdot_init, Ωdot_init, dt)

    compatibility = dynamic_compatibility_jacobians(properties)

    resultants = dynamic_element_resultant_jacobians(properties, distributed_loads, ielem)

    insert_dynamic_element_jacobians!(jacob, indices, force_scaling, assembly, ielem, 
        properties, compatibility, resultants)

    return jacob
end


"""
    dynamic_element_jacobian!(jacob, dx, x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, 
        gravity, linear_velocity, angular_velocity)

Calculate and insert the jacobian entries corresponding to a beam element for a dynamic
analysis into the system jacobian matrix.
"""
@inline function dynamic_element_jacobian!(jacob, dx, x, indices, force_scaling, 
    structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, gravity, 
    linear_velocity, angular_velocity)

    properties = dynamic_element_properties(dx, x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity, 
        linear_velocity, angular_velocity)

    properties = dynamic_element_jacobian_properties(properties, dx, x, indices, 
        force_scaling, structural_damping, assembly, ielem, prescribed_conditions, gravity)

    compatibility = dynamic_compatibility_jacobians(properties)

    resultants = dynamic_element_resultant_jacobians(properties, distributed_loads, ielem)

    insert_dynamic_element_jacobians!(jacob, indices, force_scaling, assembly, ielem, 
        properties, compatibility, resultants)

    return jacob
end

"""
    expanded_steady_element_jacobian!(jacob, x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, 
        gravity, linear_velocity, angular_velocity, linear_acceleration, 
        angular_acceleration)

Calculate and insert the jacobian entries corresponding to a beam element for a dynamic
analysis into the system jacobian matrix.
"""
@inline function expanded_steady_element_jacobian!(jacob, x, indices, force_scaling, 
    structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, gravity, 
    linear_velocity, angular_velocity, linear_acceleration, angular_acceleration)

    properties = expanded_steady_element_properties(x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity, 
        linear_velocity, angular_velocity, linear_acceleration, angular_acceleration)

    properties = expanded_steady_element_jacobian_properties(properties, x, indices, 
        force_scaling, structural_damping, assembly, ielem, prescribed_conditions, gravity)

    compatibility = expanded_compatibility_jacobians(properties)

    velocities = expanded_element_velocity_jacobians(properties)

    equilibrium = expanded_steady_element_equilibrium_jacobians(properties, distributed_loads, ielem)

    resultants = expanded_element_resultant_jacobians(properties)

    insert_expanded_steady_element_jacobians!(jacob, indices, force_scaling, assembly, ielem, 
        properties, compatibility, velocities, equilibrium, resultants)

    return jacob
end

"""
    expanded_dynamic_element_jacobian!(jacob, dx, x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, 
        gravity, linear_velocity, angular_velocity)

Calculate and insert the jacobian entries corresponding to a beam element for a dynamic
analysis into the system jacobian matrix.
"""
@inline function expanded_dynamic_element_jacobian!(jacob, dx, x, indices, force_scaling, structural_damping, 
    assembly, ielem, prescribed_conditions, distributed_loads, gravity, linear_velocity, angular_velocity)

    properties = expanded_dynamic_element_properties(dx, x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity, linear_velocity, angular_velocity)

    properties = expanded_dynamic_element_jacobian_properties(properties, dx, x, indices, 
        force_scaling, structural_damping, assembly, ielem, prescribed_conditions, gravity)

    compatibility = expanded_compatibility_jacobians(properties)

    velocities = expanded_element_velocity_jacobians(properties)

    equilibrium = expanded_dynamic_element_equilibrium_jacobians(properties, distributed_loads, ielem)

    resultants = expanded_element_resultant_jacobians(properties)

    insert_expanded_dynamic_element_jacobians!(jacob, indices, force_scaling, assembly, ielem, 
        properties, compatibility, velocities, equilibrium, resultants)

    return jacob
end

"""
    mass_matrix_element_jacobian!(jacob, gamma, x, indices, two_dimensional, force_scaling, assembly, 
        ielem, prescribed_conditions)

Calculate and insert the mass_matrix jacobian entries corresponding to a beam element into 
the system jacobian matrix.
"""
@inline function mass_matrix_element_jacobian!(jacob, gamma, x, indices, two_dimensional, force_scaling, assembly, 
    ielem, prescribed_conditions)

    properties = mass_matrix_element_jacobian_properties(x, indices, force_scaling, 
        assembly, ielem, prescribed_conditions)

    resultants = mass_matrix_element_resultant_jacobians(properties)
    
    insert_mass_matrix_element_jacobians!(jacob, gamma, indices, two_dimensional, force_scaling, assembly, ielem, 
        resultants)

    return jacob
end

"""
    expanded_mass_matrix_element_jacobian!(jacob, gamma, indices, two_dimensional, force_scaling, assembly, 
        ielem, prescribed_conditions)

Calculate and insert the mass_matrix jacobian entries corresponding to a beam element into 
the system jacobian matrix for a constant mass matrix system
"""
@inline function expanded_mass_matrix_element_jacobian!(jacob, gamma, indices, two_dimensional, force_scaling, assembly, 
    ielem, prescribed_conditions)

    properties = expanded_mass_matrix_element_jacobian_properties(assembly, ielem, prescribed_conditions)

    equilibrium = expanded_mass_matrix_element_equilibrium_jacobians(properties)
    
    insert_expanded_mass_matrix_element_jacobians!(jacob, gamma, indices, two_dimensional, force_scaling, 
        assembly, ielem, equilibrium)

    return jacob
end