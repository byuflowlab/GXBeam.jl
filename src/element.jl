"""
    Element{TF}

Composite type that defines a beam element's properties

# Fields
 - `L`: Beam element length
 - `x`: Beam element location
 - `compliance`: Beam element compliance matrix
 - `mass`: Beam element mass matrix
 - `Cab`: Transformation matrix from the undeformed beam element frame to the body frame
 - `mu`: Beam element damping coefficients
"""
struct Element{TF}
    L::TF
    x::SVector{3,TF}
    compliance::SMatrix{6,6,TF,36}
    mass::SMatrix{6,6,TF,36}
    Cab::SMatrix{3,3,TF,9}
    mu::SVector{6,TF}
end

"""
    Element(L, x, compliance, mass, Cab, mu)

Construct a beam element

# Arguments
- `L`: Length of the beam element
- `x`: Location of the beam element (the center of the beam element)
- `compliance`: Beam element compliance matrix
- `mass`: Beam element mass matrix
- `Cab`: Transformation matrix from the undeformed beam element frame to the body frame
- `mu`: Beam element damping coefficients
"""
function Element(L, x, compliance, mass, Cab, mu)
    TF = promote_type(typeof(L), eltype(x), eltype(compliance), eltype(mass), eltype(Cab), eltype(mu))
    return Element{TF}(L, x, compliance, mass, Cab, mu)
end

"""
    element_loads(x, ielem, icol_elem, force_scaling)

Extract the internal loads (`F`, `M`) for a beam element from the state variable vector.  
These loads are expressed in the deformed element frame.
"""
@inline function element_loads(x, ielem, icol_elem, force_scaling)

    icol = icol_elem[ielem]

    F = SVector(x[icol  ], x[icol+1], x[icol+2]) .* force_scaling
    M = SVector(x[icol+3], x[icol+4], x[icol+5]) .* force_scaling

    return F, M
end

"""
    expanded_element_loads(x, ielem, icol_elem, force_scaling)

Extract the internal loads of a beam element at its endpoints (`F1`, `M1`, `F2`, `M2`) 
from the state variable vector for a constant mass matrix system.
"""
@inline function expanded_element_loads(x, ielem, icol_elem, force_scaling)

    icol = icol_elem[ielem]

    F1 = SVector(x[icol  ], x[icol+1], x[icol+2]) .* force_scaling
    M1 = SVector(x[icol+3], x[icol+4], x[icol+5]) .* force_scaling
    F2 = SVector(x[icol+6], x[icol+7], x[icol+8]) .* force_scaling
    M2 = SVector(x[icol+9], x[icol+10], x[icol+11]) .* force_scaling

    return F1, M1, F2, M2
end

"""
    expanded_element_velocities(x, ielem, icol_elem, force_scaling)

Extract the velocities of a beam element from the state variable vector for a constant mass
matrix system.
"""
@inline function expanded_element_velocities(x, ielem, icol_elem)

    icol = icol_elem[ielem]

    V = SVector(x[icol+12], x[icol+13], x[icol+14])
    Ω = SVector(x[icol+15], x[icol+16], x[icol+17])

    return V, Ω
end

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

    # linear and angular acceleration
    a = -SVector{3}(gravity)
    α = @SVector zeros(3)

    return (; L, C, Cab, CtCab, Qinv, S11, S12, S21, S22, mass11, mass12, mass21, mass22, 
        u1, u2, θ1, θ2, u, θ, F, M, γ, κ, a, α)
end

"""
    steady_state_element_properties(x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity)

Calculate/extract the element properties needed to construct the residual for a steady 
state analysis
"""
@inline function steady_state_element_properties(x, indices, force_scaling, 
    structural_damping, assembly, ielem, prescribed_conditions, gravity)

    properties = static_element_properties(x, indices, force_scaling, assembly, ielem, 
        prescribed_conditions, gravity)

    @unpack L, Cab, C, CtCab, Qinv, mass11, mass12, mass21, mass22, u1, u2, θ1, θ2, 
        u, θ, γ, κ = properties

    # rotation parameter matrices
    Q = get_Q(θ)

    # linear and angular displacement of the body frame
    ub, θb = body_frame_displacement(x)

    # linear and angular velocity of the body frame
    vb, ωb = body_frame_velocity(x)

    # linear and angular acceleration of the body frame
    ab, αb = body_frame_acceleration(x)

    # distance from the rotation center
    Δx = assembly.elements[ielem].x

    # linear and angular velocity
    v = vb + cross(ωb, Δx) + cross(ωb, u)# + udot = V
    ω = ωb# + Cab'*Q*θdot = Ω

    # linear and angular acceleration
    a = ab + cross(αb, Δx) + cross(αb, u)# + cross(ω, V) + uddot = Vdot
    α = αb# + Cab'*Qdot*θdot + Cab'*Q*θddot = Ωdot

    # add gravitational acceleration
    a -= get_C(θb)*gravity

    # linear and angular velocity
    V1, Ω1 = point_velocities(x, assembly.start[ielem], indices.icol_point)
    V2, Ω2 = point_velocities(x, assembly.stop[ielem], indices.icol_point)
    V = (V1 + V2)/2
    Ω = (Ω1 + Ω2)/2

    # linear and angular momentum
    P = CtCab*mass11*CtCab'*V + CtCab*mass12*CtCab'*Ω
    H = CtCab*mass21*CtCab'*V + CtCab*mass22*CtCab'*Ω

    # save steady state properties
    properties = (; properties..., Q, V1, V2, Ω1, Ω2, V, Ω, P, H, ub, θb, vb, ωb, ab, αb, 
        Δx, v, ω, a, α)

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

        # linear and angular velocity
        v1 = vb + cross(ωb, Δx1) + cross(ωb, u1)
        v2 = vb + cross(ωb, Δx2) + cross(ωb, u2)
        ω1 = ωb
        ω2 = ωb

        # linear displacement rates 
        udot1 = V1 - v1
        udot2 = V2 - v2
        udot = (udot1 + udot2)/2

        # angular displacement rates
        θdot1 = Qinv1*C1*(Ω1 - ω1)
        θdot2 = Qinv2*C2*(Ω2 - ω2)
        θdot = (θdot1 + θdot2)/2

        # change in linear and angular displacement
        Δu = u2 - u1
        Δθ = θ2 - θ1

        # change in linear and angular displacement rates
        Δudot = udot2 - udot1
        Δθdot = θdot2 - θdot1

        # ΔQ matrix (see structural damping theory)
        ΔQ = get_ΔQ(θ, Δθ, Q)

        # strain rates
        γdot = -CtCab'*tilde(Ω - ω)*Δu + CtCab'*Δudot - L*CtCab'*tilde(Ω - ω)*Cab*e1
        κdot = Cab'*Q*Δθdot + Cab'*ΔQ*θdot

        # adjust strains to account for strain rates
        γ -= μ11*γdot
        κ -= μ22*κdot

        # save new strains and structural damping properties
        properties = (; properties..., γ, κ, γdot, κdot,
            μ11, μ22, C1, C2, Qinv1, Qinv2, Δx1, Δx2, v1, v2, ω1, ω2, 
            udot1, udot2, udot, θdot1, θdot2, θdot, Δu, Δθ, Δudot, Δθdot, ΔQ)
    end

    return properties
end

"""
    initial_condition_element_properties(x, indices, rate_vars, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity, 
        u0, θ0, V0, Ω0, Vdot0, Ωdot0)

Calculate/extract the element properties needed to construct the residual for a time-domain
analysis initialization
"""
@inline function initial_condition_element_properties(x, indices, rate_vars, 
    force_scaling, structural_damping, assembly, ielem, prescribed_conditions, gravity, 
    u0, θ0, V0, Ω0, Vdot0, Ωdot0)

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

    # linear and angular displacement of the body frame
    ub, θb = body_frame_displacement(x)

    # linear and angular velocity of the body frame
    vb, ωb = body_frame_velocity(x)

    # linear and angular acceleration of the body frame
    ab, αb = body_frame_acceleration(x)

    # distance from the rotation center
    Δx = assembly.elements[ielem].x
    
    # linear and angular velocity
    v = vb + cross(ωb, Δx) + cross(ωb, u)# + udot = V
    ω = ωb# + Cab'*Q*θdot = Ω

    # linear and angular acceleration
    a = ab + cross(αb, Δx) + cross(αb, u)# + cross(ω, V) + uddot = Vdot
    α = αb# + Cab'*Qdot*θdot + Cab'*Q*θddot = Ωdot

    # linear and angular velocity **excluding contributions from body frame motion**
    V1 = V0[assembly.start[ielem]]
    V2 = V0[assembly.stop[ielem]]
    V = (V1 + V2)/2

    Ω1 = Ω0[assembly.start[ielem]]
    Ω2 = Ω0[assembly.stop[ielem]]
    Ω = (Ω1 + Ω2)/2

    # add contributions from body frame motion to velocities
    V += v
    Ω += ω

    # linear and angular momentum
    P = CtCab*mass11*CtCab'*V + CtCab*mass12*CtCab'*Ω
    H = CtCab*mass21*CtCab'*V + CtCab*mass22*CtCab'*Ω

    # linear and angular acceleration **excluding contributions from body frame motion**
    V1dot, Ω1dot = initial_point_velocity_rates(x, assembly.start[ielem], indices.icol_point, 
        prescribed_conditions, Vdot0, Ωdot0, rate_vars)
    V2dot, Ω2dot = initial_point_velocity_rates(x, assembly.stop[ielem], indices.icol_point, 
        prescribed_conditions, Vdot0, Ωdot0, rate_vars)
    Vdot = (V1dot + V2dot)/2
    Ωdot = (Ω1dot + Ω2dot)/2

    # add contributions from body frame motion to accelerations
    Vdot += a
    Ωdot += α

    # linear and angular momentum rates
    CtCabdot = tilde(Ω-ω)*CtCab

    Pdot = CtCab*mass11*CtCab'*Vdot + CtCab*mass12*CtCab'*Ωdot +
        CtCab*mass11*CtCabdot'*V + CtCab*mass12*CtCabdot'*Ω +
        CtCabdot*mass11*CtCab'*V + CtCabdot*mass12*CtCab'*Ω
    
    Hdot = CtCab*mass21*CtCab'*Vdot + CtCab*mass22*CtCab'*Ωdot +
        CtCab*mass21*CtCabdot'*V + CtCab*mass22*CtCabdot'*Ω +
        CtCabdot*mass21*CtCab'*V + CtCabdot*mass22*CtCab'*Ω

    # overwrite acceleration terms so we don't double count them
    a = -get_C(θb)*SVector{3}(gravity)
    α = @SVector zeros(3)

    # save properties
    properties = (; L, C, Cab, CtCab, Q, Qinv, S11, S12, S21, S22, mass11, mass12, mass21, mass22, 
        u1, u2, θ1, θ2, u, θ, F, M, γ, κ, V1, V2, Ω1, Ω2, V, Ω, P, H, ub, θb, vb, ωb, ab, αb, 
        Δx, v, ω, a, α, CtCabdot, Vdot, Ωdot, Pdot, Hdot) 

    if structural_damping 
        
        # damping coefficients
        μ = assembly.elements[ielem].mu

        # damping submatrices
        μ11 = @SMatrix [μ[1] 0 0; 0 μ[2] 0; 0 0 μ[3]]
        μ22 = @SMatrix [μ[4] 0 0; 0 μ[5] 0; 0 0 μ[6]]
    
        # linear and angular displacement rates
        udot1, θdot1 = point_velocities(x, assembly.start[ielem], indices.icol_point)
        udot2, θdot2 = point_velocities(x, assembly.stop[ielem], indices.icol_point)
        udot = (udot1 + udot2)/2
        θdot = (θdot1 + θdot2)/2

        # change in linear and angular displacement
        Δu = u2 - u1
        Δθ = θ2 - θ1

        # change in linear and angular displacement rates
        Δudot = udot2 - udot1
        Δθdot = θdot2 - θdot1

        # ΔQ matrix (see structural damping theory)
        ΔQ = get_ΔQ(θ, Δθ, Q)

        # strain rates
        γdot = -CtCab'*tilde(Ω-ω)*Δu + CtCab'*Δudot - L*CtCab'*tilde(Ω-ω)*Cab*e1
        κdot = Cab'*Q*Δθdot + Cab'*ΔQ*θdot

        # adjust strains to account for strain rates
        γ -= μ11*γdot
        κ -= μ22*κdot

        # add structural damping properties
        properties = (; properties..., γ, κ, γdot, κdot,
            μ11, μ22, udot1, udot2, udot, θdot1, θdot2, θdot, Δu, Δθ, Δudot, Δθdot, ΔQ)
    end

    return properties
end

"""
    newmark_element_properties(x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, gravity, 
        Vdot_init, Ωdot_init, dt)

Calculate/extract the element properties needed to construct the residual for a newmark-
scheme time stepping analysis
"""
@inline function newmark_element_properties(x, indices, force_scaling, 
    structural_damping, assembly, ielem, prescribed_conditions, gravity, 
    Vdot_init, Ωdot_init, dt)

    properties = steady_state_element_properties(x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity)

    @unpack C, Cab, CtCab, mass11, mass12, mass21, mass22, V1, V2, Ω1, Ω2, V, Ω, ω, θb = properties

    # linear and angular acceleration
    V1dot = 2/dt*V1 - SVector{3}(Vdot_init[assembly.start[ielem]])
    Ω1dot = 2/dt*Ω1 - SVector{3}(Ωdot_init[assembly.start[ielem]])

    V2dot = 2/dt*V2 - SVector{3}(Vdot_init[assembly.stop[ielem]])
    Ω2dot = 2/dt*Ω2 - SVector{3}(Ωdot_init[assembly.stop[ielem]])

    Vdot = (V1dot + V2dot)/2
    Ωdot = (Ω1dot + Ω2dot)/2

    # linear and angular momentum rates
    CtCabdot = tilde(Ω-ω)*CtCab
    
    Pdot = CtCab*mass11*CtCab'*Vdot + CtCab*mass12*CtCab'*Ωdot +
        CtCab*mass11*CtCabdot'*V + CtCab*mass12*CtCabdot'*Ω +
        CtCabdot*mass11*CtCab'*V + CtCabdot*mass12*CtCab'*Ω
    
    Hdot = CtCab*mass21*CtCab'*Vdot + CtCab*mass22*CtCab'*Ωdot +
        CtCab*mass21*CtCabdot'*V + CtCab*mass22*CtCabdot'*Ω +
        CtCabdot*mass21*CtCab'*V + CtCabdot*mass22*CtCab'*Ω

    # overwrite acceleration terms so we don't double count them
    a = -get_C(θb)*SVector{3}(gravity)
    α = @SVector zeros(3)
    
    return (; properties..., CtCabdot, Vdot, Ωdot, Pdot, Hdot, a, α) 
end

"""
    dynamic_element_properties(dx, x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, gravity)

Calculate/extract the element properties needed to construct the residual for a dynamic 
analysis
"""
@inline function dynamic_element_properties(dx, x, indices, force_scaling, 
    structural_damping, assembly, ielem, prescribed_conditions, gravity)

    properties = steady_state_element_properties(x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity)

    @unpack C, Cab, CtCab, mass11, mass12, mass21, mass22, V1, V2, Ω1, Ω2, V, Ω, ω, θb = properties

    # linear and angular acceleration
    V1dot, Ω1dot = point_velocities(dx, assembly.start[ielem], indices.icol_point)
    V2dot, Ω2dot = point_velocities(dx, assembly.stop[ielem], indices.icol_point)
    Vdot = (V1dot + V2dot)/2
    Ωdot = (Ω1dot + Ω2dot)/2

    # linear and angular momentum rates
    CtCabdot = tilde(Ω-ω)*CtCab
    
    Pdot = CtCab*mass11*CtCab'*Vdot + CtCab*mass12*CtCab'*Ωdot +
        CtCab*mass11*CtCabdot'*V + CtCab*mass12*CtCabdot'*Ω +
        CtCabdot*mass11*CtCab'*V + CtCabdot*mass12*CtCab'*Ω
    
    Hdot = CtCab*mass21*CtCab'*Vdot + CtCab*mass22*CtCab'*Ωdot +
        CtCab*mass21*CtCabdot'*V + CtCab*mass22*CtCabdot'*Ω +
        CtCabdot*mass21*CtCab'*V + CtCabdot*mass22*CtCab'*Ω

    # overwrite acceleration terms so we don't double count them
    a = -get_C(θb)*SVector{3}(gravity)
    α = @SVector zeros(3)

    return (; properties..., CtCabdot, Vdot, Ωdot, Pdot, Hdot, a, α)
end

"""
    expanded_steady_element_properties(x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, gravity)

Calculate/extract the element properties needed to construct the residual for a constant
mass matrix system
"""
@inline function expanded_steady_element_properties(x, indices, force_scaling, 
    structural_damping, assembly, ielem, prescribed_conditions, gravity)

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
    Q = get_Q(θ)
    Qinv = get_Qinv(θ)

    C1 = get_C(θ1)
    C2 = get_C(θ2)

    # forces and moments
    F1, M1, F2, M2 = expanded_element_loads(x, ielem, indices.icol_elem, force_scaling)
    F = (F1 + F2)/2
    M = (M1 + M2)/2

    # strain and curvature
    γ = S11*F + S12*M
    κ = S21*F + S22*M

    # linear and angular displacement of the body frame
    ub, θb = body_frame_displacement(x)

    # linear and angular velocity of the body frame
    vb, ωb = body_frame_velocity(x)

    # linear and angular acceleration of the body frame
    ab, αb = body_frame_acceleration(x)

    # distance from the rotation center
    Δx = assembly.elements[ielem].x

    # linear and angular velocity
    v = vb + cross(ωb, Δx) + cross(ωb, u)# + udot = CtCab*V
    ω = ωb# + Cab'*Q*θdot = CtCab*Ω

    # linear and angular acceleration
    a = ab + cross(αb, Δx) + cross(αb, u)# + cross(ω, V) + uddot = d/dt (CtCab*V)
    α = αb# + Cab'*Qdot*θdot + Cab'*Q*θddot = d/dt (CtCab*Ω)

    # add gravitational acceleration
    a -= get_C(θb)*gravity

    # linear and angular velocity
    V1, Ω1 = point_velocities(x, assembly.start[ielem], indices.icol_point)
    V2, Ω2 = point_velocities(x, assembly.stop[ielem], indices.icol_point)
    V, Ω = expanded_element_velocities(x, ielem, indices.icol_elem)

    # linear and angular momentum
    P = mass11*V + mass12*Ω
    H = mass21*V + mass22*Ω

    # element properties
    properties = (; L, C, C1, C2, Cab, CtCab, Q, Qinv, S11, S12, S21, S22, 
        mass11, mass12, mass21, mass22, u1, u2, θ1, θ2, V1, V2, Ω1, Ω2, F1, F2, M1, M2, 
        u, θ, V, Ω, P, H, F, M, γ, κ, ub, θb, vb, ωb, ab, αb, Δx, v, ω, a, α)

    if structural_damping

        # damping coefficients
        μ = assembly.elements[ielem].mu

        # damping submatrices
        μ11 = @SMatrix [μ[1] 0 0; 0 μ[2] 0; 0 0 μ[3]]
        μ22 = @SMatrix [μ[4] 0 0; 0 μ[5] 0; 0 0 μ[6]]

        # additional rotation parameter matrices
        C1 = get_C(θ1)
        C2 = get_C(θ2)
        Qinv1 = get_Qinv(θ1)
        Qinv2 = get_Qinv(θ2)

        # distance from rotation center
        Δx1 = assembly.points[assembly.start[ielem]]
        Δx2 = assembly.points[assembly.stop[ielem]]

        # linear and angular velocity
        v1 = vb + cross(ωb, Δx1) + cross(ωb, u1)
        v2 = vb + cross(ωb, Δx2) + cross(ωb, u2)
        ω1 = ωb
        ω2 = ωb

        # linear displacement rates 
        udot1 = C1'*V1 - v1
        udot2 = C2'*V2 - v2
        udot = (udot1 + udot2)/2

        # angular displacement rates
        θdot1 = Qinv1*(Ω1 - C1*ω1)
        θdot2 = Qinv2*(Ω2 - C2*ω2)
        θdot = (θdot1 + θdot2)/2

        # change in linear and angular displacement
        Δu = u2 - u1
        Δθ = θ2 - θ1

        # change in linear and angular displacement rates
        Δudot = udot2 - udot1
        Δθdot = θdot2 - θdot1

        # ΔQ matrix (see structural damping theory)
        ΔQ = get_ΔQ(θ, Δθ, Q)

        # strain rates
        γdot = -CtCab'*tilde(CtCab*Ω - ω)*Δu + CtCab'*Δudot - L*CtCab'*tilde(CtCab*Ω - ω)*Cab*e1
        κdot = Cab'*Q*Δθdot + Cab'*ΔQ*θdot

        # adjust strains to account for strain rates
        γ -= μ11*γdot
        κ -= μ22*κdot

        # add structural damping properties
        properties = (; properties..., γ, κ, γdot, κdot,
            μ11, μ22, C1, C2, Qinv1, Qinv2, Δx1, Δx2, v1, v2, ω1, ω2, 
            udot1, udot2, udot, θdot1, θdot2, θdot, Δu, Δθ, Δudot, Δθdot, ΔQ)
    end

    return properties
end

function expanded_dynamic_element_properties(x, indices, force_scaling, structural_damping, 
    assembly, ielem, prescribed_conditions, gravity)

    properties = expanded_steady_element_properties(x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity)

    @unpack θb, a, α = properties

    # overwrite acceleration terms so we don't double count them
    a = -get_C(θb)*SVector{3}(gravity)
    α = @SVector zeros(3)

    return (; properties..., a, α)
end

"""
    compatability_residuals(properties)

Calculate the compatability residuals for the beam element
"""
@inline function compatability_residuals(properties)
   
    @unpack L, Cab, CtCab, Qinv, u1, u2, θ1, θ2, γ, κ = properties

    Δu = CtCab*(L*e1 + γ) - L*Cab*e1

    Δθ = Qinv*Cab*κ

    ru = u2 - u1 - Δu
    rθ = θ2 - θ1 - Δθ

    return (; ru, rθ)
end

"""
    expanded_element_velocity_residuals(properties)

Calculate the velocity residuals `rV` and `rΩ` of a beam element for a constant mass matrix 
system.
"""
@inline function expanded_element_velocity_residuals(properties)

    @unpack C, Cab, CtCab, Qinv, u, V, Ω, v, ω = properties
    
    rV = CtCab*V - v
    rΩ = Qinv*(Cab*Ω - C*ω)

    # @unpack CtCab, V, Ω, C1, V1, Ω1, C2, V2, Ω2 = properties
    # rV = CtCab*V - 1/2*(C1'*V1 + C2'*V2)
    # rΩ = CtCab*Ω - 1/2*(C1'*Ω1 + C2'*Ω2)

    return (; rV, rΩ)
end

"""
    expanded_element_equilibrium_residuals(properties, distributed_loads, ielem)

Calculate the equilibrium residuals of a beam element for a constant mass matrix system.
"""
@inline function expanded_element_equilibrium_residuals(properties, distributed_loads, ielem)
  
    @unpack L, C, Cab, CtCab, mass11, mass12, mass21, mass22, F1, F2, M1, M2, F, M, γ, κ, 
        V, Ω, P, H, v, ω, a, α = properties

    # initialize equilibrium residual
    rF = F2 - F1
    rM = M2 - M1

    # add loads due to internal forces/moments and stiffness
    rM += cross(L*e1 + γ, F)

    # add loads due to linear and angular acceleration (including gravity)
    rF -= mass11*CtCab'*a + mass12*CtCab'*α
    rM -= mass21*CtCab'*a + mass22*CtCab'*α

    # add loads due to linear and angular momentum
    rF -= cross(Ω, P)
    rM -= cross(Ω, H) + cross(V, P)

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
    static_element_resultants(properties, distributed_loads, ielem)

Calculate the resultant loads applied at each end of a beam element for a static analysis.
"""
@inline function static_element_resultants(properties, distributed_loads, ielem)
  
    @unpack L, C, Cab, CtCab, mass11, mass12, mass21, mass22, F, M, γ, κ, a, α = properties

    # add loads due to internal forces/moments and stiffness
    tmp = CtCab*F
    F1 = tmp
    F2 = tmp

    tmp1 = CtCab*M
    tmp2 = 1/2*CtCab*cross(L*e1 + γ, F)
    M1 = tmp1 + tmp2
    M2 = tmp1 - tmp2

    # add loads due to an accelerating reference frame
    tmp = CtCab*mass11*CtCab'*a + CtCab*mass12*CtCab'*α
    F1 -= 1/2*tmp
    F2 += 1/2*tmp

    tmp = CtCab*mass21*CtCab'*a + CtCab*mass22*CtCab'*α
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
    steady_state_element_resultants(properties, distributed_loads, ielem)

Calculate the resultant loads applied at each end of a beam element for a steady_state 
analysis.
"""
@inline function steady_state_element_resultants(properties, distributed_loads, ielem)

    resultants = static_element_resultants(properties, distributed_loads, ielem)

    @unpack F1, F2, M1, M2 = resultants

    @unpack ω, V, Ω, P, H = properties

    # add loads due to linear and angular momentum
    tmp = cross(ω, P)
    F1 -= 1/2*tmp
    F2 += 1/2*tmp

    tmp = cross(ω, H) + cross(V, P)
    M1 -= 1/2*tmp
    M2 += 1/2*tmp

    return (; resultants..., F1, F2, M1, M2)
end

"""
    dynamic_element_resultants(properties, distributed_loads, ielem)

Calculate the resultant loads applied at each end of a beam element for a dynamic 
analysis.
"""
@inline function dynamic_element_resultants(properties, distributed_loads, ielem)

    resultants = steady_state_element_resultants(properties, distributed_loads, ielem)

    @unpack F1, F2, M1, M2 = resultants

    @unpack Pdot, Hdot = properties
    
    # add loads due to linear and angular momentum rates  
    F1 -= 1/2*Pdot
    F2 += 1/2*Pdot

    M1 -= 1/2*Hdot
    M2 += 1/2*Hdot

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
    insert_element_residuals!(resid, indices, force_scaling, assembly, ielem, 
        compatability, resultants)

Insert the residual entries corresponding to a beam element into the system residual vector.
"""
@inline function insert_element_residuals!(resid, indices, force_scaling, assembly, ielem, 
    compatability, resultants)

    @unpack ru, rθ = compatability
    @unpack F1, F2, M1, M2 = resultants

    # compatability equations
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
        compatability, velocities, equilibrium, resultants)

Insert the residual entries corresponding to a beam element into the system residual vector
for a constant mass matrix system.
"""
@inline function insert_expanded_element_residuals!(resid, indices, force_scaling, assembly, 
    ielem, compatability, velocities, equilibrium, resultants)

    insert_element_residuals!(resid, indices, force_scaling, assembly, ielem, 
        compatability, resultants)

    @unpack rV, rΩ = velocities
    @unpack rF, rM = equilibrium

    irow = indices.irow_elem[ielem]

    # equilibrium residuals
    resid[irow+6:irow+8] .= rF ./ force_scaling
    resid[irow+9:irow+11] .= rM ./ force_scaling
    
    # velocity residuals
    resid[irow+12:irow+14] .= rV
    resid[irow+15:irow+17] .= rΩ

    return resid
end

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

    compatability = compatability_residuals(properties)

    resultants = static_element_resultants(properties, distributed_loads, ielem)

    insert_element_residuals!(resid, indices, force_scaling, assembly, ielem, 
        compatability, resultants)

    return resid
end

"""
    steady_state_element_residual!(resid, x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, 
        gravity)

Calculate and insert the residual entries corresponding to a beam element for a steady state 
analysis into the system residual vector.
"""
@inline function steady_state_element_residual!(resid, x, indices, force_scaling, structural_damping, 
    assembly, ielem, prescribed_conditions, distributed_loads, gravity)

    properties = steady_state_element_properties(x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity)

    compatability = compatability_residuals(properties)

    resultants = steady_state_element_resultants(properties, distributed_loads, ielem)

    insert_element_residuals!(resid, indices, force_scaling, assembly, ielem, 
        compatability, resultants)

    return resid
end

"""
    initial_condition_element_residual!(resid, x, indices, rate_vars, 
        force_scaling, structural_damping, assembly, ielem, prescribed_conditions, 
        distributed_loads, gravity, u0, θ0, V0, Ω0, Vdot0, Ωdot0)

Calculate and insert the residual entries corresponding to a beam element for the 
initialization of a time domain simulation into the system residual vector.
"""
@inline function initial_condition_element_residual!(resid, x, indices, rate_vars, 
    force_scaling, structural_damping, assembly, ielem, prescribed_conditions, 
    distributed_loads, gravity, u0, θ0, V0, Ω0, Vdot0, Ωdot0)

    properties = initial_condition_element_properties(x, indices, rate_vars, 
        force_scaling, structural_damping, assembly, ielem, prescribed_conditions, gravity, 
        u0, θ0, V0, Ω0, Vdot0, Ωdot0)

    compatability = compatability_residuals(properties)

    resultants = dynamic_element_resultants(properties, distributed_loads, ielem)

    insert_element_residuals!(resid, indices, force_scaling, assembly, ielem, 
        compatability, resultants)

    return resid
end

"""
    newmark_element_residual!(resid, x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, 
        gravity, Vdot_init, Ωdot_init, dt)

Calculate and insert the residual entries corresponding to a beam element for a 
newmark-scheme time marching analysis into the system residual vector.
"""
@inline function newmark_element_residual!(resid, x, indices, force_scaling, 
    structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, gravity, 
    Vdot_init, Ωdot_init, dt)

    properties = newmark_element_properties(x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity, 
        Vdot_init, Ωdot_init, dt)

    compatability = compatability_residuals(properties)

    resultants = dynamic_element_resultants(properties, distributed_loads, ielem)

    insert_element_residuals!(resid, indices, force_scaling, assembly, ielem, 
        compatability, resultants)

    return resid
end

"""
    dynamic_element_residual!(resid, dx, x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, 
        gravity)

Calculate and insert the residual entries corresponding to a beam element for a dynamic
analysis into the system residual vector.
"""
@inline function dynamic_element_residual!(resid, dx, x, indices, force_scaling, 
    structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, gravity)

    properties = dynamic_element_properties(dx, x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity)

    compatability = compatability_residuals(properties)

    resultants = dynamic_element_resultants(properties, distributed_loads, ielem)

    insert_element_residuals!(resid, indices, force_scaling, assembly, ielem, 
        compatability, resultants)

    return resid
end

"""
    expanded_steady_element_residual!(resid, x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, 
        gravity)

Calculate and insert the residual entries corresponding to a beam element for a constant
mass matrix system into the system residual vector.
"""
@inline function expanded_steady_element_residual!(resid, x, indices, force_scaling, 
    structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, gravity)

    properties = expanded_steady_element_properties(x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, gravity)

    compatability = compatability_residuals(properties)

    velocities = expanded_element_velocity_residuals(properties)

    equilibrium = expanded_element_equilibrium_residuals(properties, distributed_loads, ielem)

    resultants = expanded_element_resultants(properties)

    insert_expanded_element_residuals!(resid, indices, force_scaling, assembly, ielem, 
        compatability, velocities, equilibrium, resultants)

    return resid
end

"""
    expanded_dynamic_element_residual!(resid, x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, 
        gravity)

Calculate and insert the residual entries corresponding to a beam element for a constant
mass matrix system into the system residual vector.
"""
@inline function expanded_dynamic_element_residual!(resid, x, indices, force_scaling, 
    structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, gravity)

    properties = expanded_dynamic_element_properties(x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, gravity)

    compatability = compatability_residuals(properties)

    velocities = expanded_element_velocity_residuals(properties)

    equilibrium = expanded_element_equilibrium_residuals(properties, distributed_loads, ielem)

    resultants = expanded_element_resultants(properties)

    insert_expanded_element_residuals!(resid, indices, force_scaling, assembly, ielem, 
        compatability, velocities, equilibrium, resultants)

    return resid
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

    return (; properties..., u1_u1, u2_u2, θ1_θ1, θ2_θ2, C_θ1, C_θ2, C_θ3,
        Qinv_θ1, Qinv_θ2, Qinv_θ3, γ_F, γ_M, κ_F, κ_M)
end

"""
    steady_state_element_jacobian_properties(properties, x, indices, 
        force_scaling, structural_damping, assembly, ielem, prescribed_conditions, 
        gravity)

Calculate/extract the element properties needed to calculate the jacobian entries 
corresponding to a element for a steady state analysis
"""
@inline function steady_state_element_jacobian_properties(properties, x, indices, 
    force_scaling, structural_damping, assembly, ielem, prescribed_conditions, 
    gravity)

    properties = static_element_jacobian_properties(properties, x, indices, force_scaling, 
        assembly, ielem, prescribed_conditions, gravity)

    @unpack L, mass11, mass12, mass21, mass22, C, Cab, CtCab, Q, u, θ, V, Ω, θb, αb, Δx, ω, 
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

    # linear velocity
    v_vb = I3
    v_ωb = -tilde(Δx) - tilde(u)
    
    # angular velocity
    ω_ωb = I3

    # linear acceleration
    a_u = tilde(αb)
    a_ab = I3
    a_αb = -tilde(Δx) - tilde(u)
    
    # angular acceleration
    α_αb = I3

    # add rotated gravity vector
    C_θb1, C_θb2, C_θb3 = get_C_θ(θb)
    a_θb = -mul3(C_θb1, C_θb2, C_θb3, gravity)

    if structural_damping

        @unpack C1, C2, Qinv1, Qinv2, u1, u2, θ1, θ2, Ω1, Ω2, ω1, ω2, μ11, μ22, udot, θdot, 
            Δx1, Δx2, Δu, Δθ, ΔQ, Δudot, Δθdot = properties

        # rotation parameter matrices
        C1_θ1, C1_θ2, C1_θ3 = get_C_θ(C1, θ1)
        C2_θ1, C2_θ2, C2_θ3 = get_C_θ(C2, θ2)
        Qinv1_θ1, Qinv1_θ2, Qinv1_θ3 = get_Qinv_θ(θ1)
        Qinv2_θ1, Qinv2_θ2, Qinv2_θ3 = get_Qinv_θ(θ2)

        # linear displacement rates
        udot1_vb = -I3
        udot2_vb = -I3

        udot1_ωb = tilde(Δx1 + u1)
        udot2_ωb = tilde(Δx2 + u2)

        udot1_u1 = -tilde(ω1)
        udot2_u2 = -tilde(ω2)

        udot1_V1 = I3
        udot2_V2 = I3

        # angular displacement rates
        θdot1_ωb = -Qinv1*C1
        θdot2_ωb = -Qinv2*C2

        θdot1_θ1 = mul3(Qinv1_θ1, Qinv1_θ2, Qinv1_θ3, C1*(Ω1 - ω1)) + Qinv1*mul3(C1_θ1, C1_θ2, C1_θ3, Ω1 - ω1)
        θdot2_θ2 = mul3(Qinv2_θ1, Qinv2_θ2, Qinv2_θ3, C2*(Ω2 - ω2)) + Qinv2*mul3(C2_θ1, C2_θ2, C2_θ3, Ω2 - ω2)
        
        θdot1_Ω1 = Qinv1*C1
        θdot2_Ω2 = Qinv2*C2

        θdot_ωb = (θdot1_ωb + θdot2_ωb)/2

        θdot_θ1 = θdot1_θ1/2
        θdot_θ2 = θdot2_θ2/2

        θdot_Ω1 = θdot1_Ω1/2
        θdot_Ω2 = θdot2_Ω2/2

        # change in linear displacement
        Δu_u1 = -I3
        Δu_u2 =  I3

        # change in linear and angular displacement rates
        Δudot_vb = udot2_vb - udot1_vb
        Δudot_ωb = udot2_ωb - udot1_ωb

        Δudot_u1 = -udot1_u1
        Δudot_u2 =  udot2_u2

        Δudot_V1 = -udot1_V1
        Δudot_V2 =  udot2_V2

        Δθdot_ωb = θdot2_ωb - θdot1_ωb

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

        γdot_vb = CtCab'*Δudot_vb

        γdot_ωb = -CtCab'*tilde(Δu) + CtCab'*Δudot_ωb - L*CtCab'*tilde(Cab*e1)

        tmp = CtCab'*tilde(Ω - ω)
        γdot_u1 = -tmp*Δu_u1 + CtCab'*Δudot_u1 
        γdot_u2 = -tmp*Δu_u2 + CtCab'*Δudot_u2 

        tmp = -Cab'*mul3(C_θ1, C_θ2, C_θ3, tilde(Ω - ω)*Δu) + 
            Cab'*mul3(C_θ1, C_θ2, C_θ3, Δudot) - 
            L*Cab'*mul3(C_θ1, C_θ2, C_θ3, tilde(Ω - ω)*Cab*e1)
        γdot_θ1 = 1/2*tmp
        γdot_θ2 = 1/2*tmp

        γdot_V1 = CtCab'*Δudot_V1
        γdot_V2 = CtCab'*Δudot_V2

        tmp = CtCab'*tilde(Δu) + L*CtCab'*tilde(Cab*e1)
        γdot_Ω1 = 1/2*tmp
        γdot_Ω2 = 1/2*tmp

        κdot_ωb = Cab'*Q*Δθdot_ωb + Cab'*ΔQ*θdot_ωb

        tmp1 = Cab'*mul3(Q_θ1, Q_θ2, Q_θ3, Δθdot)
        tmp2 = Cab'*mul3(ΔQ_θ1, ΔQ_θ2, ΔQ_θ3, θdot) 
        tmp3 = Cab'*mul3(ΔQ_Δθ1, ΔQ_Δθ2, ΔQ_Δθ3, θdot)
        κdot_θ1 = 1/2*tmp1 + 1/2*tmp2 - tmp3 + Cab'*Q*Δθdot_θ1 + Cab'*ΔQ*θdot_θ1
        κdot_θ2 = 1/2*tmp1 + 1/2*tmp2 + tmp3 + Cab'*Q*Δθdot_θ2 + Cab'*ΔQ*θdot_θ2  

        κdot_Ω1 = Cab'*Q*Δθdot_Ω1 + Cab'*ΔQ*θdot_Ω1
        κdot_Ω2 = Cab'*Q*Δθdot_Ω2 + Cab'*ΔQ*θdot_Ω2

        # adjust strains to account for strain rates
        
        γ_vb = -μ11*γdot_vb
        γ_ωb = -μ11*γdot_ωb

        γ_u1 = -μ11*γdot_u1
        γ_u2 = -μ11*γdot_u2

        γ_θ1 = -μ11*γdot_θ1
        γ_θ2 = -μ11*γdot_θ2
        
        γ_V1 = -μ11*γdot_V1
        γ_V2 = -μ11*γdot_V2
        
        γ_Ω1 = -μ11*γdot_Ω1
        γ_Ω2 = -μ11*γdot_Ω2
        
        κ_ωb = -μ22*κdot_ωb

        κ_θ1 = -μ22*κdot_θ1
        κ_θ2 = -μ22*κdot_θ2
        
        κ_Ω1 = -μ22*κdot_Ω1
        κ_Ω2 = -μ22*κdot_Ω2

    else
        
        γ_vb = @SMatrix zeros(3,3)
        γ_ωb = @SMatrix zeros(3,3)

        γ_u1 = @SMatrix zeros(3,3)
        γ_u2 = @SMatrix zeros(3,3)

        γ_θ1 = @SMatrix zeros(3,3)
        γ_θ2 = @SMatrix zeros(3,3)
        
        γ_V1 = @SMatrix zeros(3,3)
        γ_V2 = @SMatrix zeros(3,3)
        
        γ_Ω1 = @SMatrix zeros(3,3)
        γ_Ω2 = @SMatrix zeros(3,3)
        
        κ_ωb = @SMatrix zeros(3,3)

        κ_θ1 = @SMatrix zeros(3,3)
        κ_θ2 = @SMatrix zeros(3,3)
        
        κ_Ω1 = @SMatrix zeros(3,3)
        κ_Ω2 = @SMatrix zeros(3,3)
         
    end

    return (; properties..., γ_vb, γ_ωb, γ_u1, γ_u2, γ_θ1, γ_θ2, γ_V1, γ_V2, γ_Ω1, γ_Ω2, 
        κ_ωb, κ_θ1, κ_θ2, κ_Ω1, κ_Ω2, P_θ, P_V, P_Ω, H_θ, H_V, H_Ω, v_vb, v_ωb, ω_ωb, 
        a_θb, a_ab, a_αb, a_u, α_αb)
end

"""
    initial_condition_element_jacobian_properties(properties, x, indices, rate_vars, 
        force_scaling, structural_damping, assembly, ielem, 
        prescribed_conditions, gravity, 
        u0, θ0, V0, Ω0, Vdot0, Ωdot0)

Calculate/extract the element properties needed to calculate the jacobian entries 
corresponding to a element for a Newmark scheme time marching analysis
"""
@inline function initial_condition_element_jacobian_properties(properties, x, indices, 
    rate_vars, force_scaling, structural_damping, assembly, ielem, 
    prescribed_conditions, gravity, 
    u0, θ0, V0, Ω0, Vdot0, Ωdot0)

    @unpack L, C, Cab, CtCab, CtCabdot, Q, S11, S12, S21, S22, mass11, mass12, mass21, mass22, 
        u, θ, V, Ω, Vdot, Ωdot, θb, ωb, αb, Δx, ω = properties

    # starting node linear and angular displacement
    u1_u1, θ1_θ1 = initial_point_displacement_jacobian(assembly.start[ielem], indices.icol_point, 
        prescribed_conditions, rate_vars)

    # ending node linear and angular displacement
    u2_u2, θ2_θ2 = initial_point_displacement_jacobian(assembly.stop[ielem], indices.icol_point, 
        prescribed_conditions, rate_vars)

    # rotation parameter matrices
    C_θ1, C_θ2, C_θ3 = get_C_θ(C, θ)
    Q_θ1, Q_θ2, Q_θ3 = get_Q_θ(Q, θ)
    Qinv_θ1, Qinv_θ2, Qinv_θ3 = get_Qinv_θ(θ)

    # strain and curvature
    γ_F, γ_M, κ_F, κ_M = S11, S12, S21, S22
   
    # linear and angular velocity
    v_u = tilde(ωb)
    v_vb = I3
    v_ωb = -tilde(Δx) - tilde(u)
    ω_ωb = I3   

    # linear and angular acceleration
    a_u = tilde(αb)
    a_ab = I3
    a_αb = -tilde(Δx) - tilde(u)
    α_αb = I3

    # add contributions from body frame motion to velocities
    V_u = v_u
    V_vb = v_vb
    V_ωb = v_ωb
    Ω_ωb = ω_ωb

    # linear and angular momentum
    P_u = CtCab*mass11*CtCab'*V_u
    P_θ = mul3(C_θ1', C_θ2', C_θ3', Cab*(mass11*CtCab'*V + mass12*CtCab'*Ω)) + 
        CtCab*mass11*Cab'*mul3(C_θ1, C_θ2, C_θ3, V) + 
        CtCab*mass12*Cab'*mul3(C_θ1, C_θ2, C_θ3, Ω)
    P_vb = CtCab*mass11*CtCab'*V_vb
    P_ωb = CtCab*mass11*CtCab'*V_ωb + CtCab*mass12*CtCab'*Ω_ωb

    H_u = CtCab*mass21*CtCab'*V_u
    H_θ = mul3(C_θ1', C_θ2', C_θ3', Cab*(mass21*CtCab'*V + mass22*CtCab'*Ω)) + 
        CtCab*mass21*Cab'*mul3(C_θ1, C_θ2, C_θ3, V) + 
        CtCab*mass22*Cab'*mul3(C_θ1, C_θ2, C_θ3, Ω)
    H_vb = CtCab*mass21*CtCab'*V_vb
    H_ωb = CtCab*mass21*CtCab'*V_ωb + CtCab*mass22*CtCab'*Ω_ωb

    # starting node linear and angular velocity rates
    V1dot_V1dot, Ω1dot_Ω1dot = initial_point_velocity_rate_jacobian(assembly.start[ielem], 
        indices.icol_point, prescribed_conditions, rate_vars)

    # ending node linear and angular velocity rates
    V2dot_V2dot, Ω2dot_Ω2dot = initial_point_velocity_rate_jacobian(assembly.stop[ielem], 
        indices.icol_point, prescribed_conditions, rate_vars)
 
    # add contributions from body frame motion to accelerations
    Vdot_u = a_u
    Vdot_ab = a_ab
    Vdot_αb = a_αb
    Ωdot_αb = α_αb

    # linear and angular momentum rates
    Pdot_V = CtCab*mass11*CtCabdot' + CtCabdot*mass11*CtCab'
    Pdot_Ω = CtCab*mass12*CtCabdot' + CtCabdot*mass12*CtCab'
    Pdot_Vdot = CtCab*mass11*CtCab'
    Pdot_Ωdot = CtCab*mass12*CtCab'
    Pdot_u = Pdot_V*V_u + Pdot_Vdot*Vdot_u
    Pdot_θ = mul3(C_θ1', C_θ2', C_θ3', Cab*mass11*CtCab'*Vdot) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass12*CtCab'*Ωdot) +
        CtCab*mass11*Cab'*mul3(C_θ1, C_θ2, C_θ3, Vdot) + 
        CtCab*mass12*Cab'*mul3(C_θ1, C_θ2, C_θ3, Ωdot) +
        tilde(Ω - ω)*mul3(C_θ1', C_θ2', C_θ3', Cab*mass11*CtCab'*V) + 
        tilde(Ω - ω)*mul3(C_θ1', C_θ2', C_θ3', Cab*mass12*CtCab'*Ω) +
        CtCabdot*mass11*Cab'*mul3(C_θ1, C_θ2, C_θ3, V) + 
        CtCabdot*mass12*Cab'*mul3(C_θ1, C_θ2, C_θ3, Ω) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass11*CtCabdot'*V) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass12*CtCabdot'*Ω) +
        -CtCab*mass11*Cab'*mul3(C_θ1, C_θ2, C_θ3, tilde(Ω - ω)*V) + 
        -CtCab*mass12*Cab'*mul3(C_θ1, C_θ2, C_θ3, tilde(Ω - ω)*Ω)
    Pdot_vb = Pdot_V*V_vb
    Pdot_ωb = Pdot_V*V_ωb + Pdot_Ω*Ω_ωb
    Pdot_ab = Pdot_Vdot*Vdot_ab
    Pdot_αb = Pdot_Vdot*Vdot_αb + Pdot_Ωdot*Ωdot_αb

    Hdot_V = CtCab*mass21*CtCabdot' + CtCabdot*mass21*CtCab'
    Hdot_Ω = CtCab*mass22*CtCabdot' + CtCabdot*mass22*CtCab'
    Hdot_Vdot = CtCab*mass21*CtCab'
    Hdot_Ωdot = CtCab*mass22*CtCab'
    Hdot_u = Hdot_V*V_u + Hdot_Vdot*Vdot_u
    Hdot_θ = mul3(C_θ1', C_θ2', C_θ3', Cab*mass21*CtCab'*Vdot) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass22*CtCab'*Ωdot) +
        CtCab*mass21*Cab'*mul3(C_θ1, C_θ2, C_θ3, Vdot) + 
        CtCab*mass22*Cab'*mul3(C_θ1, C_θ2, C_θ3, Ωdot) +
        tilde(Ω - ω)*mul3(C_θ1', C_θ2', C_θ3', Cab*mass21*CtCab'*V) + 
        tilde(Ω - ω)*mul3(C_θ1', C_θ2', C_θ3', Cab*mass22*CtCab'*Ω) +
        CtCabdot*mass21*Cab'*mul3(C_θ1, C_θ2, C_θ3, V) + 
        CtCabdot*mass22*Cab'*mul3(C_θ1, C_θ2, C_θ3, Ω) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass21*CtCabdot'*V) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass22*CtCabdot'*Ω) +
        -CtCab*mass21*Cab'*mul3(C_θ1, C_θ2, C_θ3, tilde(Ω - ω)*V) + 
        -CtCab*mass22*Cab'*mul3(C_θ1, C_θ2, C_θ3, tilde(Ω - ω)*Ω)
    Hdot_vb = Hdot_V*V_vb
    Hdot_ωb = Hdot_V*V_ωb + Hdot_Ω*Ω_ωb
    Hdot_ab = Hdot_Vdot*Vdot_ab
    Hdot_αb = Hdot_Vdot*Vdot_αb + Hdot_Ωdot*Ωdot_αb

    # overwrite acceleration terms so we don't double count them
    a_u = @SMatrix zeros(3,3)
    a_θb = -mul3(get_C_θ(θb)..., gravity)
    a_ab = @SMatrix zeros(3,3)
    a_αb = @SMatrix zeros(3,3)
    α_αb = @SMatrix zeros(3,3)

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
        tmp = CtCab'*tilde(Ω - ω)
        γdot_u1 = -tmp*Δu_u1
        γdot_u2 = -tmp*Δu_u2
        γdot_u1dot = CtCab'*Δudot_u1dot
        γdot_u2dot = CtCab'*Δudot_u2dot

        tmp = -Cab'*mul3(C_θ1, C_θ2, C_θ3, tilde(Ω - ω)*Δu) + 
            Cab'*mul3(C_θ1, C_θ2, C_θ3, Δudot) - 
            L*Cab'*mul3(C_θ1, C_θ2, C_θ3, tilde(Ω - ω)*Cab*e1)
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

    return (; properties..., u1_u1, u2_u2, θ1_θ1, θ2_θ2, C_θ1, C_θ2, C_θ3,
        Qinv_θ1, Qinv_θ2, Qinv_θ3, γ_u1, γ_u2, γ_θ1, γ_θ2, γ_u1dot, γ_u2dot, γ_F, γ_M,   
        κ_θ1, κ_θ2, κ_θ1dot, κ_θ2dot, κ_F, κ_M, 
        V_u, V_vb, V_ωb, Ω_ωb, P_u, P_θ, P_vb, P_ωb, H_u, H_θ, H_vb, H_ωb,
        V1dot_V1dot, V2dot_V2dot, Ω1dot_Ω1dot, Ω2dot_Ω2dot, 
        Pdot_vb, Pdot_ωb, Pdot_ab, Pdot_αb, Pdot_u, Pdot_θ, Pdot_Vdot, Pdot_Ωdot,
        Hdot_vb, Hdot_ωb, Hdot_ab, Hdot_αb, Hdot_u, Hdot_θ, Hdot_Vdot, Hdot_Ωdot, 
        v_vb, v_ωb, ω_ωb, a_θb, a_ab, a_αb, a_u, α_αb)
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

    properties = steady_state_element_jacobian_properties(properties, x, indices, 
        force_scaling, structural_damping, assembly, ielem, 
        prescribed_conditions, gravity)

    @unpack C, Cab, CtCab, CtCabdot, mass11, mass12, mass21, mass22, V, Ω, ω, Vdot, Ωdot, 
        C_θ1, C_θ2, C_θ3 = properties

    # linear and angular momentum rates

    Pdot_ωb = -CtCab*mass11*CtCab'*tilde(V) - CtCab*mass12*CtCab'*tilde(Ω) +
        tilde(CtCab*mass11*CtCab'*V) + tilde(CtCab*mass12*CtCab'*Ω)

    Pdot_θ = mul3(C_θ1', C_θ2', C_θ3', Cab*mass11*CtCab'*Vdot) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass12*CtCab'*Ωdot) +
        CtCab*mass11*Cab'*mul3(C_θ1, C_θ2, C_θ3, Vdot) + 
        CtCab*mass12*Cab'*mul3(C_θ1, C_θ2, C_θ3, Ωdot) +
        tilde(Ω - ω)*mul3(C_θ1', C_θ2', C_θ3', Cab*mass11*CtCab'*V) + 
        tilde(Ω - ω)*mul3(C_θ1', C_θ2', C_θ3', Cab*mass12*CtCab'*Ω) +
        CtCabdot*mass11*Cab'*mul3(C_θ1, C_θ2, C_θ3, V) + 
        CtCabdot*mass12*Cab'*mul3(C_θ1, C_θ2, C_θ3, Ω) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass11*CtCabdot'*V) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass12*CtCabdot'*Ω) +
        -CtCab*mass11*Cab'*mul3(C_θ1, C_θ2, C_θ3, tilde(Ω - ω)*V) + 
        -CtCab*mass12*Cab'*mul3(C_θ1, C_θ2, C_θ3, tilde(Ω - ω)*Ω)

    Pdot_V = CtCab*mass11*CtCabdot' + CtCabdot*mass11*CtCab'

    Pdot_Ω = CtCab*mass12*CtCabdot' + CtCabdot*mass12*CtCab' +
        CtCab*mass11*CtCab'*tilde(V) + CtCab*mass12*CtCab'*tilde(Ω) +
        -tilde(CtCab*mass11*CtCab'*V) - tilde(CtCab*mass12*CtCab'*Ω)
    
    Pdot_V += 2/dt*CtCab*mass11*CtCab'
    
    Pdot_Ω += 2/dt*CtCab*mass12*CtCab'

    Hdot_ωb = -CtCab*mass21*CtCab'*tilde(V) - CtCab*mass22*CtCab'*tilde(Ω) +
        tilde(CtCab*mass21*CtCab'*V) + tilde(CtCab*mass22*CtCab'*Ω)

    Hdot_θ = mul3(C_θ1', C_θ2', C_θ3', Cab*mass21*CtCab'*Vdot) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass22*CtCab'*Ωdot) +
        CtCab*mass21*Cab'*mul3(C_θ1, C_θ2, C_θ3, Vdot) + 
        CtCab*mass22*Cab'*mul3(C_θ1, C_θ2, C_θ3, Ωdot) +
        tilde(Ω - ω)*mul3(C_θ1', C_θ2', C_θ3', Cab*mass21*CtCab'*V) + 
        tilde(Ω - ω)*mul3(C_θ1', C_θ2', C_θ3', Cab*mass22*CtCab'*Ω) +
        CtCabdot*mass21*Cab'*mul3(C_θ1, C_θ2, C_θ3, V) + 
        CtCabdot*mass22*Cab'*mul3(C_θ1, C_θ2, C_θ3, Ω) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass21*CtCabdot'*V) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass22*CtCabdot'*Ω) +
        -CtCab*mass21*Cab'*mul3(C_θ1, C_θ2, C_θ3, tilde(Ω - ω)*V) + 
        -CtCab*mass22*Cab'*mul3(C_θ1, C_θ2, C_θ3, tilde(Ω - ω)*Ω)

    Hdot_V = CtCab*mass21*CtCabdot' + CtCabdot*mass21*CtCab'

    Hdot_Ω = CtCab*mass22*CtCabdot' + CtCabdot*mass22*CtCab' +
        CtCab*mass21*CtCab'*tilde(V) + CtCab*mass22*CtCab'*tilde(Ω) +
        -tilde(CtCab*mass21*CtCab'*V) - tilde(CtCab*mass22*CtCab'*Ω)

    Hdot_V += 2/dt*CtCab*mass21*CtCab'

    Hdot_Ω += 2/dt*CtCab*mass22*CtCab'

    # overwrite acceleration terms so we don't double count them
    a_u = @SMatrix zeros(3,3)
    a_ab = @SMatrix zeros(3,3)
    a_αb = @SMatrix zeros(3,3)
    α_αb = @SMatrix zeros(3,3)

    return (; properties..., Pdot_ωb, Pdot_θ, Pdot_V, Pdot_Ω, Hdot_ωb, Hdot_θ, Hdot_V, Hdot_Ω,
        a_u, a_ab, a_αb, α_αb) 
end

"""
    dynamic_element_jacobian_properties(properties, dx, x, indices, 
        force_scaling, assembly, ielem, prescribed_conditions, gravity)

Calculate/extract the element properties needed to calculate the jacobian entries 
corresponding to a element for a dynamic analysis
"""
@inline function dynamic_element_jacobian_properties(properties, dx, x, indices, 
    force_scaling, structural_damping, assembly, ielem, prescribed_conditions, gravity)

    properties = steady_state_element_jacobian_properties(properties, x, indices, 
        force_scaling, structural_damping, assembly, ielem, 
        prescribed_conditions, gravity)

    @unpack C, Cab, CtCab, CtCabdot, mass11, mass12, mass21, mass22, V, Ω, Vdot, Ωdot, ω, C_θ1, C_θ2, C_θ3 = properties

    # linear and angular momentum rates

    Pdot_ωb = -CtCab*mass11*CtCab'*tilde(V) - CtCab*mass12*CtCab'*tilde(Ω) +
        tilde(CtCab*mass11*CtCab'*V) + tilde(CtCab*mass12*CtCab'*Ω)

    Pdot_θ = mul3(C_θ1', C_θ2', C_θ3', Cab*mass11*CtCab'*Vdot) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass12*CtCab'*Ωdot) +
        CtCab*mass11*Cab'*mul3(C_θ1, C_θ2, C_θ3, Vdot) + 
        CtCab*mass12*Cab'*mul3(C_θ1, C_θ2, C_θ3, Ωdot) +
        tilde(Ω - ω)*mul3(C_θ1', C_θ2', C_θ3', Cab*mass11*CtCab'*V) + 
        tilde(Ω - ω)*mul3(C_θ1', C_θ2', C_θ3', Cab*mass12*CtCab'*Ω) +
        CtCabdot*mass11*Cab'*mul3(C_θ1, C_θ2, C_θ3, V) + 
        CtCabdot*mass12*Cab'*mul3(C_θ1, C_θ2, C_θ3, Ω) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass11*CtCabdot'*V) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass12*CtCabdot'*Ω) +
        -CtCab*mass11*Cab'*mul3(C_θ1, C_θ2, C_θ3, tilde(Ω - ω)*V) + 
        -CtCab*mass12*Cab'*mul3(C_θ1, C_θ2, C_θ3, tilde(Ω - ω)*Ω)

    Pdot_V = CtCab*mass11*CtCabdot' + CtCabdot*mass11*CtCab'

    Pdot_Ω = CtCab*mass12*CtCabdot' + CtCabdot*mass12*CtCab' +
        CtCab*mass11*CtCab'*tilde(V) + CtCab*mass12*CtCab'*tilde(Ω) +
        -tilde(CtCab*mass11*CtCab'*V) - tilde(CtCab*mass12*CtCab'*Ω)
    
    Hdot_ωb = -CtCab*mass21*CtCab'*tilde(V) - CtCab*mass22*CtCab'*tilde(Ω) +
        tilde(CtCab*mass21*CtCab'*V) + tilde(CtCab*mass22*CtCab'*Ω)

    Hdot_θ = mul3(C_θ1', C_θ2', C_θ3', Cab*mass21*CtCab'*Vdot) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass22*CtCab'*Ωdot) +
        CtCab*mass21*Cab'*mul3(C_θ1, C_θ2, C_θ3, Vdot) + 
        CtCab*mass22*Cab'*mul3(C_θ1, C_θ2, C_θ3, Ωdot) +
        tilde(Ω - ω)*mul3(C_θ1', C_θ2', C_θ3', Cab*mass21*CtCab'*V) + 
        tilde(Ω - ω)*mul3(C_θ1', C_θ2', C_θ3', Cab*mass22*CtCab'*Ω) +
        CtCabdot*mass21*Cab'*mul3(C_θ1, C_θ2, C_θ3, V) + 
        CtCabdot*mass22*Cab'*mul3(C_θ1, C_θ2, C_θ3, Ω) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass21*CtCabdot'*V) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass22*CtCabdot'*Ω) +
        -CtCab*mass21*Cab'*mul3(C_θ1, C_θ2, C_θ3, tilde(Ω - ω)*V) + 
        -CtCab*mass22*Cab'*mul3(C_θ1, C_θ2, C_θ3, tilde(Ω - ω)*Ω)

    Hdot_V = CtCab*mass21*CtCabdot' + CtCabdot*mass21*CtCab'

    Hdot_Ω = CtCab*mass22*CtCabdot' + CtCabdot*mass22*CtCab' +
        CtCab*mass21*CtCab'*tilde(V) + CtCab*mass22*CtCab'*tilde(Ω) +
        -tilde(CtCab*mass21*CtCab'*V) - tilde(CtCab*mass22*CtCab'*Ω)

    # overwrite acceleration terms so we don't double count them
    a_u = @SMatrix zeros(3,3)
    a_ab = @SMatrix zeros(3,3)
    a_αb = @SMatrix zeros(3,3)
    α_αb = @SMatrix zeros(3,3)

    return (; properties..., Pdot_ωb, Pdot_θ, Pdot_V, Pdot_Ω, Hdot_ωb, Hdot_θ, Hdot_V, Hdot_Ω,
        a_u, a_ab, a_αb, α_αb) 
end

"""
    expanded_element_jacobian_properties(properties, x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity)

Calculate/extract the element properties needed to calculate the jacobian entries 
corresponding to a element for a steady state analysis
"""
@inline function expanded_steady_element_jacobian_properties(properties, x, indices, 
    force_scaling, structural_damping, assembly, ielem, prescribed_conditions, gravity)

    @unpack L, C1, C2, C, Cab, CtCab, Q, mass11, mass12, mass21, mass22, S11, S12, S21, S22,
        u1, u2, θ1, θ2, V1, V2, Ω1, Ω2, u, θ, V, Ω, θb, ωb, αb, Δx, ω = properties

    # linear and angular displacement
    u1_u1, θ1_θ1 = point_displacement_jacobians(assembly.start[ielem], prescribed_conditions)
    u2_u2, θ2_θ2 = point_displacement_jacobians(assembly.stop[ielem], prescribed_conditions)

    # rotation parameter matrices
    C_θ1, C_θ2, C_θ3 = get_C_θ(C, θ)
    Q_θ1, Q_θ2, Q_θ3 = get_Q_θ(Q, θ)
    Qinv_θ1, Qinv_θ2, Qinv_θ3 = get_Qinv_θ(θ)

    C1_θ1, C1_θ2, C1_θ3 = get_C_θ(C1, θ1)
    C2_θ1, C2_θ2, C2_θ3 = get_C_θ(C2, θ2)

    # strain and curvature
    γ_F = S11
    γ_M = S12
    κ_F = S21
    κ_M = S22

    # linear and angular momentum
    P_V = mass11
    P_Ω = mass12
    H_V = mass21
    H_Ω = mass22

    # linear and angular velocity
    v_vb = I3
    v_ωb = -tilde(Δx) - tilde(u)

    v_u = tilde(ωb)

    ω_ωb = I3

    # linear and angular acceleration
    a_θb = -mul3(get_C_θ(θb)..., gravity)
    a_ab = I3
    a_αb = -tilde(Δx) - tilde(u)
    a_u = tilde(αb)

    α_αb = I3

    if structural_damping

        @unpack C1, C2, Qinv1, Qinv2, u1, u2, θ1, θ2, Ω1, Ω2, ω1, ω2, μ11, μ22, udot, θdot, 
            Δx1, Δx2, Δu, Δθ, ΔQ, Δudot, Δθdot = properties

        # rotation parameter matrices
        Qinv1_θ1, Qinv1_θ2, Qinv1_θ3 = get_Qinv_θ(θ1)
        Qinv2_θ1, Qinv2_θ2, Qinv2_θ3 = get_Qinv_θ(θ2)

        # linear displacement rates
        udot1_vb = -I3
        udot2_vb = -I3

        udot1_ωb = tilde(Δx1 + u1)
        udot2_ωb = tilde(Δx2 + u2)

        udot1_u1 = -tilde(ω1)
        udot2_u2 = -tilde(ω2)

        udot1_θ1 = mul3(C1_θ1', C1_θ2', C1_θ3', V1)
        udot2_θ2 = mul3(C2_θ1', C2_θ2', C2_θ3', V2)

        udot1_V1 = C1'
        udot2_V2 = C2'

        # angular displacement rates
        θdot1_ωb = -Qinv1*C1
        θdot2_ωb = -Qinv2*C2

        θdot1_θ1 = mul3(Qinv1_θ1, Qinv1_θ2, Qinv1_θ3, Ω1 - C1*ω1) - Qinv1*mul3(C1_θ1, C1_θ2, C1_θ3, ω1)
        θdot2_θ2 = mul3(Qinv2_θ1, Qinv2_θ2, Qinv2_θ3, Ω2 - C2*ω2) - Qinv2*mul3(C2_θ1, C2_θ2, C2_θ3, ω2)
        
        θdot1_Ω1 = Qinv1
        θdot2_Ω2 = Qinv2

        θdot_ωb = (θdot1_ωb + θdot2_ωb)/2

        θdot_θ1 = θdot1_θ1/2
        θdot_θ2 = θdot2_θ2/2

        θdot_Ω1 = θdot1_Ω1/2
        θdot_Ω2 = θdot2_Ω2/2

        # change in linear displacement
        Δu_u1 = -I3
        Δu_u2 =  I3

        # change in linear and angular displacement rates
        Δudot_vb = udot2_vb - udot1_vb
        Δudot_ωb = udot2_ωb - udot1_ωb

        Δudot_u1 = -udot1_u1
        Δudot_u2 =  udot2_u2

        Δudot_θ1 = -udot1_θ1
        Δudot_θ2 =  udot2_θ2

        Δudot_V1 = -udot1_V1
        Δudot_V2 =  udot2_V2

        Δθdot_ωb = θdot2_ωb - θdot1_ωb

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

        γdot_vb = CtCab'*Δudot_vb

        γdot_ωb = -CtCab'*tilde(Δu) + CtCab'*Δudot_ωb - L*CtCab'*tilde(Cab*e1)

        tmp = CtCab'*tilde(C'*Ω - ω)
        γdot_u1 = -tmp*Δu_u1 + CtCab'*Δudot_u1 
        γdot_u2 = -tmp*Δu_u2 + CtCab'*Δudot_u2 

        tmp = -Cab'*mul3(C_θ1, C_θ2, C_θ3, tilde(C'*Ω - ω)*Δu) + 
            Cab'*mul3(C_θ1, C_θ2, C_θ3, Δudot) - 
            L*Cab'*mul3(C_θ1, C_θ2, C_θ3, tilde(C'*Ω - ω)*Cab*e1) +
            CtCab'*tilde(Δu)*mul3(C_θ1', C_θ2', C_θ3', Ω) + 
            L*CtCab'*tilde(Cab*e1)*mul3(C_θ1', C_θ2', C_θ3', Ω)
        γdot_θ1 = 1/2*tmp + CtCab'*Δudot_θ1
        γdot_θ2 = 1/2*tmp + CtCab'*Δudot_θ2

        γdot_V1 = CtCab'*Δudot_V1
        γdot_V2 = CtCab'*Δudot_V2

        γdot_Ω = CtCab'*tilde(Δu)*C' + L*CtCab'*tilde(Cab*e1)*C'

        κdot_ωb = Cab'*Q*Δθdot_ωb + Cab'*ΔQ*θdot_ωb

        tmp1 = Cab'*mul3(Q_θ1, Q_θ2, Q_θ3, Δθdot)
        tmp2 = Cab'*mul3(ΔQ_θ1, ΔQ_θ2, ΔQ_θ3, θdot) 
        tmp3 = Cab'*mul3(ΔQ_Δθ1, ΔQ_Δθ2, ΔQ_Δθ3, θdot)
        κdot_θ1 = 1/2*tmp1 + 1/2*tmp2 - tmp3 + Cab'*Q*Δθdot_θ1 + Cab'*ΔQ*θdot_θ1
        κdot_θ2 = 1/2*tmp1 + 1/2*tmp2 + tmp3 + Cab'*Q*Δθdot_θ2 + Cab'*ΔQ*θdot_θ2  

        κdot_Ω1 = Cab'*Q*Δθdot_Ω1 + Cab'*ΔQ*θdot_Ω1
        κdot_Ω2 = Cab'*Q*Δθdot_Ω2 + Cab'*ΔQ*θdot_Ω2

        # adjust strains to account for strain rates

        γ_vb = -μ11*γdot_vb
        γ_ωb = -μ11*γdot_ωb

        γ_u1 = -μ11*γdot_u1
        γ_u2 = -μ11*γdot_u2

        γ_θ1 = -μ11*γdot_θ1
        γ_θ2 = -μ11*γdot_θ2
        
        γ_V1 = -μ11*γdot_V1
        γ_V2 = -μ11*γdot_V2
        
        γ_Ω = -μ11*γdot_Ω
        
        κ_ωb = -μ22*κdot_ωb

        κ_θ1 = -μ22*κdot_θ1
        κ_θ2 = -μ22*κdot_θ2
        
        κ_Ω1 = -μ22*κdot_Ω1
        κ_Ω2 = -μ22*κdot_Ω2

    else
        
        γ_vb = @SMatrix zeros(3,3)
        γ_ωb = @SMatrix zeros(3,3)

        γ_u1 = @SMatrix zeros(3,3)
        γ_u2 = @SMatrix zeros(3,3)

        γ_θ1 = @SMatrix zeros(3,3)
        γ_θ2 = @SMatrix zeros(3,3)
        
        γ_V1 = @SMatrix zeros(3,3)
        γ_V2 = @SMatrix zeros(3,3)
        
        γ_Ω = @SMatrix zeros(3,3)
        
        κ_ωb = @SMatrix zeros(3,3)

        κ_θ1 = @SMatrix zeros(3,3)
        κ_θ2 = @SMatrix zeros(3,3)
        
        κ_Ω1 = @SMatrix zeros(3,3)
        κ_Ω2 = @SMatrix zeros(3,3)
         
    end

    return (; properties..., u1_u1, u2_u2, θ1_θ1, θ2_θ2, 
        C1_θ1, C1_θ2, C1_θ3, C2_θ1, C2_θ2, C2_θ3, C_θ1, C_θ2, C_θ3, Qinv_θ1, Qinv_θ2, Qinv_θ3,
        γ_F, γ_M, γ_vb, γ_ωb, γ_u1, γ_u2, γ_θ1, γ_θ2, γ_V1, γ_V2, γ_Ω, 
        κ_F, κ_M, κ_ωb, κ_θ1, κ_θ2, κ_Ω1, κ_Ω2,
        P_V, P_Ω, H_V, H_Ω, v_vb, v_ωb, v_u, ω_ωb, a_θb, a_ab, a_αb, α_αb, a_u
        )
end

"""
    expanded_dynamic_element_jacobian_properties(properties, x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity)

Calculate/extract the element properties needed to calculate the jacobian entries 
corresponding to a element for a dynamic analysis
"""
@inline function expanded_dynamic_element_jacobian_properties(properties, x, indices, 
    force_scaling, structural_damping, assembly, ielem, prescribed_conditions, gravity)

    properties = expanded_steady_element_jacobian_properties(properties, x, indices, 
        force_scaling, structural_damping, assembly, ielem, prescribed_conditions, gravity)

    # overwrite acceleration terms so we don't double count them
    a_u = @SMatrix zeros(3,3)
    a_ab = @SMatrix zeros(3,3)
    a_αb = @SMatrix zeros(3,3)
    α_αb = @SMatrix zeros(3,3)

    return (; properties..., a_u, a_ab, a_αb, α_αb)
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

    # displacement rates
    udot1_udot1, θdot1_θdot1 = point_displacement_jacobians(assembly.start[ielem], prescribed_conditions)
    udot2_udot2, θdot2_θdot2 = point_displacement_jacobians(assembly.stop[ielem], prescribed_conditions)

    # linear and angular momentum rates
    Pdot_Vdot = mass11
    Pdot_Ωdot = mass12
    Hdot_Vdot = mass21
    Hdot_Ωdot = mass22

    return (; udot1_udot1, udot2_udot2, θdot1_θdot1, θdot2_θdot2, 
        Pdot_Vdot, Pdot_Ωdot, Hdot_Vdot, Hdot_Ωdot)
end

@inline function static_compatability_jacobians(properties)
   
    @unpack L, Cab, CtCab, Qinv, γ, κ, u1_u1, u2_u2, θ1_θ1, θ2_θ2, γ_F, γ_M, κ_F, κ_M, 
        C_θ1, C_θ2, C_θ3, Qinv_θ1, Qinv_θ2, Qinv_θ3 = properties

    ru_u1 = -u1_u1
    ru_u2 = u2_u2

    Δu_θ = mul3(C_θ1', C_θ2', C_θ3', Cab*(L*e1 + γ))
    ru_θ1 = -1/2*Δu_θ*θ1_θ1
    ru_θ2 = -1/2*Δu_θ*θ2_θ2
    
    Δu_F = CtCab*γ_F
    ru_F = -Δu_F
    
    Δu_M = CtCab*γ_M
    ru_M = -Δu_M

    Δθ_θ = mul3(Qinv_θ1, Qinv_θ2, Qinv_θ3, Cab*κ)
    rθ_θ1 = -θ1_θ1 - 1/2*Δθ_θ*θ1_θ1
    rθ_θ2 = θ2_θ2 - 1/2*Δθ_θ*θ2_θ2
    
    Δθ_F = Qinv*Cab*κ_F
    rθ_F = -Δθ_F
    
    Δθ_M = Qinv*Cab*κ_M
    rθ_M = -Δθ_M

    return (; ru_u1, ru_u2, ru_θ1, ru_θ2, ru_F, ru_M, rθ_θ1, rθ_θ2, rθ_F, rθ_M)
end

@inline function initial_condition_compatability_jacobians(properties)
   
    jacobians = static_compatability_jacobians(properties)

    @unpack ru_u1, ru_u2, ru_θ1, ru_θ2, rθ_θ1, rθ_θ2 = jacobians

    @unpack Cab, CtCab, Qinv, u1_u1, u2_u2, θ1_θ1, θ2_θ2, γ_u1, γ_u2, γ_θ1, γ_θ2, 
        γ_u1dot, γ_u2dot, κ_θ1, κ_θ2, κ_θ1dot, κ_θ2dot = properties

    ru_u1 -= CtCab*γ_u1*u1_u1
    ru_u2 -= CtCab*γ_u2*u2_u2

    ru_θ1 -= CtCab*γ_θ1*θ1_θ1
    ru_θ2 -= CtCab*γ_θ2*θ2_θ2
    
    ru_u1dot = -CtCab*γ_u1dot
    ru_u2dot = -CtCab*γ_u2dot

    rθ_θ1 -= Qinv*Cab*κ_θ1*θ1_θ1
    rθ_θ2 -= Qinv*Cab*κ_θ2*θ2_θ2
      
    rθ_θ1dot = -Qinv*Cab*κ_θ1dot
    rθ_θ2dot = -Qinv*Cab*κ_θ2dot

    return (; jacobians..., ru_u1, ru_u2, ru_θ1, ru_θ2, ru_u1dot, ru_u2dot, rθ_θ1, rθ_θ2, 
        rθ_θ1dot, rθ_θ2dot)
end

@inline function dynamic_compatability_jacobians(properties)
   
    jacobians = static_compatability_jacobians(properties)

    @unpack Cab, CtCab, Qinv, u1_u1, u2_u2, θ1_θ1, θ2_θ2, γ_vb, γ_ωb, γ_u1, γ_u2, γ_θ1, γ_θ2, 
        γ_V1, γ_V2, γ_Ω1, γ_Ω2, κ_ωb, κ_θ1, κ_θ2, κ_Ω1, κ_Ω2 = properties

    @unpack ru_u1, ru_u2, ru_θ1, ru_θ2, rθ_θ1, rθ_θ2 = jacobians

    ru_vb = -CtCab*γ_vb

    ru_ωb = -CtCab*γ_ωb

    ru_u1 -= CtCab*γ_u1*u1_u1
    ru_u2 -= CtCab*γ_u2*u2_u2

    ru_θ1 -= CtCab*γ_θ1*θ1_θ1
    ru_θ2 -= CtCab*γ_θ2*θ2_θ2
    
    ru_V1 = -CtCab*γ_V1
    ru_V2 = -CtCab*γ_V2

    ru_Ω1 = -CtCab*γ_Ω1
    ru_Ω2 = -CtCab*γ_Ω2

    rθ_ωb = -Qinv*Cab*κ_ωb

    rθ_θ1 -= Qinv*Cab*κ_θ1*θ1_θ1
    rθ_θ2 -= Qinv*Cab*κ_θ2*θ2_θ2
    
    rθ_Ω1 = -Qinv*Cab*κ_Ω1
    rθ_Ω2 = -Qinv*Cab*κ_Ω2

    return (; jacobians..., ru_vb, ru_ωb, ru_u1, ru_u2, ru_θ1, ru_θ2, ru_V1, ru_V2, ru_Ω1, ru_Ω2, 
        rθ_ωb, rθ_θ1, rθ_θ2, rθ_Ω1, rθ_Ω2)
end

@inline function expanded_compatability_jacobians(properties)
   
    @unpack L, Cab, CtCab, Qinv, γ, κ, u1_u1, u2_u2, θ1_θ1, θ2_θ2, C_θ1, C_θ2, C_θ3, 
        Qinv_θ1, Qinv_θ2, Qinv_θ3, γ_vb, γ_ωb, γ_u1, γ_u2, γ_θ1, γ_θ2, 
        γ_V1, γ_V2, γ_Ω, γ_F, γ_M, κ_ωb, κ_θ1, κ_θ2, κ_Ω1, κ_Ω2, κ_F, κ_M = properties

    ru_vb = -CtCab*γ_vb

    ru_ωb = -CtCab*γ_ωb

    ru_u1 = -u1_u1 - CtCab*γ_u1*u1_u1
    ru_u2 =  u2_u2 - CtCab*γ_u2*u2_u2

    Δu_θ = mul3(C_θ1', C_θ2', C_θ3', Cab*(L*e1 + γ))
    ru_θ1 = -(1/2*Δu_θ + CtCab*γ_θ1)*θ1_θ1
    ru_θ2 = -(1/2*Δu_θ + CtCab*γ_θ2)*θ2_θ2
    
    ru_V1 = -CtCab*γ_V1
    ru_V2 = -CtCab*γ_V2

    ru_Ω = -CtCab*γ_Ω

    Δu_F = CtCab*γ_F
    ru_F1 = -1/2*Δu_F
    ru_F2 = -1/2*Δu_F

    Δu_M = CtCab*γ_M
    ru_M1 = -1/2*Δu_M
    ru_M2 = -1/2*Δu_M

    rθ_ωb = -Qinv*Cab*κ_ωb

    Δθ_θ = mul3(Qinv_θ1, Qinv_θ2, Qinv_θ3, Cab*κ)
    rθ_θ1 = -θ1_θ1 - 1/2*Δθ_θ*θ1_θ1 - Qinv*Cab*κ_θ1*θ1_θ1
    rθ_θ2 =  θ2_θ2 - 1/2*Δθ_θ*θ2_θ2 - Qinv*Cab*κ_θ2*θ2_θ2
    
    rθ_Ω1 = -Qinv*Cab*κ_Ω1
    rθ_Ω2 = -Qinv*Cab*κ_Ω2

    Δθ_F = Qinv*Cab*κ_F
    rθ_F1 = -1/2*Δθ_F
    rθ_F2 = -1/2*Δθ_F

    Δθ_M = Qinv*Cab*κ_M
    rθ_M1 = -1/2*Δθ_M
    rθ_M2 = -1/2*Δθ_M

    return (; ru_vb, ru_ωb, ru_u1, ru_u2, ru_θ1, ru_θ2, ru_V1, ru_V2, ru_Ω, 
        ru_F1, ru_F2, ru_M1, ru_M2, 
        rθ_ωb, rθ_θ1, rθ_θ2, rθ_Ω1, rθ_Ω2, rθ_F1, rθ_F2, rθ_M1, rθ_M2)
end

"""
    expanded_element_velocity_jacobians(properties)

Calculate the jacobians of the element velocity residuals for a constant mass matrix system.
"""
@inline function expanded_element_velocity_jacobians(properties)

    @unpack C, Cab, CtCab, Qinv, V, Ω, ω, v_vb, v_ωb, v_u, ω_ωb, u1_u1, u2_u2, θ1_θ1, θ2_θ2, 
        C_θ1, C_θ2, C_θ3, Qinv_θ1, Qinv_θ2, Qinv_θ3 = properties
    
    rV_vb = -v_vb

    rV_ωb = -v_ωb

    tmp = -v_u
    rV_u1 = 1/2*tmp*u1_u1
    rV_u2 = 1/2*tmp*u2_u2

    tmp = mul3(C_θ1', C_θ2', C_θ3', Cab*V)
    rV_θ1 = 1/2*tmp*θ1_θ1
    rV_θ2 = 1/2*tmp*θ2_θ2

    rV_V = CtCab

    rΩ_ωb = -Qinv*C*ω_ωb

    tmp = mul3(Qinv_θ1, Qinv_θ2, Qinv_θ3, Cab*Ω - C*ω) - Qinv*mul3(C_θ1, C_θ2, C_θ3, ω)
    rΩ_θ1 = 1/2*tmp*θ1_θ1
    rΩ_θ2 = 1/2*tmp*θ2_θ2

    rΩ_Ω = Qinv*Cab

    # @unpack CtCab, V, Ω, C1, V1, Ω1, C2, V2, Ω2 = properties
    # rV = CtCab*V - 1/2*(C1'*V1 + C2'*V2)
    # rΩ = CtCab*Ω - 1/2*(C1'*Ω1 + C2'*Ω2)

    return (; rV_vb, rV_ωb, rV_u1, rV_u2, rV_θ1, rV_θ2, rV_V, rΩ_ωb, rΩ_θ1, rΩ_θ2, rΩ_Ω)
end

@inline function expanded_mass_matrix_velocity_jacobians(properties)

    @unpack u1dot_u1dot, u2dot_u2dot, θ1dot_θ1dot, θ2dot_θ2dot = properties

    rV_u1dot = -u1dot_u1dot ./ 2
    rV_u2dot = -u2dot_u2dot ./ 2
    rΩ_θ1dot = -θ1dot_θ1dot ./ 2
    rΩ_θ2dot = -θ2dot_θ2dot ./ 2

    return (; rV_u1dot, rV_u2dot, rΩ_θ1dot, rΩ_θ2dot)
end

"""
    expanded_element_equilibrium_jacobians(properties)

Calculate the jacobians of the element equilibrium residuals for a constant mass matrix system.
"""
@inline function expanded_element_equilibrium_jacobians(properties, distributed_loads, ielem)

    @unpack L, C, Cab, CtCab, mass11, mass12, mass21, mass22, F1, F2, M1, M2, 
        V, Ω, P, H, F, M, γ, κ, v, ω, a, α, γ_vb, γ_ωb, γ_F, γ_M, γ_u1, γ_u2, γ_θ1, γ_θ2, γ_V1, γ_V2, γ_Ω, 
        a_θb, a_ab, a_αb, α_αb, a_u, P_V, P_Ω, H_V, H_Ω, 
        C_θ1, C_θ2, C_θ3, u1_u1, u2_u2, θ1_θ1, θ2_θ2 = properties

    # initialize equilibrium residual
    rF_F1 = -I3
    rF_F2 =  I3

    rM_M1 = -I3
    rM_M2 =  I3

    # add loads due to internal loads and stiffness
    
    tmp1 =  tilde(L*e1 + γ)
    tmp2 = -tilde(F)  

    rM_vb = tmp2*γ_vb

    rM_ωb = tmp2*γ_ωb

    rM_u1 = tmp2*γ_u1*u1_u1
    rM_u2 = tmp2*γ_u2*u2_u2

    rM_θ1 = tmp2*γ_θ1*θ1_θ1
    rM_θ2 = tmp2*γ_θ2*θ2_θ2

    rM_V1 = tmp2*γ_V1
    rM_V2 = tmp2*γ_V2

    rM_Ω = tmp2*γ_Ω

    rM_F1 = 1/2*(tmp1 + tmp2*γ_F)
    rM_F2 = 1/2*(tmp1 + tmp2*γ_F)

    rM_M1 += 1/2*tmp2*γ_M
    rM_M2 += 1/2*tmp2*γ_M

    # add loads due to linear and angular acceleration (including gravity)
    
    tmp = mass11*CtCab'*a_θb
    rF_θb = -tmp

    tmp = mass11*CtCab'*a_ab
    rF_ab = -tmp

    tmp = mass11*CtCab'*a_αb + mass12*CtCab'*α_αb
    rF_αb = -tmp

    tmp = mass11*CtCab'*a_u
    rF_u1 = -1/2*tmp*u1_u1
    rF_u2 = -1/2*tmp*u2_u2

    tmp = mass21*CtCab'*a_u    
    rM_u1 -= 1/2*tmp*u1_u1
    rM_u2 -= 1/2*tmp*u2_u2

    tmp = mass11*Cab'*mul3(C_θ1, C_θ2, C_θ3, a) + mass12*Cab'*mul3(C_θ1, C_θ2, C_θ3, α)
    
    rF_θ1 = -1/2*tmp*θ1_θ1
    rF_θ2 = -1/2*tmp*θ2_θ2

    tmp = mass21*CtCab'*a_θb
    rM_θb = -tmp

    tmp = mass21*CtCab'*a_ab
    rM_ab = -tmp

    tmp = mass21*CtCab'*a_αb + mass22*CtCab'*α_αb
    rM_αb = -tmp

    tmp = mass21*Cab'*mul3(C_θ1, C_θ2, C_θ3, a) + mass22*Cab'*mul3(C_θ1, C_θ2, C_θ3, α)
    rM_θ1 -= 1/2*tmp*θ1_θ1
    rM_θ2 -= 1/2*tmp*θ2_θ2

    # add loads due to linear and angular momentum
    rF_V = -tilde(Ω)*P_V 
    rF_Ω = -tilde(Ω)*P_Ω + tilde(P)

    rM_V = -tilde(Ω)*H_V - tilde(V)*P_V + tilde(P)
    rM_Ω -= tilde(Ω)*H_Ω + tilde(V)*P_Ω - tilde(H)

    # add distributed loads
    if haskey(distributed_loads, ielem)
        dload = distributed_loads[ielem]

        tmp = Cab'*mul3(C_θ1, C_θ2, C_θ3, dload.f1 + dload.f2)
        rF_θ1 += 1/2*tmp*θ1_θ1
        rF_θ2 += 1/2*tmp*θ2_θ2

        tmp = Cab'*mul3(C_θ1, C_θ2, C_θ3, dload.m1 + dload.m2)
        rM_θ1 += 1/2*tmp*θ1_θ1
        rM_θ2 += 1/2*tmp*θ2_θ2

    end

    return (; rF_θb, rF_ab, rF_αb, rF_F1, rF_F2, rM_F1, rM_F2, rM_M1, rM_M2, 
        rF_u1, rF_u2, rF_θ1, rF_θ2, rF_V, rF_Ω,  
        rM_θb, rM_vb, rM_ωb, rM_ab, rM_αb, rM_u1, rM_u2, rM_θ1, rM_θ2, rM_V1, rM_V2, rM_V, rM_Ω)
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

"""
    static_element_resultant_jacobians(properties, distributed_loads, ielem)

Calculate the jacobians for the resultant loads applied at each end of a beam element 
for a static analysis.
"""
@inline function static_element_resultant_jacobians(properties, distributed_loads, ielem)

    @unpack L, Cab, CtCab, mass11, mass12, mass21, mass22, F, M, γ, κ, a, α, 
        C_θ1, C_θ2, C_θ3, θ1_θ1, θ2_θ2, γ_F, γ_M = properties

    # add loads due to internal loads and stiffness
    
    tmp = mul3(C_θ1', C_θ2', C_θ3', Cab*F)

    F1_θ1 = 1/2*tmp*θ1_θ1
    F2_θ1 = 1/2*tmp*θ1_θ1
    
    F1_θ2 = 1/2*tmp*θ2_θ2
    F2_θ2 = 1/2*tmp*θ2_θ2

    F1_F = CtCab
    F2_F = CtCab

    tmp1 = mul3(C_θ1', C_θ2', C_θ3', Cab*M)
    tmp2 = 1/2*mul3(C_θ1', C_θ2', C_θ3', Cab*cross(L*e1 + γ, F))
    
    M1_θ1 = 1/2*(tmp1 + tmp2)*θ1_θ1
    M2_θ1 = 1/2*(tmp1 - tmp2)*θ1_θ1

    M1_θ2 = 1/2*(tmp1 + tmp2)*θ2_θ2
    M2_θ2 = 1/2*(tmp1 - tmp2)*θ2_θ2

    tmp = 1/2*CtCab*(tilde(L*e1 + γ) - tilde(F)*γ_F)
    M1_F = tmp
    M2_F = -tmp

    tmp = 1/2*CtCab*tilde(F)*γ_M
    M1_M = CtCab - tmp
    M2_M = CtCab + tmp
    
    # add loads due to linear and angular acceleration (including gravity)
    
    tmp = 1/2*(
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass11*CtCab'*a) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass12*CtCab'*α) + 
        CtCab*mass11*Cab'*mul3(C_θ1, C_θ2, C_θ3, a) + 
        CtCab*mass12*Cab'*mul3(C_θ1, C_θ2, C_θ3, α)
    )
    
    F1_θ1 -= 1/2*tmp*θ1_θ1
    F2_θ1 += 1/2*tmp*θ1_θ1
    
    F1_θ2 -= 1/2*tmp*θ2_θ2
    F2_θ2 += 1/2*tmp*θ2_θ2

    tmp = 1/2*(
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass21*CtCab'*a) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass22*CtCab'*α) + 
        CtCab*mass21*Cab'*mul3(C_θ1, C_θ2, C_θ3, a) + 
        CtCab*mass22*Cab'*mul3(C_θ1, C_θ2, C_θ3, α)
    )
    
    M1_θ1 -= 1/2*tmp*θ1_θ1
    M2_θ1 += 1/2*tmp*θ1_θ1

    M1_θ2 -= 1/2*tmp*θ2_θ2
    M2_θ2 += 1/2*tmp*θ2_θ2

    # add distributed loads
    if haskey(distributed_loads, ielem)
        dload = distributed_loads[ielem]

        tmp = mul3(C_θ1', C_θ2', C_θ3', dload.f1_follower)
        F1_θ1 += 1/2*tmp*θ1_θ1
        F1_θ2 += 1/2*tmp*θ2_θ2

        tmp = mul3(C_θ1', C_θ2', C_θ3', dload.f2_follower)
        F2_θ1 -= 1/2*tmp*θ1_θ1
        F2_θ2 -= 1/2*tmp*θ2_θ2

        tmp = mul3(C_θ1', C_θ2', C_θ3', dload.m1_follower)
        M1_θ1 += 1/2*tmp*θ1_θ1
        M1_θ2 += 1/2*tmp*θ2_θ2

        tmp = mul3(C_θ1', C_θ2', C_θ3', dload.m2_follower)
        M2_θ1 -= 1/2*tmp*θ1_θ1
        M2_θ2 -= 1/2*tmp*θ2_θ2
    end

    return (; F1_θ1, F1_θ2, F1_F, F2_θ1, F2_θ2, F2_F,
        M1_θ1, M1_θ2, M1_F, M1_M, M2_θ1, M2_θ2, M2_F, M2_M)
end

"""
    steady_state_element_resultant_jacobians(properties, distributed_loads, ielem)

Calculate the jacobians for the resultant loads applied at each end of a beam element 
for a steady state analysis.
"""
@inline function steady_state_element_resultant_jacobians(properties, distributed_loads, ielem)

    jacobians = static_element_resultant_jacobians(properties, distributed_loads, ielem)

    @unpack L, CtCab, mass11, mass12, mass21, mass22, F, V, Ω, P, H, ω, 
        γ_vb, γ_ωb, γ_u1, γ_u2, γ_θ1, γ_θ2, γ_V1, γ_V2, γ_Ω1, γ_Ω2,
        ω_ωb, a_u, a_θb, a_ab, a_αb, α_αb, P_θ, P_V, P_Ω, H_θ, H_V, H_Ω, u1_u1, u2_u2, θ1_θ1, θ2_θ2 = properties

    @unpack F1_θ1, F1_θ2, F2_θ1, F2_θ2, M1_θ1, M1_θ2, M2_θ1, M2_θ2 = jacobians

    # add loads due to internal forces/moments and stiffness
    tmp = 1/2*CtCab*tilde(F)
   
    M1_vb = -tmp*γ_vb
    M2_vb = tmp*γ_vb

    M1_ωb = -tmp*γ_ωb
    M2_ωb = tmp*γ_ωb  
    
    M1_u1 = -tmp*γ_u1*u1_u1
    M2_u1 = tmp*γ_u1*u1_u1

    M1_u2 = -tmp*γ_u2*u2_u2
    M2_u2 =  tmp*γ_u2*u2_u2

    M1_θ1 -= tmp*γ_θ1*θ1_θ1
    M2_θ1 += tmp*γ_θ1*θ1_θ1

    M1_θ2 -= tmp*γ_θ2*θ2_θ2
    M2_θ2 += tmp*γ_θ2*θ2_θ2

    M1_V1 = -tmp*γ_V1
    M2_V1 =  tmp*γ_V1

    M1_V2 = -tmp*γ_V2
    M2_V2 =  tmp*γ_V2

    M1_Ω1 = -tmp*γ_Ω1
    M2_Ω1 =  tmp*γ_Ω1

    M1_Ω2 = -tmp*γ_Ω2
    M2_Ω2 =  tmp*γ_Ω2

    # add loads due to linear and angular acceleration (including gravity)
    
    tmp = 1/2*CtCab*mass11*CtCab'*a_θb
    F1_θb = -tmp
    F2_θb = tmp

    tmp = 1/2*CtCab*mass11*CtCab'*a_ab
    F1_ab = -tmp
    F2_ab = tmp

    tmp = 1/2*(CtCab*mass11*CtCab'*a_αb + CtCab*mass12*CtCab'*α_αb)
    F1_αb = -tmp
    F2_αb = tmp
    
    tmp = 1/2*CtCab*mass11*CtCab'*a_u
    
    F1_u1 = -1/2*tmp*u1_u1
    F2_u1 = 1/2*tmp*u1_u1

    F1_u2 = -1/2*tmp*u2_u2    
    F2_u2 = 1/2*tmp*u2_u2

    tmp = 1/2*CtCab*mass21*CtCab'*a_θb
    M1_θb = -tmp
    M2_θb = tmp

    tmp = 1/2*CtCab*mass21*CtCab'*a_ab
    M1_ab = -tmp
    M2_ab = tmp

    tmp = 1/2*(CtCab*mass21*CtCab'*a_αb + CtCab*mass22*CtCab'*α_αb)
    M1_αb = -tmp
    M2_αb = tmp

    tmp = 1/2*CtCab*mass21*CtCab'*a_u
    
    M1_u1 -= 1/2*tmp*u1_u1
    M2_u1 += 1/2*tmp*u1_u1

    M1_u2 -= 1/2*tmp*u2_u2
    M2_u2 += 1/2*tmp*u2_u2

    # add loads due to linear and angular momentum

    tmp = 1/2*tilde(P)*ω_ωb

    F1_ωb = tmp
    F2_ωb = -tmp

    tmp = 1/2*tilde(ω)*P_θ

    F1_θ1 -= 1/2*tmp*θ1_θ1
    F2_θ1 += 1/2*tmp*θ1_θ1
    
    F1_θ2 -= 1/2*tmp*θ2_θ2
    F2_θ2 += 1/2*tmp*θ2_θ2

    tmp = 1/2*tilde(ω)*P_V

    F1_V1 = -1/2*tmp
    F2_V1 = 1/2*tmp

    F1_V2 = -1/2*tmp
    F2_V2 = 1/2*tmp

    tmp = 1/2*tilde(ω)*P_Ω

    F1_Ω1 = -1/2*tmp
    F2_Ω1 = 1/2*tmp

    F1_Ω2 = -1/2*tmp
    F2_Ω2 = 1/2*tmp

    tmp = 1/2*tilde(H)*ω_ωb

    M1_ωb += tmp
    M2_ωb -= tmp

    tmp = 1/2*(tilde(ω)*H_θ + tilde(V)*P_θ)
    
    M1_θ1 -= 1/2*tmp*θ1_θ1
    M2_θ1 += 1/2*tmp*θ1_θ1

    M1_θ2 -= 1/2*tmp*θ2_θ2
    M2_θ2 += 1/2*tmp*θ2_θ2

    tmp = 1/2*(tilde(ω)*H_V + tilde(V)*P_V - tilde(P))
    
    M1_V1 -= 1/2*tmp
    M2_V1 += 1/2*tmp

    M1_V2 -= 1/2*tmp
    M2_V2 += 1/2*tmp

    tmp = 1/2*(tilde(ω)*H_Ω + tilde(V)*P_Ω)
    
    M1_Ω1 -= 1/2*tmp
    M2_Ω1 += 1/2*tmp

    M1_Ω2 -= 1/2*tmp
    M2_Ω2 += 1/2*tmp

    return (; jacobians..., 
        F1_θb, F1_ωb, F1_ab, F1_αb, F1_u1, F1_u2, F1_θ1, F1_θ2, F1_V1, F1_V2, F1_Ω1, F1_Ω2, 
        F2_θb, F2_ωb, F2_ab, F2_αb, F2_u1, F2_u2, F2_θ1, F2_θ2, F2_V1, F2_V2, F2_Ω1, F2_Ω2, 
        M1_θb, M1_vb, M1_ωb, M1_ab, M1_αb, M1_u1, M1_u2, M1_θ1, M1_θ2, M1_V1, M1_V2, M1_Ω1, M1_Ω2, 
        M2_θb, M2_vb, M2_ωb, M2_ab, M2_αb, M2_u1, M2_u2, M2_θ1, M2_θ2, M2_V1, M2_V2, M2_Ω1, M2_Ω2)

end

"""
    initial_condition_element_resultant_jacobians(properties, distributed_loads, ielem)

Calculate the jacobians for the resultant loads applied at each end of V beam element 
for the initialization of a time domain analysis.
"""
@inline function initial_condition_element_resultant_jacobians(properties, distributed_loads, ielem)

    jacobians = static_element_resultant_jacobians(properties, distributed_loads, ielem)

    @unpack L, CtCab, mass11, mass12, mass21, mass22, F, V, Ω, P, H, ω, 
        γ_u1, γ_u2, γ_θ1, γ_θ2, γ_u1dot, γ_u2dot, V_u, V_vb, V_ωb, Ω_ωb, 
        P_u, P_θ, P_vb, P_ωb, H_u, H_θ, H_vb, H_ωb, 
        u1_u1, u2_u2, θ1_θ1, θ2_θ2, 
        V1dot_V1dot, V2dot_V2dot, Ω1dot_Ω1dot, Ω2dot_Ω2dot,
        Pdot_vb, Pdot_ωb, Pdot_ab, Pdot_αb, Pdot_u, Pdot_θ, Pdot_Vdot, Pdot_Ωdot,
        Hdot_vb, Hdot_ωb, Hdot_ab, Hdot_αb, Hdot_u, Hdot_θ, Hdot_Vdot, Hdot_Ωdot, 
        v_vb, v_ωb, ω_ωb, a_u, a_θb, a_ab, a_αb, α_αb = properties

    @unpack F1_θ1, F1_θ2, F2_θ1, F2_θ2, M1_θ1, M1_θ2, M2_θ1, M2_θ2 = jacobians

    # add loads due to internal forces/moments and stiffness
    tmp = 1/2*CtCab*tilde(F)
    
    M1_u1 = -tmp*γ_u1*u1_u1
    M2_u1 = tmp*γ_u1*u1_u1

    M1_u2 = -tmp*γ_u2*u2_u2
    M2_u2 =  tmp*γ_u2*u2_u2

    M1_θ1 -= tmp*γ_θ1*θ1_θ1
    M2_θ1 += tmp*γ_θ1*θ1_θ1

    M1_θ2 -= tmp*γ_θ2*θ2_θ2
    M2_θ2 += tmp*γ_θ2*θ2_θ2

    M1_u1dot = -tmp*γ_u1dot
    M2_u1dot =  tmp*γ_u1dot

    M1_u2dot = -tmp*γ_u2dot
    M2_u2dot =  tmp*γ_u2dot

    # add loads due to linear and angular acceleration (including gravity)
    
    tmp = 1/2*CtCab*mass11*CtCab'*a_θb
    F1_θb = -tmp
    F2_θb = tmp

    tmp = 1/2*CtCab*mass11*CtCab'*a_ab
    F1_ab = -tmp
    F2_ab = tmp

    tmp = 1/2*(CtCab*mass11*CtCab'*a_αb + CtCab*mass12*CtCab'*α_αb)
    F1_αb = -tmp
    F2_αb = tmp
    
    tmp = 1/2*CtCab*mass11*CtCab'*a_u
    
    F1_u1 = -1/2*tmp*u1_u1
    F2_u1 = 1/2*tmp*u1_u1

    F1_u2 = -1/2*tmp*u2_u2    
    F2_u2 = 1/2*tmp*u2_u2

    tmp = 1/2*CtCab*mass21*CtCab'*a_θb
    M1_θb = -tmp
    M2_θb = tmp

    tmp = 1/2*CtCab*mass21*CtCab'*a_ab
    M1_ab = -tmp
    M2_ab = tmp

    tmp = 1/2*(CtCab*mass21*CtCab'*a_αb + CtCab*mass22*CtCab'*α_αb)
    M1_αb = -tmp
    M2_αb = tmp

    tmp = 1/2*CtCab*mass21*CtCab'*a_u
    
    M1_u1 -= 1/2*tmp*u1_u1
    M2_u1 += 1/2*tmp*u1_u1

    M1_u2 -= 1/2*tmp*u2_u2
    M2_u2 += 1/2*tmp*u2_u2

    # add loads due to linear and angular momentum

    tmp = 1/2*tilde(ω)*P_vb
    F1_vb = -tmp
    F2_vb = tmp

    tmp = 1/2*(tilde(ω)*P_ωb - tilde(P)*ω_ωb)
    F1_ωb = -tmp
    F2_ωb = tmp

    tmp = 1/2*tilde(ω)*P_u

    F1_u1 -= 1/2*tmp*u1_u1
    F2_u1 += 1/2*tmp*u1_u1
    
    F1_u2 -= 1/2*tmp*u2_u2
    F2_u2 += 1/2*tmp*u2_u2

    tmp = 1/2*tilde(ω)*P_θ

    F1_θ1 -= 1/2*tmp*θ1_θ1
    F2_θ1 += 1/2*tmp*θ1_θ1
    
    F1_θ2 -= 1/2*tmp*θ2_θ2
    F2_θ2 += 1/2*tmp*θ2_θ2

    tmp = 1/2*(tilde(ω)*H_vb + tilde(V)*P_vb - tilde(P)*V_vb)

    M1_vb = -tmp
    M2_vb = tmp

    tmp = 1/2*(tilde(ω)*H_ωb + tilde(V)*P_ωb - tilde(P)*V_ωb - tilde(H))

    M1_ωb = -tmp
    M2_ωb = tmp

    tmp = 1/2*(tilde(ω)*H_u + tilde(V)*P_u - tilde(P)*V_u)

    M1_u1 -= 1/2*tmp*u1_u1
    M2_u1 += 1/2*tmp*u1_u1

    M1_u2 -= 1/2*tmp*u2_u2
    M2_u2 += 1/2*tmp*u2_u2

    tmp = 1/2*(tilde(ω)*H_θ + tilde(V)*P_θ)
    
    M1_θ1 -= 1/2*tmp*θ1_θ1
    M2_θ1 += 1/2*tmp*θ1_θ1

    M1_θ2 -= 1/2*tmp*θ2_θ2
    M2_θ2 += 1/2*tmp*θ2_θ2

    # add loads due to linear and angular momentum rates

    tmp = 1/2*Pdot_vb

    F1_vb -= tmp
    F2_vb += tmp

    tmp = 1/2*Pdot_ωb

    F1_ωb -= tmp
    F2_ωb += tmp

    tmp = 1/2*Pdot_ab

    F1_ab -= tmp
    F2_ab += tmp

    tmp = 1/2*Pdot_αb

    F1_αb -= tmp
    F2_αb += tmp

    tmp = 1/2*Pdot_u

    F1_u1 -= 1/2*tmp*u1_u1
    F2_u1 += 1/2*tmp*u1_u1

    F1_u2 -= 1/2*tmp*u2_u2
    F2_u2 += 1/2*tmp*u2_u2

    tmp = 1/2*Pdot_θ

    F1_θ1 -= 1/2*tmp*θ1_θ1
    F2_θ1 += 1/2*tmp*θ1_θ1

    F1_θ1 -= 1/2*tmp*θ1_θ1
    F2_θ1 += 1/2*tmp*θ1_θ1

    tmp = 1/2*Pdot_Vdot

    F1_V1dot =  -1/2*tmp*V1dot_V1dot
    F2_V1dot =  1/2*tmp*V1dot_V1dot

    F1_V2dot = -1/2*tmp*V2dot_V2dot
    F2_V2dot =  1/2*tmp*V2dot_V2dot

    tmp = 1/2*Pdot_Ωdot

    F1_Ω1dot = -1/2*tmp*Ω1dot_Ω1dot
    F2_Ω1dot =  1/2*tmp*Ω1dot_Ω1dot

    F1_Ω2dot = -1/2*tmp*Ω2dot_Ω2dot
    F2_Ω2dot =  1/2*tmp*Ω2dot_Ω2dot

    tmp = 1/2*Hdot_vb

    M1_vb -= tmp
    M2_vb += tmp

    tmp = 1/2*Hdot_ωb

    M1_ωb -= tmp
    M2_ωb += tmp

    tmp = 1/2*Hdot_ab

    M1_ab -= tmp
    M2_ab += tmp

    tmp = 1/2*Hdot_αb

    M1_αb -= tmp
    M2_αb += tmp

    tmp = 1/2*Hdot_u

    M1_u1 -= 1/2*tmp*u1_u1
    M2_u1 += 1/2*tmp*u1_u1

    M1_u2 -= 1/2*tmp*u2_u2
    M2_u2 += 1/2*tmp*u2_u2

    tmp = 1/2*Hdot_θ

    M1_θ1 -= 1/2*tmp*θ1_θ1
    M2_θ1 += 1/2*tmp*θ1_θ1

    M1_θ1 -= 1/2*tmp*θ1_θ1
    M2_θ1 += 1/2*tmp*θ1_θ1

    tmp = 1/2*Hdot_Vdot

    M1_V1dot = -1/2*tmp*V1dot_V1dot
    M2_V1dot =  1/2*tmp*V1dot_V1dot

    M1_V2dot = -1/2*tmp*V2dot_V2dot
    M2_V2dot =  1/2*tmp*V2dot_V2dot

    tmp = 1/2*Hdot_Ωdot

    M1_Ω1dot = -1/2*tmp*Ω1dot_Ω1dot
    M2_Ω1dot =  1/2*tmp*Ω1dot_Ω1dot

    M1_Ω2dot = -1/2*tmp*Ω2dot_Ω2dot
    M2_Ω2dot =  1/2*tmp*Ω2dot_Ω2dot

    return (; jacobians..., 
        F1_θb, F1_vb, F1_ωb, F1_ab, F1_αb, F1_u1, F1_u2, F1_θ1, F1_θ2, F1_V1dot, F1_V2dot, F1_Ω1dot, F1_Ω2dot, 
        F2_θb, F2_vb, F2_ωb, F2_ab, F2_αb, F2_u1, F2_u2, F2_θ1, F2_θ2, F2_V1dot, F2_V2dot, F2_Ω1dot, F2_Ω2dot, 
        M1_θb, M1_vb, M1_ωb, M1_ab, M1_αb, M1_u1, M1_u2, M1_θ1, M1_θ2, M1_u1dot, M1_u2dot, M1_V1dot, M1_V2dot, M1_Ω1dot, M1_Ω2dot, 
        M2_θb, M2_vb, M2_ωb, M2_ab, M2_αb, M2_u1, M2_u2, M2_θ1, M2_θ2, M2_u1dot, M2_u2dot, M2_V1dot, M2_V2dot, M2_Ω1dot, M2_Ω2dot)
end

"""
    newmark_element_resultant_jacobians(properties, distributed_loads, ielem)

Calculate the jacobians for the resultant loads applied at each end of a beam element 
for a Newmark scheme time marching analysis.
"""
@inline function newmark_element_resultant_jacobians(properties, distributed_loads, ielem)

    jacobians = steady_state_element_resultant_jacobians(properties, distributed_loads, ielem)

    @unpack L, Pdot_ωb, Pdot_θ, Pdot_V, Pdot_Ω, Hdot_ωb, Hdot_θ, Hdot_V, Hdot_Ω, θ1_θ1, θ2_θ2 = properties
    
    @unpack F1_ωb, F1_θ1, F1_θ2, F1_V1, F1_V2, F1_Ω1, F1_Ω2, 
            F2_ωb, F2_θ1, F2_θ2, F2_V1, F2_V2, F2_Ω1, F2_Ω2, 
            M1_ωb, M1_θ1, M1_θ2, M1_V1, M1_V2, M1_Ω1, M1_Ω2, 
            M2_ωb, M2_θ1, M2_θ2, M2_V1, M2_V2, M2_Ω1, M2_Ω2 = jacobians

    # add loads due to linear and angular momentum rates

    tmp = 1/2*Pdot_ωb

    F1_ωb -= tmp
    F2_ωb += tmp

    tmp = 1/2*Pdot_θ

    F1_θ1 -= 1/2*tmp*θ1_θ1
    F2_θ1 += 1/2*tmp*θ1_θ1

    F1_θ2 -= 1/2*tmp*θ2_θ2
    F2_θ2 += 1/2*tmp*θ2_θ2

    tmp = 1/2*Pdot_V

    F1_V1 -= 1/2*tmp
    F2_V1 += 1/2*tmp

    F1_V2 -= 1/2*tmp
    F2_V2 += 1/2*tmp

    tmp = 1/2*Pdot_Ω

    F1_Ω1 -= 1/2*tmp
    F2_Ω1 += 1/2*tmp

    F1_Ω2 -= 1/2*tmp
    F2_Ω2 += 1/2*tmp

    tmp = 1/2*Hdot_ωb

    M1_ωb -= tmp
    M2_ωb += tmp

    tmp = 1/2*Hdot_θ

    M1_θ1 -= 1/2*tmp*θ1_θ1
    M2_θ1 += 1/2*tmp*θ1_θ1

    M1_θ2 -= 1/2*tmp*θ2_θ2
    M2_θ2 += 1/2*tmp*θ2_θ2

    tmp = 1/2*Hdot_V

    M1_V1 -= 1/2*tmp
    M2_V1 += 1/2*tmp

    M1_V2 -= 1/2*tmp
    M2_V2 += 1/2*tmp

    tmp = 1/2*Hdot_Ω

    M1_Ω1 -= 1/2*tmp
    M2_Ω1 += 1/2*tmp

    M1_Ω2 -= 1/2*tmp
    M2_Ω2 += 1/2*tmp

    return (; jacobians..., 
        F1_ωb, F1_θ1, F1_θ2, F1_V1, F1_V2, F1_Ω1, F1_Ω2, 
        F2_ωb, F2_θ1, F2_θ2, F2_V1, F2_V2, F2_Ω1, F2_Ω2, 
        M1_ωb, M1_θ1, M1_θ2, M1_V1, M1_V2, M1_Ω1, M1_Ω2, 
        M2_ωb, M2_θ1, M2_θ2, M2_V1, M2_V2, M2_Ω1, M2_Ω2)
end

"""
    dynamic_element_resultant_jacobians(properties, distributed_loads, ielem)

Calculate the jacobians for the resultant loads applied at each end of a beam element 
for a dynamic analysis.
"""
@inline function dynamic_element_resultant_jacobians(properties, distributed_loads, ielem)

    return newmark_element_resultant_jacobians(properties, distributed_loads, ielem)
end

"""
    expanded_element_resultant_jacobians(properties)

Calculate the jacobians for the resultant loads applied at each end of a beam element 
for a constant mass matrix system.
"""
@inline function expanded_element_resultant_jacobians(properties)
   
    @unpack Cab, C1, C2, CtCab, F1, F2, M1, M2, C1_θ1, C1_θ2, C1_θ3, 
        C2_θ1, C2_θ2, C2_θ3, C_θ1, C_θ2, C_θ3, θ1_θ1, θ2_θ2 = properties

    F1_θ = C1*mul3(C_θ1', C_θ2', C_θ3', Cab*F1)
    F1_θ1 = 1/2*F1_θ*θ1_θ1 + mul3(C1_θ1, C1_θ2, C1_θ3, CtCab*F1)*θ1_θ1
    F1_θ2 = 1/2*F1_θ*θ2_θ2
    F1_F1 = C1*CtCab

    F2_θ = C2*mul3(C_θ1', C_θ2', C_θ3', Cab*F2)
    F2_θ1 = 1/2*F2_θ*θ1_θ1
    F2_θ2 = 1/2*F2_θ*θ2_θ2 + mul3(C2_θ1, C2_θ2, C2_θ3, CtCab*F2)*θ2_θ2 
    F2_F2 = C2*CtCab

    M1_θ = C1*mul3(C_θ1', C_θ2', C_θ3', Cab*M1)
    M1_θ1 = 1/2*M1_θ*θ1_θ1 + mul3(C1_θ1, C1_θ2, C1_θ3, CtCab*M1)*θ1_θ1 
    M1_θ2 = 1/2*M1_θ*θ2_θ2
    M1_M1 = C1*CtCab

    M2_θ = C2*mul3(C_θ1', C_θ2', C_θ3', Cab*M2)
    M2_θ1 = 1/2*M2_θ*θ1_θ1
    M2_θ2 = 1/2*M2_θ*θ2_θ2 + mul3(C2_θ1, C2_θ2, C2_θ3, CtCab*M2)*θ2_θ2
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

"""
    expanded_mass_matrix_element_velocity_jacobians(properties)

Calculate the mass matrix jacobians of the velocity residuals `rV` and `rΩ` of an element
"""
@inline function expanded_mass_matrix_element_velocity_jacobians(properties)

    @unpack udot1_udot1, udot2_udot2, θdot1_θdot1, θdot2_θdot2 = properties

    rV_u1dot = -udot1_udot1/2
    rV_u2dot = -udot2_udot2/2
    rΩ_θ1dot = -θdot1_θdot1/2
    rΩ_θ2dot = -θdot2_θdot2/2

    return (; rV_u1dot, rV_u2dot, rΩ_θ1dot, rΩ_θ2dot)
end

@inline function insert_static_element_jacobians!(jacob, indices, force_scaling, 
    assembly, ielem, compatability, resultants)

    @unpack ru_u1, ru_u2, ru_θ1, ru_θ2, ru_F, ru_M, 
                          rθ_θ1, rθ_θ2, rθ_F, rθ_M = compatability

    @unpack F1_θ1, F1_θ2, F1_F,
            F2_θ1, F2_θ2, F2_F,
            M1_θ1, M1_θ2, M1_F, M1_M,
            M2_θ1, M2_θ2, M2_F, M2_M = resultants

    icol = indices.icol_elem[ielem]
    icol1 = indices.icol_point[assembly.start[ielem]]
    icol2 = indices.icol_point[assembly.stop[ielem]]

    # compatability equations
    irow = indices.irow_elem[ielem]

    jacob[irow:irow+2, icol1:icol1+2] .= ru_u1
    jacob[irow:irow+2, icol1+3:icol1+5] .= ru_θ1

    jacob[irow:irow+2, icol2:icol2+2] .= ru_u2
    jacob[irow:irow+2, icol2+3:icol2+5] .= ru_θ2

    jacob[irow:irow+2, icol:icol+2] .= ru_F .* force_scaling
    jacob[irow:irow+2, icol+3:icol+5] .= ru_M .* force_scaling

    jacob[irow+3:irow+5, icol1+3:icol1+5] .= rθ_θ1

    jacob[irow+3:irow+5, icol2+3:icol2+5] .= rθ_θ2

    jacob[irow+3:irow+5, icol:icol+2] .= rθ_F .* force_scaling
    jacob[irow+3:irow+5, icol+3:icol+5] .= rθ_M .* force_scaling

    # equilibrium equations for the start of the beam element
    irow1 = indices.irow_point[assembly.start[ielem]]

    @views jacob[irow1:irow1+2, icol1+3:icol1+5] .-= F1_θ1 ./ force_scaling
    @views jacob[irow1:irow1+2, icol2+3:icol2+5] .-= F1_θ2 ./ force_scaling

    jacob[irow1:irow1+2, icol:icol+2] .= -F1_F

    @views jacob[irow1+3:irow1+5, icol1+3:icol1+5] .-= M1_θ1 ./ force_scaling
    @views jacob[irow1+3:irow1+5, icol2+3:icol2+5] .-= M1_θ2 ./ force_scaling

    jacob[irow1+3:irow1+5, icol:icol+2] .= -M1_F
    jacob[irow1+3:irow1+5, icol+3:icol+5] .= -M1_M

    # equilibrium equations for the end of the beam element
    irow2 = indices.irow_point[assembly.stop[ielem]]
    
    @views jacob[irow2:irow2+2, icol1+3:icol1+5] .+= F2_θ1 ./ force_scaling
    @views jacob[irow2:irow2+2, icol2+3:icol2+5] .+= F2_θ2 ./ force_scaling

    jacob[irow2:irow2+2, icol:icol+2] .= F2_F

    @views jacob[irow2+3:irow2+5, icol1+3:icol1+5] .+= M2_θ1./ force_scaling
    @views jacob[irow2+3:irow2+5, icol2+3:icol2+5] .+= M2_θ2 ./ force_scaling

    jacob[irow2+3:irow2+5, icol:icol+2] .= M2_F
    jacob[irow2+3:irow2+5, icol+3:icol+5] .= M2_M

    return jacob
end

@inline function insert_initial_condition_element_jacobians!(jacob, indices, force_scaling, 
    assembly, ielem, compatability, resultants)

    insert_static_element_jacobians!(jacob, indices, force_scaling, assembly, ielem, 
        compatability, resultants)

    @unpack ru_u1dot, ru_u2dot, rθ_θ1dot, rθ_θ2dot = compatability

    @unpack F1_θb, F1_vb, F1_ωb, F1_ab, F1_αb, F1_u1, F1_u2, F1_V1dot, F1_V2dot, F1_Ω1dot, F1_Ω2dot,
            F2_θb, F2_vb, F2_ωb, F2_ab, F2_αb, F2_u1, F2_u2, F2_V1dot, F2_V2dot, F2_Ω1dot, F2_Ω2dot,
            M1_θb, M1_vb, M1_ωb, M1_ab, M1_αb, M1_u1, M1_u2, M1_θ1, M1_θ2, M1_u1dot, M1_u2dot, M1_V1dot, M1_V2dot, M1_Ω1dot, M1_Ω2dot, 
            M2_θb, M2_vb, M2_ωb, M2_ab, M2_αb, M2_u1, M2_u2, M2_θ1, M2_θ2, M2_u1dot, M2_u2dot, M2_V1dot, M2_V2dot, M2_Ω1dot, M2_Ω2dot = resultants

    icol1 = indices.icol_point[assembly.start[ielem]]
    icol2 = indices.icol_point[assembly.stop[ielem]]

    # compatability equations
    irow = indices.irow_elem[ielem]

    jacob[irow:irow+2, icol1+6:icol1+8] .= ru_u1dot

    jacob[irow:irow+2, icol2+6:icol2+8] .= ru_u2dot

    jacob[irow+3:irow+5, icol1+9:icol1+11] .= rθ_θ1dot

    jacob[irow+3:irow+5, icol2+9:icol2+11] .= rθ_θ2dot

    # equilibrium equations for the start of the beam element
    irow1 = indices.irow_point[assembly.start[ielem]]

    @views jacob[irow1:irow1+2, 4:6] .-= F1_θb ./ force_scaling
    @views jacob[irow1:irow1+2, 7:9] .-= F1_vb ./ force_scaling
    @views jacob[irow1:irow1+2, 10:12] .-= F1_ωb ./ force_scaling
    @views jacob[irow1:irow1+2, 13:15] .-= F1_ab ./ force_scaling
    @views jacob[irow1:irow1+2, 16:18] .-= F1_αb ./ force_scaling

    @views jacob[irow1:irow1+2, icol1:icol1+2] .-= F1_u1 ./ force_scaling
    @views jacob[irow1:irow1+2, icol1:icol1+2] .-= F1_V1dot ./ force_scaling
    @views jacob[irow1:irow1+2, icol1+3:icol1+5] .-= F1_Ω1dot ./ force_scaling

    @views jacob[irow1:irow1+2, icol2:icol2+2] .-= F1_u2 ./ force_scaling
    @views jacob[irow1:irow1+2, icol2:icol2+2] .-= F1_V2dot ./ force_scaling
    @views jacob[irow1:irow1+2, icol2+3:icol2+5] .-= F1_Ω2dot ./ force_scaling

    @views jacob[irow1+3:irow1+5, 4:6] .-= M1_θb ./ force_scaling
    @views jacob[irow1+3:irow1+5, 7:9] .-= M1_vb ./ force_scaling
    @views jacob[irow1+3:irow1+5, 10:12] .-= M1_ωb ./ force_scaling
    @views jacob[irow1+3:irow1+5, 13:15] .-= M1_ab ./ force_scaling
    @views jacob[irow1+3:irow1+5, 16:18] .-= M1_αb ./ force_scaling

    @views jacob[irow1+3:irow1+5, icol1:icol1+2] .-= M1_u1 ./ force_scaling
    @views jacob[irow1+3:irow1+5, icol1:icol1+2] .-= M1_V1dot ./ force_scaling
    @views jacob[irow1+3:irow1+5, icol1+3:icol1+5] .-= M1_Ω1dot ./ force_scaling
    @views jacob[irow1+3:irow1+5, icol1+6:icol1+8] .-= M1_u1dot ./ force_scaling

    @views jacob[irow1+3:irow1+5, icol2:icol2+2] .-= M1_u2 ./ force_scaling
    @views jacob[irow1+3:irow1+5, icol2:icol2+2] .-= M1_V2dot ./ force_scaling
    @views jacob[irow1+3:irow1+5, icol2+3:icol2+5] .-= M1_Ω2dot ./ force_scaling
    @views jacob[irow1+3:irow1+5, icol2+6:icol2+8] .-= M1_u2dot ./ force_scaling

    # equilibrium equations for the end of the beam element
    irow2 = indices.irow_point[assembly.stop[ielem]]

    @views jacob[irow2:irow2+2, 4:6] .+= F2_θb ./ force_scaling
    @views jacob[irow2:irow2+2, 7:9] .+= F2_vb ./ force_scaling
    @views jacob[irow2:irow2+2, 10:12] .+= F2_ωb ./ force_scaling
    @views jacob[irow2:irow2+2, 13:15] .+= F2_ab ./ force_scaling
    @views jacob[irow2:irow2+2, 16:18] .+= F2_αb ./ force_scaling

    @views jacob[irow2:irow2+2, icol1:icol1+2] .+= F2_u1 ./ force_scaling
    @views jacob[irow2:irow2+2, icol1:icol1+2] .+= F2_V1dot ./ force_scaling
    @views jacob[irow2:irow2+2, icol1+3:icol1+5] .+= F2_Ω1dot ./ force_scaling

    @views jacob[irow2:irow2+2, icol2:icol2+2] .+= F2_u2 ./ force_scaling
    @views jacob[irow2:irow2+2, icol2:icol2+2] .+= F2_V2dot ./ force_scaling
    @views jacob[irow2:irow2+2, icol2+3:icol2+5] .+= F2_Ω2dot ./ force_scaling

    @views jacob[irow2+3:irow2+5, 4:6] .+= M2_θb ./ force_scaling
    @views jacob[irow2+3:irow2+5, 7:9] .+= M2_vb ./ force_scaling
    @views jacob[irow2+3:irow2+5, 10:12] .+= M2_ωb ./ force_scaling
    @views jacob[irow2+3:irow2+5, 13:15] .+= M2_ab ./ force_scaling
    @views jacob[irow2+3:irow2+5, 16:18] .+= M2_αb ./ force_scaling

    @views jacob[irow2+3:irow2+5, icol1:icol1+2] .+= M2_u1 ./ force_scaling
    @views jacob[irow2+3:irow2+5, icol1:icol1+2] .+= M2_V1dot ./ force_scaling
    @views jacob[irow2+3:irow2+5, icol1+3:icol1+5] .+= M2_Ω1dot ./ force_scaling
    @views jacob[irow2+3:irow2+5, icol1+6:icol1+8] .+= M2_u1dot ./ force_scaling

    @views jacob[irow2+3:irow2+5, icol2:icol2+2] .+= M2_u2 ./ force_scaling
    @views jacob[irow2+3:irow2+5, icol2:icol2+2] .+= M2_V2dot ./ force_scaling
    @views jacob[irow2+3:irow2+5, icol2+3:icol2+5] .+= M2_Ω2dot ./ force_scaling
    @views jacob[irow2+3:irow2+5, icol2+6:icol2+8] .+= M2_u2dot ./ force_scaling

    return jacob
end

@inline function insert_dynamic_element_jacobians!(jacob, indices, force_scaling,
    assembly, ielem, compatability, resultants)

    insert_static_element_jacobians!(jacob, indices, force_scaling, assembly, ielem, 
        compatability, resultants)

    @unpack ru_vb, ru_ωb, ru_V1, ru_V2, ru_Ω1, ru_Ω2, rθ_ωb, rθ_Ω1, rθ_Ω2 = compatability

    @unpack F1_θb, F1_ωb, F1_ab, F1_αb, F1_u1, F1_u2, F1_V1, F1_V2, F1_Ω1, F1_Ω2,
            F2_θb, F2_ωb, F2_ab, F2_αb, F2_u1, F2_u2, F2_V1, F2_V2, F2_Ω1, F2_Ω2,
            M1_θb, M1_vb, M1_ωb, M1_ab, M1_αb, M1_u1, M1_u2, M1_V1, M1_V2, M1_Ω1, M1_Ω2,
            M2_θb, M2_vb, M2_ωb, M2_ab, M2_αb, M2_u1, M2_u2, M2_V1, M2_V2, M2_Ω1, M2_Ω2 = resultants

    icol1 = indices.icol_point[assembly.start[ielem]]
    icol2 = indices.icol_point[assembly.stop[ielem]]

    # compatability equations
    irow = indices.irow_elem[ielem]

    jacob[irow:irow+2, 7:9] .= ru_vb 
    jacob[irow:irow+2, 10:12] .= ru_ωb 

    jacob[irow:irow+2, icol1+6:icol1+8] .= ru_V1
    jacob[irow:irow+2, icol1+9:icol1+11] .= ru_Ω1

    jacob[irow:irow+2, icol2+6:icol2+8] .= ru_V2
    jacob[irow:irow+2, icol2+9:icol2+11] .= ru_Ω2

    jacob[irow+3:irow+5, 10:12] .= rθ_ωb 

    jacob[irow+3:irow+5, icol1+9:icol1+11] .= rθ_Ω1

    jacob[irow+3:irow+5, icol2+9:icol2+11] .= rθ_Ω2

    # equilibrium equations for the start of the beam element
    irow1 = indices.irow_point[assembly.start[ielem]]

    @views jacob[irow1:irow1+2, 4:6] .-= F1_θb ./ force_scaling
    @views jacob[irow1:irow1+2, 10:12] .-= F1_ωb ./ force_scaling
    @views jacob[irow1:irow1+2, 13:15] .-= F1_ab ./ force_scaling
    @views jacob[irow1:irow1+2, 16:18] .-= F1_αb ./ force_scaling

    @views jacob[irow1:irow1+2, icol1:icol1+2] .-= F1_u1 ./ force_scaling

    @views jacob[irow1:irow1+2, icol2:icol2+2] .-= F1_u2 ./ force_scaling

    @views jacob[irow1:irow1+2, icol1+6:icol1+8] .-= F1_V1 ./ force_scaling
    @views jacob[irow1:irow1+2, icol1+9:icol1+11] .-= F1_Ω1 ./ force_scaling

    @views jacob[irow1:irow1+2, icol2+6:icol2+8] .-= F1_V2 ./ force_scaling
    @views jacob[irow1:irow1+2, icol2+9:icol2+11] .-= F1_Ω2 ./ force_scaling

    @views jacob[irow1+3:irow1+5, 4:6] .-= M1_θb ./ force_scaling
    @views jacob[irow1+3:irow1+5, 7:9] .-= M1_vb ./ force_scaling
    @views jacob[irow1+3:irow1+5, 10:12] .-= M1_ωb ./ force_scaling
    @views jacob[irow1+3:irow1+5, 13:15] .-= M1_ab ./ force_scaling
    @views jacob[irow1+3:irow1+5, 16:18] .-= M1_αb ./ force_scaling

    @views jacob[irow1+3:irow1+5, icol1:icol1+2] .-= M1_u1 ./ force_scaling
    
    @views jacob[irow1+3:irow1+5, icol2:icol2+2] .-= M1_u2 ./ force_scaling

    @views jacob[irow1+3:irow1+5, icol1+6:icol1+8] .-= M1_V1 ./ force_scaling
    @views jacob[irow1+3:irow1+5, icol1+9:icol1+11] .-= M1_Ω1 ./ force_scaling

    @views jacob[irow1+3:irow1+5, icol2+6:icol2+8] .-= M1_V2 ./ force_scaling
    @views jacob[irow1+3:irow1+5, icol2+9:icol2+11] .-= M1_Ω2 ./ force_scaling

    # equilibrium equations for the end of the beam element
    irow2 = indices.irow_point[assembly.stop[ielem]]

    @views jacob[irow2:irow2+2, 4:6] .+= F2_θb ./ force_scaling
    @views jacob[irow2:irow2+2, 10:12] .+= F2_ωb ./ force_scaling
    @views jacob[irow2:irow2+2, 13:15] .+= F2_ab ./ force_scaling
    @views jacob[irow2:irow2+2, 16:18] .+= F2_αb ./ force_scaling

    @views jacob[irow2:irow2+2, icol1:icol1+2] .+= F2_u1 ./ force_scaling

    @views jacob[irow2:irow2+2, icol2:icol2+2] .+= F2_u2 ./ force_scaling

    @views jacob[irow2:irow2+2, icol1+6:icol1+8] .+= F2_V1 ./ force_scaling
    @views jacob[irow2:irow2+2, icol1+9:icol1+11] .+= F2_Ω1 ./ force_scaling

    @views jacob[irow2:irow2+2, icol2+6:icol2+8] .+= F2_V2 ./ force_scaling
    @views jacob[irow2:irow2+2, icol2+9:icol2+11] .+= F2_Ω2 ./ force_scaling

    @views jacob[irow2+3:irow2+5, 4:6] .+= M2_θb ./ force_scaling
    @views jacob[irow2+3:irow2+5, 7:9] .+= M2_vb ./ force_scaling
    @views jacob[irow2+3:irow2+5, 10:12] .+= M2_ωb ./ force_scaling
    @views jacob[irow2+3:irow2+5, 13:15] .+= M2_ab ./ force_scaling
    @views jacob[irow2+3:irow2+5, 16:18] .+= M2_αb ./ force_scaling

    @views jacob[irow2+3:irow2+5, icol1:icol1+2] .+= M2_u1 ./ force_scaling

    @views jacob[irow2+3:irow2+5, icol2:icol2+2] .+= M2_u2 ./ force_scaling

    @views jacob[irow2+3:irow2+5, icol1+6:icol1+8] .+= M2_V1 ./ force_scaling
    @views jacob[irow2+3:irow2+5, icol1+9:icol1+11] .+= M2_Ω1 ./ force_scaling

    @views jacob[irow2+3:irow2+5, icol2+6:icol2+8] .+= M2_V2 ./ force_scaling
    @views jacob[irow2+3:irow2+5, icol2+9:icol2+11] .+= M2_Ω2 ./ force_scaling

    return jacob
end

@inline function insert_expanded_element_jacobians!(jacob, indices, force_scaling,
    assembly, ielem, compatability, velocities, equilibrium, resultants)

    @unpack ru_vb, ru_ωb, ru_u1, ru_u2, ru_θ1, ru_θ2, ru_V1, ru_V2, ru_Ω, ru_F1, ru_F2, ru_M1, ru_M2, 
            rθ_ωb, rθ_θ1, rθ_θ2, rθ_Ω1, rθ_Ω2, rθ_F1, rθ_F2, rθ_M1, rθ_M2 = compatability
    
    @unpack rV_vb, rV_ωb, rV_u1, rV_u2, rV_θ1, rV_θ2, rV_V, 
                   rΩ_ωb,               rΩ_θ1, rΩ_θ2, rΩ_Ω = velocities

    @unpack rF_θb, rF_ab, rF_αb, rF_F1, rF_F2, rM_F1, rM_F2, rM_M1, rM_M2, 
            rF_u1, rF_u2, rF_θ1, rF_θ2, rF_V, rF_Ω,  
            rM_θb, rM_vb, rM_ωb, rM_ab, rM_αb, rM_u1, rM_u2, rM_θ1, rM_θ2, 
            rM_V1, rM_V2, rM_V, rM_Ω = equilibrium
    
    @unpack F1_θ1, F1_θ2, F1_F1, F2_θ1, F2_θ2, F2_F2, 
            M1_θ1, M1_θ2, M1_M1, M2_θ1, M2_θ2, M2_M2 = resultants

    icol = indices.icol_elem[ielem]
    icol1 = indices.icol_point[assembly.start[ielem]]
    icol2 = indices.icol_point[assembly.stop[ielem]]

    irow = indices.irow_elem[ielem]

    # element compatability residuals
    jacob[irow:irow+2, 7:9] .= ru_vb 
    jacob[irow:irow+2, 10:12] .= ru_ωb 
    jacob[irow:irow+2, icol1:icol1+2] .= ru_u1
    jacob[irow:irow+2, icol2:icol2+2] .= ru_u2
    jacob[irow:irow+2, icol1+3:icol1+5] .= ru_θ1
    jacob[irow:irow+2, icol2+3:icol2+5] .= ru_θ2
    jacob[irow:irow+2, icol1+6:icol1+8] .= ru_V1
    jacob[irow:irow+2, icol2+6:icol2+8] .= ru_V2
    jacob[irow:irow+2, icol:icol+2] .= ru_F1 .* force_scaling
    jacob[irow:irow+2, icol+3:icol+5] .= ru_M1 .* force_scaling
    jacob[irow:irow+2, icol+6:icol+8] .= ru_F2 .* force_scaling
    jacob[irow:irow+2, icol+9:icol+11] .= ru_M2 .* force_scaling
    jacob[irow:irow+2, icol+15:icol+17] .= ru_Ω

    jacob[irow+3:irow+5, 10:12] .= rθ_ωb 
    jacob[irow+3:irow+5, icol1+3:icol1+5] .= rθ_θ1
    jacob[irow+3:irow+5, icol2+3:icol2+5] .= rθ_θ2
    jacob[irow+3:irow+5, icol1+9:icol1+11] .= rθ_Ω1
    jacob[irow+3:irow+5, icol2+9:icol2+11] .= rθ_Ω2
    jacob[irow+3:irow+5, icol:icol+2] .= rθ_F1 .* force_scaling
    jacob[irow+3:irow+5, icol+3:icol+5] .= rθ_M1 .* force_scaling
    jacob[irow+3:irow+5, icol+6:icol+8] .= rθ_F2 .* force_scaling
    jacob[irow+3:irow+5, icol+9:icol+11] .= rθ_M2 .* force_scaling

    # element equilibrium residuals
    jacob[irow+6:irow+8, 4:6] .= rF_θb ./ force_scaling 
    jacob[irow+6:irow+8, 13:15] .= rF_ab ./ force_scaling 
    jacob[irow+6:irow+8, 16:18] .= rF_αb ./ force_scaling 
    jacob[irow+6:irow+8, icol1:icol1+2] .= rF_u1 ./ force_scaling
    jacob[irow+6:irow+8, icol2:icol2+2] .= rF_u2 ./ force_scaling
    jacob[irow+6:irow+8, icol1+3:icol1+5] .= rF_θ1 ./ force_scaling
    jacob[irow+6:irow+8, icol2+3:icol2+5] .= rF_θ2 ./ force_scaling
    jacob[irow+6:irow+8, icol:icol+2] .= rF_F1
    jacob[irow+6:irow+8, icol+6:icol+8] .= rF_F2
    jacob[irow+6:irow+8, icol+12:icol+14] .= rF_V ./ force_scaling
    jacob[irow+6:irow+8, icol+15:icol+17] .= rF_Ω ./ force_scaling

    jacob[irow+9:irow+11, 4:6] .= rM_θb ./ force_scaling 
    jacob[irow+9:irow+11, 7:9] .= rM_vb ./ force_scaling 
    jacob[irow+9:irow+11, 10:12] .= rM_ωb ./ force_scaling 
    jacob[irow+9:irow+11, 13:15] .= rM_ab ./ force_scaling 
    jacob[irow+9:irow+11, 16:18] .= rM_αb ./ force_scaling 
    jacob[irow+9:irow+11, icol1:icol1+2] .= rM_u1 ./ force_scaling
    jacob[irow+9:irow+11, icol2:icol2+2] .= rM_u2 ./ force_scaling
    jacob[irow+9:irow+11, icol1+3:icol1+5] .= rM_θ1 ./ force_scaling
    jacob[irow+9:irow+11, icol2+3:icol2+5] .= rM_θ2 ./ force_scaling
    jacob[irow+9:irow+11, icol1+6:icol1+8] .= rM_V1 ./ force_scaling
    jacob[irow+9:irow+11, icol2+6:icol2+8] .= rM_V2 ./ force_scaling
    jacob[irow+9:irow+11, icol:icol+2] .= rM_F1
    jacob[irow+9:irow+11, icol+3:icol+5] .= rM_M1
    jacob[irow+9:irow+11, icol+6:icol+8] .= rM_F2
    jacob[irow+9:irow+11, icol+9:icol+11] .= rM_M2
    jacob[irow+9:irow+11, icol+12:icol+14] .= rM_V ./ force_scaling
    jacob[irow+9:irow+11, icol+15:icol+17] .= rM_Ω ./ force_scaling

    # velocity residuals

    jacob[irow+12:irow+14, 7:9] .= rV_vb
    jacob[irow+12:irow+14, 10:12] .= rV_ωb
    jacob[irow+12:irow+14, icol1:icol1+2] .= rV_u1
    jacob[irow+12:irow+14, icol2:icol2+2] .= rV_u2
    jacob[irow+12:irow+14, icol1+3:icol1+5] .= rV_θ1
    jacob[irow+12:irow+14, icol2+3:icol2+5] .= rV_θ2
    jacob[irow+12:irow+14, icol+12:icol+14] .= rV_V

    jacob[irow+15:irow+17, 10:12] .= rΩ_ωb
    jacob[irow+15:irow+17, icol1+3:icol1+5] .= rΩ_θ1
    jacob[irow+15:irow+17, icol2+3:icol2+5] .= rΩ_θ2
    jacob[irow+15:irow+17, icol+15:icol+17] .= rΩ_Ω

    # equilibrium equations for the start of the beam element
    irow = indices.irow_point[assembly.start[ielem]]
    @views jacob[irow:irow+2, icol1+3:icol1+5] .-= F1_θ1 ./ force_scaling
    @views jacob[irow:irow+2, icol2+3:icol2+5] .-= F1_θ2 ./ force_scaling
    @views jacob[irow:irow+2, icol:icol+2] .-= F1_F1
    @views jacob[irow+3:irow+5, icol1+3:icol1+5] .-= M1_θ1 ./ force_scaling
    @views jacob[irow+3:irow+5, icol2+3:icol2+5] .-= M1_θ2 ./ force_scaling
    @views jacob[irow+3:irow+5, icol+3:icol+5] .-= M1_M1

    # equilibrium equations for the end of the beam element
    irow = indices.irow_point[assembly.stop[ielem]]
    @views jacob[irow:irow+2, icol1+3:icol1+5] .+= F2_θ1 ./ force_scaling
    @views jacob[irow:irow+2, icol2+3:icol2+5] .+= F2_θ2 ./ force_scaling
    @views jacob[irow:irow+2, icol+6:icol+8] .+= F2_F2
    @views jacob[irow+3:irow+5, icol1+3:icol1+5] .+= M2_θ1 ./ force_scaling
    @views jacob[irow+3:irow+5, icol2+3:icol2+5] .+= M2_θ2 ./ force_scaling
    @views jacob[irow+3:irow+5, icol+9:icol+11] .+= M2_M2

    return jacob
end

@inline function insert_mass_matrix_element_jacobians!(jacob, gamma, indices, force_scaling,
    assembly, ielem, resultants)

    @unpack F1_V1dot, F1_V2dot, F1_Ω1dot, F1_Ω2dot, 
            F2_V1dot, F2_V2dot, F2_Ω1dot, F2_Ω2dot, 
            M1_V1dot, M1_V2dot, M1_Ω1dot, M1_Ω2dot, 
            M2_V1dot, M2_V2dot, M2_Ω1dot, M2_Ω2dot = resultants

    icol1 = indices.icol_point[assembly.start[ielem]]
    icol2 = indices.icol_point[assembly.stop[ielem]]
    
    # equilibrium equations for the start of the beam element
    irow1 = indices.irow_point[assembly.start[ielem]]

    @views jacob[irow1:irow1+2, icol1+6:icol1+8] .-= F1_V1dot .* gamma ./ force_scaling
    @views jacob[irow1:irow1+2, icol1+9:icol1+11] .-= F1_Ω1dot .* gamma ./ force_scaling

    @views jacob[irow1:irow1+2, icol2+6:icol2+8] .-= F1_V2dot .* gamma ./ force_scaling
    @views jacob[irow1:irow1+2, icol2+9:icol2+11] .-= F1_Ω2dot .* gamma ./ force_scaling

    @views jacob[irow1+3:irow1+5, icol1+6:icol1+8] .-= M1_V1dot .* gamma ./ force_scaling
    @views jacob[irow1+3:irow1+5, icol1+9:icol1+11] .-= M1_Ω1dot .* gamma ./ force_scaling

    @views jacob[irow1+3:irow1+5, icol2+6:icol2+8] .-= M1_V2dot .* gamma ./ force_scaling
    @views jacob[irow1+3:irow1+5, icol2+9:icol2+11] .-= M1_Ω2dot .* gamma ./ force_scaling

    # equilibrium equations for the end of the beam element
    irow2 = indices.irow_point[assembly.stop[ielem]]

    @views jacob[irow2:irow2+2, icol1+6:icol1+8] .+= F2_V1dot .* gamma ./ force_scaling
    @views jacob[irow2:irow2+2, icol1+9:icol1+11] .+= F2_Ω1dot .* gamma ./ force_scaling

    @views jacob[irow2:irow2+2, icol2+6:icol2+8] .+= F2_V2dot .* gamma ./ force_scaling
    @views jacob[irow2:irow2+2, icol2+9:icol2+11] .+= F2_Ω2dot .* gamma ./ force_scaling

    @views jacob[irow2+3:irow2+5, icol1+6:icol1+8] .+= M2_V1dot .* gamma ./ force_scaling
    @views jacob[irow2+3:irow2+5, icol1+9:icol1+11] .+= M2_Ω1dot .* gamma ./ force_scaling

    @views jacob[irow2+3:irow2+5, icol2+6:icol2+8] .+= M2_V2dot .* gamma ./ force_scaling
    @views jacob[irow2+3:irow2+5, icol2+9:icol2+11] .+= M2_Ω2dot .* gamma ./ force_scaling

    return jacob
end

@inline function insert_expanded_mass_matrix_element_jacobians!(jacob, gamma, indices, 
    force_scaling, assembly, ielem, equilibrium, velocities)

    @unpack rF_Vdot, rF_Ωdot, rM_Vdot, rM_Ωdot = equilibrium

    rV_u1dot, rV_u2dot, rΩ_θ1dot, rΩ_θ2dot = velocities

    irow = indices.irow_elem[ielem] 
    icol = indices.icol_elem[ielem]
    icol1 = indices.icol_point[assembly.start[ielem]]
    icol2 = indices.icol_point[assembly.stop[ielem]]

    # equilibrium residuals
    @views jacob[irow+6:irow+8, icol+12:icol+14] .+= rF_Vdot .* gamma ./ force_scaling
    @views jacob[irow+6:irow+8, icol+15:icol+17] .+= rF_Ωdot .* gamma ./ force_scaling

    @views jacob[irow+9:irow+11, icol+12:icol+14] .+= rM_Vdot .* gamma ./ force_scaling
    @views jacob[irow+9:irow+11, icol+15:icol+17] .+= rM_Ωdot .* gamma ./ force_scaling

    # velocity residuals
    @views jacob[irow+12:irow+14, icol1:icol1+2] .+= rV_u1dot .* gamma
    @views jacob[irow+12:irow+14, icol2:icol2+2] .+= rV_u2dot .* gamma
    @views jacob[irow+15:irow+17, icol1+3:icol1+5] .+= rΩ_θ1dot .* gamma
    @views jacob[irow+15:irow+17, icol2+3:icol2+5] .+= rΩ_θ2dot .* gamma

    return jacob
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

    compatability = static_compatability_jacobians(properties)

    resultants = static_element_resultant_jacobians(properties, distributed_loads, ielem)

    insert_static_element_jacobians!(jacob, indices, force_scaling, assembly, ielem, 
        compatability, resultants)

    return jacob
end

"""
    steady_state_element_jacobian!(jacob, x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, 
        gravity)

Calculate and insert the jacobian entries corresponding to a beam element for a steady state 
analysis into the system jacobian matrix.
"""
@inline function steady_state_element_jacobian!(jacob, x, indices, 
    force_scaling, structural_damping, assembly, ielem, prescribed_conditions, 
    distributed_loads, gravity)

    properties = steady_state_element_properties(x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity)

    properties = steady_state_element_jacobian_properties(properties, x, indices, 
        force_scaling, structural_damping, assembly, ielem, 
        prescribed_conditions, gravity)

    compatability = dynamic_compatability_jacobians(properties)

    resultants = steady_state_element_resultant_jacobians(properties, distributed_loads, ielem)

    insert_dynamic_element_jacobians!(jacob, indices, force_scaling, assembly, ielem, 
        compatability, resultants)

    return jacob
end

"""
    initial_condition_element_jacobian!(jacob, x, indices, rate_vars, 
        force_scaling, structural_damping, assembly, ielem, prescribed_conditions, 
        distributed_loads, gravity, 
        u0, θ0, V0, Ω0, Vdot0, Ωdot0)

Calculate and insert the jacobian entries corresponding to a beam element for the 
initialization of a time domain analysis into the system jacobian matrix.
"""
@inline function initial_condition_element_jacobian!(jacob, x, indices, rate_vars, 
    force_scaling, structural_damping, assembly, ielem, prescribed_conditions, 
    distributed_loads, gravity, 
    u0, θ0, V0, Ω0, Vdot0, Ωdot0)

    properties = initial_condition_element_properties(x, indices, rate_vars, 
        force_scaling, structural_damping, assembly, ielem, prescribed_conditions, gravity, 
        u0, θ0, V0, Ω0, Vdot0, Ωdot0)

    properties = initial_condition_element_jacobian_properties(properties, x, indices, 
        rate_vars, force_scaling, structural_damping, assembly, ielem, 
        prescribed_conditions, gravity, 
        u0, θ0, V0, Ω0, Vdot0, Ωdot0)

    compatability = initial_condition_compatability_jacobians(properties)

    resultants = initial_condition_element_resultant_jacobians(properties, distributed_loads, ielem)

    insert_initial_condition_element_jacobians!(jacob, indices, force_scaling, assembly, ielem, 
        compatability, resultants)

    return jacob
end

"""
    newmark_element_jacobian!(jacob, x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, 
        gravity, Vdot_init, Ωdot_init, dt)

Calculate and insert the jacobian entries corresponding to a beam element for a Newmark-scheme
time marching analysis into the system jacobian matrix.
"""
@inline function newmark_element_jacobian!(jacob, x, indices, force_scaling, 
    structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, gravity, 
    Vdot_init, Ωdot_init, dt)

    properties = newmark_element_properties(x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity, 
        Vdot_init, Ωdot_init, dt)

    properties = newmark_element_jacobian_properties(properties, x, indices, 
        force_scaling, structural_damping, assembly, ielem, prescribed_conditions, gravity, 
        Vdot_init, Ωdot_init, dt)

    compatability = dynamic_compatability_jacobians(properties)

    resultants = newmark_element_resultant_jacobians(properties, distributed_loads, ielem)

    insert_dynamic_element_jacobians!(jacob, indices, force_scaling, assembly, ielem, 
        compatability, resultants)

    return jacob
end

"""
    dynamic_element_jacobian!(jacob, dx, x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, 
        gravity)

Calculate and insert the jacobian entries corresponding to a beam element for a dynamic
analysis into the system jacobian matrix.
"""
@inline function dynamic_element_jacobian!(jacob, dx, x, indices, force_scaling, 
    structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, gravity)

    properties = dynamic_element_properties(dx, x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity)

    properties = dynamic_element_jacobian_properties(properties, dx, x, indices, 
        force_scaling, structural_damping, assembly, ielem, prescribed_conditions, gravity)

    compatability = dynamic_compatability_jacobians(properties)

    resultants = dynamic_element_resultant_jacobians(properties, distributed_loads, ielem)

    insert_dynamic_element_jacobians!(jacob, indices, force_scaling, assembly, ielem, 
        compatability, resultants)

    return jacob
end

"""
    expanded_steady_element_jacobian!(jacob, x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, 
        gravity)

Calculate and insert the jacobian entries corresponding to a beam element for a dynamic
analysis into the system jacobian matrix.
"""
@inline function expanded_steady_element_jacobian!(jacob, x, indices, force_scaling, structural_damping, 
    assembly, ielem, prescribed_conditions, distributed_loads, gravity)

    properties = expanded_steady_element_properties(x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity)

    properties = expanded_steady_element_jacobian_properties(properties, x, indices, 
        force_scaling, structural_damping, assembly, ielem, prescribed_conditions, gravity)

    compatability = expanded_compatability_jacobians(properties)

    velocities = expanded_element_velocity_jacobians(properties)

    equilibrium = expanded_element_equilibrium_jacobians(properties, distributed_loads, ielem)

    resultants = expanded_element_resultant_jacobians(properties)

    insert_expanded_element_jacobians!(jacob, indices, force_scaling, assembly, ielem, 
        compatability, velocities, equilibrium, resultants)

    return jacob
end

"""
    expanded_dynamic_element_jacobian!(jacob, x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, 
        gravity)

Calculate and insert the jacobian entries corresponding to a beam element for a dynamic
analysis into the system jacobian matrix.
"""
@inline function expanded_dynamic_element_jacobian!(jacob, x, indices, force_scaling, structural_damping, 
    assembly, ielem, prescribed_conditions, distributed_loads, gravity)

    properties = expanded_dynamic_element_properties(x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity)

    properties = expanded_dynamic_element_jacobian_properties(properties, x, indices, 
        force_scaling, structural_damping, assembly, ielem, prescribed_conditions, gravity)

    compatability = expanded_compatability_jacobians(properties)

    velocities = expanded_element_velocity_jacobians(properties)

    equilibrium = expanded_element_equilibrium_jacobians(properties, distributed_loads, ielem)

    resultants = expanded_element_resultant_jacobians(properties)

    insert_expanded_element_jacobians!(jacob, indices, force_scaling, assembly, ielem, 
        compatability, velocities, equilibrium, resultants)

    return jacob
end

"""
    mass_matrix_element_jacobian!(jacob, gamma, x, indices, force_scaling, assembly, 
        ielem, prescribed_conditions)

Calculate and insert the mass_matrix jacobian entries corresponding to a beam element into 
the system jacobian matrix.
"""
@inline function mass_matrix_element_jacobian!(jacob, gamma, x, indices, force_scaling, assembly, 
    ielem, prescribed_conditions)

    properties = mass_matrix_element_jacobian_properties(x, indices, force_scaling, 
        assembly, ielem, prescribed_conditions)

    resultants = mass_matrix_element_resultant_jacobians(properties)
    
    insert_mass_matrix_element_jacobians!(jacob, gamma, indices, force_scaling, assembly, ielem, 
        resultants)

    return jacob
end

"""
    expanded_mass_matrix_element_jacobian!(jacob, gamma, indices, force_scaling, assembly, 
        ielem, prescribed_conditions)

Calculate and insert the mass_matrix jacobian entries corresponding to a beam element into 
the system jacobian matrix for a constant mass matrix system
"""
@inline function expanded_mass_matrix_element_jacobian!(jacob, gamma, indices, force_scaling, assembly, 
    ielem, prescribed_conditions)

    properties = expanded_mass_matrix_element_jacobian_properties(assembly, ielem, prescribed_conditions)

    equilibrium = expanded_mass_matrix_element_equilibrium_jacobians(properties)
    
    velocities = expanded_mass_matrix_element_velocity_jacobians(properties)

    insert_expanded_mass_matrix_element_jacobians!(jacob, gamma, indices, force_scaling, 
        assembly, ielem, equilibrium, velocities)

    return jacob
end
