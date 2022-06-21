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
    steady_state_element_properties(x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, gravity, x0, v0, ω0, a0, α0)

Calculate/extract the element properties needed to construct the residual for a steady 
state analysis
"""
@inline function steady_state_element_properties(x, indices, force_scaling, structural_damping, 
    assembly, ielem, prescribed_conditions, gravity, x0, v0, ω0, a0, α0)

    properties = static_element_properties(x, indices, force_scaling, assembly, ielem, 
        prescribed_conditions, gravity)

    @unpack L, Cab, C, CtCab, Qinv, mass11, mass12, mass21, mass22, u1, u2, θ1, θ2, 
        u, θ, γ, κ, a, α = properties

    # rotation parameter matrices
    Q = get_Q(θ)

    # linear and angular velocity
    V1, Ω1 = point_velocities(x, assembly.start[ielem], indices.icol_point)
    V2, Ω2 = point_velocities(x, assembly.stop[ielem], indices.icol_point)
    V = (V1 + V2)/2
    Ω = (Ω1 + Ω2)/2

    # linear and angular momentum
    P = CtCab*mass11*CtCab'*V + CtCab*mass12*CtCab'*Ω
    H = CtCab*mass21*CtCab'*V + CtCab*mass22*CtCab'*Ω

    # distance from the rotation center
    Δx = assembly.elements[ielem].x - x0

    # rigid body linear and angular velocity
    v = v0 + cross(ω0, Δx)
    ω = ω0

    # linear and angular acceleration
    a += a0 + cross(α0, Δx) + cross(α0, u)
    α += α0

    # add steady state properties
    properties = (; properties..., V1, V2, Ω1, Ω2, Q, V, Ω, P, H, γ, κ, v, ω, a, α)

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
        Δx1 = assembly.points[assembly.start[ielem]] - x0
        Δx2 = assembly.points[assembly.stop[ielem]] - x0

        # rigid body linear and angular velocity
        v1 = v0 + cross(ω0, Δx1 - x0)
        v2 = v0 + cross(ω0, Δx2 - x0)
        ω1 = ω0
        ω2 = ω0

        # linear deflection rates 
        udot1 = V1 - v1 - cross(ω1, u1)
        udot2 = V2 - v2 - cross(ω2, u2)
        udot = (udot1 + udot2)/2

        # angular deflection rates
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

        # add structural damping properties
        properties = (; properties..., γ, κ, C1, C2, Qinv1, Qinv2, v1, v2, ω1, ω2, 
            udot1, udot2, θdot1, θdot2, ΔQ, μ11, μ22, Δu, Δθ, Δudot, Δθdot, 
            udot, θdot, γdot, κdot)
    end

    return properties
end

"""
    initial_condition_element_properties(x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, gravity, x0, v0, ω0, a0, α0, 
        u0, θ0, udot0, θdot0)

Calculate/extract the element properties needed to construct the residual for a time-domain
analysis initialization
"""
@inline function initial_condition_element_properties(x, indices, force_scaling, 
    structural_damping, assembly, ielem, prescribed_conditions, gravity, x0, v0, ω0, a0, α0, 
    u0, θ0, udot0, θdot0)

    # unpack element parameters
    @unpack L, Cab, compliance, mass, mu = assembly.elements[ielem]

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

    # damping submatrices
    μ11 = @SMatrix [mu[1] 0 0; 0 mu[2] 0; 0 0 mu[3]]
    μ22 = @SMatrix [mu[4] 0 0; 0 mu[5] 0; 0 0 mu[6]]

    # linear and angular displacement
    u1 = SVector{3}(u0[assembly.start[ielem]])
    θ1 = SVector{3}(θ0[assembly.start[ielem]])

    u2 = SVector{3}(u0[assembly.stop[ielem]])
    θ2 = SVector{3}(θ0[assembly.stop[ielem]])

    u = (u1 + u2)/2
    θ = (θ1 + θ2)/2

    # rotation parameter matrices
    C = get_C(θ)
    CtCab = C'*Cab
    Q = get_Q(θ)
    Qinv = get_Qinv(θ)

    # linear and angular velocity
    V1, Ω1 = point_velocities(x, assembly.start[ielem], indices.icol_point)
    V2, Ω2 = point_velocities(x, assembly.stop[ielem], indices.icol_point)
    V = (V1 + V2)/2
    Ω = (Ω1 + Ω2)/2

    # linear and angular momentum
    P = CtCab*mass11*CtCab'*V + CtCab*mass12*CtCab'*Ω
    H = CtCab*mass21*CtCab'*V + CtCab*mass22*CtCab'*Ω

    # forces and moments
    F, M = element_loads(x, ielem, indices.icol_elem, force_scaling)

    # strain and curvature
    γ = S11*F + S12*M
    κ = S21*F + S22*M

    # distance from the rotation center
    Δx = assembly.elements[ielem].x - x0

    # rigid body linear velocity
    v = v0 + cross(ω0, Δx)

    # rigid body angular velocity
    ω = ω0

    # linear and angular acceleration
    a = a0 + cross(α0, Δx) + cross(α0, u) - gravity
    α = α0

    # linear deflection rates
    udot1 = udot0[assembly.start[ielem]]
    udot2 = udot0[assembly.stop[ielem]]
    udot = (udot1 + udot2)/2

    # angular deflection rates
    θdot1 = θdot0[assembly.start[ielem]]
    θdot2 = θdot0[assembly.stop[ielem]]
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
    γdot = CtCab'*Δudot - Cab'*tilde(Ω - ω)*C*Δu - L*Cab'*tilde(Ω - ω)*C*Cab*e1
    κdot = Cab'*Q*Δθdot + Cab'*ΔQ*θdot

    # adjust strains to account for strain rates
    if structural_damping
        γ -= μ11*γdot
        κ -= μ22*κdot
    end

    # linear and angular velocity rates
    V1dot, Ω1dot = point_displacement_rates(x, assembly.start[ielem], indices.icol_point, prescribed_conditions)
    V2dot, Ω2dot = point_displacement_rates(x, assembly.stop[ielem], indices.icol_point, prescribed_conditions)
    Vdot = (V1dot + V2dot)/2
    Ωdot = (Ω1dot + Ω2dot)/2

    # linear and angular momentum rates
    CtCabdot = tilde(Ω - ω)*CtCab

    Pdot = CtCab*mass11*CtCab'*Vdot + CtCab*mass12*CtCab'*Ωdot +
        CtCab*mass11*CtCabdot'*V + CtCab*mass12*CtCabdot'*Ω +
        CtCabdot*mass11*CtCab'*V + CtCabdot*mass12*CtCab'*Ω
    
    Hdot = CtCab*mass21*CtCab'*Vdot + CtCab*mass22*CtCab'*Ωdot +
        CtCab*mass21*CtCabdot'*V + CtCab*mass22*CtCabdot'*Ω +
        CtCabdot*mass21*CtCab'*V + CtCabdot*mass22*CtCab'*Ω

    return (; L, C, Cab, CtCab, Q, Qinv, S11, S12, S21, S22, mass11, mass12, mass21, mass22, 
        μ11, μ22, u1, u2, θ1, θ2, V1, V2, Ω1, Ω2, u, θ, V, Ω, P, H, F, M, γ, κ, v, ω, a, α, Δu,
        CtCabdot, Vdot, Ωdot, Pdot, Hdot) 
end

"""
    newmark_element_properties(x, indices, force_scaling, structural_damping, assembly, ielem,
        prescribed_conditions, gravity, x0, v0, ω0, a0, α0, Vdot_init, Ωdot_init, dt)

Calculate/extract the element properties needed to construct the residual for a newmark-
scheme time stepping analysis
"""
@inline function newmark_element_properties(x, indices, force_scaling, structural_damping, 
    assembly, ielem, prescribed_conditions, gravity, x0, v0, ω0, a0, α0, Vdot_init, 
    Ωdot_init, dt)

    properties = steady_state_element_properties(x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, gravity, x0, v0, ω0, a0, α0)

    @unpack C, Cab, CtCab, mass11, mass12, mass21, mass22, V1, V2, Ω1, Ω2, V, Ω, ω = properties

    # linear and angular velocity rates
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

    return (; properties..., CtCabdot, Vdot, Ωdot, Pdot, Hdot) 
end

"""
    dynamic_element_properties(dx, x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, gravity, x0, v0, ω0, a0, α0)

Calculate/extract the element properties needed to construct the residual for a dynamic 
analysis
"""
@inline function dynamic_element_properties(dx, x, indices, force_scaling, structural_damping, 
    assembly, ielem, prescribed_conditions, gravity, x0, v0, ω0, a0, α0)

    properties = steady_state_element_properties(x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, gravity, x0, v0, ω0, a0, α0)

    @unpack C, Cab, CtCab, mass11, mass12, mass21, mass22, V, Ω, ω = properties

    # linear and angular velocity of the beam element
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

    return (; properties..., CtCabdot, Vdot, Ωdot, Pdot, Hdot)
end

"""
    expanded_element_properties(x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, gravity, x0, v0, ω0, a0, α0)

Calculate/extract the element properties needed to construct the residual for a constant
mass matrix system
"""
@inline function expanded_element_properties(x, indices, force_scaling, structural_damping, 
    assembly, ielem, prescribed_conditions, gravity, x0, v0, ω0, a0, α0)

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

    # linear and angular velocity
    V1, Ω1 = point_velocities(x, assembly.start[ielem], indices.icol_point)
    V2, Ω2 = point_velocities(x, assembly.stop[ielem], indices.icol_point)
    V, Ω = expanded_element_velocities(x, ielem, indices.icol_elem)

    # linear and angular momentum
    P = mass11*V + mass12*Ω
    H = mass21*V + mass22*Ω

    # forces and moments
    F1, M1, F2, M2 = expanded_element_loads(x, ielem, indices.icol_elem, force_scaling)
    F = (F1 + F2)/2
    M = (M1 + M2)/2

    # strain and curvature
    γ = S11*F + S12*M
    κ = S21*F + S22*M

    # distance from the rotation center
    Δx = assembly.elements[ielem].x - x0

    # rigid body linear and angular velocity
    v = v0 + cross(ω0, Δx)
    ω = ω0

    # linear and angular acceleration
    a = a0 + cross(α0, Δx) + cross(α0, u) - gravity
    α = α0

    # element properties
    properties = (; L, C, C1, C2, Cab, CtCab, Q, Qinv, S11, S12, S21, S22, 
        mass11, mass12, mass21, mass22, u1, u2, θ1, θ2, V1, V2, Ω1, Ω2, F1, F2, M1, M2, 
        u, θ, V, Ω, P, H, F, M, γ, κ, v, ω, a, α)

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
        Δx1 = assembly.points[assembly.start[ielem]] - x0
        Δx2 = assembly.points[assembly.stop[ielem]] - x0

        # rigid body linear and angular velocity
        v1 = v0 + cross(ω0, Δx1 - x0)
        v2 = v0 + cross(ω0, Δx2 - x0)
        ω1 = ω0
        ω2 = ω0

        # linear deflection rates 
        udot1 = C1'*V1 - v1 - cross(ω1, u1)
        udot2 = C2'*V2 - v2 - cross(ω2, u2)
        udot = (udot1 + udot2)/2

        # angular deflection rates
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
        properties = (; properties..., γ, κ, C1, C2, Qinv1, Qinv2, v1, v2, ω1, ω2, 
            udot1, udot2, θdot1, θdot2, ΔQ, μ11, μ22, Δu, Δθ, Δudot, Δθdot, 
            udot, θdot, γdot, κdot)
    end

    return properties
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
    
    rV = CtCab*V - v - cross(ω, u)
    rΩ = Qinv*(Cab*Ω - C*ω)

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
        rF += CtCab'*dload.f1 + Cab'*dload.f1_follower
        rF += CtCab'*dload.f2 + Cab'*dload.f2_follower
        rM += CtCab'*dload.m1 + Cab'*dload.m1_follower
        rM += CtCab'*dload.m2 + Cab'*dload.m2_follower
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

    # add loads due to linear and angular acceleration (including gravity)
    tmp = 1/2*(CtCab*mass11*CtCab'*a + CtCab*mass12*CtCab'*α)
    F1 -= tmp
    F2 += tmp

    tmp = 1/2*(CtCab*mass21*CtCab'*a + CtCab*mass22*CtCab'*α)
    M1 -= tmp
    M2 += tmp

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

Calculate the resultant loads applied at each end of a beam element for a dynamic 
analysis.
"""
@inline function steady_state_element_resultants(properties, distributed_loads, ielem)

    resultants = static_element_resultants(properties, distributed_loads, ielem)

    @unpack F1, F2, M1, M2 = resultants

    @unpack V, Ω, P, H, v, ω = properties

    # add loads due to linear and angular momentum
    tmp = 1/2*cross(ω, P)
    F1 -= tmp
    F2 += tmp

    tmp = 1/2*(cross(ω, H) + cross(V, P))
    M1 -= tmp
    M2 += tmp

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
    tmp = 1/2*Pdot 
    F1 -= tmp
    F2 += tmp

    tmp = 1/2*Hdot
    M1 -= tmp
    M2 += tmp

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
    steady_state_element_residual!(resid, x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, distributed_loads, gravity, 
        x0, v0, ω0, a0, α0)

Calculate and insert the residual entries corresponding to a beam element for a steady state 
analysis into the system residual vector.
"""
@inline function steady_state_element_residual!(resid, x, indices, force_scaling, structural_damping, 
    assembly, ielem, prescribed_conditions, distributed_loads, gravity, x0, v0, ω0, a0, α0)

    properties = steady_state_element_properties(x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, gravity, x0, v0, ω0, a0, α0)

    compatability = compatability_residuals(properties)

    resultants = steady_state_element_resultants(properties, distributed_loads, ielem)

    insert_element_residuals!(resid, indices, force_scaling, assembly, ielem, 
        compatability, resultants)

    return resid
end

"""
    initial_condition_element_residual!(resid, x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, 
        gravity, x0, v0, ω0, a0, α0, u0, θ0, udot0, θdot0)

Calculate and insert the residual entries corresponding to a beam element for the 
initialization of a time domain simulation into the system residual vector.
"""
@inline function initial_condition_element_residual!(resid, x, indices, force_scaling, 
    structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, gravity, 
    x0, v0, ω0, a0, α0, u0, θ0, udot0, θdot0)

    properties = initial_condition_element_properties(x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity, 
        x0, v0, ω0, a0, α0, u0, θ0, udot0, θdot0)

    compatability = compatability_residuals(properties)

    resultants = dynamic_element_resultants(properties, distributed_loads, ielem)

    insert_element_residuals!(resid, indices, force_scaling, assembly, ielem, 
        compatability, resultants)

    return resid
end

"""
    newmark_element_residual!(resid, x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, distributed_loads, gravity, 
        x0, v0, ω0, a0, α0, Vdot_init, Ωdot_init, dt)

Calculate and insert the residual entries corresponding to a beam element for a 
newmark-scheme time marching analysis into the system residual vector.
"""
@inline function newmark_element_residual!(resid, x, indices, force_scaling, 
    structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, gravity, 
    x0, v0, ω0, a0, α0, Vdot_init, Ωdot_init, dt)

    properties = newmark_element_properties(x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, gravity, x0, v0, ω0, a0, α0, 
        Vdot_init, Ωdot_init, dt)

    compatability = compatability_residuals(properties)

    resultants = dynamic_element_resultants(properties, distributed_loads, ielem)

    insert_element_residuals!(resid, indices, force_scaling, assembly, ielem, 
        compatability, resultants)

    return resid
end

"""
    dynamic_element_residual!(resid, dx, x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, distributed_loads, gravity, 
        x0, v0, ω0, a0, α0)

Calculate and insert the residual entries corresponding to a beam element for a dynamic
analysis into the system residual vector.
"""
@inline function dynamic_element_residual!(resid, dx, x, indices, force_scaling, structural_damping, 
    assembly, ielem, prescribed_conditions, distributed_loads, gravity, x0, v0, ω0, a0, α0)

    properties = dynamic_element_properties(dx, x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, gravity, x0, v0, ω0, a0, α0)

    compatability = compatability_residuals(properties)

    resultants = dynamic_element_resultants(properties, distributed_loads, ielem)

    insert_element_residuals!(resid, indices, force_scaling, assembly, ielem, 
        compatability, resultants)

    return resid
end

"""
    expanded_element_residual!(resid, x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, distributed_loads, gravity, 
        x0, v0, ω0, a0, α0)

Calculate and insert the residual entries corresponding to a beam element for a constant
mass matrix system into the system residual vector.
"""
@inline function expanded_element_residual!(resid, x, indices, force_scaling, structural_damping, 
    assembly, ielem, prescribed_conditions, distributed_loads, gravity, x0, v0, ω0, a0, α0)

    properties = expanded_element_properties(x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, gravity, x0, v0, ω0, a0, α0)

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
    steady_state_element_jacobian_properties(properties, x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, gravity, x0, v0, ω0, a0, α0))

Calculate/extract the element properties needed to calculate the jacobian entries 
corresponding to a element for a steady state analysis
"""
@inline function steady_state_element_jacobian_properties(properties, x, indices, force_scaling, 
    structural_damping, assembly, ielem, prescribed_conditions, gravity, x0, v0, ω0, a0, α0)

    properties = static_element_jacobian_properties(properties, x, indices, force_scaling, 
        assembly, ielem, prescribed_conditions, gravity)

    @unpack L, mass11, mass12, mass21, mass22, C, Cab, CtCab, Q, θ, V, Ω, ω, C_θ1, C_θ2, C_θ3 = properties

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
    a_u = tilde(α0)    

    if structural_damping

        @unpack C1, C2, Qinv1, Qinv2, u1, u2, θ1, θ2, Ω1, Ω2, ω1, ω2, μ11, μ22, udot, θdot, 
            Δu, Δθ, ΔQ, Δudot, Δθdot = properties

        # rotation parameter matrices
        C1_θ1, C1_θ2, C1_θ3 = get_C_θ(C1, θ1)
        C2_θ1, C2_θ2, C2_θ3 = get_C_θ(C2, θ2)
        Qinv1_θ1, Qinv1_θ2, Qinv1_θ3 = get_Qinv_θ(θ1)
        Qinv2_θ1, Qinv2_θ2, Qinv2_θ3 = get_Qinv_θ(θ2)

        # linear deflection rates
        udot1_u1 = -tilde(ω1)
        udot2_u2 = -tilde(ω2)

        udot1_V1 = I3
        udot2_V2 = I3

        # angular deflection rates
        θdot1_θ1 = mul3(Qinv1_θ1, Qinv1_θ2, Qinv1_θ3, C1*(Ω1 - ω1)) + Qinv1*mul3(C1_θ1, C1_θ2, C1_θ3, Ω1 - ω1)
        θdot2_θ2 = mul3(Qinv2_θ1, Qinv2_θ2, Qinv2_θ3, C2*(Ω2 - ω2)) + Qinv2*mul3(C2_θ1, C2_θ2, C2_θ3, Ω2 - ω2)
        
        θdot1_Ω1 = Qinv1*C1
        θdot2_Ω2 = Qinv2*C2

        θdot_θ1 = 1/2*θdot1_θ1
        θdot_θ2 = 1/2*θdot2_θ2

        θdot_Ω1 = 1/2*θdot1_Ω1
        θdot_Ω2 = 1/2*θdot2_Ω2

        # change in linear displacement
        Δu_u1 = -I3
        Δu_u2 =  I3

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
        κ_θ1, κ_θ2, κ_Ω1, κ_Ω2, P_θ, P_V, P_Ω, H_θ, H_V, H_Ω, a_u)
end

"""
    initial_condition_element_jacobian_properties(properties, x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity)

Calculate/extract the element properties needed to calculate the jacobian entries 
corresponding to a element for a Newmark scheme time marching analysis
"""
@inline function initial_condition_element_jacobian_properties(properties, x, indices, 
    force_scaling, structural_damping, assembly, ielem, prescribed_conditions, gravity)

    @unpack L, C, Cab, CtCab, CtCabdot, S11, S12, S21, S22, mass11, mass12, mass21, mass22, 
        μ11, μ22, Δu, V, Ω = properties

    # strain and curvature
    γ_F, γ_M, κ_F, κ_M = S11, S12, S21, S22

    # strain rates
    tmp = Cab'*tilde(C*Δu) + L*Cab'*tilde(C*Cab*e1)
    γdot_Ω1 = 1/2*tmp
    γdot_Ω2 = 1/2*tmp

    if structural_damping       
        γ_Ω1 = -μ11*γdot_Ω1
        γ_Ω2 = -μ11*γdot_Ω2
    else
        γ_Ω1 = @SMatrix zeros(3,3)
        γ_Ω2 = @SMatrix zeros(3,3)    
    end

    # linear and angular velocity rates
    V1dot_V1dot, Ω1dot_Ω1dot = point_displacement_jacobians(assembly.start[ielem], prescribed_conditions)
    V2dot_V2dot, Ω2dot_Ω2dot = point_displacement_jacobians(assembly.stop[ielem], prescribed_conditions)

    # linear and angular momentum

    P_V = CtCab*mass11*CtCab'
    P_Ω = CtCab*mass12*CtCab'

    H_V = CtCab*mass21*CtCab'
    H_Ω = CtCab*mass22*CtCab'
  
    # linear and angular momentum rates

    Pdot_V = CtCab*mass11*CtCabdot' + CtCabdot*mass11*CtCab'
    
    Pdot_Ω = CtCab*mass12*CtCabdot' + CtCabdot*mass12*CtCab' +
        CtCab*mass11*CtCab'*tilde(V) + CtCab*mass12*CtCab'*tilde(Ω) +
        -tilde(CtCab*mass11*CtCab'*V) - tilde(CtCab*mass12*CtCab'*Ω)
    
    Pdot_Vdot = CtCab*mass11*CtCab'
    
    Pdot_Ωdot = CtCab*mass12*CtCab'

    Hdot_V = CtCab*mass21*CtCabdot' + CtCabdot*mass21*CtCab'

    Hdot_Ω = CtCab*mass22*CtCabdot' + CtCabdot*mass22*CtCab' +
        CtCab*mass21*CtCab'*tilde(V) + CtCab*mass22*CtCab'*tilde(Ω) +
        -tilde(CtCab*mass21*CtCab'*V) - tilde(CtCab*mass22*CtCab'*Ω)

    Hdot_Vdot = CtCab*mass21*CtCab'

    Hdot_Ωdot = CtCab*mass22*CtCab'

    return (; properties..., γ_F, γ_M, γ_Ω1, γ_Ω2, κ_F, κ_M, V1dot_V1dot, V2dot_V2dot,
        Ω1dot_Ω1dot, Ω2dot_Ω2dot, P_V, P_Ω, H_V, H_Ω, Pdot_V, Pdot_Ω, Pdot_Vdot, Pdot_Ωdot,
        Hdot_V, Hdot_Ω, Hdot_Vdot, Hdot_Ωdot) 
end

"""
    newmark_element_jacobian_properties(properties, x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity, 
        x0, v0, ω0, a0, α0, Vdot_init, Ωdot_init, dt)

Calculate/extract the element properties needed to calculate the jacobian entries 
corresponding to a element for a Newmark scheme time marching analysis
"""
@inline function newmark_element_jacobian_properties(properties, x, indices, force_scaling, 
    structural_damping, assembly, ielem, prescribed_conditions, gravity, x0, v0, ω0, a0, α0,
    Vdot_init, Ωdot_init, dt)

    properties = steady_state_element_jacobian_properties(properties, x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity, x0, v0, ω0, a0, α0)

    @unpack C, Cab, CtCab, CtCabdot, mass11, mass12, mass21, mass22, V, Ω, ω, Vdot, Ωdot, 
        C_θ1, C_θ2, C_θ3 = properties

    # linear and angular momentum rates

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

    return (; properties..., Pdot_θ, Pdot_V, Pdot_Ω, Hdot_θ, Hdot_V, Hdot_Ω) 
end

"""
    dynamic_element_jacobian_properties(properties, dx, x, indices, force_scaling, 
        assembly, ielem, prescribed_conditions, gravity, x0, v0, ω0, a0, α0)

Calculate/extract the element properties needed to calculate the jacobian entries 
corresponding to a element for a dynamic analysis
"""
@inline function dynamic_element_jacobian_properties(properties, dx, x, indices, force_scaling, 
    structural_damping, assembly, ielem, prescribed_conditions, gravity, x0, v0, ω0, a0, α0)

    properties = steady_state_element_jacobian_properties(properties, x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity, x0, v0, ω0, a0, α0)

    @unpack C, Cab, CtCab, CtCabdot, mass11, mass12, mass21, mass22, V, Ω, Vdot, Ωdot, ω, C_θ1, C_θ2, C_θ3 = properties

    # linear and angular momentum rates

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

    return (; properties..., Pdot_θ, Pdot_V, Pdot_Ω, Hdot_θ, Hdot_V, Hdot_Ω) 
end

"""
    expanded_element_jacobian_properties(properties, x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, gravity, x0, v0, ω0, a0, α0))

Calculate/extract the element properties needed to calculate the jacobian entries 
corresponding to a element for a steady state analysis
"""
@inline function expanded_element_jacobian_properties(properties, x, indices, force_scaling, 
    structural_damping, assembly, ielem, prescribed_conditions, gravity, x0, v0, ω0, a0, α0)

    @unpack L, C1, C2, C, Cab, CtCab, Q, mass11, mass12, mass21, mass22, S11, S12, S21, S22,
        u1, u2, θ1, θ2, V1, V2, Ω1, Ω2, u, θ, V, Ω, ω = properties

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

    # linear and angular acceleration
    a_u = tilde(α0)    

    if structural_damping

        @unpack μ11, μ22, Qinv1, Qinv2, ω1, ω2, udot, θdot, Δu, Δθ, ΔQ, Δudot, Δθdot = properties

        # rotation parameter matrices
        Qinv1_θ1, Qinv1_θ2, Qinv1_θ3 = get_Qinv_θ(θ1)
        Qinv2_θ1, Qinv2_θ2, Qinv2_θ3 = get_Qinv_θ(θ2)

        # linear deflection rates
        udot1_u1 = -tilde(ω1)
        udot2_u2 = -tilde(ω2)

        udot1_θ1 = mul3(C1_θ1', C1_θ2', C1_θ3', V1)
        udot2_θ2 = mul3(C2_θ1', C2_θ2', C2_θ3', V2)

        udot1_V1 = C1'
        udot2_V2 = C2'

        # angular deflection rates
        θdot1_θ1 = mul3(Qinv1_θ1, Qinv1_θ2, Qinv1_θ3, Ω1 - C1*ω1) - Qinv1*mul3(C1_θ1, C1_θ2, C1_θ3, ω1)
        θdot2_θ2 = mul3(Qinv2_θ1, Qinv2_θ2, Qinv2_θ3, Ω2 - C2*ω2) - Qinv2*mul3(C2_θ1, C2_θ2, C2_θ3, ω2)
        
        θdot1_Ω1 = Qinv1
        θdot2_Ω2 = Qinv2

        θdot_θ1 = 1/2*θdot1_θ1
        θdot_θ2 = 1/2*θdot2_θ2

        θdot_Ω1 = 1/2*θdot1_Ω1
        θdot_Ω2 = 1/2*θdot2_Ω2

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
        κ_F, κ_M, κ_θ1, κ_θ2, κ_Ω1, κ_Ω2,
        P_V, P_Ω, H_V, H_Ω, a_u,
        )
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
        C_θ1, C_θ2, C_θ3, Qinv_θ1, Qinv_θ2, Qinv_θ3, = properties

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
   
    @unpack Cab, CtCab, Qinv, γ_F, γ_M, κ_F, κ_M, γ_Ω1, γ_Ω2 = properties

    Δu_F = CtCab*γ_F
    ru_F = -Δu_F

    Δu_M = CtCab*γ_M
    ru_M = -Δu_M

    ru_Ω1 = -CtCab*γ_Ω1
    ru_Ω2 = -CtCab*γ_Ω2

    Δθ_F = Qinv*Cab*κ_F
    rθ_F = -Δθ_F

    Δθ_M = Qinv*Cab*κ_M
    rθ_M = -Δθ_M
   
    return (; ru_F, ru_M, ru_Ω1, ru_Ω2, rθ_F, rθ_M)
end

@inline function dynamic_compatability_jacobians(properties)
   
    jacobians = static_compatability_jacobians(properties)

    @unpack Cab, CtCab, Qinv, u1_u1, u2_u2, θ1_θ1, θ2_θ2, γ_u1, γ_u2, γ_θ1, γ_θ2, 
        γ_V1, γ_V2, γ_Ω1, γ_Ω2, κ_θ1, κ_θ2, κ_Ω1, κ_Ω2 = properties

    @unpack ru_u1, ru_u2, ru_θ1, ru_θ2, rθ_θ1, rθ_θ2 = jacobians

    ru_u1 -= CtCab*γ_u1*u1_u1
    ru_u2 -= CtCab*γ_u2*u2_u2

    ru_θ1 -= CtCab*γ_θ1*θ1_θ1
    ru_θ2 -= CtCab*γ_θ2*θ2_θ2
    
    ru_V1 = -CtCab*γ_V1
    ru_V2 = -CtCab*γ_V2

    ru_Ω1 = -CtCab*γ_Ω1
    ru_Ω2 = -CtCab*γ_Ω2

    rθ_θ1 -= Qinv*Cab*κ_θ1*θ1_θ1
    rθ_θ2 -= Qinv*Cab*κ_θ2*θ2_θ2
    
    rθ_Ω1 = -Qinv*Cab*κ_Ω1
    rθ_Ω2 = -Qinv*Cab*κ_Ω2

    return (; jacobians..., ru_u1, ru_u2, ru_θ1, ru_θ2, ru_V1, ru_V2, ru_Ω1, ru_Ω2, 
        rθ_θ1, rθ_θ2, rθ_Ω1, rθ_Ω2)
end

@inline function expanded_compatability_jacobians(properties)
   
    @unpack L, Cab, CtCab, Qinv, γ, κ, u1_u1, u2_u2, θ1_θ1, θ2_θ2, C_θ1, C_θ2, C_θ3, 
        Qinv_θ1, Qinv_θ2, Qinv_θ3, γ_u1, γ_u2, γ_θ1, γ_θ2, 
        γ_V1, γ_V2, γ_Ω, γ_F, γ_M, κ_θ1, κ_θ2, κ_Ω1, κ_Ω2, κ_F, κ_M = properties

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

    return (; ru_u1, ru_u2, ru_θ1, ru_θ2, ru_V1, ru_V2, ru_Ω, 
        ru_F1, ru_F2, ru_M1, ru_M2, rθ_θ1, rθ_θ2, rθ_Ω1, rθ_Ω2, rθ_F1, rθ_F2, rθ_M1, rθ_M2)
end

"""
    expanded_element_velocity_jacobians(properties)

Calculate the jacobians of the element velocity residuals for a constant mass matrix system.
"""
@inline function expanded_element_velocity_jacobians(properties)

    @unpack C, Cab, CtCab, Qinv, V, Ω, ω, u1_u1, u2_u2, θ1_θ1, θ2_θ2, 
        C_θ1, C_θ2, C_θ3, Qinv_θ1, Qinv_θ2, Qinv_θ3 = properties
    
    tmp = -tilde(ω)
    rV_u1 = 1/2*tmp*u1_u1
    rV_u2 = 1/2*tmp*u2_u2

    tmp = mul3(C_θ1', C_θ2', C_θ3', Cab*V)
    rV_θ1 = 1/2*tmp*θ1_θ1
    rV_θ2 = 1/2*tmp*θ2_θ2

    rV_V = CtCab

    tmp = mul3(Qinv_θ1, Qinv_θ2, Qinv_θ3, Cab*Ω - C*ω) - Qinv*mul3(C_θ1, C_θ2, C_θ3, ω)
    rΩ_θ1 = 1/2*tmp*θ1_θ1
    rΩ_θ2 = 1/2*tmp*θ2_θ2

    rΩ_Ω = Qinv*Cab

    return (; rV_u1, rV_u2, rV_θ1, rV_θ2, rV_V, rΩ_θ1, rΩ_θ2, rΩ_Ω)
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
        V, Ω, P, H, F, M, γ, κ, v, ω, a, α, γ_F, γ_M, γ_u1, γ_u2, γ_θ1, γ_θ2, γ_V1, γ_V2, γ_Ω, 
        a_u, P_V, P_Ω, H_V, H_Ω, 
        C_θ1, C_θ2, C_θ3, u1_u1, u2_u2, θ1_θ1, θ2_θ2 = properties

    # initialize equilibrium residual
    rF_F1 = -I3
    rF_F2 =  I3

    rM_M1 = -I3
    rM_M2 =  I3

    # add loads due to internal loads and stiffness
    
    tmp1 =  tilde(L*e1 + γ)
    tmp2 = -tilde(F)  

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
    
    tmp = mass11*CtCab'*a_u
    rF_u1 = -1/2*tmp*u1_u1
    rF_u2 = -1/2*tmp*u2_u2

    tmp = mass21*CtCab'*a_u    
    rM_u1 -= 1/2*tmp*u1_u1
    rM_u2 -= 1/2*tmp*u2_u2

    tmp = mass11*Cab'*mul3(C_θ1, C_θ2, C_θ3, a) + mass12*Cab'*mul3(C_θ1, C_θ2, C_θ3, α)
    
    rF_θ1 = -1/2*tmp*θ1_θ1
    rF_θ2 = -1/2*tmp*θ2_θ2

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

    return (; rF_F1, rF_F2, rM_F1, rM_F2, rM_M1, rM_M2, 
        rF_u1, rF_u2, rF_θ1, rF_θ2, rF_V, rF_Ω,  
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
        γ_u1, γ_u2, γ_θ1, γ_θ2, γ_V1, γ_V2, γ_Ω1, γ_Ω2,
        a_u, P_θ, P_V, P_Ω, H_θ, H_V, H_Ω, u1_u1, u2_u2, θ1_θ1, θ2_θ2 = properties

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

    M1_V1 = -tmp*γ_V1
    M2_V1 =  tmp*γ_V1

    M1_V2 = -tmp*γ_V2
    M2_V2 =  tmp*γ_V2

    M1_Ω1 = -tmp*γ_Ω1
    M2_Ω1 =  tmp*γ_Ω1

    M1_Ω2 = -tmp*γ_Ω2
    M2_Ω2 =  tmp*γ_Ω2

    # add loads due to linear and angular acceleration (including gravity)
    tmp = 1/2*CtCab*mass11*CtCab'*a_u
    
    F1_u1 = -1/2*tmp*u1_u1
    F2_u1 = 1/2*tmp*u1_u1

    F1_u2 = -1/2*tmp*u2_u2    
    F2_u2 = 1/2*tmp*u2_u2

    tmp = 1/2*CtCab*mass21*CtCab'*a_u
    
    M1_u1 -= 1/2*tmp*u1_u1
    M2_u1 += 1/2*tmp*u1_u1

    M1_u2 -= 1/2*tmp*u2_u2
    M2_u2 += 1/2*tmp*u2_u2

    # add loads due to linear and angular momentum

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
        F1_u1, F1_u2, F1_θ1, F1_θ2, F1_V1, F1_V2, F1_Ω1, F1_Ω2, 
        F2_u1, F2_u2, F2_θ1, F2_θ2, F2_V1, F2_V2, F2_Ω1, F2_Ω2, 
        M1_u1, M1_u2, M1_θ1, M1_θ2, M1_V1, M1_V2, M1_Ω1, M1_Ω2, 
        M2_u1, M2_u2, M2_θ1, M2_θ2, M2_V1, M2_V2, M2_Ω1, M2_Ω2)

end

"""
    initial_condition_element_resultant_jacobians(properties)

Calculate the jacobians for the resultant loads applied at each end of V beam element 
for the initialization of a time domain analysis.
"""
@inline function initial_condition_element_resultant_jacobians(properties)

    @unpack L, CtCab, F, M, γ, κ, V, Ω, P, H, v, ω, γ_F, γ_M, γ_Ω1, γ_Ω2, P_V, P_Ω, H_V, H_Ω, 
        Pdot_V, Pdot_Ω, Pdot_Vdot, Pdot_Ωdot, Hdot_V, Hdot_Ω, Hdot_Vdot, Hdot_Ωdot, 
        V1dot_V1dot, Ω1dot_Ω1dot, V2dot_V2dot, Ω2dot_Ω2dot = properties
    
    # loads due to internal forces/moments
    F1_F = CtCab
    F2_F = CtCab

    tmp = 1/2*CtCab*(tilde(L*e1 + γ) - tilde(F)*γ_F)
    M1_F = tmp
    M2_F = -tmp

    tmp = 1/2*CtCab*tilde(F)*γ_M
    M1_M = CtCab - tmp
    M2_M = CtCab + tmp

    tmp = 1/2*CtCab*tilde(F)
    M1_Ω1 = -tmp*γ_Ω1
    M2_Ω1 =  tmp*γ_Ω1

    M1_Ω2 = -tmp*γ_Ω2
    M2_Ω2 =  tmp*γ_Ω2

    # add loads due to linear and angular momentum

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

    tmp = 1/2*(tilde(ω)*H_V + tilde(V)*P_V - tilde(P))
    
    M1_V1 = -1/2*tmp
    M2_V1 = 1/2*tmp

    M1_V2 = -1/2*tmp
    M2_V2 = 1/2*tmp

    tmp = 1/2*(tilde(ω)*H_Ω + tilde(V)*P_Ω)
    
    M1_Ω1 -= 1/2*tmp
    M2_Ω1 += 1/2*tmp

    M1_Ω2 -= 1/2*tmp
    M2_Ω2 += 1/2*tmp

    # add loads due to linear and angular momentum rates

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

    tmp = 1/2*Pdot_Vdot

    F1_V1dot = -1/2*tmp*V1dot_V1dot
    F2_V1dot =  1/2*tmp*V1dot_V1dot

    F1_V2dot = -1/2*tmp*V2dot_V2dot
    F2_V2dot =  1/2*tmp*V2dot_V2dot

    tmp = 1/2*Pdot_Ωdot

    F1_Ω1dot = -1/2*tmp*Ω1dot_Ω1dot
    F2_Ω1dot =  1/2*tmp*Ω1dot_Ω1dot

    F1_Ω2dot = -1/2*tmp*Ω2dot_Ω2dot
    F2_Ω2dot =  1/2*tmp*Ω2dot_Ω2dot

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

    return (;
        F1_V1dot, F1_V2dot, F1_Ω1dot, F1_Ω2dot, F1_F, F1_V1, F1_V2, F1_Ω1, F1_Ω2, 
        F2_V1dot, F2_V2dot, F2_Ω1dot, F2_Ω2dot, F2_F, F2_V1, F2_V2, F2_Ω1, F2_Ω2, 
        M1_V1dot, M1_V2dot, M1_Ω1dot, M1_Ω2dot, M1_F, M1_M, M1_V1, M1_V2, M1_Ω1, M1_Ω2, 
        M2_V1dot, M2_V2dot, M2_Ω1dot, M2_Ω2dot, M2_F, M2_M, M2_V1, M2_V2, M2_Ω1, M2_Ω2)
end

"""
    newmark_element_resultant_jacobians(properties, distributed_loads, ielem)

Calculate the jacobians for the resultant loads applied at each end of a beam element 
for a Newmark scheme time marching analysis.
"""
@inline function newmark_element_resultant_jacobians(properties, distributed_loads, ielem)

    jacobians = steady_state_element_resultant_jacobians(properties, distributed_loads, ielem)

    @unpack L, Pdot_θ, Pdot_V, Pdot_Ω, Hdot_θ, Hdot_V, Hdot_Ω, θ1_θ1, θ2_θ2 = properties
    
    @unpack F1_θ1, F1_θ2, F1_V1, F1_V2, F1_Ω1, F1_Ω2, 
            F2_θ1, F2_θ2, F2_V1, F2_V2, F2_Ω1, F2_Ω2, 
            M1_θ1, M1_θ2, M1_V1, M1_V2, M1_Ω1, M1_Ω2, 
            M2_θ1, M2_θ2, M2_V1, M2_V2, M2_Ω1, M2_Ω2 = jacobians

    # add loads due to linear and angular momentum rates

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
        F1_θ1, F1_θ2, F1_V1, F1_V2, F1_Ω1, F1_Ω2, 
        F2_θ1, F2_θ2, F2_V1, F2_V2, F2_Ω1, F2_Ω2, 
        M1_θ1, M1_θ2, M1_V1, M1_V2, M1_Ω1, M1_Ω2, 
        M2_θ1, M2_θ2, M2_V1, M2_V2, M2_Ω1, M2_Ω2)
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

    @unpack ru_F, ru_M, ru_Ω1, ru_Ω2, rθ_F, rθ_M = compatability

    @unpack F1_V1dot, F1_V2dot, F1_Ω1dot, F1_Ω2dot, F1_F, F1_V1, F1_V2, F1_Ω1, F1_Ω2,
            F2_V1dot, F2_V2dot, F2_Ω1dot, F2_Ω2dot, F2_F, F2_V1, F2_V2, F2_Ω1, F2_Ω2,
            M1_V1dot, M1_V2dot, M1_Ω1dot, M1_Ω2dot, M1_F, M1_M, M1_V1, M1_V2, M1_Ω1, M1_Ω2,
            M2_V1dot, M2_V2dot, M2_Ω1dot, M2_Ω2dot, M2_F, M2_M, M2_V1, M2_V2, M2_Ω1, M2_Ω2 = resultants

    icol = indices.icol_elem[ielem]
    icol1 = indices.icol_point[assembly.start[ielem]]
    icol2 = indices.icol_point[assembly.stop[ielem]]

    # compatability equations
    irow = indices.irow_elem[ielem]

    jacob[irow:irow+2, icol:icol+2] .= ru_F .* force_scaling
    jacob[irow:irow+2, icol+3:icol+5] .= ru_M .* force_scaling

    jacob[irow:irow+2, icol1+9:icol1+11] .= ru_Ω1
    
    jacob[irow:irow+2, icol2+9:icol2+11] .= ru_Ω2

    jacob[irow+3:irow+5, icol:icol+2] .= rθ_F .* force_scaling
    jacob[irow+3:irow+5, icol+3:icol+5] .= rθ_M .* force_scaling
    
    # equilibrium equations for the start of the beam element
    irow1 = indices.irow_point[assembly.start[ielem]]

    @views jacob[irow1:irow1+2, icol1:icol1+2] .-= F1_V1dot ./ force_scaling
    @views jacob[irow1:irow1+2, icol1+3:icol1+5] .-= F1_Ω1dot ./ force_scaling

    @views jacob[irow1:irow1+2, icol2:icol2+2] .-= F1_V2dot ./ force_scaling
    @views jacob[irow1:irow1+2, icol2+3:icol2+5] .-= F1_Ω2dot ./ force_scaling

    jacob[irow1:irow1+2, icol:icol+2] .= -F1_F

    @views jacob[irow1:irow1+2, icol1+6:icol1+8] .-= F1_V1 ./ force_scaling
    @views jacob[irow1:irow1+2, icol1+9:icol1+11] .-= F1_Ω1 ./ force_scaling

    @views jacob[irow1:irow1+2, icol2+6:icol2+8] .-= F1_V2 ./ force_scaling
    @views jacob[irow1:irow1+2, icol2+9:icol2+11] .-= F1_Ω2 ./ force_scaling

    @views jacob[irow1+3:irow1+5, icol1:icol1+2] .-= M1_V1dot ./ force_scaling
    @views jacob[irow1+3:irow1+5, icol1+3:icol1+5] .-= M1_Ω1dot ./ force_scaling

    @views jacob[irow1+3:irow1+5, icol2:icol2+2] .-= M1_V2dot ./ force_scaling
    @views jacob[irow1+3:irow1+5, icol2+3:icol2+5] .-= M1_Ω2dot ./ force_scaling

    jacob[irow1+3:irow1+5, icol:icol+2] .= -M1_F
    jacob[irow1+3:irow1+5, icol+3:icol+5] .= -M1_M

    @views jacob[irow1+3:irow1+5, icol1+6:icol1+8] .-= M1_V1 ./ force_scaling
    @views jacob[irow1+3:irow1+5, icol1+9:icol1+11] .-= M1_Ω1 ./ force_scaling

    @views jacob[irow1+3:irow1+5, icol2+6:icol2+8] .-= M1_V2 ./ force_scaling
    @views jacob[irow1+3:irow1+5, icol2+9:icol2+11] .-= M1_Ω2 ./ force_scaling

    # equilibrium equations for the end of the beam element
    irow2 = indices.irow_point[assembly.stop[ielem]]

    @views jacob[irow2:irow2+2, icol1:icol1+2] .+= F2_V1dot ./ force_scaling
    @views jacob[irow2:irow2+2, icol1+3:icol1+5] .+= F2_Ω1dot ./ force_scaling

    @views jacob[irow2:irow2+2, icol2:icol2+2] .+= F2_V2dot ./ force_scaling
    @views jacob[irow2:irow2+2, icol2+3:icol2+5] .+= F2_Ω2dot ./ force_scaling

    jacob[irow2:irow2+2, icol:icol+2] .= F2_F

    @views jacob[irow2:irow2+2, icol1+6:icol1+8] .+= F2_V1 ./ force_scaling
    @views jacob[irow2:irow2+2, icol1+9:icol1+11] .+= F2_Ω1 ./ force_scaling

    @views jacob[irow2:irow2+2, icol2+6:icol2+8] .+= F2_V2 ./ force_scaling
    @views jacob[irow2:irow2+2, icol2+9:icol2+11] .+= F2_Ω2 ./ force_scaling

    @views jacob[irow2+3:irow2+5, icol1:icol1+2] .+= M2_V1dot ./ force_scaling
    @views jacob[irow2+3:irow2+5, icol1+3:icol1+5] .+= M2_Ω1dot ./ force_scaling

    @views jacob[irow2+3:irow2+5, icol2:icol2+2] .+= M2_V2dot ./ force_scaling
    @views jacob[irow2+3:irow2+5, icol2+3:icol2+5] .+= M2_Ω2dot ./ force_scaling

    jacob[irow2+3:irow2+5, icol:icol+2] .= M2_F
    jacob[irow2+3:irow2+5, icol+3:icol+5] .= M2_M

    @views jacob[irow2+3:irow2+5, icol1+6:icol1+8] .+= M2_V1 ./ force_scaling
    @views jacob[irow2+3:irow2+5, icol1+9:icol1+11] .+= M2_Ω1 ./ force_scaling

    @views jacob[irow2+3:irow2+5, icol2+6:icol2+8] .+= M2_V2 ./ force_scaling
    @views jacob[irow2+3:irow2+5, icol2+9:icol2+11] .+= M2_Ω2 ./ force_scaling

    return jacob
end

@inline function insert_dynamic_element_jacobians!(jacob, indices, force_scaling,
    assembly, ielem, compatability, resultants)

    insert_static_element_jacobians!(jacob, indices, force_scaling, assembly, ielem, 
        compatability, resultants)

    @unpack ru_V1, ru_V2, ru_Ω1, ru_Ω2, rθ_Ω1, rθ_Ω2 = compatability

    @unpack F1_u1, F1_u2, F1_V1, F1_V2, F1_Ω1, F1_Ω2,
            F2_u1, F2_u2, F2_V1, F2_V2, F2_Ω1, F2_Ω2,
            M1_u1, M1_u2, M1_V1, M1_V2, M1_Ω1, M1_Ω2,
            M2_u1, M2_u2, M2_V1, M2_V2, M2_Ω1, M2_Ω2 = resultants

    icol1 = indices.icol_point[assembly.start[ielem]]
    icol2 = indices.icol_point[assembly.stop[ielem]]

    # compatability equations
    irow = indices.irow_elem[ielem]

    jacob[irow:irow+2, icol1+6:icol1+8] .= ru_V1
    jacob[irow:irow+2, icol1+9:icol1+11] .= ru_Ω1

    jacob[irow:irow+2, icol2+6:icol2+8] .= ru_V2
    jacob[irow:irow+2, icol2+9:icol2+11] .= ru_Ω2

    jacob[irow+3:irow+5, icol1+9:icol1+11] .= rθ_Ω1

    jacob[irow+3:irow+5, icol2+9:icol2+11] .= rθ_Ω2

    # equilibrium equations for the start of the beam element
    irow1 = indices.irow_point[assembly.start[ielem]]

    @views jacob[irow1:irow1+2, icol1:icol1+2] .-= F1_u1 ./ force_scaling

    @views jacob[irow1:irow1+2, icol2:icol2+2] .-= F1_u2 ./ force_scaling

    @views jacob[irow1:irow1+2, icol1+6:icol1+8] .-= F1_V1 ./ force_scaling
    @views jacob[irow1:irow1+2, icol1+9:icol1+11] .-= F1_Ω1 ./ force_scaling

    @views jacob[irow1:irow1+2, icol2+6:icol2+8] .-= F1_V2 ./ force_scaling
    @views jacob[irow1:irow1+2, icol2+9:icol2+11] .-= F1_Ω2 ./ force_scaling

    @views jacob[irow1+3:irow1+5, icol1:icol1+2] .-= M1_u1 ./ force_scaling
    
    @views jacob[irow1+3:irow1+5, icol2:icol2+2] .-= M1_u2 ./ force_scaling

    @views jacob[irow1+3:irow1+5, icol1+6:icol1+8] .-= M1_V1 ./ force_scaling
    @views jacob[irow1+3:irow1+5, icol1+9:icol1+11] .-= M1_Ω1 ./ force_scaling

    @views jacob[irow1+3:irow1+5, icol2+6:icol2+8] .-= M1_V2 ./ force_scaling
    @views jacob[irow1+3:irow1+5, icol2+9:icol2+11] .-= M1_Ω2 ./ force_scaling

    # equilibrium equations for the end of the beam element
    irow2 = indices.irow_point[assembly.stop[ielem]]

    @views jacob[irow2:irow2+2, icol1:icol1+2] .+= F2_u1 ./ force_scaling

    @views jacob[irow2:irow2+2, icol2:icol2+2] .+= F2_u2 ./ force_scaling

    @views jacob[irow2:irow2+2, icol1+6:icol1+8] .+= F2_V1 ./ force_scaling
    @views jacob[irow2:irow2+2, icol1+9:icol1+11] .+= F2_Ω1 ./ force_scaling

    @views jacob[irow2:irow2+2, icol2+6:icol2+8] .+= F2_V2 ./ force_scaling
    @views jacob[irow2:irow2+2, icol2+9:icol2+11] .+= F2_Ω2 ./ force_scaling

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

    @unpack ru_u1, ru_u2, ru_θ1, ru_θ2, ru_V1, ru_V2, ru_Ω, ru_F1, ru_F2, ru_M1, ru_M2, 
                          rθ_θ1, rθ_θ2,                 rθ_Ω1, rθ_Ω2, rθ_F1, rθ_F2, rθ_M1, rθ_M2 = compatability
    
    @unpack rV_u1, rV_u2, rV_θ1, rV_θ2, rV_V, 
                          rΩ_θ1, rΩ_θ2, rΩ_Ω = velocities

    @unpack rF_F1, rF_F2,                 rF_u1, rF_u2, rF_θ1, rF_θ2, rF_V, rF_Ω,  
            rM_F1, rM_F2, rM_M1, rM_M2, rM_u1, rM_u2, rM_θ1, rM_θ2, rM_V, rM_Ω,  
            rM_V1, rM_V2 = equilibrium
    
    @unpack F1_θ1, F1_θ2, F1_F1, F2_θ1, F2_θ2, F2_F2, 
            M1_θ1, M1_θ2, M1_M1, M2_θ1, M2_θ2, M2_M2 = resultants

    icol = indices.icol_elem[ielem]
    icol1 = indices.icol_point[assembly.start[ielem]]
    icol2 = indices.icol_point[assembly.stop[ielem]]

    irow = indices.irow_elem[ielem]

    # element compatability residuals
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

    jacob[irow+3:irow+5, icol1+3:icol1+5] .= rθ_θ1
    jacob[irow+3:irow+5, icol2+3:icol2+5] .= rθ_θ2
    jacob[irow+3:irow+5, icol1+9:icol1+11] .= rθ_Ω1
    jacob[irow+3:irow+5, icol2+9:icol2+11] .= rθ_Ω2
    jacob[irow+3:irow+5, icol:icol+2] .= rθ_F1 .* force_scaling
    jacob[irow+3:irow+5, icol+3:icol+5] .= rθ_M1 .* force_scaling
    jacob[irow+3:irow+5, icol+6:icol+8] .= rθ_F2 .* force_scaling
    jacob[irow+3:irow+5, icol+9:icol+11] .= rθ_M2 .* force_scaling

    # element equilibrium residuals
    jacob[irow+6:irow+8, icol1:icol1+2] .= rF_u1 ./ force_scaling
    jacob[irow+6:irow+8, icol2:icol2+2] .= rF_u2 ./ force_scaling
    jacob[irow+6:irow+8, icol1+3:icol1+5] .= rF_θ1 ./ force_scaling
    jacob[irow+6:irow+8, icol2+3:icol2+5] .= rF_θ2 ./ force_scaling
    jacob[irow+6:irow+8, icol:icol+2] .= rF_F1
    jacob[irow+6:irow+8, icol+6:icol+8] .= rF_F2
    jacob[irow+6:irow+8, icol+12:icol+14] .= rF_V ./ force_scaling
    jacob[irow+6:irow+8, icol+15:icol+17] .= rF_Ω ./ force_scaling

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
    jacob[irow+12:irow+14, icol1:icol1+2] .= rV_u1
    jacob[irow+12:irow+14, icol2:icol2+2] .= rV_u2
    jacob[irow+12:irow+14, icol1+3:icol1+5] .= rV_θ1
    jacob[irow+12:irow+14, icol2+3:icol2+5] .= rV_θ2
    jacob[irow+12:irow+14, icol+12:icol+14] .= rV_V

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
    steady_state_element_jacobian!(jacob, x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, distributed_loads, gravity, x0, v0, ω0, a0, α0)

Calculate and insert the jacobian entries corresponding to a beam element for a steady state 
analysis into the system jacobian matrix.
"""
@inline function steady_state_element_jacobian!(jacob, x, indices, force_scaling, structural_damping, 
    assembly, ielem, prescribed_conditions, distributed_loads, gravity, x0, v0, ω0, a0, α0)

    properties = steady_state_element_properties(x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, gravity, x0, v0, ω0, a0, α0)

    properties = steady_state_element_jacobian_properties(properties, x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity, x0, v0, ω0, a0, α0)

    compatability = dynamic_compatability_jacobians(properties)

    resultants = steady_state_element_resultant_jacobians(properties, distributed_loads, ielem)

    insert_dynamic_element_jacobians!(jacob, indices, force_scaling, assembly, ielem, 
        compatability, resultants)

    return jacob
end

"""
    initial_condition_element_jacobian!(jacob, x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, distributed_loads, gravity, x0, v0, ω0, a0, α0,
        u0, θ0, udot0, θdot0)

Calculate and insert the jacobian entries corresponding to a beam element for the 
initialization of a time domain analysis into the system jacobian matrix.
"""
@inline function initial_condition_element_jacobian!(jacob, x, indices, force_scaling, structural_damping, 
    assembly, ielem, prescribed_conditions, distributed_loads, gravity, x0, v0, ω0, a0, α0,
    u0, θ0, udot0, θdot0)

    properties = initial_condition_element_properties(x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, gravity, x0, v0, ω0, a0, α0, u0, θ0, udot0, θdot0)

    properties = initial_condition_element_jacobian_properties(properties, x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity)

    compatability = initial_condition_compatability_jacobians(properties)

    resultants = initial_condition_element_resultant_jacobians(properties)

    insert_initial_condition_element_jacobians!(jacob, indices, force_scaling, assembly, ielem, 
        compatability, resultants)

    return jacob
end

"""
    newmark_element_jacobian!(jacob, x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, distributed_loads, gravity, x0, v0, ω0, a0, α0,
        Vdot_init, Ωdot_init, dt)

Calculate and insert the jacobian entries corresponding to a beam element for a Newmark-scheme
time marching analysis into the system jacobian matrix.
"""
@inline function newmark_element_jacobian!(jacob, x, indices, force_scaling, structural_damping, 
    assembly, ielem, prescribed_conditions, distributed_loads, gravity, x0, v0, ω0, a0, α0,
    Vdot_init, Ωdot_init, dt)

    properties = newmark_element_properties(x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, gravity, x0, v0, ω0, a0, α0, 
        Vdot_init, Ωdot_init, dt)

    properties = newmark_element_jacobian_properties(properties, x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity, 
        x0, v0, ω0, a0, α0, Vdot_init, Ωdot_init, dt)

    compatability = dynamic_compatability_jacobians(properties)

    resultants = newmark_element_resultant_jacobians(properties, distributed_loads, ielem)

    insert_dynamic_element_jacobians!(jacob, indices, force_scaling, assembly, ielem, 
        compatability, resultants)

    return jacob
end

"""
    dynamic_element_jacobian!(jacob, dx, x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, distributed_loads, gravity, 
        x0, v0, ω0, a0, α0)

Calculate and insert the jacobian entries corresponding to a beam element for a dynamic
analysis into the system jacobian matrix.
"""
@inline function dynamic_element_jacobian!(jacob, dx, x, indices, force_scaling, structural_damping, 
    assembly, ielem, prescribed_conditions, distributed_loads, gravity, x0, v0, ω0, a0, α0)

    properties = dynamic_element_properties(dx, x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, gravity, x0, v0, ω0, a0, α0)

    properties = dynamic_element_jacobian_properties(properties, dx, x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity, x0, v0, ω0, a0, α0)

    compatability = dynamic_compatability_jacobians(properties)

    resultants = dynamic_element_resultant_jacobians(properties, distributed_loads, ielem)

    insert_dynamic_element_jacobians!(jacob, indices, force_scaling, assembly, ielem, 
        compatability, resultants)

    return jacob
end

"""
    expanded_element_jacobian!(jacob, x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, distributed_loads, gravity, 
        x0, v0, ω0, a0, α0)

Calculate and insert the jacobian entries corresponding to a beam element for a dynamic
analysis into the system jacobian matrix.
"""
@inline function expanded_element_jacobian!(jacob, x, indices, force_scaling, structural_damping, 
    assembly, ielem, prescribed_conditions, distributed_loads, gravity, x0, v0, ω0, a0, α0)

    properties = expanded_element_properties(x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, gravity, x0, v0, ω0, a0, α0)

    properties = expanded_element_jacobian_properties(properties, x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity, x0, v0, ω0, a0, α0)

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
