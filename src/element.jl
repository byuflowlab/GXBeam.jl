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

    # forces and moments (in the deformed element frame)
    F, M = element_loads(x, ielem, indices.icol_elem, force_scaling)
    
    # strain and curvature (in the deformed element frame)
    γ = S11*F + S12*M
    κ = S21*F + S22*M

    # linear and angular acceleration
    a = -gravity
    α = @SVector zeros(3)

    return (; L, C, Cab, CtCab, Qinv, S11, S12, S21, S22, mass11, mass12, mass21, mass22, 
        u1, u2, θ1, θ2, u, θ, F, M, γ, κ, a, α)
end

"""
    steady_state_element_properties(x, indices, force_scaling, assembly, ielem, 
        prescribed_conditions, gravity, x0, v0, ω0, a0, α0)

Calculate/extract the element properties needed to construct the residual for a steady 
state analysis
"""
@inline function steady_state_element_properties(x, indices, force_scaling, assembly, ielem, 
    prescribed_conditions, gravity, x0, v0, ω0, a0, α0)

    properties = static_element_properties(x, indices, force_scaling, assembly, ielem, 
        prescribed_conditions, gravity)

    @unpack CtCab, mass11, mass12, mass21, mass22, u, θ, a, α = properties

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

    # undeformed linear and angular velocity
    v = v0 + cross(ω0, Δx)
    ω = ω0

    # linear and angular acceleration
    a += a0 + cross(α0, Δx) + cross(α0, u)
    α += α0

    return (; properties..., V1, V2, Ω1, Ω2, V, Ω, P, H, v, ω, a, α) 
end

"""
    initial_condition_element_properties(x, indices, force_scaling, assembly, ielem,
        prescribed_conditions, gravity, x0, v0, ω0, a0, α0, u0, θ0, udot0, θdot0)

Calculate/extract the element properties needed to construct the residual for a time-domaing
analysis initialization
"""
@inline function initial_condition_element_properties(x, indices, force_scaling, assembly, ielem,
    prescribed_conditions, gravity, x0, v0, ω0, a0, α0, u0, θ0, udot0, θdot0)

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
    u1 = SVector{3}(u0[assembly.start[ielem]])
    θ1 = SVector{3}(θ0[assembly.start[ielem]])

    u2 = SVector{3}(u0[assembly.stop[ielem]])
    θ2 = SVector{3}(θ0[assembly.stop[ielem]])

    u = (u1 + u2)/2
    θ = (θ1 + θ2)/2

    # rotation parameter matrices
    C = get_C(θ)
    CtCab = C'*Cab
    Qinv = get_Qinv(θ)

    # linear and angular velocity
    V1, Ω1 = point_velocities(x, assembly.start[ielem], indices.icol_point)
    V2, Ω2 = point_velocities(x, assembly.stop[ielem], indices.icol_point)
    V = (V1 + V2)/2
    Ω = (Ω1 + Ω2)/2

    # linear and angular momentum
    P = CtCab*mass11*CtCab'*V + CtCab*mass12*CtCab'*Ω
    H = CtCab*mass21*CtCab'*V + CtCab*mass22*CtCab'*Ω

    # forces and moments (in the deformed element frame)
    F, M = element_loads(x, ielem, indices.icol_elem, force_scaling)

    # strain and curvature (in the deformed element frame)
    γ = S11*F + S12*M
    κ = S21*F + S22*M

    # distance from the rotation center
    Δx = assembly.elements[ielem].x - x0

    # undeformed linear and angular velocity
    v = v0 + cross(ω0, Δx)
    ω = ω0

    # linear and angular acceleration
    a = a0 + cross(α0, Δx) + cross(α0, u) - gravity
    α = α0

    # linear and angular velocity rates
    V1dot, Ω1dot = point_displacement_rates(x, assembly.start[ielem], indices.icol_point, prescribed_conditions)
    V2dot, Ω2dot = point_displacement_rates(x, assembly.stop[ielem], indices.icol_point, prescribed_conditions)
    Vdot = (V1dot + V2dot)/2
    Ωdot = (Ω1dot + Ω2dot)/2

    # linear and angular momentum rates
    CtCabdot = CtCab*tilde(Ω - ω)

    Pdot = CtCab*mass11*CtCab'*Vdot + CtCab*mass12*CtCab'*Ωdot +
        CtCab*mass11*CtCabdot'*V + CtCab*mass12*CtCabdot'*Ω +
        CtCabdot*mass11*CtCab'*V + CtCabdot*mass12*CtCab'*Ω
    
    Hdot = CtCab*mass21*CtCab'*Vdot + CtCab*mass22*CtCab'*Ωdot +
        CtCab*mass21*CtCabdot'*V + CtCab*mass22*CtCabdot'*Ω +
        CtCabdot*mass21*CtCab'*V + CtCabdot*mass22*CtCab'*Ω

    return (; L, C, Cab, CtCab, Qinv, S11, S12, S21, S22, mass11, mass12, mass21, mass22, 
        u1, u2, θ1, θ2, V1, V2, Ω1, Ω2, u, θ, V, Ω, P, H, F, M, γ, κ, v, ω, a, α, 
        CtCabdot, Vdot, Ωdot, Pdot, Hdot) 
end

"""
    newmark_element_properties(x, indices, force_scaling, assembly, ielem,
        prescribed_conditions, gravity, x0, v0, ω0, a0, α0, Vdot_init, Ωdot_init, dt)

Calculate/extract the element properties needed to construct the residual for a newmark-
scheme time stepping analysis
"""
@inline function newmark_element_properties(x, indices, force_scaling, assembly, ielem,
    prescribed_conditions, gravity, x0, v0, ω0, a0, α0, Vdot_init, Ωdot_init, dt)

    properties = steady_state_element_properties(x, indices, force_scaling, assembly, ielem, 
        prescribed_conditions, gravity, x0, v0, ω0, a0, α0)

    @unpack CtCab, mass11, mass12, mass21, mass22, V1, V2, Ω1, Ω2, V, Ω, ω = properties

    # linear and angular velocity rates
    V1dot = 2/dt*V1 - SVector{3}(Vdot_init[assembly.start[ielem]])
    Ω1dot = 2/dt*Ω1 - SVector{3}(Ωdot_init[assembly.start[ielem]])

    V2dot = 2/dt*V2 - SVector{3}(Vdot_init[assembly.stop[ielem]])
    Ω2dot = 2/dt*Ω2 - SVector{3}(Ωdot_init[assembly.stop[ielem]])

    Vdot = (V1dot + V2dot)/2
    Ωdot = (Ω1dot + Ω2dot)/2

    # linear and angular momentum rates
    CtCabdot = CtCab*tilde(Ω - ω)
    
    Pdot = CtCab*mass11*CtCab'*Vdot + CtCab*mass12*CtCab'*Ωdot +
        CtCab*mass11*CtCabdot'*V + CtCab*mass12*CtCabdot'*Ω +
        CtCabdot*mass11*CtCab'*V + CtCabdot*mass12*CtCab'*Ω
    
    Hdot = CtCab*mass21*CtCab'*Vdot + CtCab*mass22*CtCab'*Ωdot +
        CtCab*mass21*CtCabdot'*V + CtCab*mass22*CtCabdot'*Ω +
        CtCabdot*mass21*CtCab'*V + CtCabdot*mass22*CtCab'*Ω

    return (; properties..., CtCabdot, Vdot, Ωdot, Pdot, Hdot) 
end

"""
    dynamic_element_properties(dx, x, indices, force_scaling, assembly, ielem,
        prescribed_conditions, gravity, x0, v0, ω0, a0, α0)

Calculate/extract the element properties needed to construct the residual for a dynamic 
analysis
"""
@inline function dynamic_element_properties(dx, x, indices, force_scaling, assembly, ielem,
    prescribed_conditions, gravity, x0, v0, ω0, a0, α0)

    properties = steady_state_element_properties(x, indices, force_scaling, assembly, ielem, 
        prescribed_conditions, gravity, x0, v0, ω0, a0, α0)

    @unpack CtCab, mass11, mass12, mass21, mass22, V, Ω, ω = properties

    # linear and angular velocity of the beam element
    V1dot, Ω1dot = point_velocities(dx, assembly.start[ielem], indices.icol_point)
    V2dot, Ω2dot = point_velocities(dx, assembly.stop[ielem], indices.icol_point)
    Vdot = (V1dot + V2dot)/2
    Ωdot = (Ω1dot + Ω2dot)/2

    # linear and angular momentum rates
    CtCabdot = CtCab*tilde(Ω - ω)
    
    Pdot = CtCab*mass11*CtCab'*Vdot + CtCab*mass12*CtCab'*Ωdot +
        CtCab*mass11*CtCabdot'*V + CtCab*mass12*CtCabdot'*Ω +
        CtCabdot*mass11*CtCab'*V + CtCabdot*mass12*CtCab'*Ω
    
    Hdot = CtCab*mass21*CtCab'*Vdot + CtCab*mass22*CtCab'*Ωdot +
        CtCab*mass21*CtCabdot'*V + CtCab*mass22*CtCabdot'*Ω +
        CtCabdot*mass21*CtCab'*V + CtCabdot*mass22*CtCab'*Ω

    return (; properties..., CtCabdot, Vdot, Ωdot, Pdot, Hdot)
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
    steady_state_element_residual!(resid, x, indices, force_scaling, assembly, ielem,  
        prescribed_conditions, distributed_loads, gravity, x0, v0, ω0, a0, α0)

Calculate and insert the residual entries corresponding to a beam element for a steady state 
analysis into the system residual vector.
"""
@inline function steady_state_element_residual!(resid, x, indices, force_scaling, assembly, ielem,  
    prescribed_conditions, distributed_loads, gravity, x0, v0, ω0, a0, α0)

    properties = steady_state_element_properties(x, indices, force_scaling, assembly, ielem,
        prescribed_conditions, gravity, x0, v0, ω0, a0, α0)

    compatability = compatability_residuals(properties)

    resultants = steady_state_element_resultants(properties, distributed_loads, ielem)

    insert_element_residuals!(resid, indices, force_scaling, assembly, ielem, 
        compatability, resultants)

    return resid
end

"""
    initial_condition_element_residual!(resid, x, indices, force_scaling, assembly, ielem,  
        prescribed_conditions, distributed_loads, gravity, x0, v0, ω0, a0, α0, 
        u0, θ0, udot0, θdot0)

Calculate and insert the residual entries corresponding to a beam element for the 
initialization of a time domain simulation into the system residual vector.
"""
@inline function initial_condition_element_residual!(resid, x, indices, force_scaling, assembly, ielem,  
    prescribed_conditions, distributed_loads, gravity, x0, v0, ω0, a0, α0,
    u0, θ0, udot0, θdot0)

    properties = initial_condition_element_properties(x, indices, force_scaling, assembly, ielem,
        prescribed_conditions, gravity, x0, v0, ω0, a0, α0, u0, θ0, udot0, θdot0)

    compatability = compatability_residuals(properties)

    resultants = dynamic_element_resultants(properties, distributed_loads, ielem)

    insert_element_residuals!(resid, indices, force_scaling, assembly, ielem, 
        compatability, resultants)

    return resid
end

"""
    newmark_element_residual!(resid, x, indices, force_scaling, assembly, ielem,  
        prescribed_conditions, distributed_loads, gravity, x0, v0, ω0, a0, α0,
        Vdot_init, Ωdot_init, dt)

Calculate and insert the residual entries corresponding to a beam element for a 
newmark-scheme time marching analysis into the system residual vector.
"""
@inline function newmark_element_residual!(resid, x, indices, force_scaling, assembly, ielem,  
    prescribed_conditions, distributed_loads, gravity, x0, v0, ω0, a0, α0, 
    Vdot_init, Ωdot_init, dt)

    properties = newmark_element_properties(x, indices, force_scaling, assembly, ielem, 
        prescribed_conditions, gravity, x0, v0, ω0, a0, α0, 
        Vdot_init, Ωdot_init, dt)

    compatability = compatability_residuals(properties)

    resultants = dynamic_element_resultants(properties, distributed_loads, ielem)

    insert_element_residuals!(resid, indices, force_scaling, assembly, ielem, 
        compatability, resultants)

    return resid
end

"""
    dynamic_element_residual!(resid, dx, x, indices, force_scaling, assembly, ielem,  
        prescribed_conditions, distributed_loads, gravity, x0, v0, ω0, a0, α0)

Calculate and insert the residual entries corresponding to a beam element for a dynamic
analysis into the system residual vector.
"""
@inline function dynamic_element_residual!(resid, dx, x, indices, force_scaling, assembly, ielem,  
    prescribed_conditions, distributed_loads, gravity, x0, v0, ω0, a0, α0)

    properties = dynamic_element_properties(dx, x, indices, force_scaling, assembly, ielem,
        prescribed_conditions, gravity, x0, v0, ω0, a0, α0)

    compatability = compatability_residuals(properties)

    resultants = dynamic_element_resultants(properties, distributed_loads, ielem)

    insert_element_residuals!(resid, indices, force_scaling, assembly, ielem, 
        compatability, resultants)

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
    steady_state_element_jacobian_properties(properties, x, indices, force_scaling, 
        assembly, ielem, prescribed_conditions, gravity, x0, v0, ω0, a0, α0))

Calculate/extract the element properties needed to calculate the jacobian entries 
corresponding to a element for a steady state analysis
"""
@inline function steady_state_element_jacobian_properties(properties, x, indices, force_scaling, 
    assembly, ielem, prescribed_conditions, gravity, x0, v0, ω0, a0, α0)

    properties = static_element_jacobian_properties(properties, x, indices, force_scaling, 
        assembly, ielem, prescribed_conditions, gravity)

    @unpack Cab, CtCab, mass11, mass12, mass21, mass22, V, Ω, C_θ1, C_θ2, C_θ3 = properties

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

    return (; properties..., P_θ, P_V, P_Ω, H_θ, H_V, H_Ω, a_u)
end

"""
    initial_condition_element_jacobian_properties(properties, x, indices, force_scaling, 
        assembly, ielem, prescribed_conditions, gravity)

Calculate/extract the element properties needed to calculate the jacobian entries 
corresponding to a element for a Newmark scheme time marching analysis
"""
@inline function initial_condition_element_jacobian_properties(properties, x, indices, force_scaling, 
    assembly, ielem, prescribed_conditions, gravity)

    @unpack CtCab, CtCabdot, S11, S12, S21, S22, mass11, mass12, mass21, mass22, V, Ω = properties

    # strain and curvature
    γ_F, γ_M, κ_F, κ_M = S11, S12, S21, S22

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
        CtCab*mass11*tilde(CtCab'*V) + CtCab*mass12*tilde(CtCab'*Ω) +
        -CtCab*tilde(mass11*CtCab'*V) - CtCab*tilde(mass12*CtCab'*Ω)
    
    Pdot_Vdot = CtCab*mass11*CtCab'
    
    Pdot_Ωdot = CtCab*mass12*CtCab'

    Hdot_V = CtCab*mass21*CtCabdot' + CtCabdot*mass21*CtCab'

    Hdot_Ω = CtCab*mass22*CtCabdot' + CtCabdot*mass22*CtCab' +
        CtCab*mass21*tilde(CtCab'*V) + CtCab*mass22*tilde(CtCab'*Ω) +
        -CtCab*tilde(mass21*CtCab'*V) - CtCab*tilde(mass22*CtCab'*Ω)

    Hdot_Vdot = CtCab*mass21*CtCab'

    Hdot_Ωdot = CtCab*mass22*CtCab'

    return (; properties..., γ_F, γ_M, κ_F, κ_M, V1dot_V1dot, V2dot_V2dot,
        Ω1dot_Ω1dot, Ω2dot_Ω2dot, P_V, P_Ω, H_V, H_Ω, Pdot_V, Pdot_Ω, Pdot_Vdot, Pdot_Ωdot,
        Hdot_V, Hdot_Ω, Hdot_Vdot, Hdot_Ωdot) 
end

"""
    newmark_element_jacobian_properties(properties, x, indices, force_scaling, 
        assembly, ielem, prescribed_conditions, gravity, x0, v0, ω0, a0, α0,
        Vdot_init, Ωdot_init, dt)

Calculate/extract the element properties needed to calculate the jacobian entries 
corresponding to a element for a Newmark scheme time marching analysis
"""
@inline function newmark_element_jacobian_properties(properties, x, indices, force_scaling, 
    assembly, ielem, prescribed_conditions, gravity, x0, v0, ω0, a0, α0,
    Vdot_init, Ωdot_init, dt)

    properties = steady_state_element_jacobian_properties(properties, x, indices, force_scaling, 
        assembly, ielem, prescribed_conditions, gravity, x0, v0, ω0, a0, α0)

    @unpack Cab, CtCab, CtCabdot, mass11, mass12, mass21, mass22, V, Ω, ω, Vdot, Ωdot, 
        C_θ1, C_θ2, C_θ3 = properties
 
    # linear and angular momentum rates

    Pdot_θ = mul3(C_θ1', C_θ2', C_θ3', Cab*mass11*CtCab'*Vdot) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass12*CtCab'*Ωdot) +
        CtCab*mass11*Cab'*mul3(C_θ1, C_θ2, C_θ3, Vdot) + 
        CtCab*mass12*Cab'*mul3(C_θ1, C_θ2, C_θ3, Ωdot) +
        mul3(C_θ1', C_θ2', C_θ3', Cab*(tilde(Ω - ω)*mass11*CtCab'*V)) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*(tilde(Ω - ω)*mass12*CtCab'*Ω)) +
        CtCabdot*mass11*Cab'*mul3(C_θ1, C_θ2, C_θ3, V) + 
        CtCabdot*mass12*Cab'*mul3(C_θ1, C_θ2, C_θ3, Ω) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass11*CtCabdot'*V) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass12*CtCabdot'*Ω) +
        -CtCab*mass11*tilde(Ω - ω)*Cab'*mul3(C_θ1, C_θ2, C_θ3, V) + 
        -CtCab*mass12*tilde(Ω - ω)*Cab'*mul3(C_θ1, C_θ2, C_θ3, Ω)

    Pdot_V = CtCab*mass11*CtCabdot' + CtCabdot*mass11*CtCab'
    
    Pdot_Ω = CtCab*mass12*CtCabdot' + CtCabdot*mass12*CtCab' +
        CtCab*mass11*tilde(CtCab'*V) + CtCab*mass12*tilde(CtCab'*Ω) +
        -CtCab*tilde(mass11*CtCab'*V) - CtCab*tilde(mass12*CtCab'*Ω)

    Pdot_V += 2/dt*CtCab*mass11*CtCab'
    
    Pdot_Ω += 2/dt*CtCab*mass12*CtCab'

    Hdot_θ = mul3(C_θ1', C_θ2', C_θ3', Cab*mass21*CtCab'*Vdot) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass22*CtCab'*Ωdot) +
        CtCab*mass21*Cab'*mul3(C_θ1, C_θ2, C_θ3, Vdot) + 
        CtCab*mass22*Cab'*mul3(C_θ1, C_θ2, C_θ3, Ωdot) +
        mul3(C_θ1', C_θ2', C_θ3', Cab*(tilde(Ω - ω)*mass21*CtCab'*V)) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*(tilde(Ω - ω)*mass22*CtCab'*Ω)) +
        CtCabdot*mass21*Cab'*mul3(C_θ1, C_θ2, C_θ3, V) + 
        CtCabdot*mass22*Cab'*mul3(C_θ1, C_θ2, C_θ3, Ω) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass21*CtCabdot'*V) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass22*CtCabdot'*Ω) +
        -CtCab*mass21*tilde(Ω - ω)*Cab'*mul3(C_θ1, C_θ2, C_θ3, V) + 
        -CtCab*mass22*tilde(Ω - ω)*Cab'*mul3(C_θ1, C_θ2, C_θ3, Ω)

    Hdot_V = CtCab*mass21*CtCabdot' + CtCabdot*mass21*CtCab'

    Hdot_Ω = CtCab*mass22*CtCabdot' + CtCabdot*mass22*CtCab' +
        CtCab*mass21*tilde(CtCab'*V) + CtCab*mass22*tilde(CtCab'*Ω) +
        -CtCab*tilde(mass21*CtCab'*V) - CtCab*tilde(mass22*CtCab'*Ω)

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
    assembly, ielem, prescribed_conditions, gravity, x0, v0, ω0, a0, α0)

    properties = steady_state_element_jacobian_properties(properties, x, indices, force_scaling, 
        assembly, ielem, prescribed_conditions, gravity, x0, v0, ω0, a0, α0)

    @unpack Cab, CtCab, CtCabdot, mass11, mass12, mass21, mass22, V, Ω, Vdot, Ωdot, ω, C_θ1, C_θ2, C_θ3 = properties

    # linear and angular momentum rates

    Pdot_θ = mul3(C_θ1', C_θ2', C_θ3', Cab*mass11*CtCab'*Vdot) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass12*CtCab'*Ωdot) +
        CtCab*mass11*Cab'*mul3(C_θ1, C_θ2, C_θ3, Vdot) + 
        CtCab*mass12*Cab'*mul3(C_θ1, C_θ2, C_θ3, Ωdot) +
        mul3(C_θ1', C_θ2', C_θ3', Cab*(tilde(Ω - ω)*mass11*CtCab'*V)) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*(tilde(Ω - ω)*mass12*CtCab'*Ω)) +
        CtCabdot*mass11*Cab'*mul3(C_θ1, C_θ2, C_θ3, V) + 
        CtCabdot*mass12*Cab'*mul3(C_θ1, C_θ2, C_θ3, Ω) +
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass11*CtCabdot'*V) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass12*CtCabdot'*Ω) +
        -CtCab*mass11*tilde(Ω - ω)*Cab'*mul3(C_θ1, C_θ2, C_θ3, V) + 
        -CtCab*mass12*tilde(Ω - ω)*Cab'*mul3(C_θ1, C_θ2, C_θ3, Ω)

    Pdot_V = CtCab*mass11*CtCabdot' + CtCabdot*mass11*CtCab'
    
    Pdot_Ω = CtCab*mass12*CtCabdot' + CtCabdot*mass12*CtCab' +
        CtCab*mass11*tilde(CtCab'*V) + CtCab*mass12*tilde(CtCab'*Ω) +
        -CtCab*tilde(mass11*CtCab'*V) - CtCab*tilde(mass12*CtCab'*Ω)

    Hdot_θ = mul3(C_θ1', C_θ2', C_θ3', Cab*mass21*CtCab'*Vdot) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass22*CtCab'*Ωdot) +
        CtCab*mass21*Cab'*mul3(C_θ1, C_θ2, C_θ3, Vdot) + 
        CtCab*mass22*Cab'*mul3(C_θ1, C_θ2, C_θ3, Ωdot) +
        mul3(C_θ1', C_θ2', C_θ3', Cab*(tilde(Ω - ω)*mass21*CtCab'*V)) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*(tilde(Ω - ω)*mass22*CtCab'*Ω)) +
        CtCabdot*mass21*Cab'*mul3(C_θ1, C_θ2, C_θ3, V) + 
        CtCabdot*mass22*Cab'*mul3(C_θ1, C_θ2, C_θ3, Ω) +
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass21*CtCabdot'*V) + 
        mul3(C_θ1', C_θ2', C_θ3', Cab*mass22*CtCabdot'*Ω) +
        -CtCab*mass21*tilde(Ω - ω)*Cab'*mul3(C_θ1, C_θ2, C_θ3, V) + 
        -CtCab*mass22*tilde(Ω - ω)*Cab'*mul3(C_θ1, C_θ2, C_θ3, Ω)

    Hdot_V = CtCab*mass21*CtCabdot' + CtCabdot*mass21*CtCab'

    Hdot_Ω = CtCab*mass22*CtCabdot' + CtCabdot*mass22*CtCab' +
        CtCab*mass21*tilde(CtCab'*V) + CtCab*mass22*tilde(CtCab'*Ω) +
        -CtCab*tilde(mass21*CtCab'*V) - CtCab*tilde(mass22*CtCab'*Ω)

    return (; properties..., Pdot_θ, Pdot_V, Pdot_Ω, Hdot_θ, Hdot_V, Hdot_Ω) 
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

@inline function compatability_jacobians(properties)
   
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
   
    @unpack Cab, CtCab, Qinv, γ_F, γ_M, κ_F, κ_M = properties

    Δu_F = CtCab*γ_F
    ru_F = -Δu_F

    Δu_M = CtCab*γ_M
    ru_M = -Δu_M

    Δθ_F = Qinv*Cab*κ_F
    rθ_F = -Δθ_F

    Δθ_M = Qinv*Cab*κ_M
    rθ_M = -Δθ_M

    return (; ru_F, ru_M, rθ_F, rθ_M)
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

    @unpack L, CtCab, mass11, mass12, mass21, mass22, V, Ω, P, H, ω, 
        a_u, P_θ, P_V, P_Ω, H_θ, H_V, H_Ω, u1_u1, u2_u2, θ1_θ1, θ2_θ2 = properties

    @unpack F1_θ1, F1_θ2, F2_θ1, F2_θ2, M1_θ1, M1_θ2, M2_θ1, M2_θ2 = jacobians

    # add loads due to linear and angular acceleration (including gravity)
    tmp = 1/2*CtCab*mass11*CtCab'*a_u
    
    F1_u1 = -1/2*tmp*u1_u1
    F2_u1 = 1/2*tmp*u1_u1

    F1_u2 = -1/2*tmp*u2_u2    
    F2_u2 = 1/2*tmp*u2_u2

    tmp = 1/2*CtCab*mass21*CtCab'*a_u
    
    M1_u1 = -1/2*tmp*u1_u1
    M2_u1 = 1/2*tmp*u1_u1

    M1_u2 = -1/2*tmp*u2_u2
    M2_u2 = 1/2*tmp*u2_u2

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
    
    M1_V1 = -1/2*tmp
    M2_V1 = 1/2*tmp

    M1_V2 = -1/2*tmp
    M2_V2 = 1/2*tmp

    tmp = 1/2*(tilde(ω)*H_Ω + tilde(V)*P_Ω)
    
    M1_Ω1 = -1/2*tmp
    M2_Ω1 = 1/2*tmp

    M1_Ω2 = -1/2*tmp
    M2_Ω2 = 1/2*tmp

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

    @unpack L, CtCab, F, M, γ, κ, V, Ω, P, H, v, ω, γ_F, γ_M, P_V, P_Ω, H_V, H_Ω, 
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
    
    M1_Ω1 = -1/2*tmp
    M2_Ω1 = 1/2*tmp

    M1_Ω2 = -1/2*tmp
    M2_Ω2 = 1/2*tmp

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
    jacob[irow2:irow2+2, icol2+3:icol2+5] .= F2_θ2 ./ force_scaling

    jacob[irow2:irow2+2, icol:icol+2] .= F2_F

    @views jacob[irow2+3:irow2+5, icol1+3:icol1+5] .+= M2_θ1./ force_scaling
    @views jacob[irow2+3:irow2+5, icol2+3:icol2+5] .+= M2_θ2 ./ force_scaling

    jacob[irow2+3:irow2+5, icol:icol+2] .= M2_F
    jacob[irow2+3:irow2+5, icol+3:icol+5] .= M2_M

    return jacob
end

@inline function insert_initial_condition_element_jacobians!(jacob, indices, force_scaling, 
    assembly, ielem, compatability, resultants)

    @unpack ru_F, ru_M, rθ_F, rθ_M = compatability

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

    @unpack F1_u1, F1_u2, F1_V1, F1_V2, F1_Ω1, F1_Ω2,
            F2_u1, F2_u2, F2_V1, F2_V2, F2_Ω1, F2_Ω2,
            M1_u1, M1_u2, M1_V1, M1_V2, M1_Ω1, M1_Ω2,
            M2_u1, M2_u2, M2_V1, M2_V2, M2_Ω1, M2_Ω2 = resultants

    icol1 = indices.icol_point[assembly.start[ielem]]
    icol2 = indices.icol_point[assembly.stop[ielem]]

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

@inline function insert_element_mass_matrix_jacobians!(jacob, gamma, indices, force_scaling,
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

    compatability = compatability_jacobians(properties)

    resultants = static_element_resultant_jacobians(properties, distributed_loads, ielem)

    insert_static_element_jacobians!(jacob, indices, force_scaling, assembly, ielem, 
        compatability, resultants)

    return jacob
end

"""
    steady_state_element_jacobian!(jacob, x, indices, force_scaling, 
        assembly, ielem, prescribed_conditions, distributed_loads, gravity, x0, v0, ω0, a0, α0)

Calculate and insert the jacobian entries corresponding to a beam element for a steady state 
analysis into the system jacobian matrix.
"""
@inline function steady_state_element_jacobian!(jacob, x, indices, force_scaling, 
    assembly, ielem, prescribed_conditions, distributed_loads, gravity, x0, v0, ω0, a0, α0)

    properties = steady_state_element_properties(x, indices, force_scaling, assembly, ielem, 
        prescribed_conditions, gravity, x0, v0, ω0, a0, α0)

    properties = steady_state_element_jacobian_properties(properties, x, indices, force_scaling, 
        assembly, ielem, prescribed_conditions, gravity, x0, v0, ω0, a0, α0)

    compatability = compatability_jacobians(properties)

    resultants = steady_state_element_resultant_jacobians(properties, distributed_loads, ielem)

    insert_dynamic_element_jacobians!(jacob, indices, force_scaling, assembly, ielem, 
        compatability, resultants)

    return jacob
end

"""
    initial_condition_element_jacobian!(jacob, x, indices, force_scaling, 
        assembly, ielem, prescribed_conditions, distributed_loads, gravity, x0, v0, ω0, a0, α0,
        u0, θ0, udot0, θdot0)

Calculate and insert the jacobian entries corresponding to a beam element for the 
initialization of a time domain analysis into the system jacobian matrix.
"""
@inline function initial_condition_element_jacobian!(jacob, x, indices, force_scaling, 
    assembly, ielem, prescribed_conditions, distributed_loads, gravity, x0, v0, ω0, a0, α0,
    u0, θ0, udot0, θdot0)

    properties = initial_condition_element_properties(x, indices, force_scaling, assembly, ielem, 
        prescribed_conditions, gravity, x0, v0, ω0, a0, α0, u0, θ0, udot0, θdot0)

    properties = initial_condition_element_jacobian_properties(properties, x, indices, force_scaling, 
        assembly, ielem, prescribed_conditions, gravity)

    compatability = initial_condition_compatability_jacobians(properties)

    resultants = initial_condition_element_resultant_jacobians(properties)

    insert_initial_condition_element_jacobians!(jacob, indices, force_scaling, assembly, ielem, 
        compatability, resultants)

    return jacob
end

"""
    newmark_element_jacobian!(jacob, x, indices, force_scaling, 
        assembly, ielem, prescribed_conditions, distributed_loads, gravity, x0, v0, ω0, a0, α0,
        Vdot_init, Ωdot_init, dt)

Calculate and insert the jacobian entries corresponding to a beam element for a Newmark-scheme
time marching analysis into the system jacobian matrix.
"""
@inline function newmark_element_jacobian!(jacob, x, indices, force_scaling, 
    assembly, ielem, prescribed_conditions, distributed_loads, gravity, x0, v0, ω0, a0, α0,
    Vdot_init, Ωdot_init, dt)

    properties = newmark_element_properties(x, indices, force_scaling, assembly, ielem, 
        prescribed_conditions, gravity, x0, v0, ω0, a0, α0, 
        Vdot_init, Ωdot_init, dt)

    properties = newmark_element_jacobian_properties(properties, x, indices, force_scaling, 
        assembly, ielem, prescribed_conditions, gravity, x0, v0, ω0, a0, α0,
        Vdot_init, Ωdot_init, dt)

    compatability = compatability_jacobians(properties)

    resultants = newmark_element_resultant_jacobians(properties, distributed_loads, ielem)

    insert_dynamic_element_jacobians!(jacob, indices, force_scaling, assembly, ielem, 
        compatability, resultants)

    return jacob
end

"""
    dynamic_element_jacobian!(jacob, dx, x, indices, force_scaling, 
        assembly, ielem, prescribed_conditions, distributed_loads, gravity, 
        x0, v0, ω0, a0, α0)

Calculate and insert the jacobian entries corresponding to a beam element for a dynamic
analysis into the system jacobian matrix.
"""
@inline function dynamic_element_jacobian!(jacob, dx, x, indices, force_scaling, 
    assembly, ielem, prescribed_conditions, distributed_loads, gravity, x0, v0, ω0, a0, α0)

    properties = dynamic_element_properties(dx, x, indices, force_scaling, assembly, ielem, 
        prescribed_conditions, gravity, x0, v0, ω0, a0, α0)

    properties = dynamic_element_jacobian_properties(properties, dx, x, indices, force_scaling, 
        assembly, ielem, prescribed_conditions, gravity, x0, v0, ω0, a0, α0)

    compatability = compatability_jacobians(properties)

    resultants = dynamic_element_resultant_jacobians(properties, distributed_loads, ielem)

    insert_dynamic_element_jacobians!(jacob, indices, force_scaling, assembly, ielem, 
        compatability, resultants)

    return jacob
end

"""
    element_mass_matrix!(jacob, gamma, x, indices, force_scaling, assembly, 
        ielem, prescribed_conditions)

Calculate and insert the mass_matrix jacobian entries corresponding to a beam element into 
the system jacobian matrix.
"""
@inline function element_mass_matrix!(jacob, gamma, x, indices, force_scaling, assembly, 
    ielem, prescribed_conditions)

    properties = mass_matrix_element_jacobian_properties(x, indices, force_scaling, 
        assembly, ielem, prescribed_conditions)

    resultants = mass_matrix_element_resultant_jacobians(properties)
    
    insert_element_mass_matrix_jacobians!(jacob, gamma, indices, force_scaling, assembly, ielem, 
        resultants)

    return jacob
end