# --- Point Properties --- #

"""
    newmark_point_jacobian_properties(properties, x, indices, force_scaling, 
        assembly, ipoint, prescribed_conditions, point_masses, gravity, 
        udot_init, θdot_init, Vdot_init, Ωdot_init, dt)

Calculate/extract the point properties needed to calculate the jacobian entries 
corresponding to a point for a Newmark scheme time marching analysis
"""
@inline function newmark_point_jacobian_properties(properties, x, indices,  
    force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity, 
    ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb, udot_init, θdot_init, Vdot_init, Ωdot_init, dt)

    properties = steady_state_point_jacobian_properties(properties, x, indices, 
        force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity,
        ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb)

    @unpack C, Cdot, mass11, mass12, mass21, mass22, V, Ω, Vdot, Ωdot, v, ω, C_θ1, C_θ2, C_θ3 = properties

    # linear and angular displacement rates
    udot_u = 2/dt*I3
    θdot_θ = 2/dt*I3

    # linear and angular momentum rates
    Pdot_ωb = -C'*mass11*C*tilde(V) - C'*mass12*C*tilde(Ω) + 
        tilde(C'*mass11*C*V) + tilde(C'*mass12*C*Ω)

    Pdot_ωb *= ωb_ωb

    Pdot_θ = mul3(C_θ1', C_θ2', C_θ3', mass11*C*Vdot) + 
        mul3(C_θ1', C_θ2', C_θ3', mass12*C*Ωdot) +
        C'*mass11*mul3(C_θ1, C_θ2, C_θ3, Vdot) + 
        C'*mass12*mul3(C_θ1, C_θ2, C_θ3, Ωdot) +
        tilde(Ω - ω)*mul3(C_θ1', C_θ2', C_θ3', mass11*C*V) + 
        tilde(Ω - ω)*mul3(C_θ1', C_θ2', C_θ3', mass12*C*Ω) +
        Cdot'*mass11*mul3(C_θ1, C_θ2, C_θ3, V) + 
        Cdot'*mass12*mul3(C_θ1, C_θ2, C_θ3, Ω) +
        mul3(C_θ1', C_θ2', C_θ3', mass11*Cdot*V) + 
        mul3(C_θ1', C_θ2', C_θ3', mass12*Cdot*Ω) +
        -C'*mass11*mul3(C_θ1, C_θ2, C_θ3, tilde(Ω - ω)*V) + 
        -C'*mass12*mul3(C_θ1, C_θ2, C_θ3, tilde(Ω - ω)*Ω)

    Pdot_V = C'*mass11*Cdot + Cdot'*mass11*C

    Pdot_Ω = C'*mass12*Cdot + Cdot'*mass12*C +
        C'*mass11*C*tilde(V) + C'*mass12*C*tilde(Ω) + 
        -tilde(C'*mass11*C*V) - tilde(C'*mass12*C*Ω)

    Pdot_V += 2/dt*C'*mass11*C
    
    Pdot_Ω += 2/dt*C'*mass12*C

    Hdot_ωb = -C'*mass21*C*tilde(V) - C'*mass22*C*tilde(Ω) +
        tilde(C'*mass21*C*V) + tilde(C'*mass22*C*Ω)

    Hdot_ωb *= ωb_ωb

    Hdot_θ = mul3(C_θ1', C_θ2', C_θ3', mass21*C*Vdot) + 
        mul3(C_θ1', C_θ2', C_θ3', mass22*C*Ωdot) +
        C'*mass21*mul3(C_θ1, C_θ2, C_θ3, Vdot) + 
        C'*mass22*mul3(C_θ1, C_θ2, C_θ3, Ωdot) +
        tilde(Ω - ω)*mul3(C_θ1', C_θ2', C_θ3', mass21*C*V) + 
        tilde(Ω - ω)*mul3(C_θ1', C_θ2', C_θ3', mass22*C*Ω) +
        Cdot'*mass21*mul3(C_θ1, C_θ2, C_θ3, V) + 
        Cdot'*mass22*mul3(C_θ1, C_θ2, C_θ3, Ω) +
        mul3(C_θ1', C_θ2', C_θ3', mass21*Cdot*V) + 
        mul3(C_θ1', C_θ2', C_θ3', mass22*Cdot*Ω) +
        -C'*mass21*mul3(C_θ1, C_θ2, C_θ3, tilde(Ω - ω)*V) + 
        -C'*mass22*mul3(C_θ1, C_θ2, C_θ3, tilde(Ω - ω)*Ω)

    Hdot_V = C'*mass21*Cdot + Cdot'*mass21*C

    Hdot_Ω = C'*mass22*Cdot + Cdot'*mass22*C +
        C'*mass21*C*tilde(V) + C'*mass22*C*tilde(Ω) +
        -tilde(C'*mass21*C*V) - tilde(C'*mass22*C*Ω)

    Hdot_V += 2/dt*C'*mass21*C

    Hdot_Ω += 2/dt*C'*mass22*C

    # overwrite acceleration terms so we don't double count them
    a_u = Z3
    a_ab = Z3
    a_αb = Z3
    α_αb = Z3

    return (; properties..., a_u, a_ab, a_αb, α_αb, udot_u, θdot_θ, 
        Pdot_ωb, Pdot_θ, Pdot_V, Pdot_Ω, 
        Hdot_ωb, Hdot_θ, Hdot_V, Hdot_Ω) 
end

# --- Point Resultants --- #

"""
    newmark_point_resultant_jacobians(properties)

Calculate the jacobians for the net loads `F` and `M` applied at a point for a Newmark 
scheme time marching analysis.
"""
@inline function newmark_point_resultant_jacobians(properties)

    jacobians = steady_state_point_resultant_jacobians(properties)

    @unpack F_ωb, F_θ, F_V, F_Ω, M_ωb, M_θ, M_V, M_Ω = jacobians

    @unpack θ_θ, Pdot_ωb, Pdot_θ, Pdot_V, Pdot_Ω, Hdot_ωb, Hdot_θ, Hdot_V, Hdot_Ω = properties

    # add loads due to linear and angular momentum rates
    F_ωb -= Pdot_ωb
    F_θ -= Pdot_θ*θ_θ
    F_V -= Pdot_V
    F_Ω -= Pdot_Ω
    
    M_ωb -= Hdot_ωb
    M_θ -= Hdot_θ*θ_θ
    M_V -= Hdot_V
    M_Ω -= Hdot_Ω

    return (; jacobians..., F_ωb, F_θ, F_V, F_Ω, M_ωb, M_θ, M_V, M_Ω)
end

# --- Velocity Residuals --- #

"""
    newmark_point_velocity_jacobians(properties)

Calculate the jacobians of the velocity residuals `rV` and `rΩ` of a point for a Newmark 
scheme time-marching analysis.
"""
@inline function newmark_point_velocity_jacobians(properties)

    jacobians = steady_state_point_velocity_jacobians(properties)

    @unpack u_u, θ_θ, udot_u, θdot_θ = properties

    @unpack rV_u, rΩ_θ = jacobians

    rV_u -= udot_u*u_u

    rΩ_θ -= θdot_θ*θ_θ

    return (; jacobians..., rV_u, rΩ_θ)
end

# --- Point Residual --- #

"""
    newmark_point_jacobian!(jacob, x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity, ub, θb, vb, ωb, ab, αb, 
        udot_init, θdot_init, Vdot_init, Ωdot_init, dt))

Calculate and insert the jacobian entries corresponding to a point for a Newmark scheme
time marching analysis into the system jacobian matrix.
"""
@inline function newmark_point_jacobian!(jacob, x, indices, force_scaling, 
    assembly, ipoint, prescribed_conditions, point_masses, gravity, ub, θb, vb, ωb, ab, αb,
    ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb, udot_init, θdot_init, Vdot_init, Ωdot_init, dt)

    properties = newmark_point_properties(x, indices, force_scaling, 
        assembly, ipoint, prescribed_conditions, point_masses, gravity, 
        ub, θb, vb, ωb, ab, αb, udot_init, θdot_init, Vdot_init, Ωdot_init, dt)

    properties = newmark_point_jacobian_properties(properties, x, indices, 
        force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity, 
        ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb, udot_init, θdot_init, Vdot_init, Ωdot_init, dt)

    resultants = newmark_point_resultant_jacobians(properties)

    velocities = newmark_point_velocity_jacobians(properties)

    insert_dynamic_point_jacobians!(jacob, indices, force_scaling, ipoint, resultants, velocities)

    return jacob
end