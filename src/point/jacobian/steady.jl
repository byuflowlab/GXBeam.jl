# --- Point Properties --- #

"""
    steady_state_point_jacobian_properties(properties, x, indices, force_scaling, 
        assembly, ipoint, prescribed_conditions, point_masses, gravity, 
        ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb)

Calculate/extract the point properties needed to calculate the jacobian entries 
corresponding to a point for a steady state analysis
"""
@inline function steady_state_point_jacobian_properties(properties, x, indices, 
    force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity,
    ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb)

    properties = static_point_jacobian_properties(properties, x, indices, force_scaling, 
        assembly, ipoint, prescribed_conditions, point_masses, gravity)

    @unpack C, mass11, mass12, mass21, mass22, u, θ, V, Ω, θb, ωb, αb, Δx, C_θ1, C_θ2, C_θ3 = properties

    # rotation parameter matrices
    Qinv_θ1, Qinv_θ2, Qinv_θ3 = get_Qinv_θ(θ)

    # linear and angular velocity
    v_u = tilde(ωb)
    v_vb = vb_vb
    v_ωb = (-tilde(Δx) - tilde(u))*ωb_ωb
    ω_ωb = ωb_ωb

    # linear and angular acceleration
    a_u = tilde(αb)
    a_ab = ab_ab
    a_αb = (-tilde(Δx) - tilde(u))*αb_αb
    α_αb = αb_αb

    # add gravitational acceleration
    C_θb1, C_θb2, C_θb3 = get_C_θ(θb)
    a_θb = -mul3(C_θb1, C_θb2, C_θb3, gravity)*θb_θb

    # linear and angular momentum
    P_θ = mul3(C_θ1', C_θ2', C_θ3', (mass11*C*V + mass12*C*Ω)) + 
        C'*mass11*mul3(C_θ1, C_θ2, C_θ3, V) + 
        C'*mass12*mul3(C_θ1, C_θ2, C_θ3, Ω)
    P_V = C'*mass11*C
    P_Ω = C'*mass12*C

    H_θ = mul3(C_θ1', C_θ2', C_θ3', (mass21*C*V + mass22*C*Ω)) + 
        C'*mass21*mul3(C_θ1, C_θ2, C_θ3, V) + 
        C'*mass22*mul3(C_θ1, C_θ2, C_θ3, Ω)
    H_V = C'*mass21*C
    H_Ω = C'*mass22*C

    return (; properties..., Qinv_θ1, Qinv_θ2, Qinv_θ3, v_u, v_vb, v_ωb, ω_ωb, 
        a_u, a_θb, a_ab, a_αb, α_αb, P_θ, P_V, P_Ω, H_θ, H_V, H_Ω) 
end

# --- Point Resultants --- #

"""
    steady_state_point_resultant_jacobians(properties)

Calculate the jacobians for the net loads `F` and `M` applied at a point for a steady 
state analysis.
"""
@inline function steady_state_point_resultant_jacobians(properties)

    jacobians = static_point_resultant_jacobians(properties)

    @unpack C, mass11, mass12, mass21, mass22, V, Ω, P, H, ω, C_θ1, C_θ2, C_θ3, 
        u_u, θ_θ, P_θ, P_V, P_Ω, H_θ, H_V, H_Ω, ω_ωb, a_u, a_θb, a_ab, a_αb, α_αb = properties

    @unpack F_θ, M_θ = jacobians

    # add loads due to linear and angular acceleration (including gravity)
    F_θb = -C'*mass11*C*a_θb
    F_ab = -C'*mass11*C*a_ab
    F_αb = -C'*mass11*C*a_αb - C'*mass12*C*α_αb
    F_u = -C'*mass11*C*a_u*u_u

    M_θb = -C'*mass21*C*a_θb
    M_ab = -C'*mass21*C*a_ab
    M_αb = -C'*mass21*C*a_αb - C'*mass22*C*α_αb
    M_u = -C'*mass21*C*a_u*u_u

    # add loads due to linear and angular momentum
    F_ωb = tilde(P)*ω_ωb
    F_θ -= tilde(ω)*P_θ*θ_θ
    F_V = -tilde(ω)*P_V
    F_Ω = -tilde(ω)*P_Ω

    M_ωb = tilde(H)*ω_ωb
    M_θ -= (tilde(ω)*H_θ + tilde(V)*P_θ)*θ_θ
    M_V = -tilde(ω)*H_V - tilde(V)*P_V + tilde(P)
    M_Ω = -tilde(ω)*H_Ω - tilde(V)*P_Ω

    return (; jacobians..., 
        F_θb, F_ωb, F_ab, F_αb, F_u, F_θ, F_V, F_Ω, 
        M_θb, M_ωb, M_ab, M_αb, M_u, M_θ, M_V, M_Ω)
end

# --- Velocity Residuals --- #

"""
    steady_state_point_velocity_jacobians(properties)

Calculate the jacobians of the velocity residuals `rV` and `rΩ` of a point for a steady 
state analysis.
"""
@inline function steady_state_point_velocity_jacobians(properties)

    @unpack C, Qinv, Ω, ω, v_vb, v_ωb, v_u, ω_ωb, u_u, θ_θ, C_θ1, C_θ2, C_θ3, Qinv_θ1, Qinv_θ2, Qinv_θ3 = properties
    
    rV_vb = -v_vb
    rV_ωb = -v_ωb
    rV_u = -v_u*u_u
    rV_V = I3

    rΩ_ωb = -Qinv*C*ω_ωb
    rΩ_θ = (mul3(Qinv_θ1, Qinv_θ2, Qinv_θ3, C*(Ω - ω)) + Qinv*mul3(C_θ1, C_θ2, C_θ3, Ω - ω))*θ_θ
    rΩ_Ω = Qinv*C

    return (; rV_vb, rV_ωb, rV_u, rV_V, rΩ_ωb, rΩ_θ, rΩ_Ω)
end

# --- Point Residual --- #

"""
    steady_state_point_jacobian!(jacob, x, indices, force_scaling, assembly, 
        ipoint, prescribed_conditions, point_masses, gravity, ub, θb, vb, ωb, ab, αb,
        ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb)

Calculate and insert the jacobian entries corresponding to a point for a steady state 
analysis into the system jacobian matrix.
"""
@inline function steady_state_point_jacobian!(jacob, x, indices, force_scaling, 
    assembly, ipoint, prescribed_conditions, point_masses, gravity, ub, θb, vb, ωb, ab, αb,
    ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb)

    properties = steady_state_point_properties(x, indices, force_scaling, assembly, ipoint, 
        prescribed_conditions, point_masses, gravity, ub, θb, vb, ωb, ab, αb)

    properties = steady_state_point_jacobian_properties(properties, x, indices,
        force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity, 
        ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb)

    resultants = steady_state_point_resultant_jacobians(properties)

    velocities = steady_state_point_velocity_jacobians(properties)

    insert_dynamic_point_jacobians!(jacob, indices, force_scaling, ipoint, resultants, velocities)

    return jacob
end