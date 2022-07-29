# --- Point Properties --- #

"""
    expanded_steady_point_jacobian_properties(properties, x, indices, force_scaling, 
        assembly, ipoint, prescribed_conditions, point_masses, gravity,
        ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb)

Calculate/extract the point properties needed to calculate the jacobian entries 
corresponding to a point for a constant mass matrix system
"""
@inline function expanded_steady_point_jacobian_properties(properties, x, indices, 
    force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity,
    ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb)

    @unpack C, mass11, mass12, mass21, mass22, u, θ, V, Ω, θb, ωb, αb, Δx = properties

    # forces and moments
    F_θ, F_F, M_θ, M_M = point_load_jacobians(x, ipoint, indices.icol_point, force_scaling, prescribed_conditions)

    # linear and angular displacement
    u_u, θ_θ = point_displacement_jacobians(ipoint, prescribed_conditions)

    # rotation parameter matrices
    C_θ1, C_θ2, C_θ3 = get_C_θ(C, θ)
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
    a_θb = -mul3(get_C_θ(θb)..., gravity)*θb_θb

    # linear and angular momentum
    P_V = mass11
    P_Ω = mass12
    H_V = mass21
    H_Ω = mass22

    return (; properties..., C_θ1, C_θ2, C_θ3, Qinv_θ1, Qinv_θ2, Qinv_θ3, u_u, θ_θ, 
        F_θ, F_F, M_θ, M_M, P_V, P_Ω, H_V, H_Ω, v_u, v_vb, v_ωb, ω_ωb, a_u, a_θb, a_ab, a_αb, α_αb) 
end

"""
    expanded_dynamic_point_jacobian_properties(properties, x, indices, 
        force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity)

Calculate/extract the point properties needed to calculate the jacobian entries 
corresponding to a point for a dynamic analysis
"""
@inline function expanded_dynamic_point_jacobian_properties(properties, x, indices, 
    force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity,
    ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb)

    properties = expanded_steady_point_jacobian_properties(properties, x, indices, 
        force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity,
        ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb)

    @unpack θb = properties

    # overwrite acceleration terms so we don't double count them
    a_u = @SMatrix zeros(3,3)
    a_θb = -mul3(get_C_θ(θb)..., gravity)*θb_θb
    a_ab = @SMatrix zeros(3,3)
    a_αb = @SMatrix zeros(3,3)
    α_αb = @SMatrix zeros(3,3)

    return (; properties..., a_u, a_θb, a_ab, a_αb, α_αb)
end

"""
    expanded_mass_matrix_point_jacobian_properties(assembly, ipoint, prescribed_conditions, point_masses)

Calculate/extract the point properties needed to calculate the mass matrix jacobian entries 
corresponding to a point for a constant mass matrix system
"""
@inline function expanded_mass_matrix_point_jacobian_properties(assembly, ipoint, prescribed_conditions, point_masses)

    # mass matrix
    mass = haskey(point_masses, ipoint) ? point_masses[ipoint].mass : @SMatrix zeros(6,6)
    
    # mass submatrices
    mass11 = mass[SVector{3}(1:3), SVector{3}(1:3)]
    mass12 = mass[SVector{3}(1:3), SVector{3}(4:6)]
    mass21 = mass[SVector{3}(4:6), SVector{3}(1:3)]
    mass22 = mass[SVector{3}(4:6), SVector{3}(4:6)]

    # linear and angular displacement rates (in the body frame)
    udot_udot, θdot_θdot = point_displacement_jacobians(ipoint, prescribed_conditions)

    # linear and angular momentum rates
    Pdot_Vdot = mass11
    Pdot_Ωdot = mass12
    Hdot_Vdot = mass21
    Hdot_Ωdot = mass22

    return (; udot_udot, θdot_θdot, Pdot_Vdot, Pdot_Ωdot, Hdot_Vdot, Hdot_Ωdot)
end

# --- Point Resultants --- #

"""
    expanded_point_resultant_jacobians(properties)

Calculate the jacobians for the net loads `F` and `M` applied at a point for a constant
mass matrix system
"""
@inline function expanded_point_resultant_jacobians(properties)

    @unpack C, mass11, mass12, mass21, mass22, F, M, V, Ω, P, H, ω, a, α, 
        C_θ1, C_θ2, C_θ3, u_u, θ_θ, P_V, P_Ω, H_V, H_Ω, 
        ω_ωb, a_u, a_θb, a_ab, a_αb, α_αb = properties

    @unpack F_θ, M_θ, F_F, M_M = properties

    # rotate loads into the appropriate frame
    F_θ = (mul3(C_θ1, C_θ2, C_θ3, F) + C*F_θ)*θ_θ
    F_F = C*F_F

    M_θ = (mul3(C_θ1, C_θ2, C_θ3, M) + C*M_θ)*θ_θ
    M_M = C*M_M

    # add loads due to linear and angular acceleration (including gravity)
    F_θb = -mass11*C*a_θb
    F_ab = -mass11*C*a_ab
    F_αb = -mass11*C*a_αb - mass12*C*α_αb
    F_u = -mass11*C*a_u*u_u
    F_θ -= (mass11*mul3(C_θ1, C_θ2, C_θ3, a) + mass12*mul3(C_θ1, C_θ2, C_θ3, α))*θ_θ

    M_θb = -mass21*C*a_θb
    M_ab = -mass21*C*a_ab
    M_αb = -mass21*C*a_αb - mass22*C*α_αb
    M_u = -mass21*C*a_u*u_u
    M_θ -= (mass21*mul3(C_θ1, C_θ2, C_θ3, a) + mass22*mul3(C_θ1, C_θ2, C_θ3, α))*θ_θ

    # add loads due to linear and angular momentum
    F_ωb = @SVector zeros(3)
    F_V = -tilde(Ω)*P_V 
    F_Ω = -tilde(Ω)*P_Ω + tilde(P)

    M_ωb = @SVector zeros(3)
    M_V = -tilde(Ω)*H_V - tilde(V)*P_V + tilde(P)
    M_Ω = -tilde(Ω)*H_Ω - tilde(V)*P_Ω + tilde(H)

    (; F_θb, F_ωb, F_ab, F_αb, F_F, F_u, F_θ, F_V, F_Ω, 
       M_θb, M_ωb, M_ab, M_αb, M_M, M_u, M_θ, M_V, M_Ω)
end

# --- Velocity Residuals --- #

"""
    expanded_point_velocity_jacobians(properties)

Calculate the jacobians of the velocity residuals `rV` and `rΩ` of a point for a constant 
mass matrix system
"""
@inline function expanded_point_velocity_jacobians(properties)

    @unpack C, Qinv, V, Ω, ω, v_u, v_vb, v_ωb, ω_ωb, u_u, θ_θ, C_θ1, C_θ2, C_θ3, Qinv_θ1, Qinv_θ2, Qinv_θ3 = properties
    
    rV_vb = -v_vb
    rV_ωb = -v_ωb
    rV_u = -v_u*u_u
    rV_θ = mul3(C_θ1', C_θ2', C_θ3', V)*θ_θ
    rV_V = C'

    rΩ_ωb = -Qinv*C*ω_ωb
    rΩ_θ = (mul3(Qinv_θ1, Qinv_θ2, Qinv_θ3, Ω - C*ω) - Qinv*mul3(C_θ1, C_θ2, C_θ3, ω))*θ_θ
    rΩ_Ω = Qinv

    return (; rV_vb, rV_ωb, rV_u, rV_θ, rV_V, rΩ_ωb, rΩ_θ, rΩ_Ω)
end

# --- Residual Placement --- #

"""
    insert_expanded_point_jacobians!(jacob, indices, force_scaling, ipoint, resultants, 
        velocities)

Insert the jacobian entries corresponding to a point for a constant mass matrix system into 
the system jacobian matrix for a constant mass matrix system.
"""
@inline function insert_expanded_point_jacobians!(jacob, indices, force_scaling, ipoint,  
    resultants, velocities)

    insert_dynamic_point_jacobians!(jacob, indices, force_scaling, ipoint, resultants, 
        velocities)

    @unpack rV_θ = velocities

    irow = indices.irow_point[ipoint]
    icol = indices.icol_point[ipoint]

    jacob[irow+6:irow+8, icol+3:icol+5] .= rV_θ

    return jacob
end

# --- Point Residual --- #

"""
    expanded_steady_point_jacobian!(jacob, x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity)

Calculate and insert the jacobian entries corresponding to a point for a constant mass 
matrix system into the system jacobian matrix.
"""
@inline function expanded_steady_point_jacobian!(jacob, x, indices, force_scaling, 
    assembly, ipoint, prescribed_conditions, point_masses, gravity,
    ub, θb, vb, ωb, ab, αb, ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb)

    properties = expanded_steady_point_properties(x, indices, force_scaling, 
        assembly, ipoint, prescribed_conditions, point_masses, gravity,
        ub, θb, vb, ωb, ab, αb)

    properties = expanded_steady_point_jacobian_properties(properties, x, indices,
        force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity,
        ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb)

    resultants = expanded_point_resultant_jacobians(properties)

    velocities = expanded_point_velocity_jacobians(properties)

    insert_expanded_point_jacobians!(jacob, indices, force_scaling, ipoint, resultants, velocities)

    return jacob
end

"""
    expanded_dynamic_point_jacobian!(jacob, x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity)

Calculate and insert the jacobian entries corresponding to a point for a constant mass 
matrix system into the system jacobian matrix.
"""
@inline function expanded_dynamic_point_jacobian!(jacob, x, indices, force_scaling, 
    assembly, ipoint, prescribed_conditions, point_masses, gravity,
    ub, θb, vb, ωb, ab, αb, ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb)

    properties = expanded_dynamic_point_properties(x, indices, force_scaling, 
        assembly, ipoint, prescribed_conditions, point_masses, gravity, ub, θb, vb, ωb, ab, αb)

    properties = expanded_dynamic_point_jacobian_properties(properties, x, indices,
        force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity,
        ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb)

    resultants = expanded_point_resultant_jacobians(properties)

    velocities = expanded_point_velocity_jacobians(properties)

    insert_expanded_point_jacobians!(jacob, indices, force_scaling, ipoint, resultants, velocities)

    return jacob
end

"""
    expanded_mass_matrix_point_jacobian!(jacob, gamma, indices, force_scaling, assembly, 
        ipoint, prescribed_conditions, point_masses)

Calculate and insert the mass_matrix jacobian entries corresponding to a point into the 
system jacobian matrix for a constant mass matrix system
"""
@inline function expanded_mass_matrix_point_jacobian!(jacob, gamma, indices, force_scaling, assembly, 
    ipoint, prescribed_conditions, point_masses)

    properties = expanded_mass_matrix_point_jacobian_properties(assembly, ipoint, prescribed_conditions, point_masses)

    resultants = mass_matrix_point_resultant_jacobians(properties)

    velocities = mass_matrix_point_velocity_jacobians(properties)
    
    insert_mass_matrix_point_jacobians!(jacob, gamma, indices, force_scaling, 
        ipoint, resultants, velocities)

    return jacob
end