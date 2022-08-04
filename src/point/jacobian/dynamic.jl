# --- Point Properties --- #

"""
    dynamic_point_jacobian_properties(properties, dx, x, indices, force_scaling, 
        assembly, ipoint, prescribed_conditions, point_masses, gravity, 
        ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb)

Calculate/extract the point properties needed to calculate the jacobian entries 
corresponding to a point for a dynamic analysis
"""
@inline function dynamic_point_jacobian_properties(properties, dx, x, indices,  
    force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity, ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb)

    properties = steady_state_point_jacobian_properties(properties, x, indices, force_scaling, 
        assembly, ipoint, prescribed_conditions, point_masses, gravity, ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb)

    @unpack C, Cdot, mass11, mass12, mass21, mass22, V, Ω, Vdot, Ωdot, v, ω, C_θ1, C_θ2, C_θ3 = properties

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

    # overwrite acceleration terms so we don't double count them
    a_u = Z3
    a_ab = Z3
    a_αb = Z3
    α_αb = Z3

    return (; properties..., a_u, a_ab, a_αb, α_αb, 
        Pdot_ωb, Pdot_θ, Pdot_V, Pdot_Ω, 
        Hdot_ωb, Hdot_θ, Hdot_V, Hdot_Ω) 
end

"""
    mass_matrix_point_jacobian_properties(x, indices, force_scaling, 
        assembly, ipoint, prescribed_conditions, point_masses)

Calculate/extract the point properties needed to calculate the mass matrix jacobian entries 
corresponding to a point
"""
@inline function mass_matrix_point_jacobian_properties(x, indices, force_scaling, 
    assembly, ipoint, prescribed_conditions, point_masses)

    # mass matrix
    mass = haskey(point_masses, ipoint) ? point_masses[ipoint].mass : @SMatrix zeros(6,6)
    
    # mass submatrices
    mass11 = mass[SVector{3}(1:3), SVector{3}(1:3)]
    mass12 = mass[SVector{3}(1:3), SVector{3}(4:6)]
    mass21 = mass[SVector{3}(4:6), SVector{3}(1:3)]
    mass22 = mass[SVector{3}(4:6), SVector{3}(4:6)]

    # linear and angular displacement
    u, θ = point_displacement(x, ipoint, indices.icol_point, prescribed_conditions)

    # linear and angular displacement rates (in the body frame)
    udot_udot, θdot_θdot = point_displacement_jacobians(ipoint, prescribed_conditions)

    # transformation matrices
    C = get_C(θ)

    # linear and angular momentum rates
    Pdot_Vdot = C'*mass11*C
    Pdot_Ωdot = C'*mass12*C
    Hdot_Vdot = C'*mass21*C
    Hdot_Ωdot = C'*mass22*C

    return (; udot_udot, θdot_θdot, Pdot_Vdot, Pdot_Ωdot, Hdot_Vdot, Hdot_Ωdot)
end

# --- Point Resultants --- #

"""
    dynamic_point_resultant_jacobians(properties)

Calculate the jacobians for the net loads `F` and `M` applied at a point for a dynamic 
analysis.
"""
@inline function dynamic_point_resultant_jacobians(properties)

    return newmark_point_resultant_jacobians(properties)
end

"""
    mass_matrix_point_resultant_jacobians(properties)

Calculate the mass matrix jacobians for the net loads `F` and `M` applied at a point
"""
@inline function mass_matrix_point_resultant_jacobians(properties)

    @unpack Pdot_Vdot, Pdot_Ωdot, Hdot_Vdot, Hdot_Ωdot = properties

    # add loads due to linear and angular momentum rates (and gravity)
    F_Vdot = -Pdot_Vdot
    F_Ωdot = -Pdot_Ωdot
    M_Vdot = -Hdot_Vdot
    M_Ωdot = -Hdot_Ωdot

    return (; F_Vdot, F_Ωdot, M_Vdot, M_Ωdot) 
end

# --- Velocity Residuals --- #

"""
    dynamic_point_velocity_jacobians(properties)

Calculate the jacobians of the velocity residuals `rV` and `rΩ` of a point for a dynamic 
analysis.
"""
@inline function dynamic_point_velocity_jacobians(properties)

    return steady_state_point_velocity_jacobians(properties)
end

"""
    mass_matrix_point_velocity_jacobians(properties)

Calculate the mass matrix jacobians of the velocity residuals `rV` and `rΩ` of a point
"""
@inline function mass_matrix_point_velocity_jacobians(properties)

    @unpack udot_udot, θdot_θdot = properties

    rV_udot = -udot_udot
    rΩ_θdot = -θdot_θdot

    return (; rV_udot, rΩ_θdot)
end

# --- Residual Placement --- #

"""
    insert_dynamic_point_jacobians!(jacob, indices, force_scaling, ipoint, resultants, 
        velocities)

Insert the jacobian entries corresponding to a point for a steady state analysis into the 
system jacobian matrix.
"""
@inline function insert_dynamic_point_jacobians!(jacob, indices, force_scaling, ipoint,  
    resultants, velocities)

    insert_static_point_jacobians!(jacob, indices, force_scaling, ipoint, resultants)

    @unpack F_θb, F_ωb, F_ab, F_αb, F_u, F_θ, F_V, F_Ω, 
            M_θb, M_ωb, M_ab, M_αb, M_u, M_θ, M_V, M_Ω = resultants

    @unpack rV_vb, rV_ωb, rV_u, rV_V, rΩ_ωb, rΩ_θ, rΩ_Ω = velocities

    irow = indices.irow_point[ipoint]
    icol = indices.icol_point[ipoint]

    jacob[irow:irow+2, 4:6] .= -F_θb ./ force_scaling
    jacob[irow:irow+2, 10:12] .= -F_ωb ./ force_scaling
    jacob[irow:irow+2, 13:15] .= -F_ab ./ force_scaling
    jacob[irow:irow+2, 16:18] .= -F_αb ./ force_scaling

    @views jacob[irow:irow+2, icol:icol+2] .-= F_u ./ force_scaling
    jacob[irow:irow+2, icol+6:icol+8] .= -F_V ./ force_scaling
    jacob[irow:irow+2, icol+9:icol+11] .= -F_Ω ./ force_scaling

    jacob[irow+3:irow+5, 4:6] .= -M_θb ./ force_scaling
    jacob[irow+3:irow+5, 10:12] .= -M_ωb ./ force_scaling
    jacob[irow+3:irow+5, 13:15] .= -M_ab ./ force_scaling
    jacob[irow+3:irow+5, 16:18] .= -M_αb ./ force_scaling

    @views jacob[irow+3:irow+5, icol:icol+2] .-= M_u ./ force_scaling
    jacob[irow+3:irow+5, icol+6:icol+8] .= -M_V ./ force_scaling
    jacob[irow+3:irow+5, icol+9:icol+11] .= -M_Ω ./ force_scaling

    jacob[irow+6:irow+8, 7:9] = rV_vb
    jacob[irow+6:irow+8, 10:12] = rV_ωb

    jacob[irow+6:irow+8, icol:icol+2] .= rV_u
    jacob[irow+6:irow+8, icol+6:icol+8] .= rV_V

    jacob[irow+9:irow+11, 10:12] = rΩ_ωb

    jacob[irow+9:irow+11, icol+3:icol+5] .= rΩ_θ
    jacob[irow+9:irow+11, icol+9:icol+11] .= rΩ_Ω

    return jacob
end

"""
    insert_mass_matrix_point_jacobians!(jacob, gamma, indices, force_scaling, ipoint, 
        resultants, velocities)

Insert the mass matrix jacobian entries corresponding to a point into the system jacobian 
matrix.
"""
@inline function insert_mass_matrix_point_jacobians!(jacob, gamma, indices, force_scaling, ipoint,  
    resultants, velocities)

    @unpack F_Vdot, F_Ωdot, M_Vdot, M_Ωdot = resultants
    @unpack rV_udot, rΩ_θdot = velocities

    irow = indices.irow_point[ipoint]
    icol = indices.icol_point[ipoint]

    @views jacob[irow:irow+2, icol+6:icol+8] .-= F_Vdot .* gamma ./ force_scaling
    @views jacob[irow:irow+2, icol+9:icol+11] .-= F_Ωdot .* gamma ./ force_scaling

    @views jacob[irow+3:irow+5, icol+6:icol+8] .-= M_Vdot .* gamma ./ force_scaling
    @views jacob[irow+3:irow+5, icol+9:icol+11] .-= M_Ωdot .* gamma ./ force_scaling

    @views jacob[irow+6:irow+8, icol:icol+2] .+= rV_udot .* gamma
    @views jacob[irow+9:irow+11, icol+3:icol+5] .+= rΩ_θdot .* gamma

    return jacob
end

# --- Point Residual --- #

"""
    dynamic_point_jacobian!(jacob, dx, x, indices, force_scaling, assembly, ipoint, 
        prescribed_conditions, point_masses, gravity, ub, θb, vb, ωb, ab, αb)

Calculate and insert the jacobian entries corresponding to a point for a dynamic 
analysis into the system jacobian matrix.
"""
@inline function dynamic_point_jacobian!(jacob, dx, x, indices, force_scaling, assembly, 
    ipoint, prescribed_conditions, point_masses, gravity, ub, θb, vb, ωb, ab, αb,
    ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb)

    properties = dynamic_point_properties(dx, x, indices, force_scaling, assembly, ipoint, 
        prescribed_conditions, point_masses, gravity, ub, θb, vb, ωb, ab, αb)

    properties = dynamic_point_jacobian_properties(properties, dx, x, indices, 
        force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity,
        ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb)

    resultants = dynamic_point_resultant_jacobians(properties)

    velocities = dynamic_point_velocity_jacobians(properties)

    insert_dynamic_point_jacobians!(jacob, indices, force_scaling, ipoint, resultants, velocities)

    return jacob
end

"""
    mass_matrix_point_jacobian!(jacob, gamma, x, indices, force_scaling, assembly, 
        ipoint, prescribed_conditions)

Calculate and insert the mass_matrix jacobian entries corresponding to a point into the 
system jacobian matrix.
"""
@inline function mass_matrix_point_jacobian!(jacob, gamma, x, indices, force_scaling, assembly, 
    ipoint, prescribed_conditions, point_masses)

    properties = mass_matrix_point_jacobian_properties(x, indices, force_scaling, assembly, ipoint, prescribed_conditions, point_masses)

    resultants = mass_matrix_point_resultant_jacobians(properties)

    velocities = mass_matrix_point_velocity_jacobians(properties)
    
    insert_mass_matrix_point_jacobians!(jacob, gamma, indices, force_scaling, ipoint,  
        resultants, velocities)

    return jacob
end