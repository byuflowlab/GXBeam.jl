# --- Helper Functions --- #

"""
    initial_point_displacement_jacobian(ipoint, icol_point, prescribed_conditions, 
        rate_vars)

Extract the displacement jacobians `u_u` and `θ_θ` of point `ipoint` from the state 
variable vector or prescribed conditions for an initial condition analysis.
"""
@inline function initial_point_displacement_jacobian(ipoint, icol_point, 
    prescribed_conditions, rate_vars)

    if haskey(prescribed_conditions, ipoint)
        u_u, θ_θ = initial_point_displacement_jacobian(icol_point[ipoint], 
            prescribed_conditions[ipoint], rate_vars)
    else
        u_u, θ_θ = initial_point_displacement_jacobian(icol_point[ipoint], 
            rate_vars)
    end

    return u_u, θ_θ
end

@inline function initial_point_displacement_jacobian(icol, 
    prescribed_conditions, rate_vars)

    # Use prescribed displacement, if applicable.  If no displacements are prescribed use 
    # a component of `u` or `θ` as state variables if the corresponding component of `Vdot` 
    # or `Ωdot` cannot be a state variable.

    # unpack prescribed conditions for the node
    @unpack pd = prescribed_conditions

    # node linear displacement
    u_u = hcat(
        ifelse(pd[1], zero(e1), ifelse(rate_vars[icol+6], zero(e1), e1)),
        ifelse(pd[2], zero(e2), ifelse(rate_vars[icol+7], zero(e2), e2)),
        ifelse(pd[3], zero(e3), ifelse(rate_vars[icol+8], zero(e3), e3))
    )

    # node angular displacement
    θ_θ = hcat(
        ifelse(pd[4], zero(e1), ifelse(rate_vars[icol+9], zero(e1), e1)),
        ifelse(pd[5], zero(e2), ifelse(rate_vars[icol+10], zero(e2), e2)),
        ifelse(pd[6], zero(e3), ifelse(rate_vars[icol+11], zero(e3), e3))
    )   

    return u_u, θ_θ
end

@inline function initial_point_displacement_jacobian(icol, rate_vars)

    # Use a component of `u` or `θ` as state variables if the corresponding component of 
    # `Vdot` or `Ωdot` cannot be a state variable.

    u_u = hcat(
        ifelse(rate_vars[icol+6], zero(e1), e1),
        ifelse(rate_vars[icol+7], zero(e2), e2),
        ifelse(rate_vars[icol+8], zero(e3), e3)
    )
    θ_θ = hcat(
        ifelse(rate_vars[icol+9], zero(e1), e1),
        ifelse(rate_vars[icol+10], zero(e2), e2),
        ifelse(rate_vars[icol+11], zero(e3), e3)
    )

    return u_u, θ_θ
end

"""
    initial_point_velocity_rate_jacobian(ipoint, icol_point, prescribed_conditions, 
        rate_vars)

Return the velocity rate jacobians `Vdot_Vdot` and `Ωdot_Ωdot` of point `ipoint`.  Note 
that `Vdot` and `Ωdot` in this case do not minclude any contributions resulting from 
body frame motion. 
"""
@inline function initial_point_velocity_rate_jacobian(ipoint, icol_point, 
    prescribed_conditions, rate_vars)

    if haskey(prescribed_conditions, ipoint)
        Vdot_Vdot, Ωdot_Ωdot = initial_point_velocity_rate_jacobian(icol_point[ipoint], 
            prescribed_conditions[ipoint], rate_vars)
    else
        Vdot_Vdot, Ωdot_Ωdot = initial_point_velocity_rate_jacobian(icol_point[ipoint], rate_vars)
    end

    return Vdot_Vdot, Ωdot_Ωdot
end

@inline function initial_point_velocity_rate_jacobian(icol, prescribed_conditions, 
    rate_vars)

    # If a displacment is prescribed, then the corresponding component of Vdot or Ωdot 
    # (relative to the body frame) is zero.  If no displacements is prescribed use the 
    # corresponding component of `Vdot` or `Ωdot` as a state variable, if possible.  
    # Otherwise, use the provided value.

    # unpack prescribed conditions for the node
    @unpack pd, u, theta = prescribed_conditions

    # node linear displacement
    Vdot_Vdot = hcat(
        ifelse(pd[1], zero(e1), ifelse(rate_vars[icol+6], e1, zero(e1))),
        ifelse(pd[2], zero(e2), ifelse(rate_vars[icol+7], e2, zero(e2))),
        ifelse(pd[3], zero(e3), ifelse(rate_vars[icol+8], e3, zero(e3)))
    )

    # node angular displacement
    Ωdot_Ωdot = hcat(
        ifelse(pd[4], zero(e1), ifelse(rate_vars[icol+9], e1, zero(e1))),
        ifelse(pd[5], zero(e2), ifelse(rate_vars[icol+10], e2, zero(e2))),
        ifelse(pd[6], zero(e3), ifelse(rate_vars[icol+11], e3, zero(e3)))
    )   

    return Vdot_Vdot, Ωdot_Ωdot
end

@inline function initial_point_velocity_rate_jacobian(icol, rate_vars)

    # Use the components of `Vdot` and `Ωdot` as state variables, if possible. Otherwise, 
    # use the provided value.

    Vdot_Vdot = hcat(
        ifelse(rate_vars[icol+6], e1, zero(e1)),
        ifelse(rate_vars[icol+7], e2, zero(e2)),
        ifelse(rate_vars[icol+8], e3, zero(e3))
    )
    Ωdot_Ωdot = hcat(
        ifelse(rate_vars[icol+9], e1, zero(e1)),
        ifelse(rate_vars[icol+10], e2, zero(e2)),
        ifelse(rate_vars[icol+11], e3, zero(e3))
    )

    return Vdot_Vdot, Ωdot_Ωdot
end

# --- Point Properties --- #

"""
    initial_condition_point_jacobian_properties(properties, x, indices, rate_vars, 
        force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity, 
        ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb, u0, θ0, V0, Ω0, Vdot0, Ωdot0)

Calculate/extract the point properties needed to calculate the jacobian entries 
corresponding to a point for a Newmark scheme time marching analysis
"""
@inline function initial_condition_point_jacobian_properties(properties, x, indices, 
    rate_vars, force_scaling, assembly, ipoint, prescribed_conditions, point_masses, 
    gravity, ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb, u0, θ0, V0, Ω0, Vdot0, Ωdot0)

    @unpack C, Cdot, mass11, mass12, mass21, mass22, u, θ, V, Ω, Cdot, Vdot, Ωdot, 
        θb, ωb, αb, Δx, v, ω, a, α = properties

    # linear and angular displacement
    u_u, θ_θ = initial_point_displacement_jacobian(ipoint, indices.icol_point, 
        prescribed_conditions, rate_vars)

    # rotation parameter matrices
    C_θ1, C_θ2, C_θ3 = get_C_θ(C, θ)
    Qinv_θ1, Qinv_θ2, Qinv_θ3 = get_Qinv_θ(θ)

    # forces and moments
    F_θ, F_F, M_θ, M_M = point_load_jacobians(x, ipoint, indices.icol_point, force_scaling, prescribed_conditions)

    # linear and angular velocity
    v_u = tilde(ωb)
    v_vb = vb_vb
    v_ωb = (-tilde(Δx)-tilde(u))*ωb_ωb
    ω_ωb = ωb_ωb

    # linear and angular acceleration
    a_u = tilde(αb)
    a_ab = ab_ab
    a_αb = (-tilde(Δx)-tilde(u))*αb_αb
    α_αb = αb_αb

    # add contributions from body frame motion to velocities
    V_u = v_u
    V_vb = v_vb
    V_ωb = v_ωb
    Ω_ωb = ω_ωb

    # linear and angular momentum
    P_u = C'*mass11*C*V_u
    P_θ = mul3(C_θ1', C_θ2', C_θ3', (mass11*C*V + mass12*C*Ω)) + 
        C'*mass11*mul3(C_θ1, C_θ2, C_θ3, V) + 
        C'*mass12*mul3(C_θ1, C_θ2, C_θ3, Ω)
    P_vb = C'*mass11*C*V_vb
    P_ωb = C'*mass11*C*V_ωb + C'*mass12*C*Ω_ωb

    H_u = C'*mass21*C*V_u
    H_θ = mul3(C_θ1', C_θ2', C_θ3', (mass21*C*V + mass22*C*Ω)) + 
        C'*mass21*mul3(C_θ1, C_θ2, C_θ3, V) + 
        C'*mass22*mul3(C_θ1, C_θ2, C_θ3, Ω)
    H_vb = C'*mass21*C*V_vb
    H_ωb = C'*mass21*C*V_ωb + C'*mass22*C*Ω_ωb

    # linear and angular acceleration
    Vdot_Vdot, Ωdot_Ωdot = initial_point_velocity_rate_jacobian(ipoint, indices.icol_point, 
        prescribed_conditions, rate_vars)

    # add contributions from body frame motion to accelerations
    Vdot_u = a_u
    Vdot_ab = a_ab
    Vdot_αb = a_αb
    Ωdot_αb = α_αb

    # linear and angular momentum rates
    Pdot_V = C'*mass11*Cdot + Cdot'*mass11*C
    Pdot_Ω = C'*mass12*Cdot + Cdot'*mass12*C
    Pdot_Vdot = C'*mass11*C
    Pdot_Ωdot = C'*mass12*C
    Pdot_u = Pdot_V*V_u + Pdot_Vdot*Vdot_u
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
    Pdot_vb = Pdot_V*V_vb
    Pdot_ωb = Pdot_V*V_ωb + Pdot_Ω*Ω_ωb
    Pdot_ab = Pdot_Vdot*Vdot_ab
    Pdot_αb = Pdot_Vdot*Vdot_αb + Pdot_Ωdot*Ωdot_αb

    Hdot_V = C'*mass21*Cdot + Cdot'*mass21*C
    Hdot_Ω = C'*mass22*Cdot + Cdot'*mass22*C
    Hdot_Vdot = C'*mass21*C
    Hdot_Ωdot = C'*mass22*C
    Hdot_u = Hdot_V*V_u + Hdot_Vdot*Vdot_u
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
    Hdot_vb = Hdot_V*V_vb
    Hdot_ωb = Hdot_V*V_ωb + Hdot_Ω*Ω_ωb
    Hdot_ab = Hdot_Vdot*Vdot_ab
    Hdot_αb = Hdot_Vdot*Vdot_αb + Hdot_Ωdot*Ωdot_αb

    # overwrite acceleration terms so we don't double count them
    a_u = @SMatrix zeros(3,3)
    a_θb = -mul3(get_C_θ(θb)..., gravity)*θb_θb
    a_ab = @SMatrix zeros(3,3)
    a_αb = @SMatrix zeros(3,3)
    α_αb = @SMatrix zeros(3,3)

    return (; properties..., C_θ1, C_θ2, C_θ3, Qinv_θ1, Qinv_θ2, Qinv_θ3, u_u, θ_θ, 
        F_θ, F_F, M_θ, M_M, V_u, V_vb, V_ωb, Ω_ωb, P_u, P_θ, P_vb, P_ωb, H_u, H_θ, H_vb, H_ωb, 
        Vdot_Vdot, Ωdot_Ωdot, 
        Pdot_vb, Pdot_ωb, Pdot_ab, Pdot_αb, Pdot_u, Pdot_θ, Pdot_Vdot, Pdot_Ωdot,
        Hdot_vb, Hdot_ωb, Hdot_ab, Hdot_αb, Hdot_u, Hdot_θ, Hdot_Vdot, Hdot_Ωdot, 
        v_vb, v_ωb, ω_ωb, a_u, a_θb, a_ab, a_αb, α_αb) 
end

# --- Point Resultants --- #

"""
    initial_condition_point_resultant_jacobians(properties)

Calculate the jacobians for the net loads `F` and `M` applied at a point for the initialization
of a time domain analysis.
"""
@inline function initial_condition_point_resultant_jacobians(properties)

    jacobians = static_point_resultant_jacobians(properties)

    @unpack  C, mass11, mass12, mass21, mass22, V, Ω, P, H, ω, u_u, θ_θ, Vdot_Vdot, Ωdot_Ωdot,
        V_u, V_vb, V_ωb, P_u, P_θ, P_vb, P_ωb, H_u, H_θ, H_vb, H_ωb,     
        Pdot_vb, Pdot_ωb, Pdot_ab, Pdot_αb, Pdot_u, Pdot_θ, Pdot_Vdot, Pdot_Ωdot,
        Hdot_vb, Hdot_ωb, Hdot_ab, Hdot_αb, Hdot_u, Hdot_θ, Hdot_Vdot, Hdot_Ωdot, 
        v_vb, v_ωb, ω_ωb, a_u, a_θb, a_ab, a_αb, α_αb = properties

    @unpack F_θ, M_θ = jacobians

    F_θ *= θ_θ
    M_θ *= θ_θ

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
    F_vb = -tilde(ω)*P_vb
    F_ωb = -tilde(ω)*P_ωb + tilde(P)*ω_ωb
    F_u -= tilde(ω)*P_u*u_u
    F_θ -= tilde(ω)*P_θ*θ_θ

    M_vb = -tilde(ω)*H_vb - tilde(V)*P_vb + tilde(P)*V_vb
    M_ωb = -tilde(ω)*H_ωb - tilde(V)*P_ωb + tilde(P)*V_ωb + tilde(H)*ω_ωb
    M_u -= (tilde(ω)*H_u + tilde(V)*P_u - tilde(P)*V_u)*u_u
    M_θ -= (tilde(ω)*H_θ + tilde(V)*P_θ)*θ_θ

    # # add loads due to linear and angular momentum rates
    F_vb -= Pdot_vb
    F_ωb -= Pdot_ωb
    F_ab -= Pdot_ab
    F_αb -= Pdot_αb
    F_u -= Pdot_u*u_u
    F_θ -= Pdot_θ*θ_θ
    F_Vdot = -Pdot_Vdot*Vdot_Vdot
    F_Ωdot = -Pdot_Ωdot*Ωdot_Ωdot

    M_vb -= Hdot_vb
    M_ωb -= Hdot_ωb
    M_ab -= Hdot_ab
    M_αb -= Hdot_αb
    M_u -= Hdot_u*u_u
    M_θ -= Hdot_θ*θ_θ
    M_Vdot = -Hdot_Vdot*Vdot_Vdot
    M_Ωdot = -Hdot_Ωdot*Ωdot_Ωdot

    return (; jacobians..., 
        F_θb, F_vb, F_ωb, F_ab, F_αb, F_u, F_θ, F_Vdot, F_Ωdot, 
        M_θb, M_vb, M_ωb, M_ab, M_αb, M_u, M_θ, M_Vdot, M_Ωdot)
end

# --- Velocity Residuals --- #

"""
    initial_condition_point_velocity_jacobians(properties)

Calculate the jacobians of the velocity residuals `rV` and `rΩ` of a point for the
initialization of a time domain analysis.
"""
@inline function initial_condition_point_velocity_jacobians(properties)
   
    @unpack C, Qinv, Ω, ω, v_vb, v_ωb, ω_ωb, u_u, θ_θ, C_θ1, C_θ2, C_θ3, Qinv_θ1, Qinv_θ2, Qinv_θ3 = properties

    rΩ_θ = (mul3(Qinv_θ1, Qinv_θ2, Qinv_θ3, C*(Ω - ω)) + Qinv*mul3(C_θ1, C_θ2, C_θ3, Ω - ω))*θ_θ

    rV_udot = -I3
    rΩ_θdot = -I3

    return (; rV_udot, rΩ_θ, rΩ_θdot)
end

# --- Residual Placement --- #

"""
    insert_initial_condition_point_jacobians!(jacob, indices, force_scaling, ipoint, resultants, 
        velocities)

Insert the jacobian entries corresponding to a point for the initialization of a 
time domain analysis into the system jacobian matrix.
"""
@inline function insert_initial_condition_point_jacobians!(jacob, indices, force_scaling, ipoint,  
    resultants, velocities)

    @unpack F_θb, F_vb, F_ωb, F_ab, F_αb, F_u, F_θ, F_F, F_Vdot, F_Ωdot, 
            M_θb, M_vb, M_ωb, M_ab, M_αb, M_u, M_θ, M_M, M_Vdot, M_Ωdot = resultants
    
    @unpack rV_udot, rΩ_θ, rΩ_θdot = velocities

    irow = indices.irow_point[ipoint]
    icol = indices.icol_point[ipoint]

    jacob[irow:irow+2, 4:6] .= -F_θb ./ force_scaling
    jacob[irow:irow+2, 7:9] .= -F_vb ./ force_scaling
    jacob[irow:irow+2, 10:12] .= -F_ωb ./ force_scaling
    jacob[irow:irow+2, 13:15] .= -F_ab ./ force_scaling
    jacob[irow:irow+2, 16:18] .= -F_αb ./ force_scaling

    jacob[irow:irow+2, icol:icol+2] .= -F_F .- F_u ./ force_scaling .- F_Vdot ./ force_scaling
    jacob[irow:irow+2, icol+3:icol+5] .= -F_θ ./ force_scaling .- F_Ωdot ./ force_scaling

    jacob[irow+3:irow+5, 4:6] .= -M_θb ./ force_scaling
    jacob[irow+3:irow+5, 7:9] .= -M_vb ./ force_scaling
    jacob[irow+3:irow+5, 10:12] .= -M_ωb ./ force_scaling
    jacob[irow+3:irow+5, 13:15] .= -M_ab ./ force_scaling
    jacob[irow+3:irow+5, 16:18] .= -M_αb ./ force_scaling

    jacob[irow+3:irow+5, icol:icol+2] .= -M_u ./ force_scaling .- M_Vdot ./ force_scaling 
    jacob[irow+3:irow+5, icol+3:icol+5] .= -M_M .- M_θ ./ force_scaling .- M_Ωdot ./ force_scaling

    jacob[irow+6:irow+8, icol+6:icol+8] .= rV_udot
    jacob[irow+9:irow+11, icol+3:icol+5] .= rΩ_θ
    jacob[irow+9:irow+11, icol+9:icol+11] .= rΩ_θdot

    return jacob
end

# --- Point Residual --- #

"""
    initial_condition_point_jacobian!(jacob, x, indices, rate_vars, 
        force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity, 
        ub, θb, vb, ωb, ab, αb, ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb, 
        u0, θ0, V0, Ω0, Vdot0, Ωdot0)

Calculate and insert the jacobian entries corresponding to a point for the initialization
of a time domain analysis into the system jacobian matrix.
"""
@inline function initial_condition_point_jacobian!(jacob, x, indices, rate_vars, 
    force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity, 
    ub, θb, vb, ωb, ab, αb, ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb, 
    u0, θ0, V0, Ω0, Vdot0, Ωdot0)

    properties = initial_condition_point_properties(x, indices, rate_vars, 
        force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity, 
        ub, θb, vb, ωb, ab, αb, u0, θ0, V0, Ω0, Vdot0, Ωdot0)
    
    properties = initial_condition_point_jacobian_properties(properties, x, indices, 
        rate_vars, force_scaling, assembly, ipoint, prescribed_conditions, 
        point_masses, gravity, ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb, u0, θ0, V0, Ω0, Vdot0, Ωdot0)

    resultants = initial_condition_point_resultant_jacobians(properties)

    velocities = initial_condition_point_velocity_jacobians(properties)

    insert_initial_condition_point_jacobians!(jacob, indices, force_scaling, ipoint, resultants, velocities)

    return jacob
end