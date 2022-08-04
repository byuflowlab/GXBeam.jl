# --- Element Properties --- #

"""
    initial_condition_element_jacobian_properties(properties, x, indices, rate_vars, 
        force_scaling, structural_damping, assembly, ielem, prescribed_conditions, gravity, 
        ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb, u0, θ0, V0, Ω0, Vdot0, Ωdot0)

Calculate/extract the element properties needed to calculate the jacobian entries 
corresponding to a element for a Newmark scheme time marching analysis
"""
@inline function initial_condition_element_jacobian_properties(properties, x, indices, rate_vars, 
    force_scaling, structural_damping, assembly, ielem, prescribed_conditions, gravity, 
    ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb, u0, θ0, V0, Ω0, Vdot0, Ωdot0)

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
    v_vb = vb_vb
    v_ωb = (-tilde(Δx) - tilde(u))*ωb_ωb
    ω_ωb = ωb_ωb   

    # linear and angular acceleration
    a_u = tilde(αb)
    a_ab = ab_ab
    a_αb = (-tilde(Δx) - tilde(u))*αb_αb
    α_αb = αb_αb

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
    a_θb = -mul3(get_C_θ(θb)..., gravity)*θb_θb
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

# --- Compatability Residual --- #

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

# --- Element Resultants --- #

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

    tmp = 1/2*(tilde(ω)*H_ωb + tilde(V)*P_ωb - tilde(P)*V_ωb - tilde(H)*ω_ωb)

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

    F1_θ2 -= 1/2*tmp*θ2_θ2
    F2_θ2 += 1/2*tmp*θ2_θ2

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

    M1_θ2 -= 1/2*tmp*θ2_θ2
    M2_θ2 += 1/2*tmp*θ2_θ2

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

# --- Element Residual --- #

"""
    initial_condition_element_jacobian!(jacob, x, indices, rate_vars, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, 
        gravity, ub, θb, vb, ωb, ab, αb, u0, θ0, V0, Ω0, Vdot0, Ωdot0)

Calculate and insert the jacobian entries corresponding to a beam element for the 
initialization of a time domain analysis into the system jacobian matrix.
"""
@inline function initial_condition_element_jacobian!(jacob, x, indices, rate_vars, force_scaling, 
    structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, gravity, 
    ub, θb, vb, ωb, ab, αb, ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb, u0, θ0, V0, Ω0, 
    Vdot0, Ωdot0)

    properties = initial_condition_element_properties(x, indices, rate_vars, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity, 
        ub, θb, vb, ωb, ab, αb, u0, θ0, V0, Ω0, Vdot0, Ωdot0)

    properties = initial_condition_element_jacobian_properties(properties, x, indices, 
        rate_vars, force_scaling, structural_damping, assembly, ielem, prescribed_conditions, 
        gravity, ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb, u0, θ0, V0, Ω0, Vdot0, Ωdot0)

    compatability = initial_condition_compatability_jacobians(properties)

    resultants = initial_condition_element_resultant_jacobians(properties, distributed_loads, ielem)

    insert_initial_condition_element_jacobians!(jacob, indices, force_scaling, assembly, ielem, 
        compatability, resultants)

    return jacob
end
