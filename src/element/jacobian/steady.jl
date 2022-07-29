# --- Element Properties --- #

"""
    steady_state_element_jacobian_properties(properties, x, indices, 
        force_scaling, structural_damping, assembly, ielem, prescribed_conditions, 
        gravity, ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb)

Calculate/extract the element properties needed to calculate the jacobian entries 
corresponding to a element for a steady state analysis
"""
@inline function steady_state_element_jacobian_properties(properties, x, indices, 
    force_scaling, structural_damping, assembly, ielem, prescribed_conditions, 
    gravity, ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb)

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
    v_vb = vb_vb
    v_ωb = (-tilde(Δx) - tilde(u))*ωb_ωb
    
    # angular velocity
    ω_ωb = ωb_ωb

    # linear acceleration
    a_u = tilde(αb)
    a_ab = ab_ab
    a_αb = (-tilde(Δx) - tilde(u))*αb_αb
    
    # angular acceleration
    α_αb = αb_αb

    # add rotated gravity vector
    C_θb1, C_θb2, C_θb3 = get_C_θ(θb)
    a_θb = -mul3(C_θb1, C_θb2, C_θb3, gravity)*θb_θb

    if structural_damping

        @unpack C1, C2, Qinv1, Qinv2, u1, u2, θ1, θ2, Ω1, Ω2, ω1, ω2, μ11, μ22, udot, θdot, 
            Δx1, Δx2, Δu, Δθ, ΔQ, Δudot, Δθdot = properties

        # rotation parameter matrices
        C1_θ1, C1_θ2, C1_θ3 = get_C_θ(C1, θ1)
        C2_θ1, C2_θ2, C2_θ3 = get_C_θ(C2, θ2)
        Qinv1_θ1, Qinv1_θ2, Qinv1_θ3 = get_Qinv_θ(θ1)
        Qinv2_θ1, Qinv2_θ2, Qinv2_θ3 = get_Qinv_θ(θ2)

        # linear displacement rates
        udot1_vb = -vb_vb
        udot2_vb = -vb_vb

        udot1_ωb = tilde(Δx1 + u1)*ωb_ωb
        udot2_ωb = tilde(Δx2 + u2)*ωb_ωb

        udot1_u1 = -tilde(ω1)
        udot2_u2 = -tilde(ω2)

        udot1_V1 = I3
        udot2_V2 = I3

        # angular displacement rates
        θdot1_ωb = -Qinv1*C1*ωb_ωb
        θdot2_ωb = -Qinv2*C2*ωb_ωb

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

        γdot_ωb = -CtCab'*tilde(Δu)*ωb_ωb + CtCab'*Δudot_ωb - L*CtCab'*tilde(Cab*e1)*ωb_ωb

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

# --- Element Resultants --- #

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

# --- Element Residual --- #

"""
    steady_state_element_jacobian!(jacob, x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, 
        gravity, ub, θb, vb, ωb, ab, αb, ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb)

Calculate and insert the jacobian entries corresponding to a beam element for a steady state 
analysis into the system jacobian matrix.
"""
@inline function steady_state_element_jacobian!(jacob, x, indices, 
    force_scaling, structural_damping, assembly, ielem, prescribed_conditions, 
    distributed_loads, gravity, ub, θb, vb, ωb, ab, αb, ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb)

    properties = steady_state_element_properties(x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity, 
        ub, θb, vb, ωb, ab, αb)

    properties = steady_state_element_jacobian_properties(properties, x, indices, 
        force_scaling, structural_damping, assembly, ielem, 
        prescribed_conditions, gravity, ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb)

    compatability = dynamic_compatability_jacobians(properties)

    resultants = steady_state_element_resultant_jacobians(properties, distributed_loads, ielem)

    insert_dynamic_element_jacobians!(jacob, indices, force_scaling, assembly, ielem, 
        compatability, resultants)

    return jacob
end
