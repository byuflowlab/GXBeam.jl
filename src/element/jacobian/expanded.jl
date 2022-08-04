# --- Element Properties --- #

"""
    expanded_element_jacobian_properties(properties, x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity)

Calculate/extract the element properties needed to calculate the jacobian entries 
corresponding to a element for a steady state analysis
"""
@inline function expanded_steady_element_jacobian_properties(properties, x, indices, 
    force_scaling, structural_damping, assembly, ielem, prescribed_conditions, gravity,
    ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb)

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

    if structural_damping

        @unpack C1, C2, Qinv1, Qinv2, u1, u2, θ1, θ2, Ω1, Ω2, ω1, ω2, μ11, μ22, udot, θdot, 
            Δx1, Δx2, Δu, Δθ, ΔQ, Δudot, Δθdot = properties

        # rotation parameter matrices
        Qinv1_θ1, Qinv1_θ2, Qinv1_θ3 = get_Qinv_θ(θ1)
        Qinv2_θ1, Qinv2_θ2, Qinv2_θ3 = get_Qinv_θ(θ2)

        # linear displacement rates
        udot1_vb = -vb_vb
        udot2_vb = -vb_vb

        udot1_ωb = tilde(Δx1 + u1)*ωb_ωb
        udot2_ωb = tilde(Δx2 + u2)*ωb_ωb

        udot1_u1 = -tilde(ω1)
        udot2_u2 = -tilde(ω2)

        udot1_θ1 = mul3(C1_θ1', C1_θ2', C1_θ3', V1)
        udot2_θ2 = mul3(C2_θ1', C2_θ2', C2_θ3', V2)

        udot1_V1 = C1'
        udot2_V2 = C2'

        # angular displacement rates
        θdot1_ωb = -Qinv1*C1*ωb_ωb
        θdot2_ωb = -Qinv2*C2*ωb_ωb

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

        γdot_ωb = -CtCab'*tilde(Δu)*ω_ωb + CtCab'*Δudot_ωb - L*CtCab'*tilde(Cab*e1)*ω_ωb

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
    expanded_dynamic_element_jacobian_properties(properties, dx, x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity,
        ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb)

Calculate/extract the element properties needed to calculate the jacobian entries 
corresponding to a element for a dynamic analysis
"""
@inline function expanded_dynamic_element_jacobian_properties(properties, dx, x, indices, 
    force_scaling, structural_damping, assembly, ielem, prescribed_conditions, gravity,
    ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb)

    properties = expanded_steady_element_jacobian_properties(properties, x, indices, 
        force_scaling, structural_damping, assembly, ielem, prescribed_conditions, gravity,
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

# --- Compatability Residaul --- #

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

# --- Velocity Residual --- #

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

# --- Equilibrium Residual --- #

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

# --- Element Resultants --- # 

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

# --- Insert Element Residual --- #

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

# --- Element Residual --- #

"""
    expanded_steady_element_jacobian!(jacob, x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, 
        gravity, ub, θb, vb, ωb, ab, αb, ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb)

Calculate and insert the jacobian entries corresponding to a beam element for a dynamic
analysis into the system jacobian matrix.
"""
@inline function expanded_steady_element_jacobian!(jacob, x, indices, force_scaling, structural_damping, 
    assembly, ielem, prescribed_conditions, distributed_loads, gravity, ub, θb, vb, ωb, ab, αb, 
    ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb)

    properties = expanded_steady_element_properties(x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity, ub, θb, vb, ωb, ab, αb)

    properties = expanded_steady_element_jacobian_properties(properties, x, indices, 
        force_scaling, structural_damping, assembly, ielem, prescribed_conditions, gravity,
        ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb)

    compatability = expanded_compatability_jacobians(properties)

    velocities = expanded_element_velocity_jacobians(properties)

    equilibrium = expanded_element_equilibrium_jacobians(properties, distributed_loads, ielem)

    resultants = expanded_element_resultant_jacobians(properties)

    insert_expanded_element_jacobians!(jacob, indices, force_scaling, assembly, ielem, 
        compatability, velocities, equilibrium, resultants)

    return jacob
end

"""
    expanded_dynamic_element_jacobian!(jacob, dx, x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, 
        gravity, ub, θb, vb, ωb, ab, αb, ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb)

Calculate and insert the jacobian entries corresponding to a beam element for a dynamic
analysis into the system jacobian matrix.
"""
@inline function expanded_dynamic_element_jacobian!(jacob, dx, x, indices, force_scaling, structural_damping, 
    assembly, ielem, prescribed_conditions, distributed_loads, gravity, ub, θb, vb, ωb, ab, αb,
    ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb)

    properties = expanded_dynamic_element_properties(dx, x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity, ub, θb, vb, ωb, ab, αb)

    properties = expanded_dynamic_element_jacobian_properties(properties, dx, x, indices, 
        force_scaling, structural_damping, assembly, ielem, prescribed_conditions, gravity,
        ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb)

    compatability = expanded_compatability_jacobians(properties)

    velocities = expanded_element_velocity_jacobians(properties)

    equilibrium = expanded_element_equilibrium_jacobians(properties, distributed_loads, ielem)

    resultants = expanded_element_resultant_jacobians(properties)

    insert_expanded_element_jacobians!(jacob, indices, force_scaling, assembly, ielem, 
        compatability, velocities, equilibrium, resultants)

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