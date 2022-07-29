# --- Element Properties --- #

"""
    newmark_element_jacobian_properties(properties, x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity, 
        ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb, Vdot_init, Ωdot_init, dt)

Calculate/extract the element properties needed to calculate the jacobian entries 
corresponding to a element for a Newmark scheme time marching analysis
"""
@inline function newmark_element_jacobian_properties(properties, x, indices, 
    force_scaling, structural_damping, assembly, ielem, prescribed_conditions, gravity, 
    ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb, Vdot_init, Ωdot_init, dt)

    properties = steady_state_element_jacobian_properties(properties, x, indices, 
        force_scaling, structural_damping, assembly, ielem, 
        prescribed_conditions, gravity, ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb)

    @unpack C, Cab, CtCab, CtCabdot, mass11, mass12, mass21, mass22, V, Ω, ω, Vdot, Ωdot, 
        C_θ1, C_θ2, C_θ3 = properties

    # linear and angular momentum rates

    Pdot_ωb = -CtCab*mass11*CtCab'*tilde(V) - CtCab*mass12*CtCab'*tilde(Ω) +
        tilde(CtCab*mass11*CtCab'*V) + tilde(CtCab*mass12*CtCab'*Ω)

    Pdot_ωb *= ωb_ωb

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

    Hdot_ωb = -CtCab*mass21*CtCab'*tilde(V) - CtCab*mass22*CtCab'*tilde(Ω) +
        tilde(CtCab*mass21*CtCab'*V) + tilde(CtCab*mass22*CtCab'*Ω)

    Hdot_ωb *= ωb_ωb

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


    return (; properties..., Pdot_ωb, Pdot_θ, Pdot_V, Pdot_Ω, Hdot_ωb, Hdot_θ, Hdot_V, Hdot_Ω) 
end

# --- Element Resultants --- #

"""
    newmark_element_resultant_jacobians(properties, distributed_loads, ielem)

Calculate the jacobians for the resultant loads applied at each end of a beam element 
for a Newmark scheme time marching analysis.
"""
@inline function newmark_element_resultant_jacobians(properties, distributed_loads, ielem)

    jacobians = steady_state_element_resultant_jacobians(properties, distributed_loads, ielem)

    @unpack L, Pdot_ωb, Pdot_θ, Pdot_V, Pdot_Ω, Hdot_ωb, Hdot_θ, Hdot_V, Hdot_Ω, θ1_θ1, θ2_θ2 = properties
    
    @unpack F1_ωb, F1_θ1, F1_θ2, F1_V1, F1_V2, F1_Ω1, F1_Ω2, 
            F2_ωb, F2_θ1, F2_θ2, F2_V1, F2_V2, F2_Ω1, F2_Ω2, 
            M1_ωb, M1_θ1, M1_θ2, M1_V1, M1_V2, M1_Ω1, M1_Ω2, 
            M2_ωb, M2_θ1, M2_θ2, M2_V1, M2_V2, M2_Ω1, M2_Ω2 = jacobians

    # add loads due to linear and angular momentum rates

    tmp = 1/2*Pdot_ωb

    F1_ωb -= tmp
    F2_ωb += tmp

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

    tmp = 1/2*Hdot_ωb

    M1_ωb -= tmp
    M2_ωb += tmp

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
        F1_ωb, F1_θ1, F1_θ2, F1_V1, F1_V2, F1_Ω1, F1_Ω2, 
        F2_ωb, F2_θ1, F2_θ2, F2_V1, F2_V2, F2_Ω1, F2_Ω2, 
        M1_ωb, M1_θ1, M1_θ2, M1_V1, M1_V2, M1_Ω1, M1_Ω2, 
        M2_ωb, M2_θ1, M2_θ2, M2_V1, M2_V2, M2_Ω1, M2_Ω2)
end

# --- Element Residual --- #

"""
    newmark_element_jacobian!(jacob, x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, 
        gravity, ub, θb, vb, ωb, ab, αb, ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb, Vdot_init, Ωdot_init, dt)

Calculate and insert the jacobian entries corresponding to a beam element for a Newmark-scheme
time marching analysis into the system jacobian matrix.
"""
@inline function newmark_element_jacobian!(jacob, x, indices, force_scaling, 
    structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, gravity, 
    ub, θb, vb, ωb, ab, αb, ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb, Vdot_init, Ωdot_init, dt)

    properties = newmark_element_properties(x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity, 
        ub, θb, vb, ωb, ab, αb, Vdot_init, Ωdot_init, dt)

    properties = newmark_element_jacobian_properties(properties, x, indices, 
        force_scaling, structural_damping, assembly, ielem, prescribed_conditions, gravity, 
        ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb, Vdot_init, Ωdot_init, dt)

    compatability = dynamic_compatability_jacobians(properties)

    resultants = newmark_element_resultant_jacobians(properties, distributed_loads, ielem)

    insert_dynamic_element_jacobians!(jacob, indices, force_scaling, assembly, ielem, 
        compatability, resultants)

    return jacob
end