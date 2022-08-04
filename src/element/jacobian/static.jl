# --- Element Properties --- #

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

# --- Compatability Residual --- #

@inline function static_compatability_jacobians(properties)
   
    @unpack L, Cab, CtCab, Qinv, γ, κ, u1_u1, u2_u2, θ1_θ1, θ2_θ2, γ_F, γ_M, κ_F, κ_M, 
        C_θ1, C_θ2, C_θ3, Qinv_θ1, Qinv_θ2, Qinv_θ3 = properties

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

# --- Element Resultants --- #

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

# --- Insert Element Residual --- #

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
    @views jacob[irow2:irow2+2, icol2+3:icol2+5] .+= F2_θ2 ./ force_scaling

    jacob[irow2:irow2+2, icol:icol+2] .= F2_F

    @views jacob[irow2+3:irow2+5, icol1+3:icol1+5] .+= M2_θ1./ force_scaling
    @views jacob[irow2+3:irow2+5, icol2+3:icol2+5] .+= M2_θ2 ./ force_scaling

    jacob[irow2+3:irow2+5, icol:icol+2] .= M2_F
    jacob[irow2+3:irow2+5, icol+3:icol+5] .= M2_M

    return jacob
end

# --- Element Residual --- #

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

    compatability = static_compatability_jacobians(properties)

    resultants = static_element_resultant_jacobians(properties, distributed_loads, ielem)

    insert_static_element_jacobians!(jacob, indices, force_scaling, assembly, ielem, 
        compatability, resultants)

    return jacob
end