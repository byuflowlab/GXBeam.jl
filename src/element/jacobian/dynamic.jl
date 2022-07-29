# --- Element Properties --- #

"""
    dynamic_element_jacobian_properties(properties, dx, x, indices, 
        force_scaling, assembly, ielem, prescribed_conditions, gravity, 
        ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb)

Calculate/extract the element properties needed to calculate the jacobian entries 
corresponding to a element for a dynamic analysis
"""
@inline function dynamic_element_jacobian_properties(properties, dx, x, indices, 
    force_scaling, structural_damping, assembly, ielem, prescribed_conditions, gravity,
    ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb)

    properties = steady_state_element_jacobian_properties(properties, x, indices, 
        force_scaling, structural_damping, assembly, ielem, 
        prescribed_conditions, gravity, ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb)

    @unpack C, Cab, CtCab, CtCabdot, mass11, mass12, mass21, mass22, V, Ω, Vdot, Ωdot, ω, C_θ1, C_θ2, C_θ3 = properties

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

    # overwrite acceleration terms so we don't double count them
    a_u = @SMatrix zeros(3,3)
    a_ab = @SMatrix zeros(3,3)
    a_αb = @SMatrix zeros(3,3)
    α_αb = @SMatrix zeros(3,3)

    return (; properties..., Pdot_ωb, Pdot_θ, Pdot_V, Pdot_Ω, Hdot_ωb, Hdot_θ, Hdot_V, Hdot_Ω,
        a_u, a_ab, a_αb, α_αb) 
end

"""
    mass_matrix_element_jacobian_properties(x, indices, force_scaling, 
    assembly, ielem, prescribed_conditions)

Calculate/extract the element properties needed to calculate the mass matrix jacobian entries 
corresponding to a beam element
"""
@inline function mass_matrix_element_jacobian_properties(x, indices, force_scaling, 
    assembly, ielem, prescribed_conditions)

    # element properties
    @unpack L, Cab, compliance, mass = assembly.elements[ielem]
   
    # scale mass matrix by the element length
    mass *= L

    # mass submatrices
    mass11 = mass[SVector{3}(1:3), SVector{3}(1:3)]
    mass12 = mass[SVector{3}(1:3), SVector{3}(4:6)]
    mass21 = mass[SVector{3}(4:6), SVector{3}(1:3)]
    mass22 = mass[SVector{3}(4:6), SVector{3}(4:6)]

    # linear and angular displacement
    u1, θ1 = point_displacement(x, assembly.start[ielem], indices.icol_point, prescribed_conditions)
    u2, θ2 = point_displacement(x, assembly.stop[ielem], indices.icol_point, prescribed_conditions)
    u = (u1 + u2)/2
    θ = (θ1 + θ2)/2

    # transformation matrices
    C = get_C(θ)
    Ct = C'
    CtCab = Ct*Cab

    # linear and angular momentum rates
    Pdot_Vdot = CtCab*mass11*CtCab'
    Pdot_Ωdot = CtCab*mass12*CtCab'
    Hdot_Vdot = CtCab*mass21*CtCab'
    Hdot_Ωdot = CtCab*mass22*CtCab'

    return (; L, Pdot_Vdot, Pdot_Ωdot, Hdot_Vdot, Hdot_Ωdot)
end

# --- Compatability Residual --- #

@inline function dynamic_compatability_jacobians(properties)
   
    jacobians = static_compatability_jacobians(properties)

    @unpack Cab, CtCab, Qinv, u1_u1, u2_u2, θ1_θ1, θ2_θ2, γ_vb, γ_ωb, γ_u1, γ_u2, γ_θ1, γ_θ2, 
        γ_V1, γ_V2, γ_Ω1, γ_Ω2, κ_ωb, κ_θ1, κ_θ2, κ_Ω1, κ_Ω2 = properties

    @unpack ru_u1, ru_u2, ru_θ1, ru_θ2, rθ_θ1, rθ_θ2 = jacobians

    ru_vb = -CtCab*γ_vb

    ru_ωb = -CtCab*γ_ωb

    ru_u1 -= CtCab*γ_u1*u1_u1
    ru_u2 -= CtCab*γ_u2*u2_u2

    ru_θ1 -= CtCab*γ_θ1*θ1_θ1
    ru_θ2 -= CtCab*γ_θ2*θ2_θ2
    
    ru_V1 = -CtCab*γ_V1
    ru_V2 = -CtCab*γ_V2

    ru_Ω1 = -CtCab*γ_Ω1
    ru_Ω2 = -CtCab*γ_Ω2

    rθ_ωb = -Qinv*Cab*κ_ωb

    rθ_θ1 -= Qinv*Cab*κ_θ1*θ1_θ1
    rθ_θ2 -= Qinv*Cab*κ_θ2*θ2_θ2
    
    rθ_Ω1 = -Qinv*Cab*κ_Ω1
    rθ_Ω2 = -Qinv*Cab*κ_Ω2

    return (; jacobians..., ru_vb, ru_ωb, ru_u1, ru_u2, ru_θ1, ru_θ2, ru_V1, ru_V2, ru_Ω1, ru_Ω2, 
        rθ_ωb, rθ_θ1, rθ_θ2, rθ_Ω1, rθ_Ω2)
end

# --- Element Resultants --- #

"""
    dynamic_element_resultant_jacobians(properties, distributed_loads, ielem)

Calculate the jacobians for the resultant loads applied at each end of a beam element 
for a dynamic analysis.
"""
@inline function dynamic_element_resultant_jacobians(properties, distributed_loads, ielem)

    return newmark_element_resultant_jacobians(properties, distributed_loads, ielem)
end

"""
    mass_matrix_element_resultant_jacobians(properties)

Calculate the mass matrix jacobians for the resultant loads applied at each end of a beam element 
"""
@inline function mass_matrix_element_resultant_jacobians(properties)

    @unpack Pdot_Vdot, Pdot_Ωdot, Hdot_Vdot, Hdot_Ωdot = properties

    tmp = 1/2*Pdot_Vdot

    F1_V1dot = -1/2*tmp
    F2_V1dot =  1/2*tmp

    F1_V2dot = -1/2*tmp
    F2_V2dot =  1/2*tmp

    tmp = 1/2*Pdot_Ωdot

    F1_Ω1dot = -1/2*tmp
    F2_Ω1dot =  1/2*tmp

    F1_Ω2dot = -1/2*tmp
    F2_Ω2dot =  1/2*tmp

    tmp = 1/2*Hdot_Vdot

    M1_V1dot = -1/2*tmp
    M2_V1dot =  1/2*tmp

    M1_V2dot = -1/2*tmp
    M2_V2dot =  1/2*tmp

    tmp = 1/2*Hdot_Ωdot

    M1_Ω1dot = -1/2*tmp
    M2_Ω1dot =  1/2*tmp

    M1_Ω2dot = -1/2*tmp
    M2_Ω2dot =  1/2*tmp

    return (; 
        F1_V1dot, F1_V2dot, F1_Ω1dot, F1_Ω2dot, 
        F2_V1dot, F2_V2dot, F2_Ω1dot, F2_Ω2dot, 
        M1_V1dot, M1_V2dot, M1_Ω1dot, M1_Ω2dot, 
        M2_V1dot, M2_V2dot, M2_Ω1dot, M2_Ω2dot)
end

# --- Insert Element Residual --- #

@inline function insert_dynamic_element_jacobians!(jacob, indices, force_scaling,
    assembly, ielem, compatability, resultants)

    insert_static_element_jacobians!(jacob, indices, force_scaling, assembly, ielem, 
        compatability, resultants)

    @unpack ru_vb, ru_ωb, ru_V1, ru_V2, ru_Ω1, ru_Ω2, rθ_ωb, rθ_Ω1, rθ_Ω2 = compatability

    @unpack F1_θb, F1_ωb, F1_ab, F1_αb, F1_u1, F1_u2, F1_V1, F1_V2, F1_Ω1, F1_Ω2,
            F2_θb, F2_ωb, F2_ab, F2_αb, F2_u1, F2_u2, F2_V1, F2_V2, F2_Ω1, F2_Ω2,
            M1_θb, M1_vb, M1_ωb, M1_ab, M1_αb, M1_u1, M1_u2, M1_V1, M1_V2, M1_Ω1, M1_Ω2,
            M2_θb, M2_vb, M2_ωb, M2_ab, M2_αb, M2_u1, M2_u2, M2_V1, M2_V2, M2_Ω1, M2_Ω2 = resultants

    icol1 = indices.icol_point[assembly.start[ielem]]
    icol2 = indices.icol_point[assembly.stop[ielem]]

    # compatability equations
    irow = indices.irow_elem[ielem]

    jacob[irow:irow+2, 7:9] .= ru_vb 
    jacob[irow:irow+2, 10:12] .= ru_ωb 

    jacob[irow:irow+2, icol1+6:icol1+8] .= ru_V1
    jacob[irow:irow+2, icol1+9:icol1+11] .= ru_Ω1

    jacob[irow:irow+2, icol2+6:icol2+8] .= ru_V2
    jacob[irow:irow+2, icol2+9:icol2+11] .= ru_Ω2

    jacob[irow+3:irow+5, 10:12] .= rθ_ωb 

    jacob[irow+3:irow+5, icol1+9:icol1+11] .= rθ_Ω1

    jacob[irow+3:irow+5, icol2+9:icol2+11] .= rθ_Ω2

    # equilibrium equations for the start of the beam element
    irow1 = indices.irow_point[assembly.start[ielem]]

    @views jacob[irow1:irow1+2, 4:6] .-= F1_θb ./ force_scaling
    @views jacob[irow1:irow1+2, 10:12] .-= F1_ωb ./ force_scaling
    @views jacob[irow1:irow1+2, 13:15] .-= F1_ab ./ force_scaling
    @views jacob[irow1:irow1+2, 16:18] .-= F1_αb ./ force_scaling

    @views jacob[irow1:irow1+2, icol1:icol1+2] .-= F1_u1 ./ force_scaling

    @views jacob[irow1:irow1+2, icol2:icol2+2] .-= F1_u2 ./ force_scaling

    @views jacob[irow1:irow1+2, icol1+6:icol1+8] .-= F1_V1 ./ force_scaling
    @views jacob[irow1:irow1+2, icol1+9:icol1+11] .-= F1_Ω1 ./ force_scaling

    @views jacob[irow1:irow1+2, icol2+6:icol2+8] .-= F1_V2 ./ force_scaling
    @views jacob[irow1:irow1+2, icol2+9:icol2+11] .-= F1_Ω2 ./ force_scaling

    @views jacob[irow1+3:irow1+5, 4:6] .-= M1_θb ./ force_scaling
    @views jacob[irow1+3:irow1+5, 7:9] .-= M1_vb ./ force_scaling
    @views jacob[irow1+3:irow1+5, 10:12] .-= M1_ωb ./ force_scaling
    @views jacob[irow1+3:irow1+5, 13:15] .-= M1_ab ./ force_scaling
    @views jacob[irow1+3:irow1+5, 16:18] .-= M1_αb ./ force_scaling

    @views jacob[irow1+3:irow1+5, icol1:icol1+2] .-= M1_u1 ./ force_scaling
    
    @views jacob[irow1+3:irow1+5, icol2:icol2+2] .-= M1_u2 ./ force_scaling

    @views jacob[irow1+3:irow1+5, icol1+6:icol1+8] .-= M1_V1 ./ force_scaling
    @views jacob[irow1+3:irow1+5, icol1+9:icol1+11] .-= M1_Ω1 ./ force_scaling

    @views jacob[irow1+3:irow1+5, icol2+6:icol2+8] .-= M1_V2 ./ force_scaling
    @views jacob[irow1+3:irow1+5, icol2+9:icol2+11] .-= M1_Ω2 ./ force_scaling

    # equilibrium equations for the end of the beam element
    irow2 = indices.irow_point[assembly.stop[ielem]]

    @views jacob[irow2:irow2+2, 4:6] .+= F2_θb ./ force_scaling
    @views jacob[irow2:irow2+2, 10:12] .+= F2_ωb ./ force_scaling
    @views jacob[irow2:irow2+2, 13:15] .+= F2_ab ./ force_scaling
    @views jacob[irow2:irow2+2, 16:18] .+= F2_αb ./ force_scaling

    @views jacob[irow2:irow2+2, icol1:icol1+2] .+= F2_u1 ./ force_scaling

    @views jacob[irow2:irow2+2, icol2:icol2+2] .+= F2_u2 ./ force_scaling

    @views jacob[irow2:irow2+2, icol1+6:icol1+8] .+= F2_V1 ./ force_scaling
    @views jacob[irow2:irow2+2, icol1+9:icol1+11] .+= F2_Ω1 ./ force_scaling

    @views jacob[irow2:irow2+2, icol2+6:icol2+8] .+= F2_V2 ./ force_scaling
    @views jacob[irow2:irow2+2, icol2+9:icol2+11] .+= F2_Ω2 ./ force_scaling

    @views jacob[irow2+3:irow2+5, 4:6] .+= M2_θb ./ force_scaling
    @views jacob[irow2+3:irow2+5, 7:9] .+= M2_vb ./ force_scaling
    @views jacob[irow2+3:irow2+5, 10:12] .+= M2_ωb ./ force_scaling
    @views jacob[irow2+3:irow2+5, 13:15] .+= M2_ab ./ force_scaling
    @views jacob[irow2+3:irow2+5, 16:18] .+= M2_αb ./ force_scaling

    @views jacob[irow2+3:irow2+5, icol1:icol1+2] .+= M2_u1 ./ force_scaling

    @views jacob[irow2+3:irow2+5, icol2:icol2+2] .+= M2_u2 ./ force_scaling

    @views jacob[irow2+3:irow2+5, icol1+6:icol1+8] .+= M2_V1 ./ force_scaling
    @views jacob[irow2+3:irow2+5, icol1+9:icol1+11] .+= M2_Ω1 ./ force_scaling

    @views jacob[irow2+3:irow2+5, icol2+6:icol2+8] .+= M2_V2 ./ force_scaling
    @views jacob[irow2+3:irow2+5, icol2+9:icol2+11] .+= M2_Ω2 ./ force_scaling

    return jacob
end

@inline function insert_mass_matrix_element_jacobians!(jacob, gamma, indices, force_scaling,
    assembly, ielem, resultants)

    @unpack F1_V1dot, F1_V2dot, F1_Ω1dot, F1_Ω2dot, 
            F2_V1dot, F2_V2dot, F2_Ω1dot, F2_Ω2dot, 
            M1_V1dot, M1_V2dot, M1_Ω1dot, M1_Ω2dot, 
            M2_V1dot, M2_V2dot, M2_Ω1dot, M2_Ω2dot = resultants

    icol1 = indices.icol_point[assembly.start[ielem]]
    icol2 = indices.icol_point[assembly.stop[ielem]]
    
    # equilibrium equations for the start of the beam element
    irow1 = indices.irow_point[assembly.start[ielem]]

    @views jacob[irow1:irow1+2, icol1+6:icol1+8] .-= F1_V1dot .* gamma ./ force_scaling
    @views jacob[irow1:irow1+2, icol1+9:icol1+11] .-= F1_Ω1dot .* gamma ./ force_scaling

    @views jacob[irow1:irow1+2, icol2+6:icol2+8] .-= F1_V2dot .* gamma ./ force_scaling
    @views jacob[irow1:irow1+2, icol2+9:icol2+11] .-= F1_Ω2dot .* gamma ./ force_scaling

    @views jacob[irow1+3:irow1+5, icol1+6:icol1+8] .-= M1_V1dot .* gamma ./ force_scaling
    @views jacob[irow1+3:irow1+5, icol1+9:icol1+11] .-= M1_Ω1dot .* gamma ./ force_scaling

    @views jacob[irow1+3:irow1+5, icol2+6:icol2+8] .-= M1_V2dot .* gamma ./ force_scaling
    @views jacob[irow1+3:irow1+5, icol2+9:icol2+11] .-= M1_Ω2dot .* gamma ./ force_scaling

    # equilibrium equations for the end of the beam element
    irow2 = indices.irow_point[assembly.stop[ielem]]

    @views jacob[irow2:irow2+2, icol1+6:icol1+8] .+= F2_V1dot .* gamma ./ force_scaling
    @views jacob[irow2:irow2+2, icol1+9:icol1+11] .+= F2_Ω1dot .* gamma ./ force_scaling

    @views jacob[irow2:irow2+2, icol2+6:icol2+8] .+= F2_V2dot .* gamma ./ force_scaling
    @views jacob[irow2:irow2+2, icol2+9:icol2+11] .+= F2_Ω2dot .* gamma ./ force_scaling

    @views jacob[irow2+3:irow2+5, icol1+6:icol1+8] .+= M2_V1dot .* gamma ./ force_scaling
    @views jacob[irow2+3:irow2+5, icol1+9:icol1+11] .+= M2_Ω1dot .* gamma ./ force_scaling

    @views jacob[irow2+3:irow2+5, icol2+6:icol2+8] .+= M2_V2dot .* gamma ./ force_scaling
    @views jacob[irow2+3:irow2+5, icol2+9:icol2+11] .+= M2_Ω2dot .* gamma ./ force_scaling

    return jacob
end

# --- Element Residual --- #

"""
    dynamic_element_jacobian!(jacob, dx, x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, 
        gravity, ub, θb, vb, ωb, ab, αb)

Calculate and insert the jacobian entries corresponding to a beam element for a dynamic
analysis into the system jacobian matrix.
"""
@inline function dynamic_element_jacobian!(jacob, dx, x, indices, force_scaling, 
    structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, gravity, 
    ub, θb, vb, ωb, ab, αb, ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb)

    properties = dynamic_element_properties(dx, x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity, ub, θb, vb, ωb, ab, αb)

    properties = dynamic_element_jacobian_properties(properties, dx, x, indices, 
        force_scaling, structural_damping, assembly, ielem, prescribed_conditions, gravity,
        ub_ub, θb_θb, vb_vb, ωb_ωb, ab_ab, αb_αb)

    compatability = dynamic_compatability_jacobians(properties)

    resultants = dynamic_element_resultant_jacobians(properties, distributed_loads, ielem)

    insert_dynamic_element_jacobians!(jacob, indices, force_scaling, assembly, ielem, 
        compatability, resultants)

    return jacob
end

"""
    mass_matrix_element_jacobian!(jacob, gamma, x, indices, force_scaling, assembly, 
        ielem, prescribed_conditions)

Calculate and insert the mass_matrix jacobian entries corresponding to a beam element into 
the system jacobian matrix.
"""
@inline function mass_matrix_element_jacobian!(jacob, gamma, x, indices, force_scaling, assembly, 
    ielem, prescribed_conditions)

    properties = mass_matrix_element_jacobian_properties(x, indices, force_scaling, 
        assembly, ielem, prescribed_conditions)

    resultants = mass_matrix_element_resultant_jacobians(properties)
    
    insert_mass_matrix_element_jacobians!(jacob, gamma, indices, force_scaling, assembly, ielem, 
        resultants)

    return jacob
end
