# --- Element Properties --- #

"""
    initial_condition_element_properties(x, indices, rate_vars, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity, 
        u0, θ0, V0, Ω0, Vdot0, Ωdot0)

Calculate/extract the element properties needed to construct the residual for a time-domain
analysis initialization
"""
@inline function initial_condition_element_properties(x, indices, rate_vars, 
    force_scaling, structural_damping, assembly, ielem, prescribed_conditions, gravity, 
    ub, θb, vb, ωb, ab, αb, u0, θ0, V0, Ω0, Vdot0, Ωdot0)

    # unpack element parameters
    @unpack L, Cab, compliance, mass = assembly.elements[ielem]

    # scale compliance and mass matrices by the element length
    compliance *= L
    mass *= L

    # compliance submatrices
    S11 = compliance[SVector{3}(1:3), SVector{3}(1:3)]
    S12 = compliance[SVector{3}(1:3), SVector{3}(4:6)]
    S21 = compliance[SVector{3}(4:6), SVector{3}(1:3)]
    S22 = compliance[SVector{3}(4:6), SVector{3}(4:6)]

    # mass submatrices
    mass11 = mass[SVector{3}(1:3), SVector{3}(1:3)]
    mass12 = mass[SVector{3}(1:3), SVector{3}(4:6)]
    mass21 = mass[SVector{3}(4:6), SVector{3}(1:3)]
    mass22 = mass[SVector{3}(4:6), SVector{3}(4:6)]

    # linear and angular displacement
    u1, θ1 = initial_point_displacement(x, assembly.start[ielem], indices.icol_point, 
        prescribed_conditions, u0, θ0, rate_vars)
    u2, θ2 = initial_point_displacement(x, assembly.stop[ielem], indices.icol_point, 
        prescribed_conditions, u0, θ0, rate_vars)
    u = (u1 + u2)/2
    θ = (θ1 + θ2)/2

    # rotation parameter matrices
    C = get_C(θ)
    CtCab = C'*Cab
    Q = get_Q(θ)
    Qinv = get_Qinv(θ)

    # forces and moments
    F, M = element_loads(x, ielem, indices.icol_elem, force_scaling)

    # strain and curvature
    γ = S11*F + S12*M
    κ = S21*F + S22*M

    # distance from the rotation center
    Δx = assembly.elements[ielem].x
    
    # (prescribed) linear and angular velocity
    v = vb + cross(ωb, Δx) + cross(ωb, u)# + udot = V
    ω = ωb# + Cab'*Q*θdot = Ω

    # (prescribed) linear and angular acceleration
    a = ab + cross(αb, Δx) + cross(αb, u)# + cross(ω, V) + uddot = Vdot
    α = αb# + Cab'*Qdot*θdot + Cab'*Q*θddot = Ωdot

    # linear and angular velocity **relative to the body frame**
    V1 = V0[assembly.start[ielem]]
    V2 = V0[assembly.stop[ielem]]
    V = (V1 + V2)/2

    Ω1 = Ω0[assembly.start[ielem]]
    Ω2 = Ω0[assembly.stop[ielem]]
    Ω = (Ω1 + Ω2)/2

    # linear and angular velocity (including body frame motion)
    V += v
    Ω += ω

    # linear and angular momentum (including body frame motion)
    P = CtCab*mass11*CtCab'*V + CtCab*mass12*CtCab'*Ω
    H = CtCab*mass21*CtCab'*V + CtCab*mass22*CtCab'*Ω

    # linear and angular acceleration **relative to the body frame**
    V1dot, Ω1dot = initial_point_velocity_rates(x, assembly.start[ielem], indices.icol_point, 
        prescribed_conditions, Vdot0, Ωdot0, rate_vars)
    V2dot, Ω2dot = initial_point_velocity_rates(x, assembly.stop[ielem], indices.icol_point, 
        prescribed_conditions, Vdot0, Ωdot0, rate_vars)
    Vdot = (V1dot + V2dot)/2
    Ωdot = (Ω1dot + Ω2dot)/2

    # linear and angular acceleration (including body frame motion)
    Vdot += a
    Ωdot += α

    # linear and angular momentum rates (including body frame motion)
    CtCabdot = tilde(Ω-ω)*CtCab

    Pdot = CtCab*mass11*CtCab'*Vdot + CtCab*mass12*CtCab'*Ωdot +
        CtCab*mass11*CtCabdot'*V + CtCab*mass12*CtCabdot'*Ω +
        CtCabdot*mass11*CtCab'*V + CtCabdot*mass12*CtCab'*Ω
    
    Hdot = CtCab*mass21*CtCab'*Vdot + CtCab*mass22*CtCab'*Ωdot +
        CtCab*mass21*CtCabdot'*V + CtCab*mass22*CtCabdot'*Ω +
        CtCabdot*mass21*CtCab'*V + CtCabdot*mass22*CtCab'*Ω

    # NOTE: All acceleration terms except gravity are included in Vdot and Ωdot.  We need to 
    # overwrite our acceleration terms so we don't double count them.
    a = -get_C(θb)*SVector{3}(gravity)
    α = @SVector zeros(3)

    # save properties
    properties = (; L, C, Cab, CtCab, Q, Qinv, S11, S12, S21, S22, mass11, mass12, mass21, mass22, 
        u1, u2, θ1, θ2, u, θ, F, M, γ, κ, V1, V2, Ω1, Ω2, V, Ω, P, H, ub, θb, vb, ωb, ab, αb, 
        Δx, v, ω, a, α, CtCabdot, Vdot, Ωdot, Pdot, Hdot) 

    if structural_damping 
        
        # damping coefficients
        μ = assembly.elements[ielem].mu

        # damping submatrices
        μ11 = @SMatrix [μ[1] 0 0; 0 μ[2] 0; 0 0 μ[3]]
        μ22 = @SMatrix [μ[4] 0 0; 0 μ[5] 0; 0 0 μ[6]]
    
        # linear and angular displacement rates
        udot1, θdot1 = point_velocities(x, assembly.start[ielem], indices.icol_point)
        udot2, θdot2 = point_velocities(x, assembly.stop[ielem], indices.icol_point)
        udot = (udot1 + udot2)/2
        θdot = (θdot1 + θdot2)/2

        # change in linear and angular displacement
        Δu = u2 - u1
        Δθ = θ2 - θ1

        # change in linear and angular displacement rates
        Δudot = udot2 - udot1
        Δθdot = θdot2 - θdot1

        # ΔQ matrix (see structural damping theory)
        ΔQ = get_ΔQ(θ, Δθ, Q)

        # strain rates
        γdot = -CtCab'*tilde(Ω-ω)*Δu + CtCab'*Δudot - L*CtCab'*tilde(Ω-ω)*Cab*e1
        κdot = Cab'*Q*Δθdot + Cab'*ΔQ*θdot

        # adjust strains to account for strain rates
        γ -= μ11*γdot
        κ -= μ22*κdot

        # add structural damping properties
        properties = (; properties..., γ, κ, γdot, κdot,
            μ11, μ22, udot1, udot2, udot, θdot1, θdot2, θdot, Δu, Δθ, Δudot, Δθdot, ΔQ)
    end

    return properties
end

# --- Element Residual --- #

"""
    initial_condition_element_residual!(resid, x, indices, rate_vars, 
        force_scaling, structural_damping, assembly, ielem, prescribed_conditions, 
        distributed_loads, gravity, ub, θb, vb, ωb, ab, αb, u0, θ0, V0, Ω0, Vdot0, Ωdot0)

Calculate and insert the residual entries corresponding to a beam element for the 
initialization of a time domain simulation into the system residual vector.
"""
@inline function initial_condition_element_residual!(resid, x, indices, rate_vars, 
    force_scaling, structural_damping, assembly, ielem, prescribed_conditions, 
    distributed_loads, gravity, ub, θb, vb, ωb, ab, αb, u0, θ0, V0, Ω0, Vdot0, Ωdot0)

    properties = initial_condition_element_properties(x, indices, rate_vars, 
        force_scaling, structural_damping, assembly, ielem, prescribed_conditions, gravity, ub, θb, vb, ωb, ab, αb, 
        u0, θ0, V0, Ω0, Vdot0, Ωdot0)

    compatability = compatability_residuals(properties)

    resultants = dynamic_element_resultants(properties, distributed_loads, ielem)

    insert_element_residuals!(resid, indices, force_scaling, assembly, ielem, 
        compatability, resultants)

    return resid
end

# --- Insert Element Residual --- #

@inline function insert_initial_condition_element_jacobians!(jacob, indices, force_scaling, 
    assembly, ielem, compatability, resultants)

    insert_static_element_jacobians!(jacob, indices, force_scaling, assembly, ielem, 
        compatability, resultants)

    @unpack ru_u1dot, ru_u2dot, rθ_θ1dot, rθ_θ2dot = compatability

    @unpack F1_θb, F1_vb, F1_ωb, F1_ab, F1_αb, F1_u1, F1_u2, F1_V1dot, F1_V2dot, F1_Ω1dot, F1_Ω2dot,
            F2_θb, F2_vb, F2_ωb, F2_ab, F2_αb, F2_u1, F2_u2, F2_V1dot, F2_V2dot, F2_Ω1dot, F2_Ω2dot,
            M1_θb, M1_vb, M1_ωb, M1_ab, M1_αb, M1_u1, M1_u2, M1_θ1, M1_θ2, M1_u1dot, M1_u2dot, M1_V1dot, M1_V2dot, M1_Ω1dot, M1_Ω2dot, 
            M2_θb, M2_vb, M2_ωb, M2_ab, M2_αb, M2_u1, M2_u2, M2_θ1, M2_θ2, M2_u1dot, M2_u2dot, M2_V1dot, M2_V2dot, M2_Ω1dot, M2_Ω2dot = resultants

    icol1 = indices.icol_point[assembly.start[ielem]]
    icol2 = indices.icol_point[assembly.stop[ielem]]

    # compatability equations
    irow = indices.irow_elem[ielem]

    jacob[irow:irow+2, icol1+6:icol1+8] .= ru_u1dot

    jacob[irow:irow+2, icol2+6:icol2+8] .= ru_u2dot

    jacob[irow+3:irow+5, icol1+9:icol1+11] .= rθ_θ1dot

    jacob[irow+3:irow+5, icol2+9:icol2+11] .= rθ_θ2dot

    # equilibrium equations for the start of the beam element
    irow1 = indices.irow_point[assembly.start[ielem]]

    @views jacob[irow1:irow1+2, 4:6] .-= F1_θb ./ force_scaling
    @views jacob[irow1:irow1+2, 7:9] .-= F1_vb ./ force_scaling
    @views jacob[irow1:irow1+2, 10:12] .-= F1_ωb ./ force_scaling
    @views jacob[irow1:irow1+2, 13:15] .-= F1_ab ./ force_scaling
    @views jacob[irow1:irow1+2, 16:18] .-= F1_αb ./ force_scaling

    @views jacob[irow1:irow1+2, icol1:icol1+2] .-= F1_u1 ./ force_scaling
    @views jacob[irow1:irow1+2, icol1:icol1+2] .-= F1_V1dot ./ force_scaling
    @views jacob[irow1:irow1+2, icol1+3:icol1+5] .-= F1_Ω1dot ./ force_scaling

    @views jacob[irow1:irow1+2, icol2:icol2+2] .-= F1_u2 ./ force_scaling
    @views jacob[irow1:irow1+2, icol2:icol2+2] .-= F1_V2dot ./ force_scaling
    @views jacob[irow1:irow1+2, icol2+3:icol2+5] .-= F1_Ω2dot ./ force_scaling

    @views jacob[irow1+3:irow1+5, 4:6] .-= M1_θb ./ force_scaling
    @views jacob[irow1+3:irow1+5, 7:9] .-= M1_vb ./ force_scaling
    @views jacob[irow1+3:irow1+5, 10:12] .-= M1_ωb ./ force_scaling
    @views jacob[irow1+3:irow1+5, 13:15] .-= M1_ab ./ force_scaling
    @views jacob[irow1+3:irow1+5, 16:18] .-= M1_αb ./ force_scaling

    @views jacob[irow1+3:irow1+5, icol1:icol1+2] .-= M1_u1 ./ force_scaling
    @views jacob[irow1+3:irow1+5, icol1:icol1+2] .-= M1_V1dot ./ force_scaling
    @views jacob[irow1+3:irow1+5, icol1+3:icol1+5] .-= M1_Ω1dot ./ force_scaling
    @views jacob[irow1+3:irow1+5, icol1+6:icol1+8] .-= M1_u1dot ./ force_scaling

    @views jacob[irow1+3:irow1+5, icol2:icol2+2] .-= M1_u2 ./ force_scaling
    @views jacob[irow1+3:irow1+5, icol2:icol2+2] .-= M1_V2dot ./ force_scaling
    @views jacob[irow1+3:irow1+5, icol2+3:icol2+5] .-= M1_Ω2dot ./ force_scaling
    @views jacob[irow1+3:irow1+5, icol2+6:icol2+8] .-= M1_u2dot ./ force_scaling

    # equilibrium equations for the end of the beam element
    irow2 = indices.irow_point[assembly.stop[ielem]]

    @views jacob[irow2:irow2+2, 4:6] .+= F2_θb ./ force_scaling
    @views jacob[irow2:irow2+2, 7:9] .+= F2_vb ./ force_scaling
    @views jacob[irow2:irow2+2, 10:12] .+= F2_ωb ./ force_scaling
    @views jacob[irow2:irow2+2, 13:15] .+= F2_ab ./ force_scaling
    @views jacob[irow2:irow2+2, 16:18] .+= F2_αb ./ force_scaling

    @views jacob[irow2:irow2+2, icol1:icol1+2] .+= F2_u1 ./ force_scaling
    @views jacob[irow2:irow2+2, icol1:icol1+2] .+= F2_V1dot ./ force_scaling
    @views jacob[irow2:irow2+2, icol1+3:icol1+5] .+= F2_Ω1dot ./ force_scaling

    @views jacob[irow2:irow2+2, icol2:icol2+2] .+= F2_u2 ./ force_scaling
    @views jacob[irow2:irow2+2, icol2:icol2+2] .+= F2_V2dot ./ force_scaling
    @views jacob[irow2:irow2+2, icol2+3:icol2+5] .+= F2_Ω2dot ./ force_scaling

    @views jacob[irow2+3:irow2+5, 4:6] .+= M2_θb ./ force_scaling
    @views jacob[irow2+3:irow2+5, 7:9] .+= M2_vb ./ force_scaling
    @views jacob[irow2+3:irow2+5, 10:12] .+= M2_ωb ./ force_scaling
    @views jacob[irow2+3:irow2+5, 13:15] .+= M2_ab ./ force_scaling
    @views jacob[irow2+3:irow2+5, 16:18] .+= M2_αb ./ force_scaling

    @views jacob[irow2+3:irow2+5, icol1:icol1+2] .+= M2_u1 ./ force_scaling
    @views jacob[irow2+3:irow2+5, icol1:icol1+2] .+= M2_V1dot ./ force_scaling
    @views jacob[irow2+3:irow2+5, icol1+3:icol1+5] .+= M2_Ω1dot ./ force_scaling
    @views jacob[irow2+3:irow2+5, icol1+6:icol1+8] .+= M2_u1dot ./ force_scaling

    @views jacob[irow2+3:irow2+5, icol2:icol2+2] .+= M2_u2 ./ force_scaling
    @views jacob[irow2+3:irow2+5, icol2:icol2+2] .+= M2_V2dot ./ force_scaling
    @views jacob[irow2+3:irow2+5, icol2+3:icol2+5] .+= M2_Ω2dot ./ force_scaling
    @views jacob[irow2+3:irow2+5, icol2+6:icol2+8] .+= M2_u2dot ./ force_scaling

    return jacob
end

