# --- Element Properties --- #

"""
    steady_state_element_properties(x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, gravity, ub, θb, vb, ωb, ab, αb)

Calculate/extract the element properties needed to construct the residual for a steady 
state analysis
"""
@inline function steady_state_element_properties(x, indices, force_scaling, structural_damping, 
    assembly, ielem, prescribed_conditions, gravity, ub, θb, vb, ωb, ab, αb)

    properties = static_element_properties(x, indices, force_scaling, assembly, ielem, 
        prescribed_conditions, gravity)

    @unpack L, Cab, C, CtCab, Qinv, mass11, mass12, mass21, mass22, u1, u2, θ1, θ2, 
        u, θ, γ, κ = properties

    # rotation parameter matrices
    Q = get_Q(θ)

    # distance from the rotation center
    Δx = assembly.elements[ielem].x

    # (prescribed) linear and angular velocity
    v = vb + cross(ωb, Δx) + cross(ωb, u)# + udot = V
    ω = ωb# + Cab'*Q*θdot = Ω

    # (prescribed) linear and angular acceleration
    a = ab + cross(αb, Δx) + cross(αb, u)# + cross(ω, V) + uddot = Vdot
    α = αb# + Cab'*Qdot*θdot + Cab'*Q*θddot = Ωdot

    # add gravitational acceleration
    a -= get_C(θb)*gravity

    # linear and angular velocity (including body frame motion)
    V1, Ω1 = point_velocities(x, assembly.start[ielem], indices.icol_point)
    V2, Ω2 = point_velocities(x, assembly.stop[ielem], indices.icol_point)
    V = (V1 + V2)/2
    Ω = (Ω1 + Ω2)/2

    # linear and angular momentum (including body frame motion)
    P = CtCab*mass11*CtCab'*V + CtCab*mass12*CtCab'*Ω
    H = CtCab*mass21*CtCab'*V + CtCab*mass22*CtCab'*Ω

    # save steady state properties
    properties = (; properties..., Q, V1, V2, Ω1, Ω2, V, Ω, P, H, ub, θb, vb, ωb, ab, αb, 
        Δx, v, ω, a, α)

    if structural_damping

        # damping coefficients
        μ = assembly.elements[ielem].mu

        # damping submatrices
        μ11 = @SMatrix [μ[1] 0 0; 0 μ[2] 0; 0 0 μ[3]]
        μ22 = @SMatrix [μ[4] 0 0; 0 μ[5] 0; 0 0 μ[6]]

        # rotation parameter matrices
        C1 = get_C(θ1)
        C2 = get_C(θ2)
        Qinv1 = get_Qinv(θ1)
        Qinv2 = get_Qinv(θ2)

        # distance from the reference location
        Δx1 = assembly.points[assembly.start[ielem]]
        Δx2 = assembly.points[assembly.stop[ielem]]

        # linear and angular velocity
        v1 = vb + cross(ωb, Δx1) + cross(ωb, u1)
        v2 = vb + cross(ωb, Δx2) + cross(ωb, u2)
        ω1 = ωb
        ω2 = ωb

        # linear displacement rates 
        udot1 = V1 - v1
        udot2 = V2 - v2
        udot = (udot1 + udot2)/2

        # angular displacement rates
        θdot1 = Qinv1*C1*(Ω1 - ω1)
        θdot2 = Qinv2*C2*(Ω2 - ω2)
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
        γdot = -CtCab'*tilde(Ω - ω)*Δu + CtCab'*Δudot - L*CtCab'*tilde(Ω - ω)*Cab*e1
        κdot = Cab'*Q*Δθdot + Cab'*ΔQ*θdot

        # adjust strains to account for strain rates
        γ -= μ11*γdot
        κ -= μ22*κdot

        # save new strains and structural damping properties
        properties = (; properties..., γ, κ, γdot, κdot,
            μ11, μ22, C1, C2, Qinv1, Qinv2, Δx1, Δx2, v1, v2, ω1, ω2, 
            udot1, udot2, udot, θdot1, θdot2, θdot, Δu, Δθ, Δudot, Δθdot, ΔQ)
    end

    return properties
end

# --- Element Resultants --- #

"""
    steady_state_element_resultants(properties, distributed_loads, ielem)

Calculate the resultant loads applied at each end of a beam element for a steady_state 
analysis.
"""
@inline function steady_state_element_resultants(properties, distributed_loads, ielem)

    resultants = static_element_resultants(properties, distributed_loads, ielem)

    @unpack F1, F2, M1, M2 = resultants

    @unpack ω, V, Ω, P, H = properties

    # add loads due to linear and angular momentum
    tmp = cross(ω, P)
    F1 -= 1/2*tmp
    F2 += 1/2*tmp

    tmp = cross(ω, H) + cross(V, P)
    M1 -= 1/2*tmp
    M2 += 1/2*tmp

    return (; resultants..., F1, F2, M1, M2)
end

# --- Element Residual --- #

"""
    steady_state_element_residual!(resid, x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, distributed_loads, gravity, 
        ub, θb, vb, ωb, ab, αb)

Calculate and insert the residual entries corresponding to a beam element for a steady state 
analysis into the system residual vector.
"""
@inline function steady_state_element_residual!(resid, x, indices, force_scaling, 
    structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, gravity, 
    ub, θb, vb, ωb, ab, αb)

    properties = steady_state_element_properties(x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity, 
        ub, θb, vb, ωb, ab, αb)

    compatability = compatability_residuals(properties)

    resultants = steady_state_element_resultants(properties, distributed_loads, ielem)

    insert_element_residuals!(resid, indices, force_scaling, assembly, ielem, 
        compatability, resultants)

    return resid
end
