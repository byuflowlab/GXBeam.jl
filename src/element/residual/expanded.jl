# --- Helper Functions --- #

"""
    expanded_element_loads(x, ielem, icol_elem, force_scaling)

Extract the internal loads of a beam element at its endpoints (`F1`, `M1`, `F2`, `M2`) 
from the state variable vector for a constant mass matrix system.
"""
@inline function expanded_element_loads(x, ielem, icol_elem, force_scaling)

    icol = icol_elem[ielem]

    F1 = SVector(x[icol  ], x[icol+1], x[icol+2]) .* force_scaling
    M1 = SVector(x[icol+3], x[icol+4], x[icol+5]) .* force_scaling
    F2 = SVector(x[icol+6], x[icol+7], x[icol+8]) .* force_scaling
    M2 = SVector(x[icol+9], x[icol+10], x[icol+11]) .* force_scaling

    return F1, M1, F2, M2
end

"""
    expanded_element_velocities(x, ielem, icol_elem, force_scaling)

Extract the velocities of a beam element from the state variable vector for a constant mass
matrix system.
"""
@inline function expanded_element_velocities(x, ielem, icol_elem)

    icol = icol_elem[ielem]

    V = SVector(x[icol+12], x[icol+13], x[icol+14])
    Ω = SVector(x[icol+15], x[icol+16], x[icol+17])

    return V, Ω
end

# --- Element Properties --- # 

"""
    expanded_steady_element_properties(x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, gravity, ub, θb, vb, ωb, ab, αb)

Calculate/extract the element properties needed to construct the residual for a constant
mass matrix system
"""
@inline function expanded_steady_element_properties(x, indices, force_scaling, 
    structural_damping, assembly, ielem, prescribed_conditions, gravity, ub, θb, vb, ωb, ab, αb)

    # unpack element parameters
    @unpack L, Cab, compliance, mass = assembly.elements[ielem]

    # scale compliance and mass matrices by the element length
    compliance *= L
    mass *= L

    # compliance submatrices (in the deformed element frame)
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
    u1, θ1 = point_displacement(x, assembly.start[ielem], indices.icol_point, prescribed_conditions)
    u2, θ2 = point_displacement(x, assembly.stop[ielem], indices.icol_point, prescribed_conditions)
    u = (u1 + u2)/2
    θ = (θ1 + θ2)/2

    # rotation parameter matrices
    C = get_C(θ)
    CtCab = C'*Cab
    Q = get_Q(θ)
    Qinv = get_Qinv(θ)

    C1 = get_C(θ1)
    C2 = get_C(θ2)

    # forces and moments
    F1, M1, F2, M2 = expanded_element_loads(x, ielem, indices.icol_elem, force_scaling)
    F = (F1 + F2)/2
    M = (M1 + M2)/2

    # strain and curvature
    γ = S11*F + S12*M
    κ = S21*F + S22*M

    # distance from the rotation center
    Δx = assembly.elements[ielem].x

    # linear and angular velocity
    v = vb + cross(ωb, Δx) + cross(ωb, u)# + udot = CtCab*V
    ω = ωb# + Cab'*Q*θdot = CtCab*Ω

    # linear and angular acceleration
    a = ab + cross(αb, Δx) + cross(αb, u)# + cross(ω, V) + uddot = d/dt (CtCab*V)
    α = αb# + Cab'*Qdot*θdot + Cab'*Q*θddot = d/dt (CtCab*Ω)

    # add gravitational acceleration
    a -= get_C(θb)*gravity

    # linear and angular velocity (including body frame motion)
    V1, Ω1 = point_velocities(x, assembly.start[ielem], indices.icol_point)
    V2, Ω2 = point_velocities(x, assembly.stop[ielem], indices.icol_point)
    V, Ω = expanded_element_velocities(x, ielem, indices.icol_elem)

    # linear and angular momentum (including body frame motion)
    P = mass11*V + mass12*Ω
    H = mass21*V + mass22*Ω

    # element properties
    properties = (; L, C, C1, C2, Cab, CtCab, Q, Qinv, S11, S12, S21, S22, 
        mass11, mass12, mass21, mass22, u1, u2, θ1, θ2, V1, V2, Ω1, Ω2, F1, F2, M1, M2, 
        u, θ, V, Ω, P, H, F, M, γ, κ, ub, θb, vb, ωb, ab, αb, Δx, v, ω, a, α)

    if structural_damping

        # damping coefficients
        μ = assembly.elements[ielem].mu

        # damping submatrices
        μ11 = @SMatrix [μ[1] 0 0; 0 μ[2] 0; 0 0 μ[3]]
        μ22 = @SMatrix [μ[4] 0 0; 0 μ[5] 0; 0 0 μ[6]]

        # additional rotation parameter matrices
        C1 = get_C(θ1)
        C2 = get_C(θ2)
        Qinv1 = get_Qinv(θ1)
        Qinv2 = get_Qinv(θ2)

        # distance from rotation center
        Δx1 = assembly.points[assembly.start[ielem]]
        Δx2 = assembly.points[assembly.stop[ielem]]

        # linear and angular velocity
        v1 = vb + cross(ωb, Δx1) + cross(ωb, u1)
        v2 = vb + cross(ωb, Δx2) + cross(ωb, u2)
        ω1 = ωb
        ω2 = ωb

        # linear displacement rates 
        udot1 = C1'*V1 - v1
        udot2 = C2'*V2 - v2
        udot = (udot1 + udot2)/2

        # angular displacement rates
        θdot1 = Qinv1*(Ω1 - C1*ω1)
        θdot2 = Qinv2*(Ω2 - C2*ω2)
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
        γdot = -CtCab'*tilde(CtCab*Ω - ω)*Δu + CtCab'*Δudot - L*CtCab'*tilde(CtCab*Ω - ω)*Cab*e1
        κdot = Cab'*Q*Δθdot + Cab'*ΔQ*θdot

        # adjust strains to account for strain rates
        γ -= μ11*γdot
        κ -= μ22*κdot

        # add structural damping properties
        properties = (; properties..., γ, κ, γdot, κdot,
            μ11, μ22, C1, C2, Qinv1, Qinv2, Δx1, Δx2, v1, v2, ω1, ω2, 
            udot1, udot2, udot, θdot1, θdot2, θdot, Δu, Δθ, Δudot, Δθdot, ΔQ)
    end

    return properties
end

function expanded_dynamic_element_properties(x, indices, force_scaling, structural_damping, 
    assembly, ielem, prescribed_conditions, gravity, ub, θb, vb, ωb, ab, αb)

    properties = expanded_steady_element_properties(x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity, ub, θb, vb, ωb, ab, αb)

    @unpack θb, a, α = properties

    # NOTE: All acceleration terms except gravity are included in Vdot and Ωdot
    
    # overwrite our acceleration terms so we don't double count them.
    a = -get_C(θb)*SVector{3}(gravity)
    α = @SVector zeros(3)

    return (; properties..., a, α)
end

# --- Velocity Residuals --- #

"""
    expanded_element_velocity_residuals(properties)

Calculate the velocity residuals `rV` and `rΩ` of a beam element for a constant mass matrix 
system.
"""
@inline function expanded_element_velocity_residuals(properties)

    @unpack C, Cab, CtCab, Qinv, u, V, Ω, v, ω = properties
    
    rV = CtCab*V - v
    rΩ = Qinv*(Cab*Ω - C*ω)

    # @unpack CtCab, V, Ω, C1, V1, Ω1, C2, V2, Ω2 = properties
    # rV = CtCab*V - 1/2*(C1'*V1 + C2'*V2)
    # rΩ = CtCab*Ω - 1/2*(C1'*Ω1 + C2'*Ω2)

    return (; rV, rΩ)
end

# --- Equilibrium Residuals --- #

"""
    expanded_element_equilibrium_residuals(properties, distributed_loads, ielem)

Calculate the equilibrium residuals of a beam element for a constant mass matrix system.
"""
@inline function expanded_element_equilibrium_residuals(properties, distributed_loads, ielem)
  
    @unpack L, C, Cab, CtCab, mass11, mass12, mass21, mass22, F1, F2, M1, M2, F, M, γ, κ, 
        V, Ω, P, H, v, ω, a, α = properties

    # initialize equilibrium residual
    rF = F2 - F1
    rM = M2 - M1

    # add loads due to internal forces/moments and stiffness
    rM += cross(L*e1 + γ, F)

    # add loads due to linear and angular acceleration (including gravity)
    rF -= mass11*CtCab'*a + mass12*CtCab'*α
    rM -= mass21*CtCab'*a + mass22*CtCab'*α

    # add loads due to linear and angular momentum
    rF -= cross(Ω, P)
    rM -= cross(Ω, H) + cross(V, P)

    # add distributed loads
    if haskey(distributed_loads, ielem)
        dload = distributed_loads[ielem]
        f1 = CtCab'*dload.f1 + Cab'*dload.f1_follower
        f2 = CtCab'*dload.f2 + Cab'*dload.f2_follower
        rF += f1 + f2
        m1 = CtCab'*dload.m1 + Cab'*dload.m1_follower
        m2 = CtCab'*dload.m2 + Cab'*dload.m2_follower
        rM += m1 + m2
    end

    return (; rF, rM)
end

# --- Element Resultants --- #

"""
    expanded_element_resultants(properties)

Calculate the resultant loads applied at each end of a beam element for a constant mass
matrix system.
"""
@inline function expanded_element_resultants(properties)
   
    @unpack CtCab, C1, C2, F1, F2, M1, M2 = properties

    F1 = C1*CtCab*F1
    F2 = C2*CtCab*F2
    M1 = C1*CtCab*M1
    M2 = C2*CtCab*M2

    return (; F1, F2, M1, M2)
end

# --- Insert Residual --- #

"""
    insert_expanded_element_residuals!(resid, indices, force_scaling, assembly, ielem, 
        compatability, velocities, equilibrium, resultants)

Insert the residual entries corresponding to a beam element into the system residual vector
for a constant mass matrix system.
"""
@inline function insert_expanded_element_residuals!(resid, indices, force_scaling, assembly, 
    ielem, compatability, velocities, equilibrium, resultants)

    insert_element_residuals!(resid, indices, force_scaling, assembly, ielem, 
        compatability, resultants)

    @unpack rV, rΩ = velocities
    @unpack rF, rM = equilibrium

    irow = indices.irow_elem[ielem]

    # equilibrium residuals
    resid[irow+6:irow+8] .= rF ./ force_scaling
    resid[irow+9:irow+11] .= rM ./ force_scaling
    
    # velocity residuals
    resid[irow+12:irow+14] .= rV
    resid[irow+15:irow+17] .= rΩ

    return resid
end

# --- Element Residual --- #

"""
    expanded_steady_element_residual!(resid, x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, 
        gravity)

Calculate and insert the residual entries corresponding to a beam element for a constant
mass matrix system into the system residual vector.
"""
@inline function expanded_steady_element_residual!(resid, x, indices, force_scaling, 
    structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, gravity, 
    ub, θb, vb, ωb, ab, αb)

    properties = expanded_steady_element_properties(x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, gravity, ub, θb, vb, ωb, ab, αb)

    compatability = compatability_residuals(properties)

    velocities = expanded_element_velocity_residuals(properties)

    equilibrium = expanded_element_equilibrium_residuals(properties, distributed_loads, ielem)

    resultants = expanded_element_resultants(properties)

    insert_expanded_element_residuals!(resid, indices, force_scaling, assembly, ielem, 
        compatability, velocities, equilibrium, resultants)

    return resid
end

"""
    expanded_dynamic_element_residual!(resid, x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, 
        gravity)

Calculate and insert the residual entries corresponding to a beam element for a constant
mass matrix system into the system residual vector.
"""
@inline function expanded_dynamic_element_residual!(resid, x, indices, force_scaling, 
    structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, gravity, 
    ub, θb, vb, ωb, ab, αb)

    properties = expanded_dynamic_element_properties(x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, gravity, ub, θb, vb, ωb, ab, αb)

    compatability = compatability_residuals(properties)

    velocities = expanded_element_velocity_residuals(properties)

    equilibrium = expanded_element_equilibrium_residuals(properties, distributed_loads, ielem)

    resultants = expanded_element_resultants(properties)

    insert_expanded_element_residuals!(resid, indices, force_scaling, assembly, ielem, 
        compatability, velocities, equilibrium, resultants)

    return resid
end