"""
    expanded_element_loads(x, ielem, icol_elem, force_scaling)

Extract the internal loads of a beam element at its endpoints (`F1`, `M1`, `F2`, `M2`) 
from the state variable vector for the expanded system.  These loads are expressed in the 
deformed element frame.
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

Extract the velocities of a beam element from the state variable vector of the expanded
system.  These velocities are expressed in the deformed frame of reference.
"""
@inline function expanded_element_velocities(x, ielem, icol_elem)

    icol = icol_elem[ielem]

    V = SVector(x[icol+12], x[icol+13], x[icol+14])
    Ω = SVector(x[icol+15], x[icol+16], x[icol+17])

    return V, Ω
end

"""
    expanded_velocity_residuals(properties)

Calculate the velocity residuals `rV` and `rΩ` for a beam element
"""
@inline function expanded_velocity_residuals(properties)

    @unpack C, C1, C2, V, V1, V2, Ω, Ω1, Ω2 = properties
    
    rV = C'*V - (C1'*V1 + C2'*V2)/2
    rΩ = C'*Ω - (C1'*Ω1 + C2'*Ω2)/2

    return (; rV, rΩ)
end

"""
    expanded_equilibrium_residuals(properties, distributed_loads, ielem)

Calculate the equilibrium residuals for a beam element for a static analysis.
"""
@inline function static_equilibrium_residuals(properties, distributed_loads, ielem)
  
    @unpack L, Cab, mass11, mass12, mass21, mass22, C, CtCab, F1, F2, M1, M2, F, M, γ, κ, a, α = properties

    P = mass11*V + mass12*Ω
    H = mass21*V + mass22*Ω

    Pdot = mass11*Vdot + mass12*Ωdot
    Hdot = mass21*Vdot + mass22*Ωdot

    # initialize equilibrium residual
    rF = Cab*F2 - Cab*F1
    rM = Cab*M2 - Cab*M1

    # add loads due to internal forces/moments and stiffness
    rM += Cab*cross(L*e1 + γ, F)

    # add loads due to linear and angular acceleration (including gravity)
    rF -= mass11*C*a + mass12*C*α
    rM -= mass21*C*a + mass22*C*α

    # add loads (in the deformed local frame) due to linear and angular momentum
    # Note that C*Cdot' = tilde(C'*Ω - ω)
    rF -= cross(C*ω, P) + tilde(C'*Ω - ω)*P
    rM -= cross(C*ω, H) + cross(V, P) + tilde(C'*Ω - ω)*H

    # add loads (in the deformed local frame) due to linear and angular momentum rates  
    rF -= Pdot
    rM -= Hdot

    # add distributed loads
    if haskey(distributed_loads, ielem)
        dload = distributed_loads[ielem]
        rF += C*dload.f1 + dload.f1_follower
        rF += C*dload.f2 + dload.f2_follower
        rM += C*dload.m1 + dload.m1_follower
        rM += C*dload.m2 + dload.m2_follower
    end

    return (; rF, rM)
end

"""
    expanded_resultants(properties)

Calculate the resultants at each end of the beam element for an expanded system
"""
@inline function expanded_element_resultants(properties)
   
    @unpack CtCab, C1, C2, F1, F2, M1, M2 = properties

    F1 = C1*CtCab*F1
    F2 = C2*CtCab*F2
    M1 = C1*CtCab*M1
    M2 = C2*CtCab*M2

    return (; F1, F2, M1, M2)
end

"""
    insert_expanded_element_residuals!(resid, indices, force_scaling, assembly, ielem, 
        compatability, equilibrium, resultants)

Insert the residual entries corresponding to a beam element into the system residual vector
for the expanded system.
"""
@inline function insert_expanded_element_residuals!(resid, indices, force_scaling, assembly, 
    ielem, compatability, equilibrium, velocities, resultants)

    @unpack ru, rθ = compatability
    @unpack rF, rM = equilibrium
    @unpack rV, rΩ = velocities
    @unpack F1, F2, M1, M2 = resultants

    # compatability and equilibrium equations
    irow = indices.irow_elem[ielem]
    resid[irow:irow+2] .= ru
    resid[irow+3:irow+5] .= rθ
    resid[irow+6:irow+8] .= rF ./ force_scaling
    resid[irow+9:irow+11] .= rM ./ force_scaling
    resid[irow+12:irow+14] .= rV
    resid[irow+15:irow+17] .= rΩ

    # equilibrium equations for the start of the beam element
    irow = indices.irow_point[assembly.start[ielem]]
    resid[irow:irow+2] .-= F1 ./ force_scaling
    resid[irow+3:irow+5] .-= M1 ./ force_scaling

    # equilibrium equations for the end of the beam element
    irow = indices.irow_point[assembly.stop[ielem]]
    resid[irow:irow+2] .+= F2 ./ force_scaling
    resid[irow+3:irow+5] .+= M2 ./ force_scaling

    return resid
end

@inline function element_resultant_jacobians(properties)
   
    @unpack Cab, C1, C2, CtCab, F1, F2, M1, M2, C1_θ1, C1_θ2, C1_θ3, 
        C2_θ1, C2_θ2, C2_θ3, C_θ1, C_θ2, C_θ3, θ1_θ1, θ2_θ2 = properties

    F1_θ = C1*mul3(C_θ1', C_θ2', C_θ3', Cab*F1)
    F1_θ1 = (mul3(C1_θ1, C1_θ2, C1_θ3, CtCab*F1) + 1/2*F1_θ)*θ1_θ1
    F1_θ2 = 1/2*F1_θ*θ2_θ2
    F1_F1 = C1*CtCab

    F2_θ = C2*mul3(C_θ1', C_θ2', C_θ3', Cab*F2)
    F2_θ1 = 1/2*F2_θ*θ1_θ1
    F2_θ2 = (mul3(C2_θ1, C2_θ2, C2_θ3, CtCab*F2) + 1/2*F2_θ)*θ2_θ2
    F2_F2 = C2*CtCab

    M1_θ = C1*mul3(C_θ1', C_θ2', C_θ3', Cab*M1)
    M1_θ1 = (mul3(C1_θ1, C1_θ2, C1_θ3, CtCab*M1) + 1/2*M1_θ)*θ1_θ1
    M1_θ2 = 1/2*M1_θ*θ2_θ2
    M1_M1 = C1*CtCab

    M2_θ = C2*mul3(C_θ1', C_θ2', C_θ3', Cab*M2)
    M2_θ1 = 1/2*M2_θ*θ1_θ1
    M2_θ2 = (mul3(C2_θ1, C2_θ2, C2_θ3, CtCab*M2) + 1/2*M2_θ)*θ2_θ2
    M2_M2 = C2*CtCab

    return (; F1_θ1, F1_θ2, F1_F1, F2_θ1, F2_θ2, F2_F2, 
        M1_θ1, M1_θ2, M1_M1, M2_θ1, M2_θ2, M2_M2)
end

@inline function element_velocity_jacobians(properties)

    @unpack C, C1, C2, V, V1, V2, Ω, Ω1, Ω2, 
        C_θ1, C_θ2, C_θ3, C1_θ1, C1_θ2, C1_θ3, C2_θ1, C2_θ2, C2_θ3, θ1_θ1, θ2_θ2 = properties
    
    rV_θ1 = 1/2*(mul3(C_θ1', C_θ2', C_θ3', V) - mul3(C1_θ1', C1_θ2', C1_θ3', V1))*θ1_θ1
    rV_θ2 = 1/2*(mul3(C_θ1', C_θ2', C_θ3', V) - mul3(C2_θ1', C2_θ2', C2_θ3', V2))*θ2_θ2
    rV_V = C'
    rV_V1 = -C1'/2
    rV_V2 = -C2'/2

    rΩ_θ1 = 1/2*(mul3(C_θ1', C_θ2', C_θ3', Ω) - mul3(C1_θ1', C1_θ2', C1_θ3', Ω1))*θ1_θ1
    rΩ_θ2 = 1/2*(mul3(C_θ1', C_θ2', C_θ3', Ω) - mul3(C2_θ1', C2_θ2', C2_θ3', Ω2))*θ2_θ2
    rΩ_Ω = C'
    rΩ_Ω1 = -C1'/2
    rΩ_Ω2 = -C2'/2

    return (; rV_θ1, rV_θ2, rV_V, rV_V1, rV_V2, rΩ_θ1, rΩ_θ2, rΩ_Ω, rΩ_Ω1, rΩ_Ω2)
end