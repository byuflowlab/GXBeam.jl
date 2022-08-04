# --- Helper Functions --- #

"""
    element_loads(x, ielem, icol_elem, force_scaling)

Extract the internal loads (`F`, `M`) for a beam element from the state variable vector.  
These loads are expressed in the deformed element frame.
"""
@inline function element_loads(x, ielem, icol_elem, force_scaling)

    icol = icol_elem[ielem]

    F = SVector(x[icol  ], x[icol+1], x[icol+2]) .* force_scaling
    M = SVector(x[icol+3], x[icol+4], x[icol+5]) .* force_scaling

    return F, M
end

# --- Element Properties --- #

"""
    static_element_properties(x, indices, force_scaling, assembly, ielem, 
        prescribed_conditions, gravity)

Calculate/extract the element properties needed to construct the residual for a static 
analysis
"""
@inline function static_element_properties(x, indices, force_scaling, assembly, ielem, 
    prescribed_conditions, gravity)

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
    Qinv = get_Qinv(θ)

    # forces and moments
    F, M = element_loads(x, ielem, indices.icol_elem, force_scaling)
    
    # strain and curvature
    γ = S11*F + S12*M
    κ = S21*F + S22*M

    # linear and angular acceleration
    a = -SVector{3}(gravity)
    α = @SVector zeros(3)

    return (; L, C, Cab, CtCab, Qinv, S11, S12, S21, S22, mass11, mass12, mass21, mass22, 
        u1, u2, θ1, θ2, u, θ, F, M, γ, κ, a, α)
end

# --- Compatability Residuals --- #

"""
    compatability_residuals(properties)

Calculate the compatability residuals for the beam element
"""
@inline function compatability_residuals(properties)
   
    @unpack L, Cab, CtCab, Qinv, u1, u2, θ1, θ2, γ, κ = properties

    Δu = CtCab*(L*e1 + γ) - L*Cab*e1

    Δθ = Qinv*Cab*κ

    ru = u2 - u1 - Δu
    rθ = θ2 - θ1 - Δθ

    return (; ru, rθ)
end

# --- Element Resultants --- #

"""
    static_element_resultants(properties, distributed_loads, ielem)

Calculate the resultant loads applied at each end of a beam element for a static analysis.
"""
@inline function static_element_resultants(properties, distributed_loads, ielem)
  
    @unpack L, C, Cab, CtCab, mass11, mass12, mass21, mass22, F, M, γ, κ, a, α = properties

    # add loads due to internal forces/moments and stiffness
    tmp = CtCab*F
    F1 = tmp
    F2 = tmp

    tmp1 = CtCab*M
    tmp2 = 1/2*CtCab*cross(L*e1 + γ, F)
    M1 = tmp1 + tmp2
    M2 = tmp1 - tmp2

    # add loads due to an accelerating reference frame
    tmp = CtCab*mass11*CtCab'*a + CtCab*mass12*CtCab'*α
    F1 -= 1/2*tmp
    F2 += 1/2*tmp

    tmp = CtCab*mass21*CtCab'*a + CtCab*mass22*CtCab'*α
    M1 -= 1/2*tmp
    M2 += 1/2*tmp

    # add distributed loads
    if haskey(distributed_loads, ielem)
        dload = distributed_loads[ielem]
        F1 += dload.f1 + C'*dload.f1_follower
        F2 -= dload.f2 + C'*dload.f2_follower
        M1 += dload.m1 + C'*dload.m1_follower
        M2 -= dload.m2 + C'*dload.m2_follower
    end  
   
    return (; F1, F2, M1, M2)
end

# --- Insert Residual --- #

"""
    insert_element_residuals!(resid, indices, force_scaling, assembly, ielem, 
        compatability, resultants)

Insert the residual entries corresponding to a beam element into the system residual vector.
"""
@inline function insert_element_residuals!(resid, indices, force_scaling, assembly, ielem, 
    compatability, resultants)

    @unpack ru, rθ = compatability
    @unpack F1, F2, M1, M2 = resultants

    # compatability equations
    irow = indices.irow_elem[ielem]
    resid[irow:irow+2] .= ru
    resid[irow+3:irow+5] .= rθ

    # equilibrium equations for the start of the beam element
    irow = indices.irow_point[assembly.start[ielem]]
    @views resid[irow:irow+2] .-= F1 ./ force_scaling
    @views resid[irow+3:irow+5] .-= M1 ./ force_scaling

    # equilibrium equations for the end of the beam element
    irow = indices.irow_point[assembly.stop[ielem]]
    @views resid[irow:irow+2] .+= F2 ./ force_scaling
    @views resid[irow+3:irow+5] .+= M2 ./ force_scaling

    return resid
end

# --- Element Residual --- #

"""
    static_element_residual!(resid, x, indices, force_scaling, assembly, ielem,  
        prescribed_conditions, distributed_loads, gravity)

Calculate and insert the residual entries corresponding to a beam element for a static 
analysis into the system residual vector.
"""
@inline function static_element_residual!(resid, x, indices, force_scaling, assembly, ielem,  
    prescribed_conditions, distributed_loads, gravity)

    properties = static_element_properties(x, indices, force_scaling, assembly, ielem, 
        prescribed_conditions, gravity)

    compatability = compatability_residuals(properties)

    resultants = static_element_resultants(properties, distributed_loads, ielem)

    insert_element_residuals!(resid, indices, force_scaling, assembly, ielem, 
        compatability, resultants)

    return resid
end
