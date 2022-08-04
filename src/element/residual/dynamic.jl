# --- Element Properties --- #

"""
    dynamic_element_properties(dx, x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, gravity, ub, θb, vb, ωb, ab, αb)

Calculate/extract the element properties needed to construct the residual for a dynamic 
analysis
"""
@inline function dynamic_element_properties(dx, x, indices, force_scaling, 
    structural_damping, assembly, ielem, prescribed_conditions, gravity, ub, θb, vb, ωb, ab, αb)

    properties = steady_state_element_properties(x, indices, force_scaling, 
        structural_damping, assembly, ielem, prescribed_conditions, gravity, ub, θb, vb, ωb, ab, αb)

    @unpack C, Cab, CtCab, mass11, mass12, mass21, mass22, V1, V2, Ω1, Ω2, V, Ω, ω, θb = properties

    # linear and angular acceleration (including body frame motion)
    V1dot, Ω1dot = point_velocities(dx, assembly.start[ielem], indices.icol_point)
    V2dot, Ω2dot = point_velocities(dx, assembly.stop[ielem], indices.icol_point)
    Vdot = (V1dot + V2dot)/2
    Ωdot = (Ω1dot + Ω2dot)/2

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

    return (; properties..., CtCabdot, Vdot, Ωdot, Pdot, Hdot, a, α)
end

# --- Element Resultants --- #

"""
    dynamic_element_resultants(properties, distributed_loads, ielem)

Calculate the resultant loads applied at each end of a beam element for a dynamic 
analysis.
"""
@inline function dynamic_element_resultants(properties, distributed_loads, ielem)

    resultants = steady_state_element_resultants(properties, distributed_loads, ielem)

    @unpack F1, F2, M1, M2 = resultants

    @unpack Pdot, Hdot = properties
    
    # add loads due to linear and angular momentum rates  
    F1 -= 1/2*Pdot
    F2 += 1/2*Pdot

    M1 -= 1/2*Hdot
    M2 += 1/2*Hdot

    return (; resultants..., F1, F2, M1, M2)
end

# --- Element Residual --- #

"""
    dynamic_element_residual!(resid, dx, x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, distributed_loads, gravity, 
        ub, θb, vb, ωb, ab, αb)

Calculate and insert the residual entries corresponding to a beam element for a dynamic
analysis into the system residual vector.
"""
@inline function dynamic_element_residual!(resid, dx, x, indices, force_scaling, structural_damping, 
    assembly, ielem, prescribed_conditions, distributed_loads, gravity, ub, θb, vb, ωb, ab, αb)

    properties = dynamic_element_properties(dx, x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, gravity, ub, θb, vb, ωb, ab, αb)

    compatability = compatability_residuals(properties)

    resultants = dynamic_element_resultants(properties, distributed_loads, ielem)

    insert_element_residuals!(resid, indices, force_scaling, assembly, ielem, 
        compatability, resultants)

    return resid
end