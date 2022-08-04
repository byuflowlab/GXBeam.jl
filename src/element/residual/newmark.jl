# --- Element Properties --- #

"""
    newmark_element_properties(x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, gravity, ub, θb, vb, ωb, ab, αb, 
        Vdot_init, Ωdot_init, dt)

Calculate/extract the element properties needed to construct the residual for a newmark-
scheme time stepping analysis
"""
@inline function newmark_element_properties(x, indices, force_scaling, structural_damping, 
    assembly, ielem, prescribed_conditions, gravity, ub, θb, vb, ωb, ab, αb, 
    Vdot_init, Ωdot_init, dt)

    properties = steady_state_element_properties(x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, gravity, ub, θb, vb, ωb, ab, αb)

    @unpack C, Cab, CtCab, mass11, mass12, mass21, mass22, V1, V2, Ω1, Ω2, V, Ω, ω, θb = properties

    # linear and angular acceleration (including body frame motion)
    V1dot = 2/dt*V1 - SVector{3}(Vdot_init[assembly.start[ielem]])
    Ω1dot = 2/dt*Ω1 - SVector{3}(Ωdot_init[assembly.start[ielem]])

    V2dot = 2/dt*V2 - SVector{3}(Vdot_init[assembly.stop[ielem]])
    Ω2dot = 2/dt*Ω2 - SVector{3}(Ωdot_init[assembly.stop[ielem]])

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

# --- Element Residual --- #

"""
    newmark_element_residual!(resid, x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, distributed_loads, gravity, 
        ub, θb, vb, ωb, ab, αb, Vdot_init, Ωdot_init, dt)

Calculate and insert the residual entries corresponding to a beam element for a 
newmark-scheme time marching analysis into the system residual vector.
"""
@inline function newmark_element_residual!(resid, x, indices, force_scaling, structural_damping, 
    assembly, ielem, prescribed_conditions, distributed_loads, gravity, 
    ub, θb, vb, ωb, ab, αb, Vdot_init, Ωdot_init, dt)

    properties = newmark_element_properties(x, indices, force_scaling, structural_damping, 
        assembly, ielem, prescribed_conditions, gravity, 
        ub, θb, vb, ωb, ab, αb, Vdot_init, Ωdot_init, dt)

    compatability = compatability_residuals(properties)

    resultants = dynamic_element_resultants(properties, distributed_loads, ielem)

    insert_element_residuals!(resid, indices, force_scaling, assembly, ielem, 
        compatability, resultants)

    return resid
end