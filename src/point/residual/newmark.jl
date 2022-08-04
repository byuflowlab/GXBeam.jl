# --- Helper Functions --- #

"""
    point_velocity_rates(x, ipoint, icol)

Extract the velocity rates `Vdot` and `Ωdot` of point `ipoint` from the state variable 
vector for the intialization of a time domain analysis.
"""
@inline function point_velocity_rates(x, ipoint, icol)

    Vdot, Ωdot = point_velocity_rates(x, icol[ipoint])

    return Vdot, Ωdot
end

@inline function point_velocity_rates(x, icol)

    Vdot = SVector(x[icol], x[icol+1], x[icol+2])
    Ωdot = SVector(x[icol+3], x[icol+4], x[icol+5])

    return Vdot, Ωdot
end

# --- Point Properties --- #

"""
    newmark_point_properties(x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity, ub, θb, vb, ωb, ab, αb, 
        udot_init, θdot_init, Vdot_init, Ωdot_init, dt)

Calculate/extract the point properties needed to construct the residual for a newmark-scheme
time stepping analysis
"""
@inline function newmark_point_properties(x, indices, force_scaling, assembly, ipoint, 
    prescribed_conditions, point_masses, gravity, ub, θb, vb, ωb, ab, αb,
    udot_init, θdot_init, Vdot_init, Ωdot_init, dt)

    properties = steady_state_point_properties(x, indices, force_scaling, assembly, ipoint, 
        prescribed_conditions, point_masses, gravity, ub, θb, vb, ωb, ab, αb)

    @unpack C, mass11, mass12, mass21, mass22, u, θ, V, Ω, ω, θb = properties

    # linear and angular displacement rates
    udot = 2/dt*u - SVector{3}(udot_init[ipoint])
    θdot = 2/dt*θ - SVector{3}(θdot_init[ipoint])

    # linear and angular velocity rates (in the deformed local frame)
    Vdot = 2/dt*V - SVector{3}(Vdot_init[ipoint])
    Ωdot = 2/dt*Ω - SVector{3}(Ωdot_init[ipoint])

    # linear and angular momentum rates
    Cdot = -C*tilde(Ω - ω)

    Pdot = C'*mass11*C*Vdot + C'*mass12*C*Ωdot +
        C'*mass11*Cdot*V + C'*mass12*Cdot*Ω +
        Cdot'*mass11*C*V + Cdot'*mass12*C*Ω
    
    Hdot = C'*mass21*C*Vdot + C'*mass22*C*Ωdot +
        C'*mass21*Cdot*V + C'*mass22*Cdot*Ω +
        Cdot'*mass21*C*V + Cdot'*mass22*C*Ω

    # overwrite acceleration terms so we don't double count them
    a = -get_C(θb)*SVector{3}(gravity)
    α = @SVector zeros(3)

    return (; properties..., udot, θdot, Cdot, Vdot, Ωdot, Pdot, Hdot, a, α) 
end

# --- Point Residual --- #

"""
    newmark_point_residual!(resid, x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity, ub, θb, vb, ωb, ab, αb,
        udot_init, θdot_init, Vdot_init, Ωdot_init, dt)

Calculate and insert the residual entries corresponding to a point for a newmark-scheme 
time marching analysis into the system residual vector.
"""
@inline function newmark_point_residual!(resid, x, indices, force_scaling, assembly, ipoint, 
    prescribed_conditions, point_masses, gravity, ub, θb, vb, ωb, ab, αb,
    udot_init, θdot_init, Vdot_init, Ωdot_init, dt)

    irow = indices.irow_point[ipoint]

    properties = newmark_point_properties(x, indices, force_scaling, assembly, ipoint, 
        prescribed_conditions, point_masses, gravity, ub, θb, vb, ωb, ab, αb,
        udot_init, θdot_init, Vdot_init, Ωdot_init, dt)

    F, M = dynamic_point_resultants(properties)

    rV, rΩ = dynamic_point_velocity_residuals(properties)

    resid[irow:irow+2] .= -F ./ force_scaling
    resid[irow+3:irow+5] .= -M ./ force_scaling
    resid[irow+6:irow+8] .= rV
    resid[irow+9:irow+11] .= rΩ

    return resid
end