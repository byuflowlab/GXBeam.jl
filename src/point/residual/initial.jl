# --- Helper Functions --- #

"""
    initial_point_displacement(x, ipoint, icol_point, prescribed_conditions, 
        rate_vars)

Extract the displacements `u` and `θ` of point `ipoint` from the state variable vector or 
prescribed conditions for an initial condition analysis.
"""
@inline function initial_point_displacement(x, ipoint, icol_point, 
    prescribed_conditions, u0, θ0, rate_vars)

    if haskey(prescribed_conditions, ipoint)
        u, θ = initial_point_displacement(x, icol_point[ipoint], 
            prescribed_conditions[ipoint], u0[ipoint], θ0[ipoint], rate_vars)
    else
        u, θ = initial_point_displacement(x, icol_point[ipoint], 
            u0[ipoint], θ0[ipoint], rate_vars)
    end

    return u, θ
end

@inline function initial_point_displacement(x, icol, prescribed_conditions, u0, θ0, rate_vars)

    # Use prescribed displacement, if applicable.  If no displacements are prescribed use 
    # a component of `u` or `θ` as state variables if the corresponding component of `Vdot` 
    # or `Ωdot` cannot be a state variable.

    # unpack prescribed conditions for the node
    @unpack pd, u, theta = prescribed_conditions

    # node linear displacement
    u = SVector(
        ifelse(pd[1], u[1], ifelse(rate_vars[icol+6], u0[1], x[icol])),
        ifelse(pd[2], u[2], ifelse(rate_vars[icol+7], u0[2], x[icol+1])),
        ifelse(pd[3], u[3], ifelse(rate_vars[icol+8], u0[3], x[icol+2]))
    )

    # node angular displacement
    θ = SVector(
        ifelse(pd[4], theta[1], ifelse(rate_vars[icol+9], θ0[1], x[icol+3])),
        ifelse(pd[5], theta[2], ifelse(rate_vars[icol+10], θ0[2], x[icol+4])),
        ifelse(pd[6], theta[3], ifelse(rate_vars[icol+11], θ0[3], x[icol+5]))
    )   

    return u, θ
end

@inline function initial_point_displacement(x, icol, u0, θ0, rate_vars)

    # Use a component of `u` or `θ` as state variables if the corresponding component of 
    # `Vdot` or `Ωdot` cannot be a state variable.

    u = SVector{3}(
        ifelse(rate_vars[icol+6], u0[1], x[icol]),
        ifelse(rate_vars[icol+7], u0[2], x[icol+1]),
        ifelse(rate_vars[icol+8], u0[3], x[icol+2])
    )
    θ = SVector{3}(
        ifelse(rate_vars[icol+9], θ0[1], x[icol+3]),
        ifelse(rate_vars[icol+10], θ0[2], x[icol+4]),
        ifelse(rate_vars[icol+11], θ0[3], x[icol+5])
    )

    return u, θ
end

"""
    initial_point_velocity_rates(x, ipoint, icol_point, prescribed_conditions, 
        Vdot0, Ωdot0, rate_vars)

Extract the velocity rates `Vdot` and `Ωdot` of point `ipoint` from the state variable 
vector or provided initial conditions.  Note that `Vdot` and `Ωdot` in this case do not 
include any contributions resulting from body frame motion. 
"""
@inline function initial_point_velocity_rates(x, ipoint, icol_point, prescribed_conditions, 
Vdot0, Ωdot0, rate_vars)

    if haskey(prescribed_conditions, ipoint)
        Vdot, Ωdot = initial_point_velocity_rates(x, icol_point[ipoint], prescribed_conditions[ipoint], 
            Vdot0[ipoint], Ωdot0[ipoint], rate_vars)
    else
        Vdot, Ωdot = initial_point_velocity_rates(x, icol_point[ipoint], 
            Vdot0[ipoint], Ωdot0[ipoint], rate_vars)
    end

    return Vdot, Ωdot
end

@inline function initial_point_velocity_rates(x, icol, prescribed_conditions, 
    Vdot0, Ωdot0, rate_vars)

    # If a displacment is prescribed, then the corresponding component of Vdot or Ωdot 
    # (relative to the body frame) is zero.  If no displacements is prescribed use the 
    # corresponding component of `Vdot` or `Ωdot` as a state variable, if possible.  
    # Otherwise, use the provided value.

    # unpack prescribed conditions for the node
    @unpack pd, u, theta = prescribed_conditions

    # node linear displacement
    Vdot = SVector(
        ifelse(pd[1], zero(eltype(x)), ifelse(rate_vars[icol+6], x[icol], Vdot0[1])),
        ifelse(pd[2], zero(eltype(x)), ifelse(rate_vars[icol+7], x[icol+1], Vdot0[2])),
        ifelse(pd[3], zero(eltype(x)), ifelse(rate_vars[icol+8], x[icol+2], Vdot0[3]))
    )

    # node angular displacement
    Ωdot = SVector(
        ifelse(pd[4], zero(eltype(x)), ifelse(rate_vars[icol+9], x[icol+3], Ωdot0[1])),
        ifelse(pd[5], zero(eltype(x)), ifelse(rate_vars[icol+10], x[icol+4], Ωdot0[2])),
        ifelse(pd[6], zero(eltype(x)), ifelse(rate_vars[icol+11], x[icol+5], Ωdot0[3]))
    )

    return Vdot, Ωdot
end

@inline function initial_point_velocity_rates(x, icol, Vdot0, Ωdot0, rate_vars)

    # Use the components of `Vdot` and `Ωdot` as state variables, if possible. Otherwise, 
    # use the provided value.

    Vdot = SVector{3}(
        ifelse(rate_vars[icol+6], x[icol], Vdot0[1]),
        ifelse(rate_vars[icol+7], x[icol+1], Vdot0[2]),
        ifelse(rate_vars[icol+8], x[icol+2], Vdot0[3])
    )
    Ωdot = SVector{3}(
        ifelse(rate_vars[icol+9], x[icol+3], Ωdot0[1]),
        ifelse(rate_vars[icol+10], x[icol+4], Ωdot0[2]),
        ifelse(rate_vars[icol+11], x[icol+5], Ωdot0[3])
    )

    return Vdot, Ωdot
end

# --- Point Properties --- #

"""
    initial_condition_point_properties(x, indices, rate_vars, force_scaling, 
        assembly, ipoint, prescribed_conditions, point_masses, gravity, 
        ub, θb, vb, ωb, ab, αb, u0, θ0, V0, Ω0, Vdot0, Ωdot0)

Calculate/extract the point properties needed to construct the residual for a time domain 
analysis initialization.
"""
@inline function initial_condition_point_properties(x, indices, rate_vars,
    force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity,
    ub, θb, vb, ωb, ab, αb, u0, θ0, V0, Ω0, Vdot0, Ωdot0)

    # mass matrix
    mass = haskey(point_masses, ipoint) ? point_masses[ipoint].mass : @SMatrix zeros(6,6)
    
    # mass submatrices
    mass11 = mass[SVector{3}(1:3), SVector{3}(1:3)]
    mass12 = mass[SVector{3}(1:3), SVector{3}(4:6)]
    mass21 = mass[SVector{3}(4:6), SVector{3}(1:3)]
    mass22 = mass[SVector{3}(4:6), SVector{3}(4:6)]

    # linear and angular displacement
    u, θ = initial_point_displacement(x, ipoint, indices.icol_point, 
        prescribed_conditions, u0, θ0, rate_vars)

    # linear and angular displacement rates
    udot, θdot = point_velocities(x, ipoint, indices.icol_point)

    # rotation parameter matrices
    C = get_C(θ)
    Qinv = get_Qinv(θ)

    # forces and moments
    F, M = point_loads(x, ipoint, indices.icol_point, force_scaling, prescribed_conditions)

    # distance from the rotation center
    Δx = assembly.points[ipoint]
    
    # linear and angular velocity
    v = vb + cross(ωb, Δx) + cross(ωb, u)# + udot = V
    ω = ωb# + Q*θdot = Ω

    # linear and angular acceleration
    a = ab + cross(αb, Δx) + cross(αb, u)# + cross(ω, V) + uddot = Vdot
    α = αb# + Qdot*θdot + Q*θddot = Ωdot

    # linear and angular velocity **excluding contributions from body frame motion**
    V = SVector{3}(V0[ipoint])
    Ω = SVector{3}(Ω0[ipoint])

    # add contributions from body frame motion to velocities
    V += v
    Ω += ω

    # linear and angular momentum
    P = C'*mass11*C*V + C'*mass12*C*Ω
    H = C'*mass21*C*V + C'*mass22*C*Ω

    # linear and angular acceleration
    Vdot, Ωdot = initial_point_velocity_rates(x, ipoint, indices.icol_point, 
        prescribed_conditions, Vdot0, Ωdot0, rate_vars)

    # add contributions from body frame motion to accelerations
    Vdot += a
    Ωdot += α

    # linear and angular acceleration
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

    return (; C, Qinv, mass11, mass12, mass21, mass22, u, θ, V, Ω, P, H, F, M, ub, θb, 
        vb, ωb, ab, αb, Δx, v, ω, a, α, udot, θdot, Cdot, Vdot, Ωdot, Pdot, Hdot) 
end

# --- Point Residual --- #

"""
    initial_condition_point_residual!(resid, x, indices, rate_vars,  
        force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity, ub, θb, vb, ωb, ab, αb, 
        u0, θ0, V0, Ω0, Vdot0, Ωdot0)

Calculate and insert the residual entries corresponding to a point for the initialization 
of a time domain analysis into the system residual vector.
"""
@inline function initial_condition_point_residual!(resid, x, indices, rate_vars,  
    force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity, ub, θb, vb, ωb, ab, αb, 
    u0, θ0, V0, Ω0, Vdot0, Ωdot0)

    irow = indices.irow_point[ipoint]

    properties = initial_condition_point_properties(x, indices, rate_vars, force_scaling, 
        assembly, ipoint, prescribed_conditions, point_masses, gravity, ub, θb, vb, ωb, ab, αb, 
        u0, θ0, V0, Ω0, Vdot0, Ωdot0)

    F, M = dynamic_point_resultants(properties)

    rV, rΩ = dynamic_point_velocity_residuals(properties)

    resid[irow:irow+2] .= -F ./ force_scaling
    resid[irow+3:irow+5] .= -M ./ force_scaling
    resid[irow+6:irow+8] .= rV
    resid[irow+9:irow+11] .= rΩ

    return resid
end