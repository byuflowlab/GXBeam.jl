# --- Helper Functions --- #

"""
    point_displacement_rates(dx, ipoint, icol, prescribed_conditions)

Extract the displacement rates `udot` and `θdot` of point `ipoint` from the rate variable vector.
"""
@inline function point_displacement_rates(dx, ipoint, icol_point, prescribed_conditions)

    if haskey(prescribed_conditions, ipoint)
        udot, θdot = point_displacement_rates(dx, icol_point[ipoint], prescribed_conditions[ipoint])
    else
        udot, θdot = point_displacement_rates(dx, icol_point[ipoint])
    end

    return udot, θdot
end

@inline function point_displacement_rates(dx, icol, prescribed_conditions)

    # unpack prescribed conditions for the node
    @unpack pd = prescribed_conditions

    # node linear displacement rate
    udot = SVector(ifelse(pd[1], zero(eltype(dx)), dx[icol  ]),
                   ifelse(pd[2], zero(eltype(dx)), dx[icol+1]),
                   ifelse(pd[3], zero(eltype(dx)), dx[icol+2]))

    # node angular displacement rate
    θdot = SVector(ifelse(pd[4], zero(eltype(dx)), dx[icol+3]),
                   ifelse(pd[5], zero(eltype(dx)), dx[icol+4]),
                   ifelse(pd[6], zero(eltype(dx)), dx[icol+5]))

    return udot, θdot
end

@inline function point_displacement_rates(dx, icol)

    udot = SVector(dx[icol], dx[icol+1], dx[icol+2])
    θdot = SVector(dx[icol+3], dx[icol+4], dx[icol+5])

    return udot, θdot
end

# --- Point Properties --- #

"""
    dynamic_point_properties(dx, x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity, ub, θb, vb, ωb, ab, αb)

Calculate/extract the point properties needed to construct the residual for a dynamic 
analysis
"""
@inline function dynamic_point_properties(dx, x, indices, force_scaling, assembly, ipoint, 
    prescribed_conditions, point_masses, gravity, ub, θb, vb, ωb, ab, αb)

    properties = steady_state_point_properties(x, indices, force_scaling, assembly, ipoint, 
        prescribed_conditions, point_masses, gravity, ub, θb, vb, ωb, ab, αb)

    @unpack C, mass11, mass12, mass21, mass22, u, θ, V, Ω, ω, θb = properties

    # displacement rates
    udot, θdot = point_displacement_rates(dx, ipoint, indices.icol_point, prescribed_conditions)

    # velocity rates
    Vdot, Ωdot = point_velocities(dx, ipoint, indices.icol_point)

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

# --- Point Resultants --- #

"""
    dynamic_point_resultants(properties)

Calculate the net loads `F` and `M` applied at a point for a dynamic analysis.
"""
@inline function dynamic_point_resultants(properties)

    F, M = steady_state_point_resultants(properties)

    @unpack Pdot, Hdot = properties

    # add loads due to linear and angular momentum rates  
    F -= Pdot
    M -= Hdot

    return F, M
end

# --- Velocity Residuals --- #

"""
    dynamic_point_velocity_residuals(properties)

Calculate the velocity residuals `rV` and `rΩ` for a point for a dynamic analysis.
"""
@inline function dynamic_point_velocity_residuals(properties)

    residuals = steady_state_point_velocity_residuals(properties)

    @unpack rV, rΩ = residuals

    @unpack udot, θdot = properties
    
    rV -= udot
    rΩ -= θdot

    return (; rV, rΩ)
end

# --- Point Residual --- #

"""
    dynamic_point_residual!(resid, dx, x, indices, force_scaling, assembly, ipoint, 
        prescribed_conditions, point_masses, gravity, ub, θb, vb, ωb, ab, αb)

Calculate and insert the residual entries corresponding to a point for a dynamic analysis 
into the system residual vector.
"""
@inline function dynamic_point_residual!(resid, dx, x, indices, force_scaling, 
    assembly, ipoint, prescribed_conditions, point_masses, gravity, ub, θb, vb, ωb, ab, αb)

    irow = indices.irow_point[ipoint]

    properties = dynamic_point_properties(dx, x, indices, force_scaling, assembly, 
        ipoint, prescribed_conditions, point_masses, gravity, ub, θb, vb, ωb, ab, αb)

    F, M = dynamic_point_resultants(properties)

    rV, rΩ = dynamic_point_velocity_residuals(properties)

    resid[irow:irow+2] .= -F ./ force_scaling
    resid[irow+3:irow+5] .= -M ./ force_scaling
    resid[irow+6:irow+8] .= rV
    resid[irow+9:irow+11] .= rΩ

    return resid
end