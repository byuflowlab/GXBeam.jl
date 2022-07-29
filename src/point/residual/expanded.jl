# --- Point Properties --- #

"""
    expanded_steady_point_properties(x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity, ub, θb, vb, ωb, ab, αb)

Calculate/extract the point properties needed to construct the residual for a constant 
mass matrix system.
"""
@inline function expanded_steady_point_properties(x, indices, force_scaling, assembly, 
    ipoint, prescribed_conditions, point_masses, gravity, ub, θb, vb, ωb, ab, αb)

    # mass matrix
    mass = haskey(point_masses, ipoint) ? point_masses[ipoint].mass : @SMatrix zeros(6,6)
    
    # mass submatrices
    mass11 = mass[SVector{3}(1:3), SVector{3}(1:3)]
    mass12 = mass[SVector{3}(1:3), SVector{3}(4:6)]
    mass21 = mass[SVector{3}(4:6), SVector{3}(1:3)]
    mass22 = mass[SVector{3}(4:6), SVector{3}(4:6)]

    # linear and angular displacement
    u, θ = point_displacement(x, ipoint, indices.icol_point, prescribed_conditions)

    # rotation parameter matrices
    C = get_C(θ)
    Qinv = get_Qinv(θ)

    # forces and moments
    F, M = point_loads(x, ipoint, indices.icol_point, force_scaling, prescribed_conditions)

    # distance from the rotation center
    Δx = assembly.points[ipoint]

    # linear and angular velocity
    v = vb + cross(ωb, Δx) + cross(ωb, u)# + udot = C'*V
    ω = ωb# + Q*θdot = C'*Ω

    # linear and angular acceleration
    a = ab + cross(αb, Δx) + cross(αb, u)# + cross(ω, V) + uddot = d/dt (C'*V)
    α = αb# + Qdot*θdot + Q*θddot = d/dt (C'*Ω)

    # add gravitational acceleration
    a -= get_C(θb)*gravity

    # linear and angular velocity
    V, Ω = point_velocities(x, ipoint, indices.icol_point)

    # linear and angular momentum
    P = mass11*V + mass12*Ω
    H = mass21*V + mass22*Ω

    return (; C, Qinv, mass11, mass12, mass21, mass22, u, θ, V, Ω, P, H, F, M, ub, θb, 
        vb, ωb, ab, αb, Δx, v, ω, a, α) 
end

@inline function expanded_dynamic_point_properties(x, indices, force_scaling, assembly, 
    ipoint, prescribed_conditions, point_masses, gravity, ub, θb, vb, ωb, ab, αb)

    properties = expanded_steady_point_properties(x, indices, force_scaling, assembly, 
        ipoint, prescribed_conditions, point_masses, gravity, ub, θb, vb, ωb, ab, αb)

    @unpack θb, a, α = properties

    # NOTE: All acceleration terms except gravity are included in Vdot and Ωdot

    # overwrite acceleration terms so we don't double count them.
    a = -get_C(θb)*SVector{3}(gravity)
    α = @SVector zeros(3)

    return (; properties..., a, α)
end

# --- Point Resultants --- #

"""
    expanded_point_resultants(properties)

Calculate the net loads `F` and `M` applied at a point for a constant mass matrix system.
"""
@inline function expanded_point_resultants(properties)

    @unpack C, mass11, mass12, mass21, mass22, F, M, V, Ω, P, H, ω, a, α = properties

    # rotate loads into the deformed frame
    F = C*F
    M = C*M

    # add loads due to linear and angular acceleration (including gravity)
    F -= mass11*C*a + mass12*C*α
    M -= mass21*C*a + mass22*C*α

    # add loads due to linear and angular momentum
    F -= cross(Ω, P)
    M -= cross(Ω, H) + cross(V, P)

    return F, M
end

# --- Velocity Residuals --- #

"""
    expanded_point_velocity_residuals(properties)

Calculate the velocity residuals `rV` and `rΩ` for a point for a steady state analysis.
"""
@inline function expanded_point_velocity_residuals(properties)

    @unpack u, C, Qinv, V, Ω, v, ω = properties
    
    rV = C'*V - v
    rΩ = Qinv*(Ω - C*ω)

    return (; rV, rΩ)
end

# --- Point Residual --- #

"""
    expanded_steady_point_residual!(resid, x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity, ub, θb, vb, ωb, ab, αb)

Calculate and insert the residual entries corresponding to a point into the system residual 
vector for a constant mass matrix system.
"""
@inline function expanded_steady_point_residual!(resid, x, indices, force_scaling, 
    assembly, ipoint, prescribed_conditions, point_masses, gravity, ub, θb, vb, ωb, ab, αb)

    irow = indices.irow_point[ipoint]

    properties = expanded_steady_point_properties(x, indices, force_scaling, assembly, ipoint, 
        prescribed_conditions, point_masses, gravity, ub, θb, vb, ωb, ab, αb)

    F, M = expanded_point_resultants(properties)

    rV, rΩ = expanded_point_velocity_residuals(properties)

    resid[irow:irow+2] .= -F ./ force_scaling
    resid[irow+3:irow+5] .= -M ./ force_scaling
    resid[irow+6:irow+8] .= rV
    resid[irow+9:irow+11] .= rΩ

    return resid
end

"""
    expanded_dynamic_point_residual!(resid, x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity, ub, θb, vb, ωb, ab, αb)

Calculate and insert the residual entries corresponding to a point into the system residual 
vector for a constant mass matrix system.
"""
@inline function expanded_dynamic_point_residual!(resid, x, indices, force_scaling, 
    assembly, ipoint, prescribed_conditions, point_masses, gravity, ub, θb, vb, ωb, ab, αb)

    irow = indices.irow_point[ipoint]

    properties = expanded_dynamic_point_properties(x, indices, force_scaling, assembly, ipoint, 
        prescribed_conditions, point_masses, gravity, ub, θb, vb, ωb, ab, αb)

    F, M = expanded_point_resultants(properties)

    rV, rΩ = expanded_point_velocity_residuals(properties)

    resid[irow:irow+2] .= -F ./ force_scaling
    resid[irow+3:irow+5] .= -M ./ force_scaling
    resid[irow+6:irow+8] .= rV
    resid[irow+9:irow+11] .= rΩ

    return resid
end