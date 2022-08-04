# --- Helper Functions --- #

"""
    point_velocities(x, ipoint, icol_point)

Extract the velocities `V` and `Ω` of point `ipoint` from the state variable vector
"""
@inline function point_velocities(x, ipoint, icol_point)

    V, Ω = point_velocities(x, icol_point[ipoint])

    return V, Ω
end

@inline function point_velocities(x, icol)

    V = SVector(x[icol+6], x[icol+7], x[icol+8])
    Ω = SVector(x[icol+9], x[icol+10], x[icol+11])

    return V, Ω
end

# --- Point Properties --- #

"""
    steady_state_point_properties(x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity, ub, θb, vb, ωb, ab, αb)

Calculate/extract the point properties needed to construct the residual for a steady state 
analysis
"""
@inline function steady_state_point_properties(x, indices, force_scaling, assembly, ipoint, 
    prescribed_conditions, point_masses, gravity, ub, θb, vb, ωb, ab, αb)

    properties = static_point_properties(x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity)

    @unpack mass11, mass12, mass21, mass22, u, θ, C = properties

    # rotation parameter matrices
    Qinv = get_Qinv(θ)

    # distance from the rotation center (in the body frame)
    Δx = assembly.points[ipoint]

    # linear and angular velocity
    v = vb + cross(ωb, Δx) + cross(ωb, u)# + udot = V
    ω = ωb# + Cab'*Q*θdot = Ω

    # linear and angular acceleration
    a = ab + cross(αb, Δx) + cross(αb, u)# + cross(ω, V) + uddot = Vdot
    α = αb# + Cab'*Qdot*θdot + Cab'*Q*θddot = Ωdot

    # add gravitational acceleration
    a -= get_C(θb)*gravity

    # linear and angular velocity
    V, Ω = point_velocities(x, ipoint, indices.icol_point)

    # linear and angular momentum
    P = C'*mass11*C*V + C'*mass12*C*Ω
    H = C'*mass21*C*V + C'*mass22*C*Ω

    return (; properties..., Qinv, V, Ω, P, H, ub, θb, vb, ωb, ab, αb, Δx, v, ω, a, α)
end

# --- Point Resultants --- #

"""
    steady_state_point_resultants(properties)

Calculate the net loads `F` and `M` applied at a point for a steady_state analysis.
"""
@inline function steady_state_point_resultants(properties)

    F, M = static_point_resultants(properties)

    @unpack ω, V, Ω, P, H = properties

    # # add loads due to linear and angular momentum rates  
    F -= cross(ω, P)
    M -= cross(ω, H) + cross(V, P)

    return F, M
end

# --- Velocity Residuals --- #

"""
    steady_state_point_velocity_residuals(properties)

Calculate the velocity residuals `rV` and `rΩ` for a point for a steady state analysis.
"""
@inline function steady_state_point_velocity_residuals(properties)

    @unpack u, C, Qinv, V, Ω, v, ω = properties
    
    rV = V - v
    rΩ = Qinv*C*(Ω - ω)

    return (; rV, rΩ)
end

# --- Point Residual --- #

"""
    steady_state_point_residual!(resid, x, indices, force_scaling, assembly, ipoint, 
        prescribed_conditions, point_masses, gravity)

Calculate and insert the residual entries corresponding to a point for a steady state analysis into the 
system residual vector.
"""
@inline function steady_state_point_residual!(resid, x, indices, force_scaling, assembly, 
    ipoint, prescribed_conditions, point_masses, gravity, ub, θb, vb, ωb, ab, αb)

    irow = indices.irow_point[ipoint]

    properties = steady_state_point_properties(x, indices, force_scaling, assembly, ipoint, 
        prescribed_conditions, point_masses, gravity, ub, θb, vb, ωb, ab, αb)

    F, M = steady_state_point_resultants(properties)

    rV, rΩ = steady_state_point_velocity_residuals(properties)
       
    resid[irow:irow+2] .= -F ./ force_scaling
    resid[irow+3:irow+5] .= -M ./ force_scaling
    resid[irow+6:irow+8] .= rV
    resid[irow+9:irow+11] .= rΩ

    return resid
end