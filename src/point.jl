"""
    point_loads(x, ipoint, icol, force_scaling, prescribed_conditions)

Extract the loads `F` and `M` of point `ipoint` from the state variable vector or 
prescribed conditions.
"""
@inline function point_loads(x, ipoint, icol, force_scaling, prescribed_conditions)

    if haskey(prescribed_conditions, ipoint)
        F, M = point_loads(x, icol[ipoint], force_scaling, prescribed_conditions[ipoint])
    else
        F, M = point_loads(x, icol[ipoint], force_scaling)
    end

    return F, M
end

@inline function point_loads(x, icol, force_scaling, prescribed_conditions)

    @unpack pl, F, M, Ff, Mf = prescribed_conditions

    # combine non-follower and follower loads
    _, θ = point_displacement(x, icol, prescribed_conditions)
    C = get_C(θ)
    F = F + C'*Ff
    M = M + C'*Mf

    # set state variable loads explicitly
    F = SVector(ifelse(pl[1], F[1], x[icol]*force_scaling),
                ifelse(pl[2], F[2], x[icol+1]*force_scaling),
                ifelse(pl[3], F[3], x[icol+2]*force_scaling))
    M = SVector(ifelse(pl[4], M[1], x[icol+3]*force_scaling),
                ifelse(pl[5], M[2], x[icol+4]*force_scaling),
                ifelse(pl[6], M[3], x[icol+5]*force_scaling))
    
    return F, M
end

@inline function point_loads(x, icol, force_scaling)

    F = @SVector zeros(eltype(x), 3)
    M = @SVector zeros(eltype(x), 3)

    return F, M
end

"""
    point_load_jacobians(x, ipoint, icol, force_scaling, prescribed_conditions)

Calculate the load jacobians `F_θ`, `F_F`, `M_θ`, and `M_M` of point `ipoint`.
"""
@inline function point_load_jacobians(x, ipoint, icol, force_scaling, prescribed_conditions)

    if haskey(prescribed_conditions, ipoint)
        F_θ, F_F, M_θ, M_M = point_load_jacobians(x, icol[ipoint], force_scaling, prescribed_conditions[ipoint])
    else
        F_θ, F_F, M_θ, M_M = point_load_jacobians(x, icol[ipoint], force_scaling)
    end

    return F_θ, F_F, M_θ, M_M
end

@inline function point_load_jacobians(x, icol, force_scaling, prescribed_conditions)

    @unpack pl, Ff, Mf = prescribed_conditions

    _, θ = point_displacement(x, icol, prescribed_conditions)
    C = get_C(θ)

    _, θ_θ = point_displacement_jacobians(prescribed_conditions)
    C_θ1, C_θ2, C_θ3 = get_C_θ(C, θ)

    # solve for the jacobian wrt theta of the follower forces
    F_θ = @SMatrix zeros(eltype(x), 3, 3)
    for i = 1:3
        rot_θ = @SMatrix [
            C_θ1[i,1] C_θ2[i,1] C_θ3[i,1];
            C_θ1[i,2] C_θ2[i,2] C_θ3[i,2];
            C_θ1[i,3] C_θ2[i,3] C_θ3[i,3]
            ]
        F_θ += rot_θ*Ff[i]
    end

    F1_θ = ifelse(pl[1], SVector(F_θ[1,1], F_θ[1,2], F_θ[1,3]), @SVector zeros(eltype(x), 3))
    F2_θ = ifelse(pl[2], SVector(F_θ[2,1], F_θ[2,2], F_θ[2,3]), @SVector zeros(eltype(x), 3))
    F3_θ = ifelse(pl[3], SVector(F_θ[3,1], F_θ[3,2], F_θ[3,3]), @SVector zeros(eltype(x), 3))
    F_θ = vcat(F1_θ', F2_θ', F3_θ')

    # solve for the jacobian wrt theta of the follower moments
    M_θ = @SMatrix zeros(eltype(x), 3, 3)
    for i = 1:3
        rot_θ = @SMatrix [
            C_θ1[i,1] C_θ2[i,1] C_θ3[i,1];
            C_θ1[i,2] C_θ2[i,2] C_θ3[i,2];
            C_θ1[i,3] C_θ2[i,3] C_θ3[i,3]
            ]
        M_θ += rot_θ*Mf[i]
    end
    
    M1_θ = ifelse(pl[4], SVector(M_θ[1,1], M_θ[1,2], M_θ[1,3]), @SVector zeros(eltype(x), 3))
    M2_θ = ifelse(pl[5], SVector(M_θ[2,1], M_θ[2,2], M_θ[2,3]), @SVector zeros(eltype(x), 3))
    M3_θ = ifelse(pl[6], SVector(M_θ[3,1], M_θ[3,2], M_θ[3,3]), @SVector zeros(eltype(x), 3))
    M_θ = vcat(M1_θ', M2_θ', M3_θ')

    F_F = hcat(ifelse(pl[1], zero(e1), e1),
               ifelse(pl[2], zero(e2), e2),
               ifelse(pl[3], zero(e3), e3))

    M_M = hcat(ifelse(pl[4], zero(e1), e1),
               ifelse(pl[5], zero(e2), e2),
               ifelse(pl[6], zero(e3), e3))

    return F_θ, F_F, M_θ, M_M
end

@inline function point_load_jacobians(x, icol, force_scaling)

    F_θ = @SMatrix zeros(eltype(x), 3, 3)
    F_F = @SMatrix zeros(eltype(x), 3, 3)
    M_θ = @SMatrix zeros(eltype(x), 3, 3)
    M_M = @SMatrix zeros(eltype(x), 3, 3)

    return F_θ, F_F, M_θ, M_M
end

"""
    point_displacement(x, ipoint, icol_point, prescribed_conditions)

Extract the displacements `u` and `θ` of point `ipoint` from the state variable vector or 
prescribed conditions.
"""
@inline function point_displacement(x, ipoint, icol_point, prescribed_conditions)

    if haskey(prescribed_conditions, ipoint)
        u, θ = point_displacement(x, icol_point[ipoint], prescribed_conditions[ipoint])
    else
        u, θ = point_displacement(x, icol_point[ipoint])
    end

    return u, θ
end

@inline function point_displacement(x, icol, prescribed_conditions)

    # unpack prescribed conditions for the node
    @unpack pd, u, theta = prescribed_conditions

    # node linear displacement
    u = SVector(
        ifelse(pd[1], u[1], x[icol]),
        ifelse(pd[2], u[2], x[icol+1]),
        ifelse(pd[3], u[3], x[icol+2])
    )

    # node angular displacement
    θ = SVector(
        ifelse(pd[4], theta[1], x[icol+3]),
        ifelse(pd[5], theta[2], x[icol+4]),
        ifelse(pd[6], theta[3], x[icol+5])
    )

    return u, θ
end

@inline function point_displacement(x, icol)

    u = SVector(x[icol], x[icol+1], x[icol+2])
    θ = SVector(x[icol+3], x[icol+4], x[icol+5])

    return u, θ
end

"""
    point_displacement_jacobians(ipoint, prescribed_conditions)

Calculate the displacement jacobians `u_u` and `θ_θ` of point `ipoint`.
"""
@inline function point_displacement_jacobians(ipoint, prescribed_conditions)

    if haskey(prescribed_conditions, ipoint)
        u_u, θ_θ = point_displacement_jacobians(prescribed_conditions[ipoint])
    else
        u_u, θ_θ = point_displacement_jacobians()
    end

    return u_u, θ_θ
end

@inline function point_displacement_jacobians(prescribed_conditions)

    @unpack pd = prescribed_conditions

    u_u = hcat(ifelse(pd[1], zero(e1), e1),
               ifelse(pd[2], zero(e2), e2),
               ifelse(pd[3], zero(e3), e3))

    θ_θ = hcat(ifelse(pd[4], zero(e1), e1),
               ifelse(pd[5], zero(e2), e2),
               ifelse(pd[6], zero(e3), e3))

    return u_u, θ_θ
end

@inline function point_displacement_jacobians()

    u_u = I3
    θ_θ = I3

    return u_u, θ_θ
end

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
    initial_point_displacement_jacobian(ipoint, icol_point, prescribed_conditions, 
        rate_vars)

Extract the displacement jacobians `u_u` and `θ_θ` of point `ipoint` from the state 
variable vector or prescribed conditions for an initial condition analysis.
"""
@inline function initial_point_displacement_jacobian(ipoint, icol_point, 
    prescribed_conditions, rate_vars)

    if haskey(prescribed_conditions, ipoint)
        u_u, θ_θ = initial_point_displacement_jacobian(icol_point[ipoint], 
            prescribed_conditions[ipoint], rate_vars)
    else
        u_u, θ_θ = initial_point_displacement_jacobian(icol_point[ipoint], 
            rate_vars)
    end

    return u_u, θ_θ
end

@inline function initial_point_displacement_jacobian(icol, 
    prescribed_conditions, rate_vars)

    # unpack prescribed conditions for the node
    @unpack pd = prescribed_conditions

    # node linear displacement
    u_u = hcat(
        ifelse(pd[1], zero(e1), ifelse(rate_vars[icol+6], zero(e1), e1)),
        ifelse(pd[2], zero(e2), ifelse(rate_vars[icol+7], zero(e2), e2)),
        ifelse(pd[3], zero(e3), ifelse(rate_vars[icol+8], zero(e3), e3))
    )

    # node angular displacement
    θ_θ = hcat(
        ifelse(pd[4], zero(e1), ifelse(rate_vars[icol+9], zero(e1), e1)),
        ifelse(pd[5], zero(e2), ifelse(rate_vars[icol+10], zero(e2), e2)),
        ifelse(pd[6], zero(e3), ifelse(rate_vars[icol+11], zero(e3), e3))
    )   

    return u_u, θ_θ
end

@inline function initial_point_displacement_jacobian(icol, rate_vars)

    u_u = hcat(
        ifelse(rate_vars[icol+6], zero(e1), e1),
        ifelse(rate_vars[icol+7], zero(e2), e2),
        ifelse(rate_vars[icol+8], zero(e3), e3)
    )
    θ_θ = hcat(
        ifelse(rate_vars[icol+9], zero(e1), e1),
        ifelse(rate_vars[icol+10], zero(e2), e2),
        ifelse(rate_vars[icol+11], zero(e3), e3)
    )

    return u_u, θ_θ
end

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

"""
    initial_point_velocity_rate_jacobian(ipoint, icol_point, prescribed_conditions, 
        rate_vars)

Return the velocity rate jacobians `Vdot_Vdot` and `Ωdot_Ωdot` of point `ipoint`.  Note 
that `Vdot` and `Ωdot` in this case do not minclude any contributions resulting from 
body frame motion. 
"""
@inline function initial_point_velocity_rate_jacobian(ipoint, icol_point, 
    prescribed_conditions, rate_vars)

    if haskey(prescribed_conditions, ipoint)
        Vdot_Vdot, Ωdot_Ωdot = initial_point_velocity_rate_jacobian(icol_point[ipoint], 
            prescribed_conditions[ipoint], rate_vars)
    else
        Vdot_Vdot, Ωdot_Ωdot = initial_point_velocity_rate_jacobian(icol_point[ipoint], rate_vars)
    end

    return Vdot_Vdot, Ωdot_Ωdot
end

@inline function initial_point_velocity_rate_jacobian(icol, prescribed_conditions, 
    rate_vars)

    # unpack prescribed conditions for the node
    @unpack pd, u, theta = prescribed_conditions

    # node linear displacement
    Vdot_Vdot = hcat(
        ifelse(pd[1], zero(e1), ifelse(rate_vars[icol+6], e1, zero(e1))),
        ifelse(pd[2], zero(e2), ifelse(rate_vars[icol+7], e2, zero(e2))),
        ifelse(pd[3], zero(e3), ifelse(rate_vars[icol+8], e3, zero(e3)))
    )

    # node angular displacement
    Ωdot_Ωdot = hcat(
        ifelse(pd[4], zero(e1), ifelse(rate_vars[icol+9], e1, zero(e1))),
        ifelse(pd[5], zero(e2), ifelse(rate_vars[icol+10], e2, zero(e2))),
        ifelse(pd[6], zero(e3), ifelse(rate_vars[icol+11], e3, zero(e3)))
    )   

    return Vdot_Vdot, Ωdot_Ωdot
end

@inline function initial_point_velocity_rate_jacobian(icol, rate_vars)

    Vdot_Vdot = hcat(
        ifelse(rate_vars[icol+6], e1, zero(e1)),
        ifelse(rate_vars[icol+7], e2, zero(e2)),
        ifelse(rate_vars[icol+8], e3, zero(e3))
    )
    Ωdot_Ωdot = hcat(
        ifelse(rate_vars[icol+9], e1, zero(e1)),
        ifelse(rate_vars[icol+10], e2, zero(e2)),
        ifelse(rate_vars[icol+11], e3, zero(e3))
    )

    return Vdot_Vdot, Ωdot_Ωdot
end

"""
    static_point_properties(x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity)

Calculate/extract the point properties needed to construct the residual for a static 
analysis
"""
@inline function static_point_properties(x, indices, force_scaling, assembly, ipoint,  
    prescribed_conditions, point_masses, gravity)

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

    # forces and moments
    F, M = point_loads(x, ipoint, indices.icol_point, force_scaling, prescribed_conditions)

    # linear and angular acceleration
    a = -SVector{3}(gravity)
    α = @SVector zeros(3)

    return (; mass11, mass12, mass21, mass22, u, θ, C, F, M, a, α)
end

"""
    steady_state_point_properties(x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity)

Calculate/extract the point properties needed to construct the residual for a steady state 
analysis
"""
@inline function steady_state_point_properties(x, indices, force_scaling, 
    assembly, ipoint, prescribed_conditions, point_masses, gravity)

    properties = static_point_properties(x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity)

    @unpack mass11, mass12, mass21, mass22, u, θ, C = properties

    # rotation parameter matrices
    Qinv = get_Qinv(θ)

    # linear and angular displacement of the body frame
    ub, θb = body_frame_displacement(x)

    # linear and angular velocity of the body frame
    vb, ωb = body_frame_velocity(x)

    # linear and angular acceleration of the body frame
    ab, αb = body_frame_acceleration(x)

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

"""
    initial_condition_point_properties(x, indices, rate_vars, force_scaling, 
        assembly, ipoint, prescribed_conditions, point_masses, gravity, 
        u0, θ0, V0, Ω0, Vdot0, Ωdot0)

Calculate/extract the point properties needed to construct the residual for a time domain 
analysis initialization.
"""
@inline function initial_condition_point_properties(x, indices, rate_vars,
    force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity,
    u0, θ0, V0, Ω0, Vdot0, Ωdot0)

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

    # linear and angular displacement of the body frame
    ub, θb = body_frame_displacement(x)

    # linear and angular velocity of the body frame
    vb, ωb = body_frame_velocity(x)

    # linear and angular acceleration of the body frame
    ab, αb = body_frame_acceleration(x)

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

"""
    newmark_point_properties(x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity, udot_init, θdot_init, 
        Vdot_init, Ωdot_init, dt)

Calculate/extract the point properties needed to construct the residual for a newmark-scheme
time stepping analysis
"""
@inline function newmark_point_properties(x, indices, force_scaling, 
    assembly, ipoint, prescribed_conditions, point_masses, gravity, 
    udot_init, θdot_init, Vdot_init, Ωdot_init, dt)

    properties = steady_state_point_properties(x, indices, force_scaling, 
        assembly, ipoint, prescribed_conditions, point_masses, gravity)

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

"""
    dynamic_point_properties(dx, x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity)

Calculate/extract the point properties needed to construct the residual for a dynamic 
analysis
"""
@inline function dynamic_point_properties(dx, x, indices, force_scaling, 
    assembly, ipoint, prescribed_conditions, point_masses, gravity)

    properties = steady_state_point_properties(x, indices, force_scaling, 
        assembly, ipoint, prescribed_conditions, point_masses, gravity)

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

"""
    expanded_steady_point_properties(x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity)

Calculate/extract the point properties needed to construct the residual for a constant 
mass matrix system.
"""
@inline function expanded_steady_point_properties(x, indices, force_scaling, assembly, 
    ipoint, prescribed_conditions, point_masses, gravity)

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

    # linear and angular displacement of the body frame
    ub, θb = body_frame_displacement(x)

    # linear and angular velocity of the body frame
    vb, ωb = body_frame_velocity(x)

    # linear and angular acceleration of the body frame
    ab, αb = body_frame_acceleration(x)

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

function expanded_dynamic_point_properties(x, indices, force_scaling, assembly, 
    ipoint, prescribed_conditions, point_masses, gravity)

    properties = expanded_steady_point_properties(x, indices, force_scaling, assembly, 
        ipoint, prescribed_conditions, point_masses, gravity)

    @unpack θb, a, α = properties

    # overwrite acceleration terms so we don't double count them
    a = -get_C(θb)*SVector{3}(gravity)
    α = @SVector zeros(3)

    return (; properties..., a, α)
end

"""
    static_point_resultants(properties)

Calculate the net loads `F` and `M` applied at a point for a static analysis.
"""
@inline function static_point_resultants(properties)

    @unpack C, mass11, mass12, mass21, mass22, F, M, a, α = properties

    # add acceleration loads (including gravity)
    F -= C'*mass11*C*a + C'*mass12*C*α
    M -= C'*mass21*C*a + C'*mass22*C*α

    return F, M
end

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

"""
    static_point_residual!(resid, x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity)

Calculate and insert the residual entries corresponding to a point for a static analysis 
into the system residual vector.
"""
@inline function static_point_residual!(resid, x, indices, force_scaling, assembly, ipoint,  
    prescribed_conditions, point_masses, gravity)

    irow = indices.irow_point[ipoint]

    properties = static_point_properties(x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity)

    F, M = static_point_resultants(properties)

    resid[irow:irow+2] .= -F ./ force_scaling
    resid[irow+3:irow+5] .= -M ./ force_scaling

    return resid
end

"""
    steady_state_point_residual!(resid, x, indices, force_scaling, 
        assembly, ipoint, prescribed_conditions, point_masses, gravity)

Calculate and insert the residual entries corresponding to a point for a steady state analysis into the 
system residual vector.
"""
@inline function steady_state_point_residual!(resid, x, indices, force_scaling, 
    assembly, ipoint, prescribed_conditions, point_masses, gravity)

    irow = indices.irow_point[ipoint]

    properties = steady_state_point_properties(x, indices, force_scaling, 
        assembly, ipoint, prescribed_conditions, point_masses, gravity)

    F, M = steady_state_point_resultants(properties)

    rV, rΩ = steady_state_point_velocity_residuals(properties)
       
    resid[irow:irow+2] .= -F ./ force_scaling
    resid[irow+3:irow+5] .= -M ./ force_scaling
    resid[irow+6:irow+8] .= rV
    resid[irow+9:irow+11] .= rΩ

    return resid
end

"""
    initial_condition_point_residual!(resid, x, indices, rate_vars,  
        force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity, 
        u0, θ0, V0, Ω0, Vdot0, Ωdot0)

Calculate and insert the residual entries corresponding to a point for the initialization 
of a time domain analysis into the system residual vector.
"""
@inline function initial_condition_point_residual!(resid, x, indices, rate_vars,  
    force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity, 
    u0, θ0, V0, Ω0, Vdot0, Ωdot0)

    irow = indices.irow_point[ipoint]

    properties = initial_condition_point_properties(x, indices, rate_vars, force_scaling, 
        assembly, ipoint, prescribed_conditions, point_masses, gravity, 
        u0, θ0, V0, Ω0, Vdot0, Ωdot0)

    F, M = dynamic_point_resultants(properties)

    rV, rΩ = dynamic_point_velocity_residuals(properties)

    resid[irow:irow+2] .= -F ./ force_scaling
    resid[irow+3:irow+5] .= -M ./ force_scaling
    resid[irow+6:irow+8] .= rV
    resid[irow+9:irow+11] .= rΩ

    return resid
end

"""
    newmark_point_residual!(resid, x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity, udot_init, θdot_init, 
        Vdot_init, Ωdot_init, dt)

Calculate and insert the residual entries corresponding to a point for a newmark-scheme 
time marching analysis into the system residual vector.
"""
@inline function newmark_point_residual!(resid, x, indices, force_scaling, 
    assembly, ipoint, prescribed_conditions, point_masses, gravity, 
    udot_init, θdot_init, Vdot_init, Ωdot_init, dt)

    irow = indices.irow_point[ipoint]

    properties = newmark_point_properties(x, indices, force_scaling, 
        assembly, ipoint, prescribed_conditions, point_masses, gravity, 
        udot_init, θdot_init, Vdot_init, Ωdot_init, dt)

    F, M = dynamic_point_resultants(properties)

    rV, rΩ = dynamic_point_velocity_residuals(properties)

    resid[irow:irow+2] .= -F ./ force_scaling
    resid[irow+3:irow+5] .= -M ./ force_scaling
    resid[irow+6:irow+8] .= rV
    resid[irow+9:irow+11] .= rΩ

    return resid
end

"""
    dynamic_point_residual!(resid, dx, x, indices, force_scaling, assembly, ipoint, 
        prescribed_conditions, point_masses, gravity)

Calculate and insert the residual entries corresponding to a point for a dynamic analysis 
into the system residual vector.
"""
@inline function dynamic_point_residual!(resid, dx, x, indices, force_scaling, 
    assembly, ipoint, prescribed_conditions, point_masses, gravity)

    irow = indices.irow_point[ipoint]

    properties = dynamic_point_properties(dx, x, indices, force_scaling, assembly, 
        ipoint, prescribed_conditions, point_masses, gravity)

    F, M = dynamic_point_resultants(properties)

    rV, rΩ = dynamic_point_velocity_residuals(properties)

    resid[irow:irow+2] .= -F ./ force_scaling
    resid[irow+3:irow+5] .= -M ./ force_scaling
    resid[irow+6:irow+8] .= rV
    resid[irow+9:irow+11] .= rΩ

    return resid
end

"""
    expanded_steady_point_residual!(resid, x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity)

Calculate and insert the residual entries corresponding to a point into the system residual 
vector for a constant mass matrix system.
"""
@inline function expanded_steady_point_residual!(resid, x, indices, force_scaling, 
    assembly, ipoint, prescribed_conditions, point_masses, gravity)

    irow = indices.irow_point[ipoint]

    properties = expanded_steady_point_properties(x, indices, force_scaling, assembly, ipoint, 
        prescribed_conditions, point_masses, gravity)

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
        prescribed_conditions, point_masses, gravity)

Calculate and insert the residual entries corresponding to a point into the system residual 
vector for a constant mass matrix system.
"""
@inline function expanded_dynamic_point_residual!(resid, x, indices, force_scaling, 
    assembly, ipoint, prescribed_conditions, point_masses, gravity)

    irow = indices.irow_point[ipoint]

    properties = expanded_dynamic_point_properties(x, indices, force_scaling, assembly, ipoint, 
        prescribed_conditions, point_masses, gravity)

    F, M = expanded_point_resultants(properties)

    rV, rΩ = expanded_point_velocity_residuals(properties)

    resid[irow:irow+2] .= -F ./ force_scaling
    resid[irow+3:irow+5] .= -M ./ force_scaling
    resid[irow+6:irow+8] .= rV
    resid[irow+9:irow+11] .= rΩ

    return resid
end

"""
    static_point_jacobian_properties(properties, x, indices, force_scaling, assembly, 
        ipoint, prescribed_conditions, point_masses, gravity)

Calculate/extract the point properties needed to calculate the jacobian entries 
corresponding to a point for a static analysis
"""
@inline function static_point_jacobian_properties(properties, x, indices, 
    force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity)

    @unpack C, θ = properties

    # forces and moments
    F_θ, F_F, M_θ, M_M = point_load_jacobians(x, ipoint, indices.icol_point, force_scaling, 
        prescribed_conditions)

    # linear and angular displacement
    u_u, θ_θ = point_displacement_jacobians(ipoint, prescribed_conditions)

    # rotation parameter matrices
    C_θ1, C_θ2, C_θ3 = get_C_θ(C, θ)

    return (; properties..., C_θ1, C_θ2, C_θ3, u_u, θ_θ, F_θ, M_θ, F_F, M_M)
end

"""
    steady_state_point_jacobian_properties(properties, x, indices, force_scaling, 
        assembly, ipoint, prescribed_conditions, point_masses, gravity)

Calculate/extract the point properties needed to calculate the jacobian entries 
corresponding to a point for a steady state analysis
"""
@inline function steady_state_point_jacobian_properties(properties, x, indices, 
    force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity)

    properties = static_point_jacobian_properties(properties, x, indices, force_scaling, 
        assembly, ipoint, prescribed_conditions, point_masses, gravity)

    @unpack C, mass11, mass12, mass21, mass22, u, θ, V, Ω, θb, ωb, αb, Δx, C_θ1, C_θ2, C_θ3 = properties

    # rotation parameter matrices
    Qinv_θ1, Qinv_θ2, Qinv_θ3 = get_Qinv_θ(θ)

    # linear and angular momentum
    P_θ = mul3(C_θ1', C_θ2', C_θ3', (mass11*C*V + mass12*C*Ω)) + 
        C'*mass11*mul3(C_θ1, C_θ2, C_θ3, V) + 
        C'*mass12*mul3(C_θ1, C_θ2, C_θ3, Ω)
    P_V = C'*mass11*C
    P_Ω = C'*mass12*C

    H_θ = mul3(C_θ1', C_θ2', C_θ3', (mass21*C*V + mass22*C*Ω)) + 
        C'*mass21*mul3(C_θ1, C_θ2, C_θ3, V) + 
        C'*mass22*mul3(C_θ1, C_θ2, C_θ3, Ω)
    H_V = C'*mass21*C
    H_Ω = C'*mass22*C

    # linear and angular velocity
    v_u = tilde(ωb)
    v_vb = I3
    v_ωb = -tilde(Δx) - tilde(u)
    ω_ωb = I3

    # linear and angular acceleration
    a_u = tilde(αb)
    a_ab = I3
    a_αb = -tilde(Δx) - tilde(u)
    α_αb = I3

    # add rotated gravity vector
    C_θb1, C_θb2, C_θb3 = get_C_θ(θb)
    a_θb = -mul3(C_θb1, C_θb2, C_θb3, gravity)

    return (; properties..., Qinv_θ1, Qinv_θ2, Qinv_θ3, P_θ, P_V, P_Ω, H_θ, H_V, H_Ω, 
        v_u, v_vb, v_ωb, ω_ωb, a_u, a_θb, a_ab, a_αb, α_αb) 
end

"""
    initial_condition_point_jacobian_properties(properties, x, indices, rate_vars, 
        force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity, 
        u0, θ0, V0, Ω0, Vdot0, Ωdot0)

Calculate/extract the point properties needed to calculate the jacobian entries 
corresponding to a point for a Newmark scheme time marching analysis
"""
@inline function initial_condition_point_jacobian_properties(properties, x, indices, 
    rate_vars, force_scaling, assembly, ipoint, prescribed_conditions, point_masses, 
    gravity, u0, θ0, V0, Ω0, Vdot0, Ωdot0)

    @unpack C, Cdot, mass11, mass12, mass21, mass22, u, θ, V, Ω, Cdot, Vdot, Ωdot, 
        θb, ωb, αb, Δx, v, ω, a, α = properties

    # linear and angular displacement
    u_u, θ_θ = initial_point_displacement_jacobian(ipoint, indices.icol_point, 
        prescribed_conditions, rate_vars)

    # rotation parameter matrices
    C_θ1, C_θ2, C_θ3 = get_C_θ(C, θ)
    Qinv_θ1, Qinv_θ2, Qinv_θ3 = get_Qinv_θ(θ)

    # forces and moments
    F_θ, F_F, M_θ, M_M = point_load_jacobians(x, ipoint, indices.icol_point, force_scaling, prescribed_conditions)

    # update force and moment jacobians
    F_θ = F_θ * θ_θ
    M_θ = M_θ * θ_θ

    # linear and angular velocity
    v_u = tilde(ωb)*u_u
    v_vb = I3
    v_ωb = -tilde(Δx)-tilde(u)
    ω_ωb = I3   

    # linear and angular acceleration
    a_u = tilde(αb)*u_u
    a_ab = I3
    a_αb = -tilde(Δx)-tilde(u)
    α_αb = I3

    # add contributions from body frame motion to velocities
    V_u = v_u
    V_vb = v_vb
    V_ωb = v_ωb
    Ω_ωb = ω_ωb

    # linear and angular momentum
    P_u = C'*mass11*C*V_u*u_u
    P_θ = (mul3(C_θ1', C_θ2', C_θ3', (mass11*C*V + mass12*C*Ω)) + 
        C'*mass11*mul3(C_θ1, C_θ2, C_θ3, V) + 
        C'*mass12*mul3(C_θ1, C_θ2, C_θ3, Ω))*θ_θ
    P_vb = C'*mass11*C*V_vb
    P_ωb = C'*mass11*C*V_ωb + C'*mass12*C*Ω_ωb

    H_u = C'*mass21*C*V_u*u_u
    H_θ = (mul3(C_θ1', C_θ2', C_θ3', (mass21*C*V + mass22*C*Ω)) + 
        C'*mass21*mul3(C_θ1, C_θ2, C_θ3, V) + 
        C'*mass22*mul3(C_θ1, C_θ2, C_θ3, Ω))*θ_θ
    H_vb = C'*mass21*C*V_vb
    H_ωb = C'*mass21*C*V_ωb + C'*mass22*C*Ω_ωb

    # linear and angular acceleration
    Vdot_Vdot, Ωdot_Ωdot = initial_point_velocity_rate_jacobian(ipoint, indices.icol_point, 
        prescribed_conditions, rate_vars)

    # add contributions from body frame motion to accelerations
    Vdot_u = a_u
    Vdot_ab = a_ab
    Vdot_αb = a_αb
    Ωdot_αb = α_αb

    # linear and angular momentum rates
    Pdot_V = C'*mass11*Cdot + Cdot'*mass11*C
    Pdot_Ω = C'*mass12*Cdot + Cdot'*mass12*C
    Pdot_Vdot = C'*mass11*C
    Pdot_Ωdot = C'*mass12*C
    Pdot_u = Pdot_V*V_u + Pdot_Vdot*Vdot_u
        Pdot_θ = (mul3(C_θ1', C_θ2', C_θ3', mass11*C*Vdot) + 
        mul3(C_θ1', C_θ2', C_θ3', mass12*C*Ωdot) +
        C'*mass11*mul3(C_θ1, C_θ2, C_θ3, Vdot) + 
        C'*mass12*mul3(C_θ1, C_θ2, C_θ3, Ωdot) +
        tilde(Ω - ω)*mul3(C_θ1', C_θ2', C_θ3', mass11*C*V) + 
        tilde(Ω - ω)*mul3(C_θ1', C_θ2', C_θ3', mass12*C*Ω) +
        Cdot'*mass11*mul3(C_θ1, C_θ2, C_θ3, V) + 
        Cdot'*mass12*mul3(C_θ1, C_θ2, C_θ3, Ω) +
        mul3(C_θ1', C_θ2', C_θ3', mass11*Cdot*V) + 
        mul3(C_θ1', C_θ2', C_θ3', mass12*Cdot*Ω) +
        -C'*mass11*mul3(C_θ1, C_θ2, C_θ3, tilde(Ω - ω)*V) + 
        -C'*mass12*mul3(C_θ1, C_θ2, C_θ3, tilde(Ω - ω)*Ω))*θ_θ
    Pdot_vb = Pdot_V*V_vb
    Pdot_ωb = Pdot_V*V_ωb + Pdot_Ω*Ω_ωb
    Pdot_ab = Pdot_Vdot*Vdot_ab
    Pdot_αb = Pdot_Vdot*Vdot_αb + Pdot_Ωdot*Ωdot_αb

    Hdot_V = C'*mass21*Cdot + Cdot'*mass21*C
    Hdot_Ω = C'*mass22*Cdot + Cdot'*mass22*C
    Hdot_Vdot = C'*mass21*C
    Hdot_Ωdot = C'*mass22*C
    Hdot_u = Hdot_V*V_u + Hdot_Vdot*Vdot_u
    Hdot_θ = (mul3(C_θ1', C_θ2', C_θ3', mass21*C*Vdot) + 
        mul3(C_θ1', C_θ2', C_θ3', mass22*C*Ωdot) +
        C'*mass21*mul3(C_θ1, C_θ2, C_θ3, Vdot) + 
        C'*mass22*mul3(C_θ1, C_θ2, C_θ3, Ωdot) +
        tilde(Ω - ω)*mul3(C_θ1', C_θ2', C_θ3', mass21*C*V) + 
        tilde(Ω - ω)*mul3(C_θ1', C_θ2', C_θ3', mass22*C*Ω) +
        Cdot'*mass21*mul3(C_θ1, C_θ2, C_θ3, V) + 
        Cdot'*mass22*mul3(C_θ1, C_θ2, C_θ3, Ω) +
        mul3(C_θ1', C_θ2', C_θ3', mass21*Cdot*V) + 
        mul3(C_θ1', C_θ2', C_θ3', mass22*Cdot*Ω) +
        -C'*mass21*mul3(C_θ1, C_θ2, C_θ3, tilde(Ω - ω)*V) + 
        -C'*mass22*mul3(C_θ1, C_θ2, C_θ3, tilde(Ω - ω)*Ω))*θ_θ
    Hdot_vb = Hdot_V*V_vb
    Hdot_ωb = Hdot_V*V_ωb + Hdot_Ω*Ω_ωb
    Hdot_ab = Hdot_Vdot*Vdot_ab
    Hdot_αb = Hdot_Vdot*Vdot_αb + Hdot_Ωdot*Ωdot_αb

    # overwrite acceleration terms so we don't double count them
    a_u = @SMatrix zeros(3,3)
    a_θb = -mul3(get_C_θ(θb)..., gravity)
    a_ab = @SMatrix zeros(3,3)
    a_αb = @SMatrix zeros(3,3)
    α_αb = @SMatrix zeros(3,3)

    return (; properties..., C_θ1, C_θ2, C_θ3, Qinv_θ1, Qinv_θ2, Qinv_θ3, u_u, θ_θ, 
        F_θ, F_F, M_θ, M_M, V_u, V_vb, V_ωb, Ω_ωb, P_u, P_θ, P_vb, P_ωb, H_u, H_θ, H_vb, H_ωb, 
        Vdot_Vdot, Ωdot_Ωdot, 
        Pdot_vb, Pdot_ωb, Pdot_ab, Pdot_αb, Pdot_u, Pdot_θ, Pdot_Vdot, Pdot_Ωdot,
        Hdot_vb, Hdot_ωb, Hdot_ab, Hdot_αb, Hdot_u, Hdot_θ, Hdot_Vdot, Hdot_Ωdot, 
        v_vb, v_ωb, ω_ωb, a_u, a_θb, a_ab, a_αb, α_αb) 
end

"""
    newmark_point_jacobian_properties(properties, x, indices, force_scaling, 
        assembly, ipoint, prescribed_conditions, point_masses, gravity, 
        udot_init, θdot_init, Vdot_init, Ωdot_init, dt)

Calculate/extract the point properties needed to calculate the jacobian entries 
corresponding to a point for a Newmark scheme time marching analysis
"""
@inline function newmark_point_jacobian_properties(properties, x, indices,  
    force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity, 
    udot_init, θdot_init, Vdot_init, Ωdot_init, dt)

    properties = steady_state_point_jacobian_properties(properties, x, indices, 
        force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity)

    @unpack C, Cdot, mass11, mass12, mass21, mass22, V, Ω, Vdot, Ωdot, v, ω, C_θ1, C_θ2, C_θ3 = properties

    # linear and angular displacement rates
    udot_u = 2/dt*I3
    θdot_θ = 2/dt*I3

    # linear and angular momentum rates
    Pdot_ωb = -C'*mass11*C*tilde(V) - C'*mass12*C*tilde(Ω) + 
        tilde(C'*mass11*C*V) + tilde(C'*mass12*C*Ω)

    Pdot_θ = mul3(C_θ1', C_θ2', C_θ3', mass11*C*Vdot) + 
        mul3(C_θ1', C_θ2', C_θ3', mass12*C*Ωdot) +
        C'*mass11*mul3(C_θ1, C_θ2, C_θ3, Vdot) + 
        C'*mass12*mul3(C_θ1, C_θ2, C_θ3, Ωdot) +
        tilde(Ω - ω)*mul3(C_θ1', C_θ2', C_θ3', mass11*C*V) + 
        tilde(Ω - ω)*mul3(C_θ1', C_θ2', C_θ3', mass12*C*Ω) +
        Cdot'*mass11*mul3(C_θ1, C_θ2, C_θ3, V) + 
        Cdot'*mass12*mul3(C_θ1, C_θ2, C_θ3, Ω) +
        mul3(C_θ1', C_θ2', C_θ3', mass11*Cdot*V) + 
        mul3(C_θ1', C_θ2', C_θ3', mass12*Cdot*Ω) +
        -C'*mass11*mul3(C_θ1, C_θ2, C_θ3, tilde(Ω - ω)*V) + 
        -C'*mass12*mul3(C_θ1, C_θ2, C_θ3, tilde(Ω - ω)*Ω)

    Pdot_V = C'*mass11*Cdot + Cdot'*mass11*C

    Pdot_Ω = C'*mass12*Cdot + Cdot'*mass12*C +
        C'*mass11*C*tilde(V) + C'*mass12*C*tilde(Ω) + 
        -tilde(C'*mass11*C*V) - tilde(C'*mass12*C*Ω)

    Pdot_V += 2/dt*C'*mass11*C
    
    Pdot_Ω += 2/dt*C'*mass12*C

    Hdot_ωb = -C'*mass21*C*tilde(V) - C'*mass22*C*tilde(Ω) +
        tilde(C'*mass21*C*V) + tilde(C'*mass22*C*Ω)

    Hdot_θ = mul3(C_θ1', C_θ2', C_θ3', mass21*C*Vdot) + 
        mul3(C_θ1', C_θ2', C_θ3', mass22*C*Ωdot) +
        C'*mass21*mul3(C_θ1, C_θ2, C_θ3, Vdot) + 
        C'*mass22*mul3(C_θ1, C_θ2, C_θ3, Ωdot) +
        tilde(Ω - ω)*mul3(C_θ1', C_θ2', C_θ3', mass21*C*V) + 
        tilde(Ω - ω)*mul3(C_θ1', C_θ2', C_θ3', mass22*C*Ω) +
        Cdot'*mass21*mul3(C_θ1, C_θ2, C_θ3, V) + 
        Cdot'*mass22*mul3(C_θ1, C_θ2, C_θ3, Ω) +
        mul3(C_θ1', C_θ2', C_θ3', mass21*Cdot*V) + 
        mul3(C_θ1', C_θ2', C_θ3', mass22*Cdot*Ω) +
        -C'*mass21*mul3(C_θ1, C_θ2, C_θ3, tilde(Ω - ω)*V) + 
        -C'*mass22*mul3(C_θ1, C_θ2, C_θ3, tilde(Ω - ω)*Ω)

    Hdot_V = C'*mass21*Cdot + Cdot'*mass21*C

    Hdot_Ω = C'*mass22*Cdot + Cdot'*mass22*C +
        C'*mass21*C*tilde(V) + C'*mass22*C*tilde(Ω) +
        -tilde(C'*mass21*C*V) - tilde(C'*mass22*C*Ω)

    Hdot_V += 2/dt*C'*mass21*C

    Hdot_Ω += 2/dt*C'*mass22*C

    # overwrite acceleration terms so we don't double count them
    a_u = Z3
    a_ab = Z3
    a_αb = Z3
    α_αb = Z3

    return (; properties..., a_u, a_ab, a_αb, α_αb, udot_u, θdot_θ, 
        Pdot_ωb, Pdot_θ, Pdot_V, Pdot_Ω, 
        Hdot_ωb, Hdot_θ, Hdot_V, Hdot_Ω) 
end

"""
    dynamic_point_jacobian_properties(properties, dx, x, indices, force_scaling, 
        assembly, ipoint, prescribed_conditions, point_masses, gravity)

Calculate/extract the point properties needed to calculate the jacobian entries 
corresponding to a point for a dynamic analysis
"""
@inline function dynamic_point_jacobian_properties(properties, dx, x, indices,  
    force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity)

    properties = steady_state_point_jacobian_properties(properties, x, indices, force_scaling, 
        assembly, ipoint, prescribed_conditions, point_masses, gravity)

    @unpack C, Cdot, mass11, mass12, mass21, mass22, V, Ω, Vdot, Ωdot, v, ω, C_θ1, C_θ2, C_θ3 = properties

    # redefine linear and angular acceleration (accelerations are included in Vdot and Ωdot)
    a_u = Z3
    a_ab = Z3
    a_αb = Z3
    α_αb = Z3

    # linear and angular momentum rates
    Pdot_ωb = -C'*mass11*C*tilde(V) - C'*mass12*C*tilde(Ω) + 
        tilde(C'*mass11*C*V) + tilde(C'*mass12*C*Ω)

    Pdot_θ = mul3(C_θ1', C_θ2', C_θ3', mass11*C*Vdot) + 
        mul3(C_θ1', C_θ2', C_θ3', mass12*C*Ωdot) +
        C'*mass11*mul3(C_θ1, C_θ2, C_θ3, Vdot) + 
        C'*mass12*mul3(C_θ1, C_θ2, C_θ3, Ωdot) +
        tilde(Ω - ω)*mul3(C_θ1', C_θ2', C_θ3', mass11*C*V) + 
        tilde(Ω - ω)*mul3(C_θ1', C_θ2', C_θ3', mass12*C*Ω) +
        Cdot'*mass11*mul3(C_θ1, C_θ2, C_θ3, V) + 
        Cdot'*mass12*mul3(C_θ1, C_θ2, C_θ3, Ω) +
        mul3(C_θ1', C_θ2', C_θ3', mass11*Cdot*V) + 
        mul3(C_θ1', C_θ2', C_θ3', mass12*Cdot*Ω) +
        -C'*mass11*mul3(C_θ1, C_θ2, C_θ3, tilde(Ω - ω)*V) + 
        -C'*mass12*mul3(C_θ1, C_θ2, C_θ3, tilde(Ω - ω)*Ω)

    Pdot_V = C'*mass11*Cdot + Cdot'*mass11*C

    Pdot_Ω = C'*mass12*Cdot + Cdot'*mass12*C +
        C'*mass11*C*tilde(V) + C'*mass12*C*tilde(Ω) + 
        -tilde(C'*mass11*C*V) - tilde(C'*mass12*C*Ω)
    
    Hdot_ωb = -C'*mass21*C*tilde(V) - C'*mass22*C*tilde(Ω) +
        tilde(C'*mass21*C*V) + tilde(C'*mass22*C*Ω)

    Hdot_θ = mul3(C_θ1', C_θ2', C_θ3', mass21*C*Vdot) + 
        mul3(C_θ1', C_θ2', C_θ3', mass22*C*Ωdot) +
        C'*mass21*mul3(C_θ1, C_θ2, C_θ3, Vdot) + 
        C'*mass22*mul3(C_θ1, C_θ2, C_θ3, Ωdot) +
        tilde(Ω - ω)*mul3(C_θ1', C_θ2', C_θ3', mass21*C*V) + 
        tilde(Ω - ω)*mul3(C_θ1', C_θ2', C_θ3', mass22*C*Ω) +
        Cdot'*mass21*mul3(C_θ1, C_θ2, C_θ3, V) + 
        Cdot'*mass22*mul3(C_θ1, C_θ2, C_θ3, Ω) +
        mul3(C_θ1', C_θ2', C_θ3', mass21*Cdot*V) + 
        mul3(C_θ1', C_θ2', C_θ3', mass22*Cdot*Ω) +
        -C'*mass21*mul3(C_θ1, C_θ2, C_θ3, tilde(Ω - ω)*V) + 
        -C'*mass22*mul3(C_θ1, C_θ2, C_θ3, tilde(Ω - ω)*Ω)

    Hdot_V = C'*mass21*Cdot + Cdot'*mass21*C

    Hdot_Ω = C'*mass22*Cdot + Cdot'*mass22*C +
        C'*mass21*C*tilde(V) + C'*mass22*C*tilde(Ω) +
        -tilde(C'*mass21*C*V) - tilde(C'*mass22*C*Ω)

    # overwrite acceleration terms so we don't double count them
    a_u = Z3
    a_ab = Z3
    a_αb = Z3
    α_αb = Z3

    return (; properties..., a_u, a_ab, a_αb, α_αb, 
        Pdot_ωb, Pdot_θ, Pdot_V, Pdot_Ω, 
        Hdot_ωb, Hdot_θ, Hdot_V, Hdot_Ω) 
end

"""
    expanded_steady_point_jacobian_properties(properties, x, indices, force_scaling, 
        assembly, ipoint, prescribed_conditions, point_masses, gravity)

Calculate/extract the point properties needed to calculate the jacobian entries 
corresponding to a point for a constant mass matrix system
"""
@inline function expanded_steady_point_jacobian_properties(properties, x, indices, 
    force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity)

    @unpack C, mass11, mass12, mass21, mass22, u, θ, V, Ω, θb, ωb, αb, Δx = properties

    # forces and moments
    F_θ, F_F, M_θ, M_M = point_load_jacobians(x, ipoint, indices.icol_point, force_scaling, prescribed_conditions)

    # linear and angular displacement
    u_u, θ_θ = point_displacement_jacobians(ipoint, prescribed_conditions)

    # rotation parameter matrices
    C_θ1, C_θ2, C_θ3 = get_C_θ(C, θ)
    Qinv_θ1, Qinv_θ2, Qinv_θ3 = get_Qinv_θ(θ)

    # linear and angular momentum
    P_V = mass11
    P_Ω = mass12
    H_V = mass21
    H_Ω = mass22

    # linear and angular velocity
    v_u = tilde(ωb)
    v_vb = I3
    v_ωb = -tilde(Δx) - tilde(u)
    ω_ωb = I3

    # linear and angular acceleration
    a_θb = -mul3(get_C_θ(θb)..., gravity)
    a_ab = I3
    a_αb = -tilde(Δx) - tilde(u)
    a_u = tilde(αb)
    α_αb = I3

    return (; properties..., C_θ1, C_θ2, C_θ3, Qinv_θ1, Qinv_θ2, Qinv_θ3, u_u, θ_θ, 
        F_θ, F_F, M_θ, M_M, P_V, P_Ω, H_V, H_Ω, v_u, v_vb, v_ωb, ω_ωb, a_u, a_θb, a_ab, a_αb, α_αb) 
end

"""
    expanded_dynamic_point_jacobian_properties(properties, x, indices, 
        force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity)

Calculate/extract the point properties needed to calculate the jacobian entries 
corresponding to a point for a dynamic analysis
"""
@inline function expanded_dynamic_point_jacobian_properties(properties, x, indices, 
    force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity)

    properties = expanded_steady_point_jacobian_properties(properties, x, indices, 
        force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity)

    # overwrite acceleration terms so we don't double count them
    a_u = @SMatrix zeros(3,3)
    a_ab = @SMatrix zeros(3,3)
    a_αb = @SMatrix zeros(3,3)
    α_αb = @SMatrix zeros(3,3)

    return (; properties..., a_u, a_ab, a_αb, α_αb)
end

"""
    mass_matrix_point_jacobian_properties(x, indices, force_scaling, 
        assembly, ipoint, prescribed_conditions, point_masses)

Calculate/extract the point properties needed to calculate the mass matrix jacobian entries 
corresponding to a point
"""
@inline function mass_matrix_point_jacobian_properties(x, indices, force_scaling, 
    assembly, ipoint, prescribed_conditions, point_masses)

    # mass matrix
    mass = haskey(point_masses, ipoint) ? point_masses[ipoint].mass : @SMatrix zeros(6,6)
    
    # mass submatrices
    mass11 = mass[SVector{3}(1:3), SVector{3}(1:3)]
    mass12 = mass[SVector{3}(1:3), SVector{3}(4:6)]
    mass21 = mass[SVector{3}(4:6), SVector{3}(1:3)]
    mass22 = mass[SVector{3}(4:6), SVector{3}(4:6)]

    # linear and angular displacement
    u, θ = point_displacement(x, ipoint, indices.icol_point, prescribed_conditions)

    # linear and angular displacement rates (in the body frame)
    udot_udot, θdot_θdot = point_displacement_jacobians(ipoint, prescribed_conditions)

    # transformation matrices
    C = get_C(θ)

    # linear and angular momentum rates
    Pdot_Vdot = C'*mass11*C
    Pdot_Ωdot = C'*mass12*C
    Hdot_Vdot = C'*mass21*C
    Hdot_Ωdot = C'*mass22*C

    return (; udot_udot, θdot_θdot, Pdot_Vdot, Pdot_Ωdot, Hdot_Vdot, Hdot_Ωdot)
end

"""
    expanded_mass_matrix_point_jacobian_properties(assembly, ipoint, prescribed_conditions, point_masses)

Calculate/extract the point properties needed to calculate the mass matrix jacobian entries 
corresponding to a point for a constant mass matrix system
"""
@inline function expanded_mass_matrix_point_jacobian_properties(assembly, ipoint, prescribed_conditions, point_masses)

    # mass matrix
    mass = haskey(point_masses, ipoint) ? point_masses[ipoint].mass : @SMatrix zeros(6,6)
    
    # mass submatrices
    mass11 = mass[SVector{3}(1:3), SVector{3}(1:3)]
    mass12 = mass[SVector{3}(1:3), SVector{3}(4:6)]
    mass21 = mass[SVector{3}(4:6), SVector{3}(1:3)]
    mass22 = mass[SVector{3}(4:6), SVector{3}(4:6)]

    # linear and angular displacement rates (in the body frame)
    udot_udot, θdot_θdot = point_displacement_jacobians(ipoint, prescribed_conditions)

    # linear and angular momentum rates
    Pdot_Vdot = mass11
    Pdot_Ωdot = mass12
    Hdot_Vdot = mass21
    Hdot_Ωdot = mass22

    return (; udot_udot, θdot_θdot, Pdot_Vdot, Pdot_Ωdot, Hdot_Vdot, Hdot_Ωdot)
end

"""
    static_point_resultant_jacobians(properties)

Calculate the jacobians for the net loads `F` and `M` applied at a point for a static analysis.
"""
@inline function static_point_resultant_jacobians(properties)

    @unpack C, mass11, mass12, mass21, mass22, a, α, θ_θ, C_θ1, C_θ2, C_θ3 = properties

    @unpack F_θ, M_θ, F_F, M_M = properties

    # add loads due to linear and angular acceleration (including gravity)
    F_θ -= mul3(C_θ1', C_θ2', C_θ3', mass11*C*a) + 
        mul3(C_θ1', C_θ2', C_θ3', mass12*C*α) + 
        C'*mass11*mul3(C_θ1, C_θ2, C_θ3, a) + 
        C'*mass12*mul3(C_θ1, C_θ2, C_θ3, α)

    M_θ -= mul3(C_θ1', C_θ2', C_θ3', mass21*C*a) + 
        mul3(C_θ1', C_θ2', C_θ3', mass22*C*α) + 
        C'*mass21*mul3(C_θ1, C_θ2, C_θ3, a) + 
        C'*mass22*mul3(C_θ1, C_θ2, C_θ3, α)

    F_θ *= θ_θ
    M_θ *= θ_θ

    return (; F_θ, F_F, M_θ, M_M)
end

"""
    steady_state_point_resultant_jacobians(properties)

Calculate the jacobians for the net loads `F` and `M` applied at a point for a steady 
state analysis.
"""
@inline function steady_state_point_resultant_jacobians(properties)

    jacobians = static_point_resultant_jacobians(properties)

    @unpack C, mass11, mass12, mass21, mass22, V, Ω, P, H, ω, C_θ1, C_θ2, C_θ3, 
        u_u, θ_θ, P_θ, P_V, P_Ω, H_θ, H_V, H_Ω, 
        ω_ωb, a_u, a_θb, a_ab, a_αb, α_αb = properties

    @unpack F_θ, M_θ = jacobians

    # add loads due to linear and angular acceleration (including gravity)
    F_θb = -C'*mass11*C*a_θb
    F_ab = -C'*mass11*C*a_ab
    F_αb = -C'*mass11*C*a_αb - C'*mass12*C*α_αb
    F_u = -C'*mass11*C*a_u*u_u

    M_θb = -C'*mass21*C*a_θb
    M_ab = -C'*mass21*C*a_ab
    M_αb = -C'*mass21*C*a_αb - C'*mass22*C*α_αb
    M_u = -C'*mass21*C*a_u*u_u

    # add loads due to linear and angular momentum
    F_ωb = tilde(P)*ω_ωb
    F_θ -= tilde(ω)*P_θ*θ_θ
    F_V = -tilde(ω)*P_V
    F_Ω = -tilde(ω)*P_Ω

    M_ωb = tilde(H)*ω_ωb
    M_θ -= (tilde(ω)*H_θ + tilde(V)*P_θ)*θ_θ
    M_V = -tilde(ω)*H_V - tilde(V)*P_V + tilde(P)
    M_Ω = -tilde(ω)*H_Ω - tilde(V)*P_Ω

    return (; jacobians..., 
        F_θb, F_ωb, F_ab, F_αb, F_u, F_θ, F_V, F_Ω, 
        M_θb, M_ωb, M_ab, M_αb, M_u, M_θ, M_V, M_Ω)
end

"""
    initial_condition_point_resultant_jacobians(properties)

Calculate the jacobians for the net loads `F` and `M` applied at a point for the initialization
of a time domain analysis.
"""
@inline function initial_condition_point_resultant_jacobians(properties)

    jacobians = static_point_resultant_jacobians(properties)

    @unpack  C, mass11, mass12, mass21, mass22, V, Ω, P, H, ω, u_u, θ_θ,
        V_u, V_vb, V_ωb, P_u, P_θ, P_vb, P_ωb, H_u, H_θ, H_vb, H_ωb,     
        Pdot_vb, Pdot_ωb, Pdot_ab, Pdot_αb, Pdot_u, Pdot_θ, Pdot_Vdot, Pdot_Ωdot,
        Hdot_vb, Hdot_ωb, Hdot_ab, Hdot_αb, Hdot_u, Hdot_θ, Hdot_Vdot, Hdot_Ωdot, 
        v_vb, v_ωb, ω_ωb, a_u, a_θb, a_ab, a_αb, α_αb = properties

    @unpack F_θ, M_θ = jacobians

    # add loads due to linear and angular acceleration (including gravity)
    F_θb = -C'*mass11*C*a_θb
    F_ab = -C'*mass11*C*a_ab
    F_αb = -C'*mass11*C*a_αb - C'*mass12*C*α_αb
    F_u = -C'*mass11*C*a_u*u_u

    M_θb = -C'*mass21*C*a_θb
    M_ab = -C'*mass21*C*a_ab
    M_αb = -C'*mass21*C*a_αb - C'*mass22*C*α_αb
    M_u = -C'*mass21*C*a_u*u_u

    # add loads due to linear and angular momentum
    F_vb = -tilde(ω)*P_vb
    F_ωb = -tilde(ω)*P_ωb + tilde(P)*ω_ωb
    F_u -= tilde(ω)*P_u
    F_θ -= tilde(ω)*P_θ

    M_vb = -tilde(ω)*H_vb - tilde(V)*P_vb + tilde(P)*V_vb
    M_ωb = -tilde(ω)*H_ωb - tilde(V)*P_ωb + tilde(P)*V_ωb + tilde(H)
    M_u -= tilde(ω)*H_u + tilde(V)*P_u - tilde(P)*V_u
    M_θ -= tilde(ω)*H_θ + tilde(V)*P_θ

    # # add loads due to linear and angular momentum rates
    F_vb -= Pdot_vb
    F_ωb -= Pdot_ωb
    F_ab -= Pdot_ab
    F_αb -= Pdot_αb
    F_u -= Pdot_u
    F_θ -= Pdot_θ
    F_Vdot = -Pdot_Vdot
    F_Ωdot = -Pdot_Ωdot

    M_vb -= Hdot_vb
    M_ωb -= Hdot_ωb
    M_ab -= Hdot_ab
    M_αb -= Hdot_αb
    M_u -= Hdot_u
    M_θ -= Hdot_θ
    M_Vdot = -Hdot_Vdot
    M_Ωdot = -Hdot_Ωdot

    return (; jacobians..., 
        F_θb, F_vb, F_ωb, F_ab, F_αb, F_u, F_θ, F_Vdot, F_Ωdot, 
        M_θb, M_vb, M_ωb, M_ab, M_αb, M_u, M_θ, M_Vdot, M_Ωdot)
end

"""
    newmark_point_resultant_jacobians(properties)

Calculate the jacobians for the net loads `F` and `M` applied at a point for a Newmark 
scheme time marching analysis.
"""
@inline function newmark_point_resultant_jacobians(properties)

    jacobians = steady_state_point_resultant_jacobians(properties)

    @unpack F_ωb, F_θ, F_V, F_Ω, M_ωb, M_θ, M_V, M_Ω = jacobians

    @unpack θ_θ, Pdot_ωb, Pdot_θ, Pdot_V, Pdot_Ω, Hdot_ωb, Hdot_θ, Hdot_V, Hdot_Ω = properties

    # add loads due to linear and angular momentum rates
    F_ωb -= Pdot_ωb
    F_θ -= Pdot_θ*θ_θ
    F_V -= Pdot_V
    F_Ω -= Pdot_Ω
    
    M_ωb -= Hdot_ωb
    M_θ -= Hdot_θ*θ_θ
    M_V -= Hdot_V
    M_Ω -= Hdot_Ω

    return (; jacobians..., F_ωb, F_θ, F_V, F_Ω, M_ωb, M_θ, M_V, M_Ω)
end

"""
    dynamic_point_resultant_jacobians(properties)

Calculate the jacobians for the net loads `F` and `M` applied at a point for a dynamic 
analysis.
"""
@inline function dynamic_point_resultant_jacobians(properties)

    return newmark_point_resultant_jacobians(properties)
end

"""
    expanded_point_resultant_jacobians(properties)

Calculate the jacobians for the net loads `F` and `M` applied at a point for a constant
mass matrix system
"""
@inline function expanded_point_resultant_jacobians(properties)

    @unpack C, mass11, mass12, mass21, mass22, F, M, V, Ω, P, H, ω, a, α, 
        C_θ1, C_θ2, C_θ3, u_u, θ_θ, P_V, P_Ω, H_V, H_Ω, 
        ω_ωb, a_u, a_θb, a_ab, a_αb, α_αb = properties

    @unpack F_θ, M_θ, F_F, M_M = properties

    # rotate loads into the appropriate frame
    F_θ = (mul3(C_θ1, C_θ2, C_θ3, F) + C*F_θ)*θ_θ
    F_F = C*F_F

    M_θ = (mul3(C_θ1, C_θ2, C_θ3, M) + C*M_θ)*θ_θ
    M_M = C*M_M

    # add loads due to linear and angular acceleration (including gravity)
    F_θb = -mass11*C*a_θb
    F_ab = -mass11*C*a_ab
    F_αb = -mass11*C*a_αb - mass12*C*α_αb
    F_u = -mass11*C*a_u*u_u
    F_θ -= (mass11*mul3(C_θ1, C_θ2, C_θ3, a) + mass12*mul3(C_θ1, C_θ2, C_θ3, α))*θ_θ

    M_θb = -mass21*C*a_θb
    M_ab = -mass21*C*a_ab
    M_αb = -mass21*C*a_αb - mass22*C*α_αb
    M_u = -mass21*C*a_u*u_u
    M_θ -= (mass21*mul3(C_θ1, C_θ2, C_θ3, a) + mass22*mul3(C_θ1, C_θ2, C_θ3, α))*θ_θ

    # add loads due to linear and angular momentum
    F_ωb = @SVector zeros(3)
    F_V = -tilde(Ω)*P_V 
    F_Ω = -tilde(Ω)*P_Ω + tilde(P)

    M_ωb = @SVector zeros(3)
    M_V = -tilde(Ω)*H_V - tilde(V)*P_V + tilde(P)
    M_Ω = -tilde(Ω)*H_Ω - tilde(V)*P_Ω + tilde(H)

    (; F_θb, F_ωb, F_ab, F_αb, F_F, F_u, F_θ, F_V, F_Ω, 
       M_θb, M_ωb, M_ab, M_αb, M_M, M_u, M_θ, M_V, M_Ω)
end

"""
    mass_matrix_point_resultant_jacobians(properties)

Calculate the mass matrix jacobians for the net loads `F` and `M` applied at a point
"""
@inline function mass_matrix_point_resultant_jacobians(properties)

    @unpack Pdot_Vdot, Pdot_Ωdot, Hdot_Vdot, Hdot_Ωdot = properties

    # add loads due to linear and angular momentum rates (and gravity)
    F_Vdot = -Pdot_Vdot
    F_Ωdot = -Pdot_Ωdot
    M_Vdot = -Hdot_Vdot
    M_Ωdot = -Hdot_Ωdot

    return (; F_Vdot, F_Ωdot, M_Vdot, M_Ωdot) 
end

"""
    steady_state_point_velocity_jacobians(properties)

Calculate the jacobians of the velocity residuals `rV` and `rΩ` of a point for a steady 
state analysis.
"""
@inline function steady_state_point_velocity_jacobians(properties)

    @unpack C, Qinv, Ω, ω, v_vb, v_ωb, v_u, ω_ωb, u_u, θ_θ, C_θ1, C_θ2, C_θ3, Qinv_θ1, Qinv_θ2, Qinv_θ3 = properties
    
    rV_vb = -v_vb
    rV_ωb = -v_ωb
    rV_u = -v_u*u_u
    rV_V = I3

    rΩ_ωb = -Qinv*C*ω_ωb
    rΩ_θ = (mul3(Qinv_θ1, Qinv_θ2, Qinv_θ3, C*(Ω - ω)) + Qinv*mul3(C_θ1, C_θ2, C_θ3, Ω - ω))*θ_θ
    rΩ_Ω = Qinv*C

    return (; rV_vb, rV_ωb, rV_u, rV_V, rΩ_ωb, rΩ_θ, rΩ_Ω)
end

"""
    initial_condition_point_velocity_jacobians(properties)

Calculate the jacobians of the velocity residuals `rV` and `rΩ` of a point for the
initialization of a time domain analysis.
"""
@inline function initial_condition_point_velocity_jacobians(properties)
   
    @unpack C, Qinv, Ω, ω, v_vb, v_ωb, ω_ωb, u_u, θ_θ, C_θ1, C_θ2, C_θ3, Qinv_θ1, Qinv_θ2, Qinv_θ3 = properties

    rΩ_θ = (mul3(Qinv_θ1, Qinv_θ2, Qinv_θ3, C*(Ω - ω)) + Qinv*mul3(C_θ1, C_θ2, C_θ3, Ω - ω))*θ_θ

    rV_udot = -I3
    rΩ_θdot = -I3

    return (; rV_udot, rΩ_θ, rΩ_θdot)
end

"""
    newmark_point_velocity_jacobians(properties)

Calculate the jacobians of the velocity residuals `rV` and `rΩ` of a point for a Newmark 
scheme time-marching analysis.
"""
@inline function newmark_point_velocity_jacobians(properties)

    jacobians = steady_state_point_velocity_jacobians(properties)

    @unpack u_u, θ_θ, udot_u, θdot_θ = properties

    @unpack rV_u, rΩ_θ = jacobians

    rV_u -= udot_u*u_u

    rΩ_θ -= θdot_θ*θ_θ

    return (; jacobians..., rV_u, rΩ_θ)
end

"""
    dynamic_point_velocity_jacobians(properties)

Calculate the jacobians of the velocity residuals `rV` and `rΩ` of a point for a dynamic 
analysis.
"""
@inline function dynamic_point_velocity_jacobians(properties)

    return steady_state_point_velocity_jacobians(properties)
end

"""
    expanded_point_velocity_jacobians(properties)

Calculate the jacobians of the velocity residuals `rV` and `rΩ` of a point for a constant 
mass matrix system
"""
@inline function expanded_point_velocity_jacobians(properties)

    @unpack C, Qinv, V, Ω, ω, v_u, v_vb, v_ωb, ω_ωb, u_u, θ_θ, C_θ1, C_θ2, C_θ3, Qinv_θ1, Qinv_θ2, Qinv_θ3 = properties
    
    rV_vb = -v_vb
    rV_ωb = -v_ωb
    rV_u = -v_u*u_u
    rV_θ = mul3(C_θ1', C_θ2', C_θ3', V)*θ_θ
    rV_V = C'

    rΩ_ωb = -Qinv*C*ω_ωb
    rΩ_θ = (mul3(Qinv_θ1, Qinv_θ2, Qinv_θ3, Ω - C*ω) - Qinv*mul3(C_θ1, C_θ2, C_θ3, ω))*θ_θ
    rΩ_Ω = Qinv

    return (; rV_vb, rV_ωb, rV_u, rV_θ, rV_V, rΩ_ωb, rΩ_θ, rΩ_Ω)
end

"""
    mass_matrix_point_velocity_jacobians(properties)

Calculate the mass matrix jacobians of the velocity residuals `rV` and `rΩ` of a point
"""
@inline function mass_matrix_point_velocity_jacobians(properties)

    @unpack udot_udot, θdot_θdot = properties

    rV_udot = -udot_udot
    rΩ_θdot = -θdot_θdot

    return (; rV_udot, rΩ_θdot)
end

"""
    insert_static_point_jacobians!(jacob, indices, force_scaling, ipoint, resultants)

Insert the jacobian entries corresponding to a point for a static analysis 
into the system jacobian matrix.
"""
@inline function insert_static_point_jacobians!(jacob, indices, force_scaling, ipoint, resultants)

    @unpack F_θ, F_F, M_θ, M_M = resultants

    irow = indices.irow_point[ipoint]
    icol = indices.icol_point[ipoint]

    jacob[irow:irow+2, icol:icol+2] .= -F_F
    jacob[irow:irow+2, icol+3:icol+5] .= -F_θ ./ force_scaling
    jacob[irow+3:irow+5, icol+3:icol+5] .= -M_M - M_θ ./ force_scaling

    return jacob
end

"""
    insert_initial_condition_point_jacobians!(jacob, indices, force_scaling, ipoint, resultants, 
        velocities)

Insert the jacobian entries corresponding to a point for the initialization of a 
time domain analysis into the system jacobian matrix.
"""
@inline function insert_initial_condition_point_jacobians!(jacob, indices, force_scaling, ipoint,  
    resultants, velocities)

    @unpack F_θb, F_vb, F_ωb, F_ab, F_αb, F_u, F_θ, F_F, F_Vdot, F_Ωdot, 
            M_θb, M_vb, M_ωb, M_ab, M_αb, M_u, M_θ, M_M, M_Vdot, M_Ωdot = resultants
    
    @unpack rV_udot, rΩ_θ, rΩ_θdot = velocities

    irow = indices.irow_point[ipoint]
    icol = indices.icol_point[ipoint]

    jacob[irow:irow+2, 4:6] .= -F_θb ./ force_scaling
    jacob[irow:irow+2, 7:9] .= -F_vb ./ force_scaling
    jacob[irow:irow+2, 10:12] .= -F_ωb ./ force_scaling
    jacob[irow:irow+2, 13:15] .= -F_ab ./ force_scaling
    jacob[irow:irow+2, 16:18] .= -F_αb ./ force_scaling

    jacob[irow:irow+2, icol:icol+2] .= -F_F .- F_u ./ force_scaling .- F_Vdot ./ force_scaling
    jacob[irow:irow+2, icol+3:icol+5] .= -F_θ ./ force_scaling .- F_Ωdot ./ force_scaling

    jacob[irow+3:irow+5, 4:6] .= -M_θb ./ force_scaling
    jacob[irow+3:irow+5, 7:9] .= -M_vb ./ force_scaling
    jacob[irow+3:irow+5, 10:12] .= -M_ωb ./ force_scaling
    jacob[irow+3:irow+5, 13:15] .= -M_ab ./ force_scaling
    jacob[irow+3:irow+5, 16:18] .= -M_αb ./ force_scaling

    jacob[irow+3:irow+5, icol:icol+2] .= -M_u ./ force_scaling .- M_Vdot ./ force_scaling 
    jacob[irow+3:irow+5, icol+3:icol+5] .= -M_M .- M_θ ./ force_scaling .- M_Ωdot ./ force_scaling

    jacob[irow+6:irow+8, icol+6:icol+8] .= rV_udot
    jacob[irow+9:irow+11, icol+3:icol+5] .= rΩ_θ
    jacob[irow+9:irow+11, icol+9:icol+11] .= rΩ_θdot

    return jacob
end

"""
    insert_dynamic_point_jacobians!(jacob, indices, force_scaling, ipoint, resultants, 
        velocities)

Insert the jacobian entries corresponding to a point for a steady state analysis into the 
system jacobian matrix.
"""
@inline function insert_dynamic_point_jacobians!(jacob, indices, force_scaling, ipoint,  
    resultants, velocities)

    insert_static_point_jacobians!(jacob, indices, force_scaling, ipoint, resultants)

    @unpack F_θb, F_ωb, F_ab, F_αb, F_u, F_θ, F_V, F_Ω, 
            M_θb, M_ωb, M_ab, M_αb, M_u, M_θ, M_V, M_Ω = resultants

    @unpack rV_vb, rV_ωb, rV_u, rV_V, rΩ_ωb, rΩ_θ, rΩ_Ω = velocities

    irow = indices.irow_point[ipoint]
    icol = indices.icol_point[ipoint]

    jacob[irow:irow+2, 4:6] .= -F_θb ./ force_scaling
    jacob[irow:irow+2, 10:12] .= -F_ωb ./ force_scaling
    jacob[irow:irow+2, 13:15] .= -F_ab ./ force_scaling
    jacob[irow:irow+2, 16:18] .= -F_αb ./ force_scaling

    @views jacob[irow:irow+2, icol:icol+2] .-= F_u ./ force_scaling
    jacob[irow:irow+2, icol+6:icol+8] .= -F_V ./ force_scaling
    jacob[irow:irow+2, icol+9:icol+11] .= -F_Ω ./ force_scaling

    jacob[irow+3:irow+5, 4:6] .= -M_θb ./ force_scaling
    jacob[irow+3:irow+5, 10:12] .= -M_ωb ./ force_scaling
    jacob[irow+3:irow+5, 13:15] .= -M_ab ./ force_scaling
    jacob[irow+3:irow+5, 16:18] .= -M_αb ./ force_scaling

    jacob[irow+3:irow+5, icol:icol+2] .-= M_u ./ force_scaling
    jacob[irow+3:irow+5, icol+6:icol+8] .= -M_V ./ force_scaling
    jacob[irow+3:irow+5, icol+9:icol+11] .= -M_Ω ./ force_scaling

    jacob[irow+6:irow+8, 7:9] = rV_vb
    jacob[irow+6:irow+8, 10:12] = rV_ωb

    jacob[irow+6:irow+8, icol:icol+2] .= rV_u
    jacob[irow+6:irow+8, icol+6:icol+8] .= rV_V

    jacob[irow+9:irow+11, 10:12] = rΩ_ωb

    jacob[irow+9:irow+11, icol+3:icol+5] .= rΩ_θ
    jacob[irow+9:irow+11, icol+9:icol+11] .= rΩ_Ω

    return jacob
end

"""
    insert_expanded_point_jacobians!(jacob, indices, force_scaling, ipoint, resultants, 
        velocities)

Insert the jacobian entries corresponding to a point for a constant mass matrix system into 
the system jacobian matrix for a constant mass matrix system.
"""
@inline function insert_expanded_point_jacobians!(jacob, indices, force_scaling, ipoint,  
    resultants, velocities)

    insert_dynamic_point_jacobians!(jacob, indices, force_scaling, ipoint, resultants, 
        velocities)

    @unpack rV_θ = velocities

    irow = indices.irow_point[ipoint]
    icol = indices.icol_point[ipoint]

    jacob[irow+6:irow+8, icol+3:icol+5] .= rV_θ

    return jacob
end

"""
    insert_mass_matrix_point_jacobians!(jacob, gamma, indices, force_scaling, ipoint, 
        resultants, velocities)

Insert the mass matrix jacobian entries corresponding to a point into the system jacobian 
matrix.
"""
@inline function insert_mass_matrix_point_jacobians!(jacob, gamma, indices, force_scaling, ipoint,  
    resultants, velocities)

    @unpack F_Vdot, F_Ωdot, M_Vdot, M_Ωdot = resultants
    @unpack rV_udot, rΩ_θdot = velocities

    irow = indices.irow_point[ipoint]
    icol = indices.icol_point[ipoint]

    @views jacob[irow:irow+2, icol+6:icol+8] .-= F_Vdot .* gamma ./ force_scaling
    @views jacob[irow:irow+2, icol+9:icol+11] .-= F_Ωdot .* gamma ./ force_scaling

    @views jacob[irow+3:irow+5, icol+6:icol+8] .-= M_Vdot .* gamma ./ force_scaling
    @views jacob[irow+3:irow+5, icol+9:icol+11] .-= M_Ωdot .* gamma ./ force_scaling

    @views jacob[irow+6:irow+8, icol:icol+2] .+= rV_udot .* gamma
    @views jacob[irow+9:irow+11, icol+3:icol+5] .+= rΩ_θdot .* gamma

    return jacob
end

"""
    static_point_jacobian!(jacob, x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity)

Calculate and insert the jacobian entries corresponding to a point for a static analysis 
into the system jacobian matrix.
"""
@inline function static_point_jacobian!(jacob, x, indices, force_scaling, assembly, ipoint,  
    prescribed_conditions, point_masses, gravity)

    properties = static_point_properties(x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity)

    properties = static_point_jacobian_properties(properties, x, indices, 
        force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity)

    resultants = static_point_resultant_jacobians(properties)

    insert_static_point_jacobians!(jacob, indices, force_scaling, ipoint, resultants)

    return jacob
end

"""
    steady_state_point_jacobian!(jacob, x, indices, force_scaling, assembly, 
        ipoint, prescribed_conditions, point_masses, gravity)

Calculate and insert the jacobian entries corresponding to a point for a steady state 
analysis into the system jacobian matrix.
"""
@inline function steady_state_point_jacobian!(jacob, x, indices, force_scaling, 
    assembly, ipoint, prescribed_conditions, point_masses, gravity)

    properties = steady_state_point_properties(x, indices, force_scaling, 
        assembly, ipoint, prescribed_conditions, point_masses, gravity)

    properties = steady_state_point_jacobian_properties(properties, x, indices,
        force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity)

    resultants = steady_state_point_resultant_jacobians(properties)

    velocities = steady_state_point_velocity_jacobians(properties)

    insert_dynamic_point_jacobians!(jacob, indices, force_scaling, ipoint, resultants, velocities)

    return jacob
end

"""
    initial_condition_point_jacobian!(jacob, x, indices, rate_vars, 
        force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity, 
        u0, θ0, V0, Ω0, Vdot0, Ωdot0)

Calculate and insert the jacobian entries corresponding to a point for the initialization
of a time domain analysis into the system jacobian matrix.
"""
@inline function initial_condition_point_jacobian!(jacob, x, indices, rate_vars, 
    force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity, 
    u0, θ0, V0, Ω0, Vdot0, Ωdot0)

    properties = initial_condition_point_properties(x, indices, rate_vars, 
        force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity, 
        u0, θ0, V0, Ω0, Vdot0, Ωdot0)
    
    properties = initial_condition_point_jacobian_properties(properties, x, indices, 
        rate_vars, force_scaling, assembly, ipoint, prescribed_conditions, 
        point_masses, gravity, u0, θ0, V0, Ω0, Vdot0, Ωdot0)

    resultants = initial_condition_point_resultant_jacobians(properties)

    velocities = initial_condition_point_velocity_jacobians(properties)

    insert_initial_condition_point_jacobians!(jacob, indices, force_scaling, ipoint, resultants, velocities)

    return jacob
end

"""
    newmark_point_jacobian!(jacob, x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity)

Calculate and insert the jacobian entries corresponding to a point for a Newmark scheme
time marching analysis into the system jacobian matrix.
"""
@inline function newmark_point_jacobian!(jacob, x, indices, force_scaling, 
    assembly, ipoint, prescribed_conditions, point_masses, gravity, 
    udot_init, θdot_init, Vdot_init, Ωdot_init, dt)

    properties = newmark_point_properties(x, indices, force_scaling, 
        assembly, ipoint, prescribed_conditions, point_masses, gravity, 
        udot_init, θdot_init, Vdot_init, Ωdot_init, dt)

    properties = newmark_point_jacobian_properties(properties, x, indices, 
        force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity, 
        udot_init, θdot_init, Vdot_init, Ωdot_init, dt)

    resultants = newmark_point_resultant_jacobians(properties)

    velocities = newmark_point_velocity_jacobians(properties)

    insert_dynamic_point_jacobians!(jacob, indices, force_scaling, ipoint, resultants, velocities)

    return jacob
end

"""
    dynamic_point_jacobian!(jacob, dx, x, indices, force_scaling, 
        assembly, ipoint, prescribed_conditions, point_masses, gravity)

Calculate and insert the jacobian entries corresponding to a point for a dynamic 
analysis into the system jacobian matrix.
"""
@inline function dynamic_point_jacobian!(jacob, dx, x, indices, force_scaling, 
    assembly, ipoint, prescribed_conditions, point_masses, gravity)

    properties = dynamic_point_properties(dx, x, indices, force_scaling, 
        assembly, ipoint, prescribed_conditions, point_masses, gravity)

    properties = dynamic_point_jacobian_properties(properties, dx, x, indices, 
        force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity)

    resultants = dynamic_point_resultant_jacobians(properties)

    velocities = dynamic_point_velocity_jacobians(properties)

    insert_dynamic_point_jacobians!(jacob, indices, force_scaling, ipoint, resultants, velocities)

    return jacob
end

"""
    expanded_steady_point_jacobian!(jacob, x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity)

Calculate and insert the jacobian entries corresponding to a point for a constant mass 
matrix system into the system jacobian matrix.
"""
@inline function expanded_steady_point_jacobian!(jacob, x, indices, force_scaling, 
    assembly, ipoint, prescribed_conditions, point_masses, gravity)

    properties = expanded_steady_point_properties(x, indices, force_scaling, 
        assembly, ipoint, prescribed_conditions, point_masses, gravity)

    properties = expanded_steady_point_jacobian_properties(properties, x, indices,
        force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity)

    resultants = expanded_point_resultant_jacobians(properties)

    velocities = expanded_point_velocity_jacobians(properties)

    insert_expanded_point_jacobians!(jacob, indices, force_scaling, ipoint, resultants, velocities)

    return jacob
end

"""
    expanded_dynamic_point_jacobian!(jacob, x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity)

Calculate and insert the jacobian entries corresponding to a point for a constant mass 
matrix system into the system jacobian matrix.
"""
@inline function expanded_dynamic_point_jacobian!(jacob, x, indices, force_scaling, 
    assembly, ipoint, prescribed_conditions, point_masses, gravity)

    properties = expanded_dynamic_point_properties(x, indices, force_scaling, 
        assembly, ipoint, prescribed_conditions, point_masses, gravity)

    properties = expanded_dynamic_point_jacobian_properties(properties, x, indices,
        force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity)

    resultants = expanded_point_resultant_jacobians(properties)

    velocities = expanded_point_velocity_jacobians(properties)

    insert_expanded_point_jacobians!(jacob, indices, force_scaling, ipoint, resultants, velocities)

    return jacob
end

"""
    mass_matrix_point_jacobian!(jacob, gamma, x, indices, force_scaling, assembly, 
        ipoint, prescribed_conditions)

Calculate and insert the mass_matrix jacobian entries corresponding to a point into the 
system jacobian matrix.
"""
@inline function mass_matrix_point_jacobian!(jacob, gamma, x, indices, force_scaling, assembly, 
    ipoint, prescribed_conditions, point_masses)

    properties = mass_matrix_point_jacobian_properties(x, indices, force_scaling, assembly, ipoint, prescribed_conditions, point_masses)

    resultants = mass_matrix_point_resultant_jacobians(properties)

    velocities = mass_matrix_point_velocity_jacobians(properties)
    
    insert_mass_matrix_point_jacobians!(jacob, gamma, indices, force_scaling, ipoint,  
        resultants, velocities)

    return jacob
end

"""
    expanded_mass_matrix_point_jacobian!(jacob, gamma, indices, force_scaling, assembly, 
        ipoint, prescribed_conditions, point_masses)

Calculate and insert the mass_matrix jacobian entries corresponding to a point into the 
system jacobian matrix for a constant mass matrix system
"""
@inline function expanded_mass_matrix_point_jacobian!(jacob, gamma, indices, force_scaling, assembly, 
    ipoint, prescribed_conditions, point_masses)

    properties = expanded_mass_matrix_point_jacobian_properties(assembly, ipoint, prescribed_conditions, point_masses)

    resultants = mass_matrix_point_resultant_jacobians(properties)

    velocities = mass_matrix_point_velocity_jacobians(properties)
    
    insert_mass_matrix_point_jacobians!(jacob, gamma, indices, force_scaling, 
        ipoint, resultants, velocities)

    return jacob
end