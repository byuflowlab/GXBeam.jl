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

    # unpack prescribed conditions for the node
    @unpack value, follower, isforce = prescribed_conditions

    # displacements
    u, θ = point_displacement(x, icol, prescribed_conditions)

    # transformation matrix
    C = get_C(θ)

    # calculate prescribed external forces
    F = zero(u)
    for i = 1:3
        if isforce[i]
            # add dead force
            F += SVector(I3[i,1], I3[i,2], I3[i,3])*value[i]
            # add follower force
            F += SVector(C[i,1], C[i,2], C[i,3])*follower[i]
        end
    end

    # overwrite prescribed external forces if linear displacements are prescribed instead
    F = SVector(ifelse(isforce[1], F[1], x[icol  ]*force_scaling),
                ifelse(isforce[2], F[2], x[icol+1]*force_scaling),
                ifelse(isforce[3], F[3], x[icol+2]*force_scaling))

    # calculate prescribed external moments
    M = zero(θ)
    for i = 4:6
        if isforce[i]
            # add dead moment
            M += SVector(I3[i-3,1], I3[i-3,2], I3[i-3,3])*value[i]
            # add follower moment
            M += SVector(C[i-3,1], C[i-3,2], C[i-3,3])*follower[i]
        end
    end

    # overwrite prescribed external moments if angular displacements are prescribed instead
    M = SVector(ifelse(isforce[4], M[1], x[icol+3]*force_scaling),
                ifelse(isforce[5], M[2], x[icol+4]*force_scaling),
                ifelse(isforce[6], M[3], x[icol+5]*force_scaling))

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

    @unpack value, follower, isforce = prescribed_conditions

    u, θ = point_displacement(x, icol, prescribed_conditions)

    F_F = @SMatrix zeros(eltype(x), 3, 3)
    if !isforce[1] 
        F_F += @SMatrix [1 0 0; 0 0 0; 0 0 0]
    end
    if !isforce[2] 
        F_F += @SMatrix [0 0 0; 0 1 0; 0 0 0]
    end
    if !isforce[3] 
        F_F += @SMatrix [0 0 0; 0 0 0; 0 0 1]
    end

    M_M = @SMatrix zeros(eltype(x), 3, 3)
    if !isforce[4] 
        M_M += @SMatrix [1 0 0; 0 0 0; 0 0 0]
    end
    if !isforce[5] 
        M_M += @SMatrix [0 0 0; 0 1 0; 0 0 0]
    end
    if !isforce[6] 
        M_M += @SMatrix [0 0 0; 0 0 0; 0 0 1]
    end

    C = get_C(θ)
    C_θ1, C_θ2, C_θ3 = get_C_θ(C, θ)           

    # solve for the jacobian wrt theta of the follower forces
    Fp_θ = @SMatrix zeros(eltype(x), 3, 3)
    for i = 1:3
        if isforce[i]
            rot_θ = @SMatrix [
                C_θ1[i,1] C_θ2[i,1] C_θ3[i,1];
                C_θ1[i,2] C_θ2[i,2] C_θ3[i,2];
                C_θ1[i,3] C_θ2[i,3] C_θ3[i,3]
                ]
            Fp_θ += rot_θ*follower[i]
        end
    end

    # solve for the jacobian wrt theta of the follower moments
    Mp_θ = @SMatrix zeros(eltype(x), 3, 3)
    for i = 1:3
        if isforce[i+3]
            rot_θ = @SMatrix [
                C_θ1[i,1] C_θ2[i,1] C_θ3[i,1];
                C_θ1[i,2] C_θ2[i,2] C_θ3[i,2];
                C_θ1[i,3] C_θ2[i,3] C_θ3[i,3]
                ]
            Mp_θ += rot_θ*follower[i+3]
        end
    end

    # if displacement is specified, corresponding component of follower jacobian is zero
    F1_θ = ifelse(isforce[1], SVector(Fp_θ[1,1], Fp_θ[1,2], Fp_θ[1,3]), (@SVector zeros(eltype(Fp_θ), 3)))
    F2_θ = ifelse(isforce[2], SVector(Fp_θ[2,1], Fp_θ[2,2], Fp_θ[2,3]), (@SVector zeros(eltype(Fp_θ), 3)))
    F3_θ = ifelse(isforce[3], SVector(Fp_θ[3,1], Fp_θ[3,2], Fp_θ[3,3]), (@SVector zeros(eltype(Fp_θ), 3)))
    F_θ = vcat(F1_θ', F2_θ', F3_θ')

    M1_θ = ifelse(isforce[4], SVector(Mp_θ[1,1], Mp_θ[1,2], Mp_θ[1,3]), (@SVector zeros(eltype(Fp_θ), 3)))
    M2_θ = ifelse(isforce[5], SVector(Mp_θ[2,1], Mp_θ[2,2], Mp_θ[2,3]), (@SVector zeros(eltype(Fp_θ), 3)))
    M3_θ = ifelse(isforce[6], SVector(Mp_θ[3,1], Mp_θ[3,2], Mp_θ[3,3]), (@SVector zeros(eltype(Fp_θ), 3)))
    M_θ = vcat(M1_θ', M2_θ', M3_θ')

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
    @unpack value, follower, isforce = prescribed_conditions

    # node linear displacement
    u = SVector(ifelse(isforce[1], x[icol  ], value[1]),
                ifelse(isforce[2], x[icol+1], value[2]),
                ifelse(isforce[3], x[icol+2], value[3]))

    # node angular displacement
    θ = SVector(ifelse(isforce[4], x[icol+3], value[4]),
               ifelse(isforce[5], x[icol+4], value[5]),
               ifelse(isforce[6], x[icol+5], value[6]))

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

    @unpack isforce = prescribed_conditions

    u_u = hcat(ifelse(isforce[1], e1, zero(e1)),
               ifelse(isforce[2], e2, zero(e2)),
               ifelse(isforce[3], e3, zero(e3)))

    θ_θ = hcat(ifelse(isforce[4], e1, zero(e1)),
               ifelse(isforce[5], e2, zero(e2)),
               ifelse(isforce[6], e3, zero(e3)))

    return u_u, θ_θ
end

@inline function point_displacement_jacobians()

    u_u = I3
    θ_θ = I3

    return u_u, θ_θ
end

"""
    point_displacement_rates(dx, ipoint, icol, prescribed_conditions)

Extract the displacement rates `udot` and `θdot` of point `ipoint` from the rate variable vector.
"""
@inline function point_displacement_rates(dx, ipoint, icol, prescribed_conditions)

    if haskey(prescribed_conditions, ipoint)
        udot, θdot = point_displacement_rates(dx, icol[ipoint], prescribed_conditions[ipoint])
    else
        udot, θdot = point_displacement_rates(dx, icol[ipoint])
    end

    return udot, θdot
end

@inline function point_displacement_rates(dx, icol, prescribed_conditions)

    # unpack prescribed conditions for the node
    @unpack value, follower, isforce = prescribed_conditions

    # node linear displacement
    udot = SVector(ifelse(isforce[1], dx[icol  ], zero(eltype(dx))),
                   ifelse(isforce[2], dx[icol+1], zero(eltype(dx))),
                   ifelse(isforce[3], dx[icol+2], zero(eltype(dx))))

    # node angular displacement
    θdot = SVector(ifelse(isforce[4], dx[icol+3], zero(eltype(dx))),
                   ifelse(isforce[5], dx[icol+4], zero(eltype(dx))),
                   ifelse(isforce[6], dx[icol+5], zero(eltype(dx))))

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
    point_velocity_rates(dx, ipoint, icol, prescribed_conditions)

Extract the displacement rates `udot` and `θdot` of point `ipoint` from the rate variable vector.
"""
@inline function point_velocity_rates(dx, ipoint, icol, prescribed_conditions)

    if haskey(prescribed_conditions, ipoint)
        Vdot, Ωdot = point_velocity_rates(dx, icol[ipoint], prescribed_conditions[ipoint])
    else
        Vdot, Ωdot = point_velocity_rates(dx, icol[ipoint])
    end

    return Vdot, Ωdot
end

@inline function point_velocity_rates(dx, icol, prescribed_conditions)

    # unpack prescribed conditions for the node
    @unpack value, follower, isforce = prescribed_conditions

    # node linear displacement
    Vdot = SVector(ifelse(isforce[1], dx[icol+6], zero(eltype(dx))),
                   ifelse(isforce[2], dx[icol+7], zero(eltype(dx))),
                   ifelse(isforce[3], dx[icol+8], zero(eltype(dx))))

    # node angular displacement
    Ωdot = SVector(ifelse(isforce[4], dx[icol+9], zero(eltype(dx))),
                   ifelse(isforce[5], dx[icol+10], zero(eltype(dx))),
                   ifelse(isforce[6], dx[icol+11], zero(eltype(dx))))

    return Vdot, Ωdot
end

@inline function point_velocity_rates(dx, icol)

    Vdot = SVector(dx[icol+6], dx[icol+7], dx[icol+8])
    Ωdot = SVector(dx[icol+9], dx[icol+10], dx[icol+11])

    return Vdot, Ωdot
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

    return (; C, mass11, mass12, mass21, mass22, F, M, u, θ, a, α)
end

"""
    steady_state_point_properties(x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)

Calculate/extract the point properties needed to construct the residual for a steady state 
analysis
"""
@inline function steady_state_point_properties(x, indices, force_scaling, assembly, ipoint,  
    prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)

    properties = static_point_properties(x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity)

    @unpack C, mass11, mass12, mass21, mass22, u, θ, a, α = properties

    Qinv = get_Qinv(θ)

    # linear and angular velocity
    V, Ω = point_velocities(x, ipoint, indices.icol_point)

    # linear and angular momentum
    P = C'*mass11*C*V + C'*mass12*C*Ω
    H = C'*mass21*C*V + C'*mass22*C*Ω

    # distance from the rotation center (in the body frame)
    Δx = assembly.points[ipoint] - x0

    # undeformed linear and angular velocity (in the body frame)
    v = v0 + cross(ω0, Δx)
    ω = ω0

    # linear and angular acceleration (in the body frame)
    a += a0 + cross(α0, Δx) + cross(α0, u)
    α += α0

    return (; properties..., Qinv, V, Ω, P, H, v, ω, a, α) 
end

"""
    initial_condition_point_properties(x, indices, rate_vars, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0,
        u0, θ0, V0, Ω0, Vdot0, Ωdot0)

Calculate/extract the point properties needed to construct the residual for a time domain 
analysis initialization.
"""
@inline function initial_condition_point_properties(x, indices, rate_vars, force_scaling, assembly, ipoint,  
    prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0,
    u0, θ0, V0, Ω0, Vdot0, Ωdot0)

    # mass matrix
    mass = haskey(point_masses, ipoint) ? point_masses[ipoint].mass : @SMatrix zeros(6,6)
    
    # mass submatrices
    mass11 = mass[SVector{3}(1:3), SVector{3}(1:3)]
    mass12 = mass[SVector{3}(1:3), SVector{3}(4:6)]
    mass21 = mass[SVector{3}(4:6), SVector{3}(1:3)]
    mass22 = mass[SVector{3}(4:6), SVector{3}(4:6)]

    # linear and angular displacement
    u, θ = point_displacement(x, ipoint, indices.icol_point, prescribed_conditions)
    icol = indices.icol_point[ipoint]
    if haskey(prescribed_conditions, ipoint)
        # linear displacement
        u = SVector{3}(
            prescribed_conditions[ipoint].isforce[1] && rate_vars[icol+6] ? u0[ipoint][1] : u[1], 
            prescribed_conditions[ipoint].isforce[2] && rate_vars[icol+7] ? u0[ipoint][2] : u[2],
            prescribed_conditions[ipoint].isforce[3] && rate_vars[icol+8] ? u0[ipoint][3] : u[3],
        )
        # angular displacement
        θ = SVector{3}(
            prescribed_conditions[ipoint].isforce[4] && rate_vars[icol+9] ? θ0[ipoint][1] : θ[1],
            prescribed_conditions[ipoint].isforce[5] && rate_vars[icol+10] ? θ0[ipoint][2] : θ[2],
            prescribed_conditions[ipoint].isforce[6] && rate_vars[icol+11] ? θ0[ipoint][3] : θ[3], 
        )
    else
        # linear displacement
        u = SVector{3}(
            rate_vars[icol+6] ? u0[ipoint][1] : u[1], 
            rate_vars[icol+7] ? u0[ipoint][2] : u[2],
            rate_vars[icol+8] ? u0[ipoint][3] : u[3],
        )
        # angular displacement
        θ = SVector{3}(
            rate_vars[icol+9] ? θ0[ipoint][1] : θ[1],
            rate_vars[icol+10] ? θ0[ipoint][2] : θ[2],
            rate_vars[icol+11] ? θ0[ipoint][3] : θ[3], 
        )
    end

    # rotation parameter matrices
    C = get_C(θ)
    Qinv = get_Qinv(θ)

    # linear and angular velocity
    V = SVector{3}(V0[ipoint])
    Ω = SVector{3}(Ω0[ipoint])

    # linear and angular momentum
    P = C'*mass11*C*V + C'*mass12*C*Ω
    H = C'*mass21*C*V + C'*mass22*C*Ω

    # forces and moments
    F, M = point_loads(x, ipoint, indices.icol_point, force_scaling, prescribed_conditions)

    # distance from the rotation center
    Δx = assembly.points[ipoint] - x0

    # rigid body linear and angular velocity
    v = v0 + cross(ω0, Δx)
    ω = ω0

    # linear and angular acceleration
    a = a0 + cross(α0, Δx) + cross(α0, u) - gravity
    α = α0

    # linear and angular displacement rates
    udot, θdot = point_velocities(x, ipoint, indices.icol_point)

    # linear and angular velocity rates
    Vdot, Ωdot = point_displacement_rates(x, ipoint, indices.icol_point, prescribed_conditions)
    if haskey(prescribed_conditions, ipoint)
        Vdot = SVector{3}(
            prescribed_conditions[ipoint].isforce[1] && rate_vars[icol+6] ? Vdot[1] : Vdot0[ipoint][1], 
            prescribed_conditions[ipoint].isforce[2] && rate_vars[icol+7] ? Vdot[2] : Vdot0[ipoint][2],
            prescribed_conditions[ipoint].isforce[3] && rate_vars[icol+8] ? Vdot[3] : Vdot0[ipoint][3],
        )
        Ωdot = SVector{3}(
            prescribed_conditions[ipoint].isforce[4] && rate_vars[icol+9] ? Ωdot[1] : Ωdot0[ipoint][1],
            prescribed_conditions[ipoint].isforce[5] && rate_vars[icol+10] ? Ωdot[2] : Ωdot0[ipoint][2],
            prescribed_conditions[ipoint].isforce[6] && rate_vars[icol+11] ? Ωdot[3] : Ωdot0[ipoint][3], 
        )
    else
        Vdot = SVector{3}(
            rate_vars[icol+6] ? Vdot[1] : Vdot0[ipoint][1], 
            rate_vars[icol+7] ? Vdot[2] : Vdot0[ipoint][2],
            rate_vars[icol+8] ? Vdot[3] : Vdot0[ipoint][3],
        )
        Ωdot = SVector{3}(
            rate_vars[icol+9] ? Ωdot[1] : Ωdot0[ipoint][1],
            rate_vars[icol+10] ? Ωdot[2] : Ωdot0[ipoint][2],
            rate_vars[icol+11] ? Ωdot[3] : Ωdot0[ipoint][3], 
        )
    end

    # linear and angular momentum rates
    Cdot = -C*tilde(Ω - ω)

    Pdot = C'*mass11*C*Vdot + C'*mass12*C*Ωdot +
        C'*mass11*Cdot*V + C'*mass12*Cdot*Ω +
        Cdot'*mass11*C*V + Cdot'*mass12*C*Ω
    
    Hdot = C'*mass21*C*Vdot + C'*mass22*C*Ωdot +
        C'*mass21*Cdot*V + C'*mass22*Cdot*Ω +
        Cdot'*mass21*C*V + Cdot'*mass22*C*Ω

    return (; C, Qinv, mass11, mass12, mass21, mass22, u, θ, V, Ω, P, H, F, M, v, ω, a, α,
        udot, θdot, Cdot, Vdot, Ωdot, Pdot, Hdot) 
end

"""
    newmark_point_properties(x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0,
        udot_init, θdot_init, Vdot_init, Ωdot_init, dt)

Calculate/extract the point properties needed to construct the residual for a newmark-scheme
time stepping analysis
"""
@inline function newmark_point_properties(x, indices, force_scaling, assembly, ipoint,  
    prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0, 
    udot_init, θdot_init, Vdot_init, Ωdot_init, dt)

    properties = steady_state_point_properties(x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)

    @unpack C, mass11, mass12, mass21, mass22, u, θ, V, Ω, ω = properties

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

    Pdot = C'*mass11*C*Vdot + C'*mass12*C*Ωdot +
        - C'*mass11*C*tilde(Ω - ω)*V - C'*mass12*C*tilde(Ω - ω)*Ω +
        tilde(Ω - ω)*C'*mass11*C*V + tilde(Ω - ω)*C'*mass12*C*Ω
    
    Hdot = C'*mass21*C*Vdot + C'*mass22*C*Ωdot +
        - C'*mass21*C*tilde(Ω - ω)*V - C'*mass22*C*tilde(Ω - ω)*Ω +
        tilde(Ω - ω)*C'*mass21*C*V + tilde(Ω - ω)*C'*mass22*C*Ω

    return (; properties..., udot, θdot, Cdot, Vdot, Ωdot, Pdot, Hdot) 
end

"""
    dynamic_point_properties(dx, x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)

Calculate/extract the point properties needed to construct the residual for a dynamic 
analysis
"""
@inline function dynamic_point_properties(dx, x, indices, force_scaling, assembly, ipoint,  
    prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)

    properties = steady_state_point_properties(x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)

    @unpack C, mass11, mass12, mass21, mass22, V, Ω, ω = properties 

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

    return (; properties..., udot, θdot, Cdot, Vdot, Ωdot, Pdot, Hdot) 
end

"""
    expanded_point_properties(x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)

Calculate/extract the point properties needed to construct the residual for a constant 
mass matrix system.
"""
@inline function expanded_point_properties(x, indices, force_scaling, assembly, ipoint,  
    prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)

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

    # linear and angular velocity
    V, Ω = point_velocities(x, ipoint, indices.icol_point)

    # linear and angular momentum
    P = mass11*V + mass12*Ω
    H = mass21*V + mass22*Ω

    # forces and moments
    F, M = point_loads(x, ipoint, indices.icol_point, force_scaling, prescribed_conditions)

    # distance from the rotation center (in the body frame)
    Δx = assembly.points[ipoint] - x0

    # undeformed linear and angular velocity (in the body frame)
    v = v0 + cross(ω0, Δx)
    ω = ω0

    # linear and angular acceleration (in the body frame)
    a = a0 + cross(α0, Δx) + cross(α0, u) - gravity
    α = α0

    return (; C, Qinv, mass11, mass12, mass21, mass22, u, θ, V, Ω, P, H, F, M, v, ω, a, α) 
end

"""
    steady_state_point_velocity_residuals(properties)

Calculate the velocity residuals `rV` and `rΩ` for a point for a steady state analysis.
"""
@inline function steady_state_point_velocity_residuals(properties)

    @unpack u, C, Qinv, V, Ω, v, ω = properties
    
    rV = V - v - cross(ω, u)
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
    
    rV = C'*V - v - cross(ω, u)
    rΩ = Qinv*(Ω - C*ω)

    return (; rV, rΩ)
end

"""
    static_point_resultants(properties)

Calculate the net loads `F` and `M` applied at a point for a static analysis.
"""
@inline function static_point_resultants(properties)

    @unpack C, mass11, mass12, mass21, mass22, F, M, a, α = properties

    # add loads due to linear and angular acceleration (including gravity)
    F -= C'*mass11*C*a + C'*mass12*C*α
    M -= C'*mass21*C*a + C'*mass22*C*α

    return F, M
end

"""
    steady_state_point_resultants(properties)

Calculate the net loads `F` and `M` applied at a point for a steady state analysis.
"""
@inline function steady_state_point_resultants(properties)

    F, M = static_point_resultants(properties)

    @unpack V, Ω, P, H, ω = properties  

    # add loads due to linear and angular momentum
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
    steady_state_point_residual!(resid, x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)

Calculate and insert the residual entries corresponding to a point for a steady state analysis into the 
system residual vector.
"""
@inline function steady_state_point_residual!(resid, x, indices, force_scaling, assembly, ipoint,  
    prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)

    irow = indices.irow_point[ipoint]

    properties = steady_state_point_properties(x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)

    F, M = steady_state_point_resultants(properties)

    rV, rΩ = steady_state_point_velocity_residuals(properties)

    resid[irow:irow+2] .= -F ./ force_scaling
    resid[irow+3:irow+5] .= -M ./ force_scaling
    resid[irow+6:irow+8] .= rV
    resid[irow+9:irow+11] .= rΩ

    return resid
end

"""
    initial_condition_point_residual!(resid, x, indices, rate_vars, force_scaling, 
        assembly, ipoint, prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0, 
        u0, θ0, V0, Ω0, Vdot0, Ωdot0)

Calculate and insert the residual entries corresponding to a point for the initialization 
of a time domain analysis into the system residual vector.
"""
@inline function initial_condition_point_residual!(resid, x, indices, rate_vars, force_scaling, 
    assembly, ipoint, prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0, 
    u0, θ0, V0, Ω0, Vdot0, Ωdot0)

    irow = indices.irow_point[ipoint]

    properties = initial_condition_point_properties(x, indices, rate_vars, force_scaling, 
        assembly, ipoint, prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0, 
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
        prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0, 
        udot_init, θdot_init, Vdot_init, Ωdot_init, dt)

Calculate and insert the residual entries corresponding to a point for a newmark-scheme time marching 
analysis into the system residual vector.
"""
@inline function newmark_point_residual!(resid, x, indices, force_scaling, assembly, ipoint,  
    prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0, 
    udot_init, θdot_init, Vdot_init, Ωdot_init, dt)

    irow = indices.irow_point[ipoint]

    properties = newmark_point_properties(x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0, 
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
        prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)

Calculate and insert the residual entries corresponding to a point for a dynamic analysis into the 
system residual vector.
"""
@inline function dynamic_point_residual!(resid, dx, x, indices, force_scaling, assembly, ipoint,  
    prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)

    irow = indices.irow_point[ipoint]

    properties = dynamic_point_properties(dx, x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)

    F, M = dynamic_point_resultants(properties)

    rV, rΩ = dynamic_point_velocity_residuals(properties)

    resid[irow:irow+2] .= -F ./ force_scaling
    resid[irow+3:irow+5] .= -M ./ force_scaling
    resid[irow+6:irow+8] .= rV
    resid[irow+9:irow+11] .= rΩ

    return resid
end

"""
    expanded_point_residual!(resid, x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)

Calculate and insert the residual entries corresponding to a point into the system residual 
vector for a constant mass matrix system.
"""
@inline function expanded_point_residual!(resid, x, indices, force_scaling, assembly, ipoint,  
    prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)

    irow = indices.irow_point[ipoint]

    properties = expanded_point_properties(x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)

    F, M = expanded_point_resultants(properties)

    rV, rΩ = expanded_point_velocity_residuals(properties)

    resid[irow:irow+2] .= -F ./ force_scaling
    resid[irow+3:irow+5] .= -M ./ force_scaling
    resid[irow+6:irow+8] .= rV
    resid[irow+9:irow+11] .= rΩ

    return resid
end

"""
    static_point_jacobian_properties(properties, x, indices, 
        force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity)

Calculate/extract the point properties needed to calculate the jacobian entries 
corresponding to a point for a static analysis
"""
@inline function static_point_jacobian_properties(properties, x, indices, 
    force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity)

    @unpack C, θ = properties

    # forces and moments
    F_θ, F_F, M_θ, M_M = point_load_jacobians(x, ipoint, indices.icol_point, force_scaling, prescribed_conditions)

    # linear and angular displacement
    u_u, θ_θ = point_displacement_jacobians(ipoint, prescribed_conditions)

    # rotation parameter matrices
    C_θ1, C_θ2, C_θ3 = get_C_θ(C, θ)

    return (; properties..., C_θ1, C_θ2, C_θ3, u_u, θ_θ, F_θ, M_θ, F_F, M_M)
end

"""
    steady_state_point_jacobian_properties(properties, x, indices, force_scaling, 
        assembly, ipoint, prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)

Calculate/extract the point properties needed to calculate the jacobian entries 
corresponding to a point for a steady state analysis
"""
@inline function steady_state_point_jacobian_properties(properties, x, indices, force_scaling, 
    assembly, ipoint, prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)

    properties = static_point_jacobian_properties(properties, x, indices, 
        force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity)

    @unpack C, mass11, mass12, mass21, mass22, θ, V, Ω, C_θ1, C_θ2, C_θ3 = properties

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

    # linear and angular acceleration
    a_u = tilde(α0)

    return (; properties..., Qinv_θ1, Qinv_θ2, Qinv_θ3, P_θ, P_V, P_Ω, H_θ, H_V, H_Ω, a_u) 
end

"""
    initial_condition_point_jacobian_properties(properties, x, indices, rate_vars, 
        force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity,
        x0, v0, ω0, a0, α0, u0, θ0, V0, Ω0, Vdot0, Ωdot0)

Calculate/extract the point properties needed to calculate the jacobian entries 
corresponding to a point for a Newmark scheme time marching analysis
"""
@inline function initial_condition_point_jacobian_properties(properties, x, indices, 
    rate_vars, force_scaling, assembly, ipoint, prescribed_conditions, 
    point_masses, gravity, x0, v0, ω0, a0, α0, u0, θ0, V0, Ω0, Vdot0, Ωdot0)

    @unpack C, Cdot, mass11, mass12, mass21, mass22, θ, V, Ω = properties

    # rotation parameter matrices
    C_θ1, C_θ2, C_θ3 = get_C_θ(C, θ)
    Qinv_θ1, Qinv_θ2, Qinv_θ3 = get_Qinv_θ(θ)

    # linear and angular displacement
    icol = indices.icol_point[ipoint]
    u_u, θ_θ = point_displacement_jacobians(ipoint, prescribed_conditions)
    if haskey(prescribed_conditions, ipoint)
        # linear displacement
        u_u = hcat(
            prescribed_conditions[ipoint].isforce[1] && rate_vars[icol+6] ? zero(e1) : u_u[:,1], 
            prescribed_conditions[ipoint].isforce[2] && rate_vars[icol+7] ? zero(e2) : u_u[:,2],
            prescribed_conditions[ipoint].isforce[3] && rate_vars[icol+8] ? zero(e3) : u_u[:,3],
        )
        # angular displacement
        θ_θ = hcat(
            prescribed_conditions[ipoint].isforce[4] && rate_vars[icol+9] ? zero(e1) : θ_θ[:,1],
            prescribed_conditions[ipoint].isforce[5] && rate_vars[icol+10] ? zero(e2) : θ_θ[:,2],
            prescribed_conditions[ipoint].isforce[6] && rate_vars[icol+11] ? zero(e3) : θ_θ[:,3], 
        )
    else
        # linear displacement
        u_u = hcat(
            rate_vars[icol+6] ? zero(e1) : u_u[:,1], 
            rate_vars[icol+7] ? zero(e2) : u_u[:,2],
            rate_vars[icol+8] ? zero(e3) : u_u[:,3],
        )
        # angular displacement
        θ_θ = hcat(
            rate_vars[icol+9] ? zero(e1) : θ_θ[:,1],
            rate_vars[icol+10] ? zero(e2) : θ_θ[:,2],
            rate_vars[icol+11] ? zero(e3) : θ_θ[:,3], 
        )
    end

    # forces and moments
    F_θ, F_F, M_θ, M_M = point_load_jacobians(x, ipoint, indices.icol_point, force_scaling, prescribed_conditions)

    # update force and moment jacobians
    F_θ = F_θ * θ_θ
    M_θ = M_θ * θ_θ

    # linear and angular velocity rates
    if haskey(prescribed_conditions, ipoint)
        Vdot_Vdot = hcat(
            prescribed_conditions[ipoint].isforce[1] && rate_vars[icol+6] ? e1 : zero(e1), 
            prescribed_conditions[ipoint].isforce[2] && rate_vars[icol+7] ? e2 : zero(e2),
            prescribed_conditions[ipoint].isforce[3] && rate_vars[icol+8] ? e3 : zero(e3),
        )
        Ωdot_Ωdot = hcat(
            prescribed_conditions[ipoint].isforce[4] && rate_vars[icol+9] ? e1 : zero(e1),
            prescribed_conditions[ipoint].isforce[5] && rate_vars[icol+10] ? e2 : zero(e2),
            prescribed_conditions[ipoint].isforce[6] && rate_vars[icol+11] ? e3 : zero(e3), 
        )
    else
        Vdot_Vdot = hcat(
            rate_vars[icol+6] ? e1 : zero(e1), 
            rate_vars[icol+7] ? e2 : zero(e2),
            rate_vars[icol+8] ? e3 : zero(e3),
        )
        Ωdot_Ωdot = hcat(
            rate_vars[icol+9] ? e1 : zero(e1),
            rate_vars[icol+10] ? e2 : zero(e2),
            rate_vars[icol+11] ? e3 : zero(e3), 
        )
    end

    # linear and angular momentum
    P_θ = mul3(C_θ1', C_θ2', C_θ3', (mass11*C*V + mass12*C*Ω)) + 
        C'*mass11*mul3(C_θ1, C_θ2, C_θ3, V) + 
        C'*mass12*mul3(C_θ1, C_θ2, C_θ3, Ω)

    H_θ = mul3(C_θ1', C_θ2', C_θ3', (mass21*C*V + mass22*C*Ω)) + 
        C'*mass21*mul3(C_θ1, C_θ2, C_θ3, V) + 
        C'*mass22*mul3(C_θ1, C_θ2, C_θ3, Ω)

    # linear and angular momentum rates
    Pdot_Vdot = C'*mass11*C
    Pdot_Ωdot = C'*mass12*C
    Hdot_Vdot = C'*mass21*C
    Hdot_Ωdot = C'*mass22*C

    # linear and angular acceleration
    a_u = tilde(α0)

    return (; properties..., C_θ1, C_θ2, C_θ3, Qinv_θ1, Qinv_θ2, Qinv_θ3, u_u, θ_θ, 
        F_θ, F_F, M_θ, M_M, P_θ, H_θ, Vdot_Vdot, Ωdot_Ωdot, Pdot_Vdot, Pdot_Ωdot, 
        Hdot_Vdot, Hdot_Ωdot, a_u) 
end

"""
    newmark_point_jacobian_properties(properties, x, indices, force_scaling, 
        assembly, ipoint, prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0,
        udot_init, θdot_init, Vdot_init, Ωdot_init, dt)

Calculate/extract the point properties needed to calculate the jacobian entries 
corresponding to a point for a Newmark scheme time marching analysis
"""
@inline function newmark_point_jacobian_properties(properties, x, indices, force_scaling, 
    assembly, ipoint, prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0,
    udot_init, θdot_init, Vdot_init, Ωdot_init, dt)

    properties = steady_state_point_jacobian_properties(properties, x, indices, 
        force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity, 
        x0, v0, ω0, a0, α0)

    @unpack C, Cdot, mass11, mass12, mass21, mass22, V, Ω, Vdot, Ωdot, v, ω, C_θ1, C_θ2, C_θ3 = properties

    # linear and angular deflection rates
    udot_u = 2/dt*I3
    θdot_θ = 2/dt*I3

    # linear and angular momentum rates
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

    return (; properties..., udot_u, θdot_θ, Pdot_θ, Pdot_V, Pdot_Ω, Hdot_θ, Hdot_V, Hdot_Ω) 
end

"""
    dynamic_point_jacobian_properties(properties, dx, x, indices, force_scaling, 
        assembly, ipoint, prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)

Calculate/extract the point properties needed to calculate the jacobian entries 
corresponding to a point for a dynamic analysis
"""
@inline function dynamic_point_jacobian_properties(properties, dx, x, indices, force_scaling, 
    assembly, ipoint, prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)

    properties = steady_state_point_jacobian_properties(properties, x, indices, 
        force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity, 
        x0, v0, ω0, a0, α0)

    @unpack C, Cdot, mass11, mass12, mass21, mass22, V, Ω, Vdot, Ωdot, v, ω, C_θ1, C_θ2, C_θ3 = properties

    # linear and angular momentum rates

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

    return (; properties..., Pdot_θ, Pdot_V, Pdot_Ω, Hdot_θ, Hdot_V, Hdot_Ω) 
end

"""
    expanded_point_jacobian_properties(properties, x, indices, force_scaling, 
        assembly, ipoint, prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)

Calculate/extract the point properties needed to calculate the jacobian entries 
corresponding to a point for a constant mass matrix system
"""
@inline function expanded_point_jacobian_properties(properties, x, indices, force_scaling, 
    assembly, ipoint, prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)

    @unpack C, mass11, mass12, mass21, mass22, θ, V, Ω = properties

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

    # linear and angular acceleration
    a_u = tilde(α0)

    return (; properties..., C_θ1, C_θ2, C_θ3, Qinv_θ1, Qinv_θ2, Qinv_θ3, u_u, θ_θ, 
        F_θ, F_F, M_θ, M_M, P_V, P_Ω, H_V, H_Ω, a_u) 
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
        u_u, θ_θ, P_θ, P_V, P_Ω, H_θ, H_V, H_Ω, a_u = properties

    @unpack F_θ, M_θ = jacobians

    # add loads due to linear and angular acceleration (including gravity)
    F_u = -C'*mass11*C*a_u*u_u
    M_u = -C'*mass21*C*a_u*u_u

    # add loads due to linear and angular momentum
    F_θ -= tilde(ω)*P_θ*θ_θ
    F_V = -tilde(ω)*P_V
    F_Ω = -tilde(ω)*P_Ω

    M_θ -= (tilde(ω)*H_θ + tilde(V)*P_θ)*θ_θ
    M_V = -tilde(ω)*H_V - tilde(V)*P_V + tilde(P)
    M_Ω = -tilde(ω)*H_Ω - tilde(V)*P_Ω

    return (; jacobians..., F_u, F_θ, F_V, F_Ω, M_u, M_θ, M_V, M_Ω)
end

"""
    initial_condition_point_resultant_jacobians(properties)

Calculate the jacobians for the net loads `F` and `M` applied at a point for the initialization
of a time domain analysis.
"""
@inline function initial_condition_point_resultant_jacobians(properties)

    jacobians = static_point_resultant_jacobians(properties)

    @unpack  C, mass11, mass21, V, Ω, P, H, ω, u_u, θ_θ, P_θ, H_θ, 
        Pdot_Vdot, Pdot_Ωdot, Hdot_Vdot, Hdot_Ωdot, Vdot_Vdot, Ωdot_Ωdot, a_u = properties

    @unpack F_θ, M_θ = jacobians

    # add loads due to linear and angular acceleration (including gravity)
    F_u = -C'*mass11*C*a_u*u_u
    M_u = -C'*mass21*C*a_u*u_u

    # add loads due to linear and angular momentum
    F_θ -= tilde(ω)*P_θ*θ_θ
    M_θ -= (tilde(ω)*H_θ + tilde(V)*P_θ)*θ_θ

    # add loads due to linear and angular momentum rates
    F_Vdot = -Pdot_Vdot*Vdot_Vdot
    F_Ωdot = -Pdot_Ωdot*Ωdot_Ωdot

    M_Vdot = -Hdot_Vdot*Vdot_Vdot
    M_Ωdot = -Hdot_Ωdot*Ωdot_Ωdot

    return (; jacobians..., F_u, F_θ, F_Vdot, F_Ωdot, M_u, M_θ, M_Vdot, M_Ωdot)
end

"""
    newmark_point_resultant_jacobians(properties)

Calculate the jacobians for the net loads `F` and `M` applied at a point for a Newmark 
scheme time marching analysis.
"""
@inline function newmark_point_resultant_jacobians(properties)

    jacobians = steady_state_point_resultant_jacobians(properties)

    @unpack F_θ, F_V, F_Ω, M_θ, M_V, M_Ω = jacobians

    @unpack θ_θ, Pdot_θ, Pdot_V, Pdot_Ω, Hdot_θ, Hdot_V, Hdot_Ω = properties

    # add loads due to linear and angular momentum rates
    F_θ -= Pdot_θ*θ_θ
    F_V -= Pdot_V
    F_Ω -= Pdot_Ω
    
    M_θ -= Hdot_θ*θ_θ
    M_V -= Hdot_V
    M_Ω -= Hdot_Ω

    return (; jacobians..., F_θ, F_V, F_Ω, M_θ, M_V, M_Ω)
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
        C_θ1, C_θ2, C_θ3, u_u, θ_θ, P_V, P_Ω, H_V, H_Ω, a_u = properties

    @unpack F_θ, M_θ, F_F, M_M = properties

    # rotate loads into the appropriate frame
    F_θ = (mul3(C_θ1, C_θ2, C_θ3, F) + C*F_θ)*θ_θ
    F_F = C*F_F

    M_θ = (mul3(C_θ1, C_θ2, C_θ3, M) + C*M_θ)*θ_θ
    M_M = C*M_M

    # add loads due to linear and angular acceleration (including gravity)
    F_u = -mass11*C*a_u*u_u
    M_u = -mass21*C*a_u*u_u

    F_θ -= (mass11*mul3(C_θ1, C_θ2, C_θ3, a) + mass12*mul3(C_θ1, C_θ2, C_θ3, α))*θ_θ
    M_θ -= (mass21*mul3(C_θ1, C_θ2, C_θ3, a) + mass22*mul3(C_θ1, C_θ2, C_θ3, α))*θ_θ

    # add loads due to linear and angular momentum
    F_V = -tilde(Ω)*P_V 
    F_Ω = -tilde(Ω)*P_Ω + tilde(P)

    M_V = -tilde(Ω)*H_V - tilde(V)*P_V + tilde(P)
    M_Ω = -tilde(Ω)*H_Ω - tilde(V)*P_Ω + tilde(H)

    (; F_F, F_u, F_θ, F_V, F_Ω, M_M, M_u, M_θ, M_V, M_Ω)
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

    @unpack C, Qinv, Ω, ω, u_u, θ_θ, C_θ1, C_θ2, C_θ3, Qinv_θ1, Qinv_θ2, Qinv_θ3 = properties
    
    rV_u = -tilde(ω)*u_u
    rV_V = I3

    rΩ_θ = (mul3(Qinv_θ1, Qinv_θ2, Qinv_θ3, C*(Ω - ω)) + Qinv*mul3(C_θ1, C_θ2, C_θ3, Ω - ω))*θ_θ
    rΩ_Ω = Qinv*C


    return (; rV_u, rV_V, rΩ_θ, rΩ_Ω)
end

"""
    initial_condition_point_velocity_jacobians(properties)

Calculate the jacobians of the velocity residuals `rV` and `rΩ` of a point for the
initialization of a time domain analysis.
"""
@inline function initial_condition_point_velocity_jacobians(properties)
   
    @unpack C, Qinv, Ω, ω, u_u, θ_θ, C_θ1, C_θ2, C_θ3, Qinv_θ1, Qinv_θ2, Qinv_θ3 = properties

    rV_u = -tilde(ω)*u_u
    rΩ_θ = (mul3(Qinv_θ1, Qinv_θ2, Qinv_θ3, C*(Ω - ω)) + Qinv*mul3(C_θ1, C_θ2, C_θ3, Ω - ω))*θ_θ

    rV_udot = -I3
    rΩ_θdot = -I3

    return (; rV_u, rV_udot, rΩ_θ, rΩ_θdot)
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

    @unpack C, Qinv, V, Ω, ω, u_u, θ_θ, C_θ1, C_θ2, C_θ3, Qinv_θ1, Qinv_θ2, Qinv_θ3 = properties
    
    rV_u = -tilde(ω)*u_u
    rV_θ = mul3(C_θ1', C_θ2', C_θ3', V)*θ_θ
    rV_V = C'

    rΩ_θ = (mul3(Qinv_θ1, Qinv_θ2, Qinv_θ3, Ω - C*ω) - Qinv*mul3(C_θ1, C_θ2, C_θ3, ω))*θ_θ
    rΩ_Ω = Qinv

    return (; rV_u, rV_θ, rV_V, rΩ_θ, rΩ_Ω)
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

    @unpack F_u, F_θ, F_F, F_Vdot, F_Ωdot, M_u, M_θ, M_M, M_Vdot, M_Ωdot = resultants
    @unpack rV_u, rV_udot, rΩ_θ, rΩ_θdot = velocities

    irow = indices.irow_point[ipoint]
    icol = indices.icol_point[ipoint]

    jacob[irow:irow+2, icol:icol+2] .= -F_F .- F_u ./ force_scaling .- F_Vdot ./ force_scaling
    jacob[irow:irow+2, icol+3:icol+5] .= -F_θ ./ force_scaling .- F_Ωdot ./ force_scaling

    jacob[irow+3:irow+5, icol:icol+2] .= -M_u ./ force_scaling .- M_Vdot ./ force_scaling 
    jacob[irow+3:irow+5, icol+3:icol+5] .= -M_M .- M_θ ./ force_scaling .- M_Ωdot ./ force_scaling

    jacob[irow+6:irow+8, icol:icol+2] .= rV_u
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

    @unpack F_u, F_V, F_Ω, M_u, M_θ, M_V, M_Ω = resultants
    @unpack rV_u, rV_V, rΩ_θ, rΩ_Ω = velocities

    irow = indices.irow_point[ipoint]
    icol = indices.icol_point[ipoint]

    @views jacob[irow:irow+2, icol:icol+2] .-= F_u ./ force_scaling
    jacob[irow:irow+2, icol+6:icol+8] .= -F_V ./ force_scaling
    jacob[irow:irow+2, icol+9:icol+11] .= -F_Ω ./ force_scaling

    @views jacob[irow+3:irow+5, icol:icol+2] .-= M_u ./ force_scaling
    jacob[irow+3:irow+5, icol+6:icol+8] .= -M_V ./ force_scaling
    jacob[irow+3:irow+5, icol+9:icol+11] .= -M_Ω ./ force_scaling

    jacob[irow+6:irow+8, icol:icol+2] .= rV_u
    jacob[irow+6:irow+8, icol+6:icol+8] .= rV_V

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
    steady_state_point_jacobian!(jacob, x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)

Calculate and insert the jacobian entries corresponding to a point for a steady state 
analysis into the system jacobian matrix.
"""
@inline function steady_state_point_jacobian!(jacob, x, indices, force_scaling, assembly, ipoint,  
    prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)

    properties = steady_state_point_properties(x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)

    properties = steady_state_point_jacobian_properties(properties, x, indices, 
        force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity, 
        x0, v0, ω0, a0, α0)

    resultants = steady_state_point_resultant_jacobians(properties)

    velocities = steady_state_point_velocity_jacobians(properties)

    insert_dynamic_point_jacobians!(jacob, indices, force_scaling, ipoint, resultants, velocities)

    return jacob
end

"""
    initial_condition_point_jacobian!(jacob, x, indices, rate_vars, force_scaling, 
        assembly, ipoint, prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0,
        u0, θ0, V0, Ω0, Vdot0, Ωdot0)

Calculate and insert the jacobian entries corresponding to a point for the initialization
of a time domain analysis into the system jacobian matrix.
"""
@inline function initial_condition_point_jacobian!(jacob, x, indices, rate_vars, force_scaling, 
    assembly, ipoint, prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0,
    u0, θ0, V0, Ω0, Vdot0, Ωdot0)

    properties = initial_condition_point_properties(x, indices, rate_vars, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0,
        u0, θ0, V0, Ω0, Vdot0, Ωdot0)
    
    properties = initial_condition_point_jacobian_properties(properties, x, indices, 
        rate_vars, force_scaling, assembly, ipoint, prescribed_conditions, point_masses, 
        gravity, x0, v0, ω0, a0, α0, u0, θ0, V0, Ω0, Vdot0, Ωdot0)

    resultants = initial_condition_point_resultant_jacobians(properties)

    velocities = initial_condition_point_velocity_jacobians(properties)

    insert_initial_condition_point_jacobians!(jacob, indices, force_scaling, ipoint, resultants, velocities)

    return jacob
end

"""
    newmark_point_jacobian!(jacob, x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)

Calculate and insert the jacobian entries corresponding to a point for a Newmark scheme
time marching analysis into the system jacobian matrix.
"""
@inline function newmark_point_jacobian!(jacob, x, indices, force_scaling, assembly, ipoint,  
    prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0,
    udot_init, θdot_init, Vdot_init, Ωdot_init, dt)

    properties = newmark_point_properties(x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0,
        udot_init, θdot_init, Vdot_init, Ωdot_init, dt)

    properties = newmark_point_jacobian_properties(properties, x, indices, force_scaling, 
        assembly, ipoint, prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0,
        udot_init, θdot_init, Vdot_init, Ωdot_init, dt)

    resultants = newmark_point_resultant_jacobians(properties)

    velocities = newmark_point_velocity_jacobians(properties)

    insert_dynamic_point_jacobians!(jacob, indices, force_scaling, ipoint, resultants, velocities)

    return jacob
end

"""
    dynamic_point_jacobian!(jacob, dx, x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)

Calculate and insert the jacobian entries corresponding to a point for a dynamic 
analysis into the system jacobian matrix.
"""
@inline function dynamic_point_jacobian!(jacob, dx, x, indices, force_scaling, assembly, ipoint,  
    prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)

    properties = dynamic_point_properties(dx, x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)

    properties = dynamic_point_jacobian_properties(properties, dx, x, indices, force_scaling, 
        assembly, ipoint, prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)

    resultants = dynamic_point_resultant_jacobians(properties)

    velocities = dynamic_point_velocity_jacobians(properties)

    insert_dynamic_point_jacobians!(jacob, indices, force_scaling, ipoint, resultants, velocities)

    return jacob
end

"""
    expanded_point_jacobian!(jacob, x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)

Calculate and insert the jacobian entries corresponding to a point for a constant mass 
matrix system into the system jacobian matrix.
"""
@inline function expanded_point_jacobian!(jacob, x, indices, force_scaling, assembly, ipoint,  
    prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)

    properties = expanded_point_properties(x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)

    properties = expanded_point_jacobian_properties(properties, x, indices, 
        force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity, 
        x0, v0, ω0, a0, α0)

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