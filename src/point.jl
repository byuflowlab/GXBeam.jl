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
    point_velocity_rates(dx, ipoint, icol, prescribed_conditions)

Extract the velocity rates `Vdot` and `Ωdot` of point `ipoint` from the state variable vector
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
    point_displacement_jacobians(x, ipoint, icol, prescribed_conditions)

Calculate the displacement jacobians `u_u` and `θ_θ` of point `ipoint`.
"""
@inline function point_displacement_jacobians(x, ipoint, icol, prescribed_conditions)

    if haskey(prescribed_conditions, ipoint)
        u_u, θ_θ = point_displacement_jacobians(x, icol[ipoint], prescribed_conditions[ipoint])
    else
        u_u, θ_θ = point_displacement_jacobians(x, icol[ipoint])
    end

    return u_u, θ_θ
end

@inline function point_displacement_jacobians(x, icol, prescribed_conditions)

    @unpack isforce = prescribed_conditions

    u_u = hcat(ifelse(isforce[1], e1, zero(e1)),
               ifelse(isforce[2], e2, zero(e2)),
               ifelse(isforce[3], e3, zero(e3)))

    θ_θ = hcat(ifelse(isforce[4], e1, zero(e1)),
               ifelse(isforce[5], e2, zero(e2)),
               ifelse(isforce[6], e3, zero(e3)))

    return u_u, θ_θ
end

@inline function point_displacement_jacobians(x, icol)

    u_u = I3
    θ_θ = I3

    return u_u, θ_θ
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

    F_F = hcat(ifelse(isforce[1], zero(e1), e1),
               ifelse(isforce[2], zero(e2), e2),
               ifelse(isforce[3], zero(e3), e3))

    M_M = hcat(ifelse(isforce[4], zero(e1), e1),
               ifelse(isforce[5], zero(e2), e2),
               ifelse(isforce[6], zero(e3), e3))

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
    F1_θ = ifelse(isforce[1], SVector(Fp_θ[1,1], Fp_θ[1,2], Fp_θ[1,3]), (@SVector zeros(3)))
    F2_θ = ifelse(isforce[2], SVector(Fp_θ[2,1], Fp_θ[2,2], Fp_θ[2,3]), (@SVector zeros(3)))
    F3_θ = ifelse(isforce[3], SVector(Fp_θ[3,1], Fp_θ[3,2], Fp_θ[3,3]), (@SVector zeros(3)))
    F_θ = vcat(F1_θ', F2_θ', F3_θ')

    M1_θ = ifelse(isforce[4], SVector(Mp_θ[1,1], Mp_θ[1,2], Mp_θ[1,3]), (@SVector zeros(3)))
    M2_θ = ifelse(isforce[5], SVector(Mp_θ[2,1], Mp_θ[2,2], Mp_θ[2,3]), (@SVector zeros(3)))
    M3_θ = ifelse(isforce[6], SVector(Mp_θ[3,1], Mp_θ[3,2], Mp_θ[3,3]), (@SVector zeros(3)))
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

    # point displacement
    u, θ = point_displacement(x, ipoint, indices.icol_point, prescribed_conditions)

    # transformation matrices
    C = get_C(θ)
    Ct = C'

    # externally applied forces and moments
    F, M = point_loads(x, ipoint, indices.icol_point, force_scaling, prescribed_conditions)

    # effective linear and angular acceleration
    a = -gravity
    α = zero(a)

    return (; C, Ct, mass11, mass12, mass21, mass22, u, θ, F, M, a, α)
end

"""
    steady_state_point_properties(x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)

Calculate/extract the point properties needed to construct the residual for a steady state analysis
"""
@inline function steady_state_point_properties(x, indices, force_scaling, assembly, ipoint,  
    prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)

    properties = static_point_properties(x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity)

    @unpack u, θ, C, Ct, mass11, mass12, mass21, mass22, a, α = properties

    # point velocities
    V, Ω = point_velocities(x, ipoint, indices.icol_point)

    # linear and angular momentum
    P = Ct*(mass11*C*V + mass12*C*Ω)
    H = Ct*(mass21*C*V + mass22*C*Ω)

    # undeformed velocities
    v = v0 + cross(ω0, assembly.points[ipoint] - x0)
    ω = ω0

    # effective linear and angular acceleration
    a += a0 + cross(α0, assembly.points[ipoint] - x0) + cross(α0, u)
    α += α0 

    return (; properties..., V, Ω, P, H, v, ω, a, α) 
end

"""
    initial_condition_point_properties(x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0,
        u0, θ0, udot0, θdot0)

Calculate/extract the point properties needed to construct the residual for a time domain 
analysis initialization.
"""
@inline function initial_condition_point_properties(x, indices, force_scaling, assembly, ipoint,  
    prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0,
    u0, θ0, udot0, θdot0)

    # mass matrix
    mass = haskey(point_masses, ipoint) ? point_masses[ipoint].mass : @SMatrix zeros(6,6)
    
    # mass submatrices
    mass11 = mass[SVector{3}(1:3), SVector{3}(1:3)]
    mass12 = mass[SVector{3}(1:3), SVector{3}(4:6)]
    mass21 = mass[SVector{3}(4:6), SVector{3}(1:3)]
    mass22 = mass[SVector{3}(4:6), SVector{3}(4:6)]

    # point displacement
    u = SVector{3}(u0[ipoint])
    θ = SVector{3}(θ0[ipoint])

    # transformation matrices
    C = get_C(θ)
    Ct = C'
    Q = get_Q(θ)

    # externally applied forces and moments
    F, M = point_loads(x, ipoint, indices.icol_point, force_scaling, prescribed_conditions)
    
    # point velocities
    V, Ω = point_velocities(x, ipoint, indices.icol_point)
    
    # linear and angular momentum
    P = Ct*(mass11*C*V + mass12*C*Ω)
    H = Ct*(mass21*C*V + mass22*C*Ω)

    # undeformed body frame velocities
    v = v0 + cross(ω0, assembly.points[ipoint] - x0)
    ω = ω0

    # effective linear and angular acceleration
    a = a0 - gravity + cross(α0, assembly.points[ipoint] - x0) + cross(α0, u)
    α = α0 

    # point displacement rates
    udot = SVector{3}(udot0[ipoint])
    θdot = SVector{3}(θdot0[ipoint])

    # point velocity rates
    Vdot, Ωdot = point_displacement_rates(x, ipoint, indices.icol_point, prescribed_conditions)

    # linear and angular momentum rates
    Pdot = Ct*mass11*C*Vdot + Ct*mass12*C*Ωdot
    Hdot = Ct*mass21*C*Vdot + Ct*mass22*C*Ωdot

    return (; C, Ct, Q, mass11, mass12, mass21, mass22, u, θ, F, M, V, Ω, P, H, v, ω, a, α, udot, θdot, 
        Vdot, Ωdot, Pdot, Hdot) 
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

    @unpack C, Ct, mass11, mass12, mass21, mass22, u, θ, V, Ω, P, H, ω = properties

    # transformation matrices
    Q = get_Q(θ)

    # point displacement rates
    udot = 2/dt*u - SVector{3}(udot_init[ipoint])
    θdot = 2/dt*θ - SVector{3}(θdot_init[ipoint])

    # point velocity rates
    Vdot = 2/dt*V - SVector{3}(Vdot_init[ipoint])
    Ωdot = 2/dt*Ω - SVector{3}(Ωdot_init[ipoint])

    # linear and angular momentum rates
    Pdot = Ct*mass11*C*Vdot + Ct*mass12*C*Ωdot
    Hdot = Ct*mass21*C*Vdot + Ct*mass22*C*Ωdot

    return (; properties..., Q, udot, θdot, Vdot, Ωdot, Pdot, Hdot) 
end

"""
    dynamic_point_properties(dx, x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)

Calculate/extract the point properties needed to construct the residual for a dynamic analysis
"""
@inline function dynamic_point_properties(dx, x, indices, force_scaling, assembly, ipoint,  
    prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)

    properties = steady_state_point_properties(x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)

    @unpack C, Ct, mass11, mass12, mass21, mass22, θ, V, Ω, P, H, ω = properties

    # transformation matrices
    Q = get_Q(θ)

    # point displacement rates
    udot, θdot = point_displacement_rates(dx, ipoint, indices.icol_point, prescribed_conditions)

    # point velocity rates
    Vdot, Ωdot = point_velocity_rates(dx, ipoint, indices.icol_point, prescribed_conditions)

    # linear and angular momentum rates
    Pdot = Ct*mass11*C*Vdot + Ct*mass12*C*Ωdot
    Hdot = Ct*mass21*C*Vdot + Ct*mass22*C*Ωdot

    return (; properties..., Q, udot, θdot, Vdot, Ωdot, Pdot, Hdot) 
end

"""
    static_point_resultants(properties)

Calculate the net loads `F` and `M` applied at a point for a static analysis.
"""
@inline function static_point_resultants(properties)

    @unpack Ct, mass11, mass12, mass21, mass22, F, M, a, α = properties

    # externally applied loads
    Fe, Me = F, M

    # loads due to acceleration (including gravity)
    Fa, Ma = acceleration_loads(Ct, mass11, mass12, mass21, mass22, a, α)

    # load resultants
    F = Fe + Fa
    M = Me + Ma

    return F, M
end

"""
    steady_state_point_resultants(properties)

Calculate the net loads `F` and `M` applied at a point for a steady state analysis.
"""
@inline function steady_state_point_resultants(properties)

    F, M = static_point_resultants(properties)

    @unpack Ct, V, Ω, P, H, ω = properties

    # loads due to linear and angular momentum
    F += cross(ω, P)
    M += cross(ω, H) + cross(V, P)

    return F, M
end

"""
    dynamic_point_resultants(properties)

Calculate the net loads `F` and `M` applied at a point for a dynamic analysis.
"""
@inline function dynamic_point_resultants(properties)

    F, M = steady_state_point_resultants(properties)

    @unpack Pdot, Hdot = properties

    # loads due to linear and angular momentum rates
    F += Pdot
    M += Hdot

    return F, M
end

"""
    steady_state_velocity_residuals(properties)

Calculate the velocity residuals `rV` and `rΩ` for a point for a steady state analysis.
"""
@inline function steady_state_velocity_residuals(properties)

    @unpack u, θ, V, Ω, v, ω = properties
    
    rV = V - v - cross(ω, u)
    rΩ = Ω - ω

    return rV, rΩ
end

"""
    dynamic_velocity_residuals(properties)

Calculate the velocity residuals `rV` and `rΩ` for a point for a dynamic analysis.
"""
@inline function dynamic_velocity_residuals(properties)

    rV, rΩ = steady_state_velocity_residuals(properties)

    @unpack Ct, Q, udot, θdot = properties
    
    rV -= udot
    rΩ -= Ct*Q*θdot

    return rV, rΩ
end

"""
    static_point_residual!(resid, x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity)

Calculate and insert the residual entries corresponding to a point for a static analysis into the system 
residual vector.
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

    rV, rΩ = steady_state_velocity_residuals(properties)

    resid[irow:irow+2] .= -F ./ force_scaling
    resid[irow+3:irow+5] .= -M ./ force_scaling
    resid[irow+6:irow+8] .= rV
    resid[irow+9:irow+11] .= rΩ

    return resid
end

"""
    initial_condition_point_residual!(resid, x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0, u0, θ0, 
        udot0, θdot0)

Calculate and insert the residual entries corresponding to a point for the initialization 
of a time domain analysis into the system residual vector.
"""
@inline function initial_condition_point_residual!(resid, x, indices, force_scaling, assembly, ipoint,  
    prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0, u0, θ0, 
    udot0, θdot0)

    irow = indices.irow_point[ipoint]

    properties = initial_condition_point_properties(x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0, u0, θ0, 
        udot0, θdot0)

    F, M = dynamic_point_resultants(properties)

    rV, rΩ = dynamic_velocity_residuals(properties)

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

    rV, rΩ = dynamic_velocity_residuals(properties)

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

    rV, rΩ = dynamic_velocity_residuals(properties)

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

    @unpack θ, C = properties

    # point displacement
    u_u, θ_θ = point_displacement_jacobians(x, ipoint, indices.icol_point, prescribed_conditions)

    # transformation matrices
    C_θ1, C_θ2, C_θ3 = get_C_θ(C, θ)
    Ct_θ1, Ct_θ2, Ct_θ3 = C_θ1', C_θ2', C_θ3'

    # externally applied forces and moments
    F_θ, F_F, M_θ, M_M = point_load_jacobians(x, ipoint, indices.icol_point, force_scaling, prescribed_conditions)

    return (; properties..., C_θ1, C_θ2, C_θ3, Ct_θ1, Ct_θ2, Ct_θ3, u_u, θ_θ, F_θ, F_F, M_θ, M_M)
end

"""
    steady_state_point_jacobian_properties(properties, x, indices, force_scaling, 
        assembly, ipoint, prescribed_conditions, point_masses, gravity)

Calculate/extract the point properties needed to calculate the jacobian entries 
corresponding to a point for a steady state analysis
"""
@inline function steady_state_point_jacobian_properties(properties, x, indices, force_scaling, 
    assembly, ipoint, prescribed_conditions, point_masses, gravity)

    properties = static_point_jacobian_properties(properties, x, indices, 
        force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity)

    @unpack C, Ct, mass11, mass12, mass21, mass22, V, Ω, C_θ1, C_θ2, C_θ3, Ct_θ1, Ct_θ2, Ct_θ3 = properties

    # linear momentum
    P_θ = mul3(Ct_θ1, Ct_θ2, Ct_θ3, mass11*C*V + mass12*C*Ω) + 
        Ct*mass11*mul3(C_θ1, C_θ2, C_θ3, V) +
        Ct*mass12*mul3(C_θ1, C_θ2, C_θ3, Ω)
    P_V = Ct*mass11*C
    P_Ω = Ct*mass12*C

    # angular momentum
    H_θ = mul3(Ct_θ1, Ct_θ2, Ct_θ3, mass21*C*V + mass22*C*Ω) + 
        Ct*mass21*mul3(C_θ1, C_θ2, C_θ3, V) +
        Ct*mass22*mul3(C_θ1, C_θ2, C_θ3, Ω)
    H_V = Ct*mass21*C
    H_Ω = Ct*mass22*C

    return (; properties..., P_θ, P_V, P_Ω, H_θ, H_V, H_Ω) 
end

"""
    initial_condition_point_jacobian_properties(properties, x, indices, force_scaling, 
        assembly, ipoint, prescribed_conditions, point_masses, gravity)

Calculate/extract the point properties needed to calculate the jacobian entries 
corresponding to a point for a Newmark scheme time marching analysis
"""
@inline function initial_condition_point_jacobian_properties(properties, x, indices, force_scaling, 
    assembly, ipoint, prescribed_conditions, point_masses, gravity)

    @unpack C, Ct, mass11, mass12, mass21, mass22 = properties

    # externally applied forces and moments
    _, F_F, _, M_M = point_load_jacobians(x, ipoint, indices.icol_point, force_scaling, prescribed_conditions)

    # linear momentum
    P_V = Ct*mass11*C
    P_Ω = Ct*mass12*C

    # angular momentum
    H_V = Ct*mass21*C
    H_Ω = Ct*mass22*C

    # point displacement rates
    Vdot_Vdot, Ωdot_Ωdot = point_displacement_jacobians(x, ipoint, indices.icol_point, prescribed_conditions)

    # linear and angular momentum rates
    Pdot_Vdot = Ct*mass11*C*Vdot_Vdot
    Pdot_Ωdot = Ct*mass12*C*Ωdot_Ωdot

    Hdot_Vdot = Ct*mass21*C*Vdot_Vdot
    Hdot_Ωdot = Ct*mass22*C*Ωdot_Ωdot

    return (; properties..., F_F, M_M, P_V, P_Ω, H_V, H_Ω, Pdot_Vdot, Pdot_Ωdot, Hdot_Vdot, Hdot_Ωdot)
end

"""
    newmark_point_jacobian_properties(properties, x, indices, force_scaling, 
        assembly, ipoint, prescribed_conditions, point_masses, gravity,
        udot_init, θdot_init, Vdot_init, Ωdot_init, dt)

Calculate/extract the point properties needed to calculate the jacobian entries 
corresponding to a point for a Newmark scheme time marching analysis
"""
@inline function newmark_point_jacobian_properties(properties, x, indices, force_scaling, 
    assembly, ipoint, prescribed_conditions, point_masses, gravity,
    udot_init, θdot_init, Vdot_init, Ωdot_init, dt)

    properties = steady_state_point_jacobian_properties(properties, x, indices, 
        force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity)

    @unpack C, Ct, Q, mass11, mass12, mass21, mass22, θ, V, Ω, P, H, 
        θdot, Vdot, Ωdot, P_θ, H_θ, P_V, P_Ω, H_V, H_Ω, 
        C_θ1, C_θ2, C_θ3, Ct_θ1, Ct_θ2, Ct_θ3 = properties

    # transformation matrices
    Q_θ1, Q_θ2, Q_θ3 = get_Q_θ(Q, θ)

    # point displacement rates
    udot_u = 2/dt*I
    θdot_θ = 2/dt*I

    # rotation matrix derivatives
    # TODO

    # point velocity rates
    Vdot_V = 2/dt*I
    Ωdot_Ω = 2/dt*I

    # linear and angular momentum rates
    Pdot_θ = mul3(Ct_θ1, Ct_θ2, Ct_θ3, mass11*C*Vdot) + Ct*mass11*mul3(C_θ1, C_θ2, C_θ3, Vdot) +
        mul3(Ct_θ1, Ct_θ2, Ct_θ3, mass12*C*Ωdot) + Ct*mass12*mul3(C_θ1, C_θ2, C_θ3, Ωdot)
    Pdot_V = Ct*mass11*C*Vdot_V
    Pdot_Ω = Ct*mass12*C*Ωdot_Ω

    Hdot_θ = mul3(Ct_θ1, Ct_θ2, Ct_θ3, mass21*(C*Vdot)) + Ct*mass21*mul3(C_θ1, C_θ2, C_θ3, Vdot) +
        mul3(Ct_θ1, Ct_θ2, Ct_θ3, mass22*(C*Ωdot)) + Ct*mass22*mul3(C_θ1, C_θ2, C_θ3, Ωdot)
    Hdot_V = Ct*mass21*C*Vdot_V
    Hdot_Ω = Ct*mass22*C*Ωdot_Ω

    return (; properties..., Q_θ1, Q_θ2, Q_θ3, udot_u, θdot_θ, 
        Pdot_θ, Pdot_V, Pdot_Ω, Hdot_θ, Hdot_V, Hdot_Ω) 
end

"""
    dynamic_point_jacobian_properties(properties, dx, x, indices, force_scaling, 
        assembly, ipoint, prescribed_conditions, point_masses, gravity)

Calculate/extract the point properties needed to calculate the jacobian entries 
corresponding to a point for a dynamic analysis
"""
@inline function dynamic_point_jacobian_properties(properties, dx, x, indices, force_scaling, 
    assembly, ipoint, prescribed_conditions, point_masses, gravity)

    properties = steady_state_point_jacobian_properties(properties, x, indices, 
        force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity)

    @unpack C, Ct, Q, mass11, mass12, mass21, mass22, θ, V, Ω, P, H, 
        θdot, Vdot, Ωdot, P_θ, H_θ, P_V, P_Ω, H_V, H_Ω, 
        C_θ1, C_θ2, C_θ3, Ct_θ1, Ct_θ2, Ct_θ3 = properties

    # transformation matrices
    Q_θ1, Q_θ2, Q_θ3 = get_Q_θ(Q, θ)

    # linear and angular momentum rates
    Pdot_θ = mul3(Ct_θ1, Ct_θ2, Ct_θ3, mass11*C*Vdot) + Ct*mass11*mul3(C_θ1, C_θ2, C_θ3, Vdot) +
        mul3(Ct_θ1, Ct_θ2, Ct_θ3, mass12*C*Ωdot) + Ct*mass12*mul3(C_θ1, C_θ2, C_θ3, Ωdot)
    
    Hdot_θ = mul3(Ct_θ1, Ct_θ2, Ct_θ3, mass21*(C*Vdot)) + Ct*mass21*mul3(C_θ1, C_θ2, C_θ3, Vdot) +
        mul3(Ct_θ1, Ct_θ2, Ct_θ3, mass22*C*Ωdot) + Ct*mass22*mul3(C_θ1, C_θ2, C_θ3, Ωdot)

    return (; properties..., Q_θ1, Q_θ2, Q_θ3, Pdot_θ, Hdot_θ) 
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

    # point displacement
    u, θ = point_displacement(x, ipoint, indices.icol_point, prescribed_conditions)
    
    # transformation matrices
    C = get_C(θ)
    Ct = C'
    Q = get_Q(θ)

    # point velocities
    V, Ω = point_velocities(x, ipoint, indices.icol_point)

    # linear and angular momentum
    P = Ct*(mass11*C*V + mass12*C*Ω)
    H = Ct*(mass21*C*V + mass22*C*Ω)

    # point displacement rates
    udot_udot, θdot_θdot = point_displacement_jacobians(x, ipoint, indices.icol_point, prescribed_conditions)

    # rotation matrix derivatives
    # TODO

    # point velocity rates
    Vdot_Vdot, Ωdot_Ωdot = point_displacement_jacobians(x, ipoint, indices.icol_point, prescribed_conditions)

    # linear and angular momentum rates
    Pdot_Vdot = Ct*mass11*C
    Pdot_Ωdot = Ct*mass12*C

    Hdot_Vdot = Ct*mass21*C
    Hdot_Ωdot = Ct*mass22*C

    return (; C, Ct, Q, mass11, mass12, mass21, mass22, u, θ, udot_udot, θdot_θdot, 
        Vdot_Vdot, Ωdot_Ωdot, Pdot_Vdot, Pdot_Ωdot, Hdot_Vdot, Hdot_Ωdot)
end

"""
    static_point_resultant_jacobians(properties)

Calculate the jacobians for the net loads `F` and `M` applied at a point for a static analysis.
"""
@inline function static_point_resultant_jacobians(properties)

    jacobians = _static_point_resultant_jacobians(properties)

    return finalize_point_resultant_jacobians(properties, jacobians)
end

@inline function _static_point_resultant_jacobians(properties)

    @unpack Ct, mass11, mass12, mass21, mass22, a, α, Ct_θ1, Ct_θ2, Ct_θ3 = properties

    @unpack F_θ, M_θ, F_F, M_M = properties

    # externally applied loads
    Fe_θ, Me_θ = F_θ, M_θ

    # loads due to acceleration (including gravity)
    Fa_u, Fa_θ, Ma_u, Ma_θ = acceleration_load_jacobians(Ct, mass11, mass12, mass21, mass22, a, α, Ct_θ1, Ct_θ2, Ct_θ3)

    # load resultants
    F_u = Fa_u
    F_θ = Fe_θ + Fa_θ
    M_u = Ma_u
    M_θ = Me_θ + Ma_θ

    return (; F_u, F_θ, F_F, M_u, M_θ, M_M)
end

@inline function finalize_point_resultant_jacobians(properties, jacobians)

    @unpack u_u, θ_θ = properties

    @unpack F_u, F_θ, M_u, M_θ = jacobians

    F_u = F_u*u_u
    F_θ = F_θ*θ_θ

    M_u = M_u*u_u
    M_θ = M_θ*θ_θ

    return (; jacobians..., F_u, F_θ, M_u, M_θ)
end

"""
    steady_state_point_resultant_jacobians(properties)

Calculate the jacobians for the net loads `F` and `M` applied at a point for a steady 
state analysis.
"""
@inline function steady_state_point_resultant_jacobians(properties)

    jacobians = _steady_state_point_resultant_jacobians(properties)

    return finalize_point_resultant_jacobians(properties, jacobians)
end

@inline function _steady_state_point_resultant_jacobians(properties)

    jacobians = _static_point_resultant_jacobians(properties)

    @unpack V, Ω, P, H, ω, P_θ, P_V, P_Ω, H_θ, H_V, H_Ω = properties

    @unpack F_θ, M_θ = jacobians

    F_θ += tilde(ω)*P_θ
    F_V = tilde(ω)*P_V
    F_Ω = tilde(ω)*P_Ω

    M_θ += tilde(ω)*H_θ + tilde(V)*P_θ
    M_V = tilde(ω)*H_V + tilde(V)*P_V - tilde(P)
    M_Ω = tilde(ω)*H_Ω + tilde(V)*P_Ω

    return (; jacobians..., F_θ, F_V, F_Ω, M_θ, M_V, M_Ω)
end

"""
    initial_condition_point_resultant_jacobians(properties)

Calculate the jacobians for the net loads `F` and `M` applied at a point for the initialization
of a time domain analysis.
"""
@inline function initial_condition_point_resultant_jacobians(properties)

    @unpack V, Ω, P, H, v, ω, F_F, M_M, P_V, P_Ω, H_V, H_Ω, Pdot_V, Pdot_Ω, Hdot_V, Hdot_Ω, 
        Pdot_Vdot, Pdot_Ωdot, Hdot_Vdot, Hdot_Ωdot = properties

    F_Vdot = Pdot_Vdot
    F_Ωdot = Pdot_Ωdot
    F_V = tilde(ω)*P_V + Pdot_V
    F_Ω = tilde(ω)*P_Ω + Pdot_Ω

    M_Vdot = Hdot_Vdot
    M_Ωdot = Hdot_Ωdot
    M_V = tilde(ω)*H_V + Hdot_V + tilde(V)*P_V - tilde(P)
    M_Ω = tilde(ω)*H_Ω + Hdot_Ω + tilde(V)*P_Ω

    return (; F_Vdot, F_Ωdot, F_F, F_V, F_Ω, M_Vdot, M_Ωdot, M_M, M_V, M_Ω)
end

"""
    dynamic_point_resultant_jacobians(properties)

Calculate the jacobians for the net loads `F` and `M` applied at a point for a dynamic analysis.
"""
@inline function dynamic_point_resultant_jacobians(properties)

    jacobians = _dynamic_point_resultant_jacobians(properties)

    return finalize_point_resultant_jacobians(properties, jacobians)
end

@inline function _dynamic_point_resultant_jacobians(properties)

    jacobians = _steady_state_point_resultant_jacobians(properties)

    @unpack Pdot_θ, Pdot_V, Pdot_Ω, Hdot_θ, Hdot_V, Hdot_Ω = properties
    @unpack F_θ, F_V, F_Ω, M_θ, M_V, M_Ω = jacobians

    F_θ += Pdot_θ
    F_V += Pdot_V
    F_Ω += Pdot_Ω

    M_θ += Hdot_θ
    M_V += Hdot_V
    M_Ω += Hdot_Ω

    return (; jacobians..., F_θ, F_V, F_Ω, M_θ, M_V, M_Ω)
end

"""
    mass_matrix_point_resultant_jacobians(properties)

Calculate the mass matrix jacobians for the net loads `F` and `M` applied at a point
"""
@inline function mass_matrix_point_resultant_jacobians(properties)

    @unpack Pdot_Vdot, Pdot_Ωdot, Hdot_Vdot, Hdot_Ωdot = properties

    F_Vdot = Pdot_Vdot
    F_Ωdot = Pdot_Ωdot
    M_Vdot = Hdot_Vdot
    M_Ωdot = Hdot_Ωdot

    @unpack Vdot_Vdot, Ωdot_Ωdot = properties

    F_Vdot *= Vdot_Vdot
    F_Ωdot *= Ωdot_Ωdot
    M_Vdot *= Vdot_Vdot
    M_Ωdot *= Ωdot_Ωdot

    return (; F_Vdot, F_Ωdot, M_Vdot, M_Ωdot) 
end

"""
    steady_state_velocity_jacobians(properties)

Calculate the jacobians of the velocity residuals `rV` and `rΩ` of a point for a steady 
state analysis.
"""
@inline function steady_state_velocity_jacobians(properties)

    @unpack ω, u_u = properties
    
    rV_u = -tilde(ω)*u_u
    rV_V = I3

    rΩ_Ω = I3

    return (; rV_u, rV_V, rΩ_Ω)
end

"""
    initial_condition_velocity_jacobians(properties)

Calculate the jacobians of the velocity residuals `rV` and `rΩ` of a point for the
initialization of a time domain analysis.
"""
@inline function initial_condition_velocity_jacobians(properties)
   
    rV_V = I3

    rΩ_Ω = I3

    return (; rV_V, rΩ_Ω)
end

"""
    newmark_velocity_jacobians(properties)

Calculate the jacobians of the velocity residuals `rV` and `rΩ` of a point for a Newmark 
scheme time-marching analysis.
"""
@inline function newmark_velocity_jacobians(properties)

    jacobians = steady_state_velocity_jacobians(properties)

    @unpack Ct, Q, θdot, Ct_θ1, Ct_θ2, Ct_θ3, Q_θ1, Q_θ2, Q_θ3, u_u, θ_θ, udot_u, θdot_θ = properties

    @unpack rV_u, rV_V, rΩ_Ω = jacobians

    rV_u -= udot_u*u_u

    rΩ_θ = (-mul3(Ct_θ1, Ct_θ2, Ct_θ3, Q*θdot) - Ct*mul3(Q_θ1, Q_θ2, Q_θ3, θdot) - Ct*Q*θdot_θ)*θ_θ

    return (; jacobians..., rV_u, rV_V, rΩ_θ, rΩ_Ω)
end

"""
    dynamic_velocity_jacobians(properties)

Calculate the jacobians of the velocity residuals `rV` and `rΩ` of a point for a dynamic 
analysis.
"""
@inline function dynamic_velocity_jacobians(properties)

    jacobians = steady_state_velocity_jacobians(properties)

    @unpack Ct, Q, θdot, Ct_θ1, Ct_θ2, Ct_θ3, Q_θ1, Q_θ2, Q_θ3 = properties

    @unpack rV_u, rV_V, rΩ_Ω = jacobians

    rΩ_θ = -mul3(Ct_θ1, Ct_θ2, Ct_θ3, Q*θdot) - Ct*mul3(Q_θ1, Q_θ2, Q_θ3, θdot)

    return (; jacobians..., rV_u, rV_V, rΩ_θ, rΩ_Ω)
end

"""
    mass_matrix_velocity_jacobians(properties)

Calculate the mass matrix jacobians of the velocity residuals `rV` and `rΩ` of a point
"""
@inline function mass_matrix_velocity_jacobians(properties)

    @unpack Ct, Q, udot_udot, θdot_θdot = properties

    rV_udot = -udot_udot
    rΩ_θdot = -Ct*Q*θdot_θdot

    return (; rV_udot, rΩ_θdot)
end

"""
    insert_static_point_jacobians!(jacob, indices, force_scaling, ipoint, resultants)

Insert the jacobian entries corresponding to a point for a static analysis 
into the system jacobian matrix.
"""
@inline function insert_static_point_jacobians!(jacob, indices, force_scaling, ipoint, resultants)

    @unpack F_u, F_θ, F_F, M_u, M_θ, M_M = resultants

    irow = indices.irow_point[ipoint]
    icol = indices.icol_point[ipoint]

    jacob[irow:irow+2, icol:icol+2] .= -F_F - F_u ./ force_scaling
    jacob[irow:irow+2, icol+3:icol+5] .= -F_θ ./ force_scaling
    jacob[irow+3:irow+5, icol:icol+2] .= -M_u ./ force_scaling 
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

    @unpack F_Vdot, F_Ωdot, F_F, F_V, F_Ω, M_Vdot, M_Ωdot, M_M, M_V, M_Ω = resultants
    @unpack rV_V, rΩ_Ω = velocities

    irow = indices.irow_point[ipoint]
    icol = indices.icol_point[ipoint]

    jacob[irow:irow+2, icol:icol+2] .= -F_F .- F_Vdot ./ force_scaling
    jacob[irow:irow+2, icol+3:icol+5] .= -F_Ωdot ./ force_scaling
    jacob[irow:irow+2, icol+6:icol+8] .= -F_V ./ force_scaling
    jacob[irow:irow+2, icol+9:icol+11] .= -F_Ω ./ force_scaling

    jacob[irow+3:irow+5, icol:icol+2] .= -M_Vdot ./ force_scaling 
    jacob[irow+3:irow+5, icol+3:icol+5] .= -M_M .- M_Ωdot ./ force_scaling
    jacob[irow+3:irow+5, icol+6:icol+8] .= -M_V ./ force_scaling
    jacob[irow+3:irow+5, icol+9:icol+11] .= -M_Ω ./ force_scaling

    jacob[irow+6:irow+8, icol+6:icol+8] .= rV_V

    jacob[irow+9:irow+11, icol+9:icol+11] .= rΩ_Ω

    return jacob
end

"""
    insert_steady_state_point_jacobians!(jacob, indices, force_scaling, ipoint, resultants, 
        velocities)

Insert the jacobian entries corresponding to a point for a steady state analysis into the 
system jacobian matrix.
"""
@inline function insert_steady_state_point_jacobians!(jacob, indices, force_scaling, ipoint,  
    resultants, velocities)

    insert_static_point_jacobians!(jacob, indices, force_scaling, ipoint, resultants)

    @unpack F_V, F_Ω, M_V, M_Ω = resultants
    @unpack rV_u, rV_V, rΩ_Ω = velocities

    irow = indices.irow_point[ipoint]
    icol = indices.icol_point[ipoint]

    jacob[irow:irow+2, icol+6:icol+8] .= -F_V ./ force_scaling
    jacob[irow:irow+2, icol+9:icol+11] .= -F_Ω ./ force_scaling

    jacob[irow+3:irow+5, icol+6:icol+8] .= -M_V ./ force_scaling
    jacob[irow+3:irow+5, icol+9:icol+11] .= -M_Ω ./ force_scaling

    jacob[irow+6:irow+8, icol:icol+2] .= rV_u
    jacob[irow+6:irow+8, icol+6:icol+8] .= rV_V

    jacob[irow+9:irow+11, icol+9:icol+11] .= rΩ_Ω

    return jacob
end

"""
    insert_newmark_point_jacobians!(jacob, indices, force_scaling, ipoint, resultants, 
        velocities)

Insert the jacobian entries corresponding to a point for a Newmark scheme time marching 
analysis into the system jacobian matrix.
"""
@inline function insert_newmark_point_jacobians!(jacob, indices, force_scaling, ipoint,  
    resultants, velocities)

    insert_steady_state_point_jacobians!(jacob, indices, force_scaling, ipoint, resultants, velocities)

    @unpack rΩ_θ = velocities

    irow = indices.irow_point[ipoint]
    icol = indices.icol_point[ipoint]

    jacob[irow+9:irow+11, icol+3:icol+5] .= rΩ_θ

    return jacob
end

"""
    insert_dynamic_point_jacobians!(jacob, indices, force_scaling, ipoint, resultants, 
        velocities)

Insert the jacobian entries corresponding to a point for a dynamic analysis into the 
system jacobian matrix.
"""
@inline function insert_dynamic_point_jacobians!(jacob, indices, force_scaling, ipoint,  
    resultants, velocities)

    insert_newmark_point_jacobians!(jacob, indices, force_scaling, ipoint,  
        resultants, velocities)

    return jacob
end

"""
    insert_point_mass_matrix_jacobians!(jacob, gamma, indices, force_scaling, ipoint, 
        resultants, velocities)

Insert the mass matrix jacobian entries corresponding to a point into the system jacobian 
matrix.
"""
@inline function insert_point_mass_matrix_jacobians!(jacob, gamma, indices, force_scaling, ipoint,  
    resultants, velocities)

    @unpack F_Vdot, F_Ωdot, M_Vdot, M_Ωdot = resultants
    @unpack rV_udot, rΩ_θdot = velocities

    irow = indices.irow_point[ipoint]
    icol = indices.icol_point[ipoint]

    jacob[irow:irow+2, icol+6:icol+8] .-= F_Vdot .* gamma ./ force_scaling
    jacob[irow:irow+2, icol+9:icol+11] .-= F_Ωdot .* gamma ./ force_scaling

    jacob[irow+3:irow+5, icol+6:icol+8] .-= M_Vdot .* gamma ./ force_scaling
    jacob[irow+3:irow+5, icol+9:icol+11] .-= M_Ωdot .* gamma ./ force_scaling

    jacob[irow+6:irow+8, icol:icol+2] .+= rV_udot .* gamma
    jacob[irow+9:irow+11, icol+3:icol+5] .+= rΩ_θdot .* gamma

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
        force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity)

    resultants = steady_state_point_resultant_jacobians(properties)

    velocities = steady_state_velocity_jacobians(properties)

    insert_steady_state_point_jacobians!(jacob, indices, force_scaling, ipoint, resultants, velocities)

    return jacob
end

"""
    initial_condition_point_jacobian!(jacob, x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)

Calculate and insert the jacobian entries corresponding to a point for the initialization
of a time domain analysis into the system jacobian matrix.
"""
@inline function initial_condition_point_jacobian!(jacob, x, indices, force_scaling, assembly, ipoint,  
    prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0,
    u0, θ0, udot0, θdot0)

    properties = initial_condition_point_properties(x, indices, force_scaling, assembly, ipoint,  
        prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0, u0, θ0, udot0, θdot0)

    properties = initial_condition_point_jacobian_properties(properties, x, indices, force_scaling, 
        assembly, ipoint, prescribed_conditions, point_masses, gravity)

    resultants = initial_condition_point_resultant_jacobians(properties)

    velocities = initial_condition_velocity_jacobians(properties)

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
        assembly, ipoint, prescribed_conditions, point_masses, gravity,
        udot_init, θdot_init, Vdot_init, Ωdot_init, dt)

    resultants = dynamic_point_resultant_jacobians(properties)

    velocities = newmark_velocity_jacobians(properties)

    insert_newmark_point_jacobians!(jacob, indices, force_scaling, ipoint, resultants, velocities)

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
        assembly, ipoint, prescribed_conditions, point_masses, gravity)

    resultants = dynamic_point_resultant_jacobians(properties)

    velocities = dynamic_velocity_jacobians(properties)

    insert_dynamic_point_jacobians!(jacob, indices, force_scaling, ipoint, resultants, velocities)

    return jacob
end

"""
    point_mass_matrix!(jacob, gamma, x, indices, force_scaling, assembly, 
        ipoint, prescribed_conditions)

Calculate and insert the mass_matrix jacobian entries corresponding to a point into the 
system jacobian matrix.
"""
@inline function point_mass_matrix!(jacob, gamma, x, indices, force_scaling, assembly, 
    ipoint, prescribed_conditions, point_masses)

    properties = mass_matrix_point_jacobian_properties(x, indices, force_scaling, 
        assembly, ipoint, prescribed_conditions, point_masses)

    resultants = mass_matrix_point_resultant_jacobians(properties)

    velocities = mass_matrix_velocity_jacobians(properties)
    
    insert_point_mass_matrix_jacobians!(jacob, gamma, indices, force_scaling, ipoint,  
        resultants, velocities)

    return jacob
end