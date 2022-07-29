# --- Helper Functions --- #

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

# --- Point Properties --- #

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

# --- Point Resultants --- #

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

# --- Point Residual --- #

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