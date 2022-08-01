# --- Helper Functions --- #

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

    F_F1 = ifelse(pl[1], SVector{3,eltype(x)}(0,0,0), SVector{3,eltype(x)}(1,0,0))
    F_F2 = ifelse(pl[2], SVector{3,eltype(x)}(0,0,0), SVector{3,eltype(x)}(0,1,0))
    F_F3 = ifelse(pl[3], SVector{3,eltype(x)}(0,0,0), SVector{3,eltype(x)}(0,0,1))
    F_F = hcat(F_F1, F_F2, F_F3)

    M_M1 = ifelse(pl[4], SVector{3,eltype(x)}(0,0,0), SVector{3,eltype(x)}(1,0,0))
    M_M2 = ifelse(pl[5], SVector{3,eltype(x)}(0,0,0), SVector{3,eltype(x)}(0,1,0))
    M_M3 = ifelse(pl[6], SVector{3,eltype(x)}(0,0,0), SVector{3,eltype(x)}(0,0,1))
    M_M = hcat(M_M1, M_M2, M_M3)

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

# --- Point Properties --- #

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

# --- Point Resultants --- #

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

# --- Residual Placement --- #

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

# --- Point Residual --- #

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