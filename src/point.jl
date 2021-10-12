"""
    point_variables(x, icol)
    point_variables(x, icol, prescribed_conditions)

Extract u, θ, F, M for the point described by the point state variables at `icol` in
x after incorporating the prescribed conditions in `prescribed_conditions`

Note that the degrees of freedom that are not specified in `prescribed_conditions`
are used as state variables (e.g. prescribing F[2] would mean u[2] = x[icol+1])
"""
@inline function point_variables(x, icol)
    u = SVector(x[icol], x[icol+1], x[icol+2])
    θ = SVector(x[icol+3], x[icol+4], x[icol+5])
    F = zero(u)
    M = zero(θ)

    return u, θ, F, M
end

@inline function point_variables(x, icol, prescribed_conditions, force_scaling)

    value = prescribed_conditions.value
    follower = prescribed_conditions.follower

    # get degrees of freedom to fix/set
    force_dof = prescribed_conditions.force

    # get the displacement and rotations of the point
    u = SVector(ifelse(force_dof[1], x[icol  ], value[1]),
                ifelse(force_dof[2], x[icol+1], value[2]),
                ifelse(force_dof[3], x[icol+2], value[3]))
    θ = SVector(ifelse(force_dof[4], x[icol+3], value[4]),
                ifelse(force_dof[5], x[icol+4], value[5]),
                ifelse(force_dof[6], x[icol+5], value[6]))

    # get the transformation matrix for the point
    C = get_C(θ)

    # solve for the force applied at the point due to the prescribed loads
    Fp = zero(u)
    for i = 1:3
        if force_dof[i]
            # add dead force_dof
            Fp += SVector(I3[i,1], I3[i,2], I3[i,3])*value[i]
            # add follower force_dof
            Fp += SVector(C[i,1], C[i,2], C[i,3])*follower[i]
        end
    end

    # solve for the moment applied at the point due to the prescribed loads
    Mp = zero(θ)
    for i = 4:6
        if force_dof[i]
            # add dead moment
            Mp += SVector(I3[i-3,1], I3[i-3,2], I3[i-3,3])*value[i]
            # add follower moment
            Mp += SVector(C[i-3,1], C[i-3,2], C[i-3,3])*follower[i]
        end
    end

    # overwrite external forces/moments with solved for forces/moments
    F = SVector(ifelse(force_dof[1], Fp[1], x[icol  ] * force_scaling),
                ifelse(force_dof[2], Fp[2], x[icol+1] * force_scaling),
                ifelse(force_dof[3], Fp[3], x[icol+2] * force_scaling))
    M = SVector(ifelse(force_dof[4], Mp[1], x[icol+3] * force_scaling),
                ifelse(force_dof[5], Mp[2], x[icol+4] * force_scaling),
                ifelse(force_dof[6], Mp[3], x[icol+5] * force_scaling))

    return u, θ, F, M
end

"""
    insert_point_residual!(resid, irow_elem, irow_point, u, θ, F, M, side)

Modify the equilibrium and constitutive equations to account for the point
variables given by u, θ, F, M

If `irow_elem != irow_point`, assume that the equilibrium equations have already
been modified

# Arguments
 - resid: Residual vector
 - irow_elem: Row index of the first equilibrium/compatability equation for one
       side of the beam element
 - irow_point: Row index of the first equilibrium equation for the point
 - u: Displacement of the point
 - θ: Rotation of the point
 - F: External forces imposed on the point
 - M: External moments imposed on the point
 - side: Side of beam element (-1 (left) or 1 (right))
"""
@inline function insert_point_residual!(resid, irow_elem, irow_point, u, θ, F, M, side, force_scaling)

    if irow_elem == irow_point
        # add to equilibrium and compatability equations
        for i = 1:3
            resid[irow_elem+i-1] -= F[i] ./ force_scaling
            resid[irow_elem+i+2] -= M[i] ./ force_scaling
            resid[irow_elem+i+5] += side*u[i]
            resid[irow_elem+i+8] += side*θ[i]
        end
    else
        # add to compatability equations
        for i = 1:3
            resid[irow_elem+i-1] += side*u[i]
            resid[irow_elem+i+2] += side*θ[i]
        end
    end

    return resid
end

"""
    point_residual!(resid, x, ipoint, assembly, prescribed_conditions,
        force_scaling, icol, irow_point, irow_elem1, irow_elem2)

Adds a points contributions to the residual vector

# Arguments
 - `resid`: residual vector
 - `x`: current state vector
 - `ipoint`: index of point
 - `assembly`: assembly of interconnected beam elements
 - `prescribed_conditions`: dictionary of prescribed conditions
 - `force_scaling`: scaling parameter for forces
 - `icol`: starting index for the point's state variables
 - `irow_point`: Row index of the first equilibrium equation for the point
 - `irow_elem1`: Row index of first equation for the left side of each beam element
 - `irow_elem2`: Row index of first equation for the right side of each beam element
"""
@inline function point_residual!(resid, x, ipoint, assembly, prescribed_conditions, 
    force_scaling, icol, irow_point, irow_elem1, irow_elem2)

    nelem = length(assembly.elements)

    # incorporate prescribed condition if applicable
    has_prescribed = haskey(prescribed_conditions, ipoint)
    if has_prescribed
        u, θ, F, M = point_variables(x, icol, prescribed_conditions[ipoint], force_scaling)
    else
        u, θ, F, M = point_variables(x, icol)
    end

    # search for beam elements that are connected to the specified point
    for ielem = 1:nelem
        # check left side of beam element
        if ipoint == assembly.start[ielem]
            # add to residual equations for the beam element endpoint
            side = -1
            irow_elem = irow_elem1[ielem]
            insert_point_residual!(resid, irow_elem, irow_point, u, θ, F, M, side, force_scaling)
        end
        # check right side of beam element
        if ipoint == assembly.stop[ielem]
            # add to residual equations for the beam element endpoint
            side = 1
            irow_elem = irow_elem2[ielem]
            insert_point_residual!(resid, irow_elem, irow_point, u, θ, F, M, side, force_scaling)
        end
    end

    return resid
end

"""
    point_follower_jacobians(x, icol, prescribed_conditions)

Calculate the jacobians of the follower forces/moments with respect to θ

# Arguments
 - x: Current state variable vector
 - icol: Row/Column index of the first state variable for the point
 - prescribed_conditions: Prescribed conditions for the point
"""
@inline function point_follower_jacobians(x, icol, prescribed_conditions)

    TF = eltype(x)

    # Apply time functions
    value = prescribed_conditions.value
    follower = prescribed_conditions.follower

    force_dof = prescribed_conditions.force

    θ = SVector(ifelse(force_dof[4], x[icol+3], value[4]),
                ifelse(force_dof[5], x[icol+4], value[5]),
                ifelse(force_dof[6], x[icol+5], value[6]))

    C = get_C(θ)
    C_θ1, C_θ2, C_θ3 = get_C_θ(C, θ)

    # solve for the jacobian wrt theta of the follower forces
    Fp_θ = @SMatrix zeros(TF, 3, 3)
    for i = 1:3
        if force_dof[i]
            rot_θ = @SMatrix [
                C_θ1[i,1] C_θ2[i,1] C_θ3[i,1];
                C_θ1[i,2] C_θ2[i,2] C_θ3[i,2];
                C_θ1[i,3] C_θ2[i,3] C_θ3[i,3]
                ]
            Fp_θ += rot_θ*follower[i]
        end
    end

    # solve for the jacobian wrt theta of the follower moments
    Mp_θ = @SMatrix zeros(TF, 3, 3)
    for i = 1:3
        if force_dof[i+3]
            rot_θ = @SMatrix [
                C_θ1[i,1] C_θ2[i,1] C_θ3[i,1];
                C_θ1[i,2] C_θ2[i,2] C_θ3[i,2];
                C_θ1[i,3] C_θ2[i,3] C_θ3[i,3]
                ]
            Mp_θ += rot_θ*follower[i+3]
        end
    end

    # if displacement is specified, corresponding component of follower jacobian is zero
    F1_θ = ifelse(force_dof[1], SVector(Fp_θ[1,1], Fp_θ[1,2], Fp_θ[1,3]), @SVector zeros(TF, 3))
    F2_θ = ifelse(force_dof[2], SVector(Fp_θ[2,1], Fp_θ[2,2], Fp_θ[2,3]), @SVector zeros(TF, 3))
    F3_θ = ifelse(force_dof[3], SVector(Fp_θ[3,1], Fp_θ[3,2], Fp_θ[3,3]), @SVector zeros(TF, 3))
    F_θ = vcat(F1_θ', F2_θ', F3_θ')

    M1_θ = ifelse(force_dof[4], SVector(Mp_θ[1,1], Mp_θ[1,2], Mp_θ[1,3]), @SVector zeros(TF, 3))
    M2_θ = ifelse(force_dof[5], SVector(Mp_θ[2,1], Mp_θ[2,2], Mp_θ[2,3]), @SVector zeros(TF, 3))
    M3_θ = ifelse(force_dof[6], SVector(Mp_θ[3,1], Mp_θ[3,2], Mp_θ[3,3]), @SVector zeros(TF, 3))
    M_θ = vcat(M1_θ', M2_θ', M3_θ')

    return F_θ, M_θ
end

"""
    insert_point_jacobian!(jacob, irow_elem, irow_point, icol,
        prescribed_conditions, side, F_θ, M_θ)

Modify the jacobian entries for the equilibrium and constitutive equations to
account for the point variables at icol

If irow_elem != irow_point, assume that the equilibrium equations have already been modified

# Arguments
 - jacob: Jacobian of residual vector with respect to state vectors
 - irow_elem: Row index of the first equilibrium/compatability equation for one side of the beam element
 - irow_point: Row index of the first equilibrium equation for the point
 - icol: Row/Column index of the first state variable for the point
 - prescribed_conditions: Prescribed force/displacement and moment/rotation on point
 - side: Side of beam element (-1 (left) or 1 (right))
"""
@inline function insert_point_jacobian!(jacob, irow_elem, irow_point, icol, side, prescribed_force, F_θ, M_θ, force_scaling)

    if irow_elem == irow_point
        # add jacobian entries for equilibrium and compatability equations
        for i = 1:3 # forces/displacements
            if prescribed_force[i]
                # F[i] is prescribed, u[i] is a state variable
                for j = 4:6
                    # check if θ[j-3] is a state variable
                    if prescribed_force[j]
                        # add equilibrium equation jacobian component
                        jacob[irow_elem+i-1, icol+j-1] = -F_θ[i,j-3] ./ force_scaling
                    end
                end
                # add compatability equation jacobian component
                jacob[irow_elem+i+5, icol+i-1] = side # u_u = I
            else
                # u[i] is prescribed, F[i] is a state variable

                # add equilibrium equation jacobian component
                jacob[irow_elem+i-1, icol+i-1] = -1 # F_F = I
            end
        end
        for i = 4:6 # moments/rotations
            if prescribed_force[i]
                # M[i-3] is prescribed, θ[i-3] is a state variable
                for j = 4:6
                    # check if θ[j-3] is a state variable
                    if prescribed_force[j]
                        # add equilibrium equation jacobian component
                        jacob[irow_elem+i-1, icol+j-1] = -M_θ[i-3,j-3] ./ force_scaling
                    end
                end
                # add compatability equation jacobian component
                jacob[irow_elem+i+5, icol+i-1] = side # θ_θ = I
            else
                # θ[i-3] is prescribed, M[i-3] is a state variable

                # add equilibrium equation jacobian component
                jacob[irow_elem+i-1, icol+i-1] = -1 # M_M = I
            end
        end
    else
        # add jacobian entries for compatability equations only
        for i = 1:3
            if prescribed_force[i]
                # F[i] is prescribed, u[i] is a state variable
                jacob[irow_elem+i-1, icol+i-1] = side # u_u = I
            end
        end
        for i = 4:6
            if prescribed_force[i]
                # M[i-3] is prescribed, θ[i-3] is a state variable
                jacob[irow_elem+i-1, icol+i-1] = side # θ_θ = I
            end
        end
    end

    return jacob
end


"""
    point_jacobian!(jacob, x, ipoint, assembly, prescribed_conditions, icol,
        irow_point, irow_elem1, irow_elem2)

Adds a points contributions to the residual vector

# Arguments
 - `jacob`: residual vector
 - `x`: current state vector
 - `ipoint`: index of point
 - `assembly`: assembly of interconnected beam elements
 - `prescribed_conditions`: dictionary of prescribed conditions
 - `force_scaling`: scaling parameter for forces
 - `icol`: starting index for the point's state variables
 - `irow_point`: Row index of the first equilibrium equation for the point
 - `irow_elem1`: Row index of first equation for the left side of each beam element
 - `irow_elem2`: Row index of first equation for the right side of each beam element
"""
@inline function point_jacobian!(jacob, x, ipoint, assembly, prescribed_conditions, 
    force_scaling, icol, irow_point, irow_elem1, irow_elem2)

    nelem = length(assembly.elements)

    # incorporate prescribed condition if applicable
    prescribed = haskey(prescribed_conditions, ipoint)
    if prescribed
        F_θ, M_θ = point_follower_jacobians(x, icol, prescribed_conditions[ipoint])
        force_dof = prescribed_conditions[ipoint].force
    else
        F_θ = @SMatrix zeros(3,3)
        M_θ = @SMatrix zeros(3,3)
        force_dof = @SVector ones(Bool, 6)
    end

    # search for beam elements that are connected to the specified point
    for ielem = 1:nelem
        # check left side of beam element
        if ipoint == assembly.start[ielem]
            # add to residual equations for the beam element endpoint
            side = -1
            irow_elem = irow_elem1[ielem]
            jacob = insert_point_jacobian!(jacob, irow_elem, irow_point, icol, side, force_dof, F_θ, M_θ, force_scaling)
        end
        # check right side of beam element
        if ipoint == assembly.stop[ielem]
            # add to residual equations for the beam element endpoint
            side = 1
            irow_elem = irow_elem2[ielem]
            jacob = insert_point_jacobian!(jacob, irow_elem, irow_point, icol, side, force_dof, F_θ, M_θ, force_scaling)
        end
    end

    return jacob
end
