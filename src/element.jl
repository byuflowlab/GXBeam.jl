"""
    Element{TF}

Composite type that defines a beam element's properties

# Fields
 - `L`: Length of the beam element
 - `x`: Location of the beam element (the center of the beam element)
 - `compliance`: Compliance matrix for the beam element
 - `mass`: Mass matrix for the beam element
 - `Cab`: Transformation matrix from the undeformed beam element frame to the body frame
 - `mu`: Beam element damping coefficients
"""
struct Element{TF}
    L::TF
    x::SVector{3,TF}
    compliance::SMatrix{6,6,TF,36}
    mass::SMatrix{6,6,TF,36}
    Cab::SMatrix{3,3,TF,9}
    mu::SVector{6,TF}
end

"""
    Element(L, x, compliance, mass, Cab, mu)

Construct a beam element

# Arguments
- `L`: Length of the beam element
- `x`: Location of the beam element (the center of the beam element)
- `compliance`: Beam element compliance matrix
- `mass`: Beam element mass matrix
- `Cab`: Transformation matrix from the undeformed beam element frame to the body frame
- `mu`: Beam element damping coefficients
"""
function Element(L, x, compliance, mass, Cab, mu)
    TF = promote_type(typeof(L), eltype(x), eltype(compliance), eltype(mass), eltype(Cab), eltype(mu))
    return Element{TF}(L, x, compliance, mass, Cab, mu)
end

"""
    element_strain(element, F, M)

Calculate the strain of a beam element given the resultant forces and moments applied on
the element expressed in the deformed beam element frame
"""
@inline function element_strain(element, F, M)
    C = element.compliance
    S11 = C[SVector{3}(1:3), SVector{3}(1:3)]
    S12 = C[SVector{3}(1:3), SVector{3}(4:6)]
    return S11*F + S12*M
end

"""
    element_curvature(element, F, M)

Calculate the curvature of a beam element given the resultant force and moments applied on
the element expressed in the deformed beam element frame
"""
@inline function element_curvature(element, F, M)
    C = element.compliance
    S21 = C[SVector{3}(4:6), SVector{3}(1:3)]
    S22 = C[SVector{3}(4:6), SVector{3}(4:6)]
    return S21*F + S22*M
end

"""
    element_linear_momentum(element, V, Ω)

Calculate the linear momentum of a beam element given the linear and angular velocity of
the element expressed in the deformed beam element frame
"""
@inline function element_linear_momentum(element, V, Ω)
    M = element.mass
    mass11 = M[SVector{3}(1:3), SVector{3}(1:3)]
    mass12 = M[SVector{3}(1:3), SVector{3}(4:6)]
    return mass11*V + mass12*Ω
end

"""
    element_angular_momentum(element, V, Ω)

Calculate the angular momentum of a beam element given the linear and angular velocity of
the element expressed in the deformed beam element frame
"""
@inline function element_angular_momentum(element, V, Ω)
    M = element.mass
    mass21 = M[SVector{3}(4:6), SVector{3}(1:3)]
    mass22 = M[SVector{3}(4:6), SVector{3}(4:6)]
    return mass21*V + mass22*Ω
end

"""
    static_element_state_variables(x, icol, force_scaling)

Return element state variables for a static analysis
"""
@inline function static_element_state_variables(x, icol, force_scaling)
    u = SVector(x[icol   ], x[icol+1 ], x[icol+2 ])
    θ = SVector(x[icol+3 ], x[icol+4 ], x[icol+5 ])
    F = SVector(x[icol+6 ], x[icol+7 ], x[icol+8 ]) .* force_scaling
    M = SVector(x[icol+9 ], x[icol+10], x[icol+11]) .* force_scaling
    return u, θ, F, M
end

"""
    dynamic_element_state_variables(x, icol, force_scaling)

Return element state variables for a dynamic analysis
"""
@inline function dynamic_element_state_variables(x, icol, force_scaling)
    u = SVector(x[icol   ], x[icol+1 ], x[icol+2 ]) # or Vdot
    θ = SVector(x[icol+3 ], x[icol+4 ], x[icol+5 ]) # or Ωdot
    F = SVector(x[icol+6 ], x[icol+7 ], x[icol+8 ]) .* force_scaling
    M = SVector(x[icol+9 ], x[icol+10], x[icol+11]) .* force_scaling
    V = SVector(x[icol+12], x[icol+13], x[icol+14])
    Ω = SVector(x[icol+15], x[icol+16], x[icol+17])
    return u, θ, F, M, V, Ω
end

"""
    static_element_equations(ΔL, Cab, CtCab, u, θ, F, M, γ, κ, f1, f2, m1, m2)

Calculate the element resultants for a static analysis.

# Arguments:
 - `ΔL`: length
 - `Cab`: transformation matrix from the undeformed beam element frame to the body frame
 - `CtCab`: transformation matrix from the deformed beam element frame to the body frame
 - `u`: linear displacement
 - `θ`: angular displacement
 - `F`: internal force
 - `M`: internal moment
 - `γ`: engineering strain
 - `κ`: curvature
 - `f1`: integrated force applied to the start of the beam element
 - `f2`: integrated force applied to the end of the beam element
 - `m1`: integrated moment applied to the start of the beam element
 - `m2`: integrated moment applied to the end of the beam element
"""
@inline function static_element_equations(ΔL, Cab, CtCab, u, θ, F, M, γ, κ, f1, f2, m1, m2)

    # note that a ΔL term has been incorporated into the stiffness and mass matrix 

    tmp = CtCab*F
    f_u1 = -tmp - f1
    f_u2 =  tmp - f2

    tmp1 = CtCab*M
    tmp2 = 1/2*CtCab*cross(ΔL*e1 + γ, F)
    f_ψ1 = -tmp1 - m1 - tmp2
    f_ψ2 =  tmp1 - m2 - tmp2

    tmp = 1/2*(CtCab*(ΔL*e1 + γ) - Cab*ΔL*e1)
    f_F1 =  u - tmp
    f_F2 = -u - tmp

    Qinv = get_Qinv(θ)
    tmp = 1/2*Qinv*Cab*κ
    f_M1 =  θ - tmp
    f_M2 = -θ - tmp

    return f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2
end

"""
    steady_state_element_equations(ΔL, Cab, CtCab, u, θ, F, M, γ, κ, V, Ω, P, H, v, ω, 
        f1, f2, m1, m2)

Calculate the element resultants for a steady state analysis.

# Arguments:
 - `ΔL`: length
 - `Cab`: transformation matrix from the undeformed beam element frame to the body frame
 - `CtCab`: transformation matrix from the deformed beam element frame to the body frame
 - `u`: linear displacement
 - `θ`: angular displacement
 - `F`: internal force
 - `M`: internal moment
 - `γ`: engineering strain
 - `κ`: curvature
 - `V`: linear velocity
 - `Ω`: angular velocity
 - `P`: linear momentum
 - `H`: angular momentum
 - `v`: linear velocity caused by body frame velocity
 - `ω`: angular velocity caused by body frame velocity
 - `f1`: integrated force applied to the start of the beam element
 - `f2`: integrated force applied to the end of the beam element
 - `m1`: integrated moment applied to the start of the beam element
 - `m2`: integrated moment applied to the end of the beam element
"""
@inline function steady_state_element_equations(ΔL, Cab, CtCab, u, θ, F, M, γ, κ, V, Ω, 
    P, H, v, ω, f1, f2, m1, m2)

    f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2 = static_element_equations(ΔL, Cab, 
        CtCab, u, θ, F, M, γ, κ, f1, f2, m1, m2)
    
    # note that a ΔL term has been incorporated into the stiffness and mass matrix 

    tmp = 1/2*cross(ω, CtCab*P) 
    f_u1 += tmp
    f_u2 += tmp

    tmp = 1/2*(cross(ω, CtCab*H) + CtCab*cross(V, P))
    f_ψ1 += tmp
    f_ψ2 += tmp

    f_V = CtCab*V - v - cross(ω, u)
    f_Ω = Ω - CtCab'*ω

    return f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2, f_V, f_Ω
end

"""
    dynamic_element_equations(ΔL, Cab, CtCab, u, θ, F, M, γ, κ, V, Ω, P, H, v, ω, 
        f1, f2, m1, m2, udot, θdot, CtCabPdot, CtCabHdot)

Calculate the element resultants for a dynamic analysis.

# Arguments:
 - `ΔL`: length
 - `Cab`: transformation matrix from the undeformed beam element frame to the body frame
 - `CtCab`: transformation matrix from the deformed beam element frame to the body frame
 - `u`: linear displacement
 - `θ`: angular displacement
 - `F`: internal force
 - `M`: internal moment
 - `γ`: engineering strain
 - `κ`: curvature
 - `V`: linear velocity
 - `Ω`: angular velocity
 - `P`: linear momentum
 - `H`: angular momentum
 - `v`: linear velocity from body frame velocity
 - `ω`: angular velocity from body frame velocity
 - `f1`: integrated force applied to the start of the beam element
 - `f2`: integrated force applied to the end of the beam element
 - `m1`: integrated moment applied to the start of the beam element
 - `m2`: integrated moment applied to the end of the beam element
 - `udot`: linear displacement rate, expressed in the body frame
 - `θdot`: angular displacement rate, expressed in the body frame
 - `CtCabPdot`: linear momentum rate, expressed in the body frame
 - `CtCabHdot`: angular momentum rate, expressed in the body frame
"""
@inline function dynamic_element_equations(ΔL, Cab, CtCab, u, θ, F, M, γ, κ, V, Ω, P, H, 
    v, ω, f1, f2, m1, m2, udot, θdot, CtCabPdot, CtCabHdot)

    f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2, f_V, f_Ω = 
        steady_state_element_equations(ΔL, Cab, CtCab, u, θ, F, M, γ, κ, V, Ω, P, H, v, ω, 
        f1, f2, m1, m2)

    # note that a ΔL term has been incorporated into the stiffness and mass matrix 

    tmp = 1/2*CtCabPdot
    f_u1 += tmp 
    f_u2 += tmp

    tmp = 1/2*CtCabHdot
    f_ψ1 += tmp
    f_ψ2 += tmp

    f_V -= udot

    Q = get_Q(θ)
    f_Ω -= Cab'*Q*θdot

    return f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2, f_V, f_Ω
end

"""
    static_insert_element_residual!(resid, force_scaling, irow_e1, 
        irow_p1, irow_e2, irow_p2, f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2)

Add beam element resultants to the residuals for a static analysis.  Initialize 
equilibrium and constitutive equation residuals if they are not yet initialized.

If `irow_e1 != irow_p1` and/or `irow_e2 != irow_p2`, assume the equilibrium equation 
residuals for the start and/or end of the beam element are already initialized

# Arguments
 - `resid`: system residual vector
 - `force_scaling`: scaling parameter for forces
 - `irow_e1`: row index of the first residual equation for the start of the beam element
    (a value <= 0 indicates the equations have been eliminated from the system of equations)
 - `irow_p1`: row index of the first residual equation for the point at the start of the 
    beam element
 - `irow_e2`: row index of the first residual equation for the end of the beam element
    (a value <= 0 indicates the equations have been eliminated from the system of equations)
 - `irow_p2`: row index of the first residual equation for the point at the end of the 
    beam element
 - `f_u1`, `f_u2`: beam element resultant forces
 - `f_ψ1`, `f_ψ2`: beam element resultant moments
 - `f_F1`, `f_F2`: beam element resultant linear displacements
 - `f_M1`, `f_M2`: beam element resultant angular displacements
"""
@inline function static_insert_element_residual!(resid, force_scaling, irow_e1, 
    irow_p1, irow_e2, irow_p2, f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2)

    # create/add to residual equations for the start of the beam
    if irow_e1 == irow_p1
        # add equilibrium equations
        resid[irow_p1:irow_p1+2] .= f_u1 ./ force_scaling
        resid[irow_p1+3:irow_p1+5] .= f_ψ1 ./ force_scaling
        # add compatability equations
        resid[irow_p1+6:irow_p1+8] .= f_F1
        resid[irow_p1+9:irow_p1+11] .= f_M1
    else
        # add to existing equilibrium equations
        resid[irow_p1] += f_u1[1] / force_scaling
        resid[irow_p1+1] += f_u1[2] / force_scaling
        resid[irow_p1+2] += f_u1[3] / force_scaling
        resid[irow_p1+3] += f_ψ1[1] / force_scaling
        resid[irow_p1+4] += f_ψ1[2] / force_scaling
        resid[irow_p1+5] += f_ψ1[3] / force_scaling
        if irow_e1 <= 0
            # compatability equations have been combined with other beam
            resid[irow_p1+6] += f_F1[1]
            resid[irow_p1+7] += f_F1[2]
            resid[irow_p1+8] += f_F1[3]
            resid[irow_p1+9] += f_M1[1]
            resid[irow_p1+10] += f_M1[2]
            resid[irow_p1+11] += f_M1[3]
        else
            # create compatability equations for this beam
            resid[irow_e1:irow_e1+2] .= f_F1
            resid[irow_e1+3:irow_e1+5] .= f_M1
        end
    end

    # create/add to residual equations for the end of the beam
    if irow_e2 == irow_p2
        # add equilibrium equations
        resid[irow_p2:irow_p2+2] .= f_u2 ./ force_scaling
        resid[irow_p2+3:irow_p2+5] .= f_ψ2 ./ force_scaling
        # add compatability equations
        resid[irow_p2+6:irow_p2+8] .= f_F2
        resid[irow_p2+9:irow_p2+11] .= f_M2
    else
        # add to existing equilibrium equations
        resid[irow_p2] += f_u2[1] / force_scaling
        resid[irow_p2+1] += f_u2[2] / force_scaling
        resid[irow_p2+2] += f_u2[3] / force_scaling
        resid[irow_p2+3] += f_ψ2[1] / force_scaling
        resid[irow_p2+4] += f_ψ2[2] / force_scaling
        resid[irow_p2+5] += f_ψ2[3] / force_scaling
        if irow_e2 <= 0
            # compatability equations have been combined with other beam
            resid[irow_p2+6] += f_F2[1]
            resid[irow_p2+7] += f_F2[2]
            resid[irow_p2+8] += f_F2[3]
            resid[irow_p2+9] += f_M2[1]
            resid[irow_p2+10] += f_M2[2]
            resid[irow_p2+11] += f_M2[3]
        else
            # create compatability equations for this beam
            resid[irow_e2:irow_e2+2] .= f_F2
            resid[irow_e2+3:irow_e2+5] .= f_M2
        end
    end

    return resid
end

"""
    dynamic_insert_element_residual!(resid, force_scaling, irow_e, irow_e1, irow_p1, 
        irow_e2, irow_p2, f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2, f_V, f_Ω)

Add beam element resultants into the residuals for a general dynamic analysis.  Initialize 
equilibrium and constitutive equation residuals if they are not yet initialized.

If `irow_e1 != irow_p1` and/or `irow_e2 != irow_p2`, assume the equilibrium equation 
residuals for the start and/or end of the beam element are already initialized

# Arguments
 - `resid`: system residual vector
 - `force_scaling`: scaling parameter for forces
 - `irow_e`: row index of the first linear/angular velocity residual for this element
 - `irow_e1`: row index of the first residual equation for the start of the beam element
   (a value <= 0 indicates the equations have been eliminated from the system of equations)
 - `irow_p1`: row index of the first residual equation for the point at the start of the 
   beam element
 - `irow_e2`: row index of the first residual equation for the end of the beam element
   (a value <= 0 indicates the equations have been eliminated from the system of equations)
 - `irow_p2`: row index of the first residual equation for the point at the end of the 
   beam element
 - `f_u1`, `f_u2`: beam element resultant forces
 - `f_ψ1`, `f_ψ2`: beam element resultant moments
 - `f_F1`, `f_F2`: beam element resultant linear displacements
 - `f_M1`, `f_M2`: beam element resultant angular displacements
 - `f_V`: element linear velocity residual
 - `f_Ω`: element angular velocity residual
"""
@inline function dynamic_insert_element_residual!(resid, force_scaling, irow_e, irow_e1, 
    irow_p1, irow_e2, irow_p2, f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2, f_V, f_Ω)

    static_insert_element_residual!(resid, force_scaling, irow_e1, irow_p1, irow_e2, 
        irow_p2, f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2)

    # residual equations for element
    resid[irow_e:irow_e+2] .= f_V
    resid[irow_e+3:irow_e+5] .= f_Ω

    return resid
end

"""
    static_element_residual!(resid, x, ielem, elem, distributed_loads, 
        point_masses, gvec, force_scaling, icol, irow_e1, irow_p1, irow_e2, irow_p2)

Compute and add a beam element's contributions to the residual vector for a static analysis.

# Arguments
 - `resid`: system residual vector
 - `x`: system state vector
 - `ielem`: beam element index
 - `elem`: beam element
 - `distributed_loads`: dictionary with all distributed loads
 - `point_masses`: dictionary with all point masses
 - `gvec`: gravity vector
 - `force_scaling`: scaling parameter for forces
 - `icol`: starting index for the beam element's state variables
 - `irow_e1`: row index of the first residual equation for the start of the beam element
    (a value <= 0 indicates the equations have been eliminated from the system of equations)
 - `irow_p1`: row index of the first residual equation for the point at the start of the 
    beam element
 - `irow_e2`: row index of the first residual equation for the end of the beam element
    (a value <= 0 indicates the equations have been eliminated from the system of equations)
 - `irow_p2`: row index of the first residual equation for the point at the end of the 
    beam element
"""
@inline function static_element_residual!(resid, x, ielem, elem, distributed_loads, 
    point_masses, gvec, force_scaling, icol, irow_e1, irow_p1, irow_e2, irow_p2)

    # element properties
    ΔL = elem.L
    Cab = elem.Cab
    compliance = elem.compliance
    mass = elem.mass

    # scale compliance and mass matrix by the element length (to allow zero length elements)
    compliance *= ΔL
    mass *= ΔL

    # modify mass matrix to account for point masses, if present
    if haskey(point_masses, ielem)
        mass += transform_properties(point_masses[ielem].mass, Cab)
    end

    # element state variables
    u, θ, F, M = static_element_state_variables(x, icol, force_scaling)

    # compliance submatrices
    S11 = compliance[SVector{3}(1:3), SVector{3}(1:3)]
    S12 = compliance[SVector{3}(1:3), SVector{3}(4:6)]
    S21 = compliance[SVector{3}(4:6), SVector{3}(1:3)]
    S22 = compliance[SVector{3}(4:6), SVector{3}(4:6)]

    # mass submatrices
    mass11 = mass[SVector{3}(1:3), SVector{3}(1:3)]
    mass12 = mass[SVector{3}(1:3), SVector{3}(4:6)]
    mass21 = mass[SVector{3}(4:6), SVector{3}(1:3)]
    mass22 = mass[SVector{3}(4:6), SVector{3}(4:6)]

    # rotation matrices
    C = get_C(θ)
    Ct = C'
    CtCab = C'*Cab

    # element strain and curvature
    γ = S11*F + S12*M
    κ = S21*F + S22*M

    # effective linear and angular acceleration
    ae = -gvec
    αe = zero(ae)

    # element acceleration loads (including gravitational loads)
    f1, f2, m1, m2 = acceleration_loads(mass11, mass12, mass21, mass22, CtCab, ae, αe)

    # distributed loads
    if haskey(distributed_loads, ielem)
        dload = distributed_loads[ielem]
        f1 += dload.f1 + Ct*dload.f1_follower
        f2 += dload.f2 + Ct*dload.f2_follower
        m1 += dload.m1 + Ct*dload.m1_follower
        m2 += dload.m2 + Ct*dload.m2_follower
    end

    # element resultants
    f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2 = static_element_equations(ΔL, Cab, 
        CtCab, u, θ, F, M, γ, κ, f1, f2, m1, m2)

    # insert element resultants into the residual vector
    static_insert_element_residual!(resid, force_scaling, irow_e1, irow_p1, irow_e2, 
        irow_p2, f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2)

    return resid
end

"""
    steady_state_element_residual!(resid, x, ielem, elem, distributed_loads, 
        point_masses, gvec, force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2, 
        x0, v0, ω0, a0, α0)

Compute and add a beam element's contributions to the residual vector for a steady state 
analysis.

# Arguments
 - `resid`: system residual vector
 - `x`: current state vector
 - `ielem`: beam element index
 - `elem`: beam element
 - `distributed_loads`: dictionary with all distributed loads
 - `point_masses`: dictionary with all point masses
 - `gvec`: gravity vector
 - `force_scaling`: scaling parameter for forces
 - `icol`: starting index for the beam's state variables
 - `irow_e1`: row index of the first residual equation for the start of the beam element
    (a value <= 0 indicates the equations have been eliminated from the system of equations)
 - `irow_p1`: row index of the first residual equation for the point at the start of the 
    beam element
 - `irow_e2`: row index of the first residual equation for the end of the beam element
    (a value <= 0 indicates the equations have been eliminated from the system of equations)
 - `irow_p2`: row index of the first residual equation for the point at the end of the 
    beam element
 - `x0`: body frame origin (for the current time step)
 - `v0`: body frame linear velocity (for the current time step)
 - `ω0`: body frame angular velocity (for the current time step)
 - `a0`: body frame linear acceleration (for the current time step)
 - `α0`: body frame angular acceleration (for the current time step)
"""
@inline function steady_state_element_residual!(resid, x, ielem, elem, distributed_loads, 
    point_masses, gvec, force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2, 
    x0, v0, ω0, a0, α0)

    # element properties
    ΔL = elem.L
    Cab = elem.Cab
    compliance = elem.compliance
    mass = elem.mass

    # scale mass matrix by the element length (to allow zero length elements)
    compliance *= ΔL
    mass *= ΔL

    # modify mass matrix to account for point masses, if present
    if haskey(point_masses, ielem)
        mass += transform_properties(point_masses[ielem].mass, Cab)
    end

    # element state variables
    u, θ, F, M, V, Ω = dynamic_element_state_variables(x, icol, force_scaling)

    # compliance submatrices
    S11 = compliance[SVector{3}(1:3), SVector{3}(1:3)]
    S12 = compliance[SVector{3}(1:3), SVector{3}(4:6)]
    S21 = compliance[SVector{3}(4:6), SVector{3}(1:3)]
    S22 = compliance[SVector{3}(4:6), SVector{3}(4:6)]

    # mass submatrices
    mass11 = mass[SVector{3}(1:3), SVector{3}(1:3)]
    mass12 = mass[SVector{3}(1:3), SVector{3}(4:6)]
    mass21 = mass[SVector{3}(4:6), SVector{3}(1:3)]
    mass22 = mass[SVector{3}(4:6), SVector{3}(4:6)]

    # rotation matrices
    C = get_C(θ)
    Ct = C'
    CtCab = Ct*Cab

    # element strain and curvature
    γ = S11*F + S12*M
    κ = S21*F + S22*M

    # element linear and angular momentum
    P = mass11*V + mass12*Ω
    H = mass21*V + mass22*Ω

    # body frame velocities
    v = v0 + cross(ω0, elem.x - x0)
    ω = ω0

    # body frame accelerations
    a = a0 + cross(α0, elem.x - x0)
    α = α0 

    # effective linear and angular acceleration
    ae = a - gvec + cross(α, u)
    αe = α

    # element acceleration loads (including gravitational loads)
    f1, f2, m1, m2 = acceleration_loads(mass11, mass12, mass21, mass22, CtCab, ae, αe)

    # distributed loads
    if haskey(distributed_loads, ielem)
        dload = distributed_loads[ielem]
        f1 += dload.f1 + Ct*dload.f1_follower
        f2 += dload.f2 + Ct*dload.f2_follower
        m1 += dload.m1 + Ct*dload.m1_follower
        m2 += dload.m2 + Ct*dload.m2_follower
    end

    # solve for element resultants
    f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2, f_V, f_Ω = 
        steady_state_element_equations(ΔL, Cab, CtCab, u, θ, F, M, γ, κ, V, Ω, P, H, v, ω, 
        f1, f2, m1, m2)

    # insert element resultants into the residual vector
    dynamic_insert_element_residual!(resid, force_scaling, irow_e, irow_e1, irow_p1, 
        irow_e2, irow_p2, f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2, f_V, f_Ω)

    return resid
end

"""
    initial_condition_element_residual!(resid, x, ielem, elem, distributed_loads, 
        point_masses, gvec, force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, 
        irow_p2, x0, v0, ω0, a0, α0, u0, θ0, Fdot0, Mdot0, udot0, θdot0)

Compute and add a beam element's contributions to the residual vector for an initial 
condition analysis.

# Arguments
 - `resid`: system residual vector
 - `x`: current state vector
 - `ielem`: beam element index
 - `elem`: beam element
 - `distributed_loads`: dictionary with all distributed loads
 - `point_masses`: dictionary with all point masses
 - `gvec`: gravity vector
 - `force_scaling`: scaling parameter for forces
 - `icol`: starting index for the beam's state variables
 - `irow_e1`: row index of the first residual equation for the start of the beam element
    (a value <= 0 indicates the equations have been eliminated from the system of equations)
 - `irow_p1`: row index of the first residual equation for the point at the start of the 
    beam element
 - `irow_e2`: row index of the first residual equation for the end of the beam element
    (a value <= 0 indicates the equations have been eliminated from the system of equations)
 - `irow_p2`: row index of the first residual equation for the point at the end of the 
    beam element
 - `x0`: body frame origin (for the current time step)
 - `v0`: body frame linear velocity (for the current time step)
 - `ω0`: body frame angular velocity (for the current time step)
 - `a0`: body frame linear acceleration (for the current time step)
 - `α0`: body frame angular acceleration (for the current time step)
 - `u0`: initial linear deflections for the beam element
 - `θ0`: initial angular deflections for the beam element
 - `udot0`: initial linear deflection rates for the beam element
 - `θdot0`: initial angular deflection rates for the beam element
 - `Fdot0`: initial resultant force for the beam element
 - `Mdot0`: initial resultant moment for the beam element 
"""
@inline function initial_condition_element_residual!(resid, x, ielem, elem, 
    distributed_loads, point_masses, gvec, force_scaling, icol, irow_e, irow_e1, irow_p1, 
    irow_e2, irow_p2, x0, v0, ω0, a0, α0, u, θ, udot, θdot, Fdot, Mdot)

    # element properties
    ΔL = elem.L
    Cab = elem.Cab
    compliance = elem.compliance
    mass = elem.mass
    μ = elem.mu

    # scale mass matrix by the element length (to allow zero length elements)
    compliance *= ΔL
    mass *= ΔL

    # modify mass matrix to account for point masses, if present
    if haskey(point_masses, ielem)
        mass += transform_properties(point_masses[ielem].mass, Cab)
    end

    # element state variables
    Vdot, Ωdot, F, M, V, Ω = dynamic_element_state_variables(x, icol, force_scaling)

    # compliance submatrices
    S11 = compliance[SVector{3}(1:3), SVector{3}(1:3)]
    S12 = compliance[SVector{3}(1:3), SVector{3}(4:6)]
    S21 = compliance[SVector{3}(4:6), SVector{3}(1:3)]
    S22 = compliance[SVector{3}(4:6), SVector{3}(4:6)]

    # mass submatrices
    mass11 = mass[SVector{3}(1:3), SVector{3}(1:3)]
    mass12 = mass[SVector{3}(1:3), SVector{3}(4:6)]
    mass21 = mass[SVector{3}(4:6), SVector{3}(1:3)]
    mass22 = mass[SVector{3}(4:6), SVector{3}(4:6)]

    # rotation matrices
    C = get_C(θ)
    Ct = C'
    CtCab = Ct*Cab

    # rotation matrix time derivatives
    Cdot = get_C_t(C, θ, θdot)
    Ctdot = Cdot'
    CtCabdot = Ctdot*Cab

    # element strain and curvature
    γ = S11*F + S12*M
    κ = S21*F + S22*M

    # element linear and angular momentum
    P = mass11*V + mass12*Ω
    H = mass21*V + mass22*Ω

    # element linear and angular momentum rates
    Pdot = mass11*Vdot + mass12*Ωdot
    Hdot = mass21*Vdot + mass22*Ωdot

    # linear and angular momentum (in the body frame)
    CtCabPdot = CtCabdot*P + CtCab*Pdot
    CtCabHdot = CtCabdot*H + CtCab*Hdot

    # body frame velocities
    v = v0 + cross(ω0, elem.x - x0)
    ω = ω0

    # body frame accelerations
    a = a0 + cross(α0, elem.x - x0)
    α = α0 

    # effective linear and angular acceleration
    ae = a - gvec + cross(α, u)
    αe = α

    # element acceleration loads (including gravitational loads)
    f1, f2, m1, m2 = acceleration_loads(mass11, mass12, mass21, mass22, CtCab, ae, αe)

    # distributed loads
    if haskey(distributed_loads, ielem)
        dload = distributed_loads[ielem]
        f1 += dload.f1 + Ct*dload.f1_follower
        f2 += dload.f2 + Ct*dload.f2_follower
        m1 += dload.m1 + Ct*dload.m1_follower
        m2 += dload.m2 + Ct*dload.m2_follower
    end

    # damping loads
    Fd = SVector(μ[1]*Fdot[1], μ[2]*Fdot[2], μ[3]*Fdot[3])
    Md = SVector(μ[4]*Mdot[1], μ[5]*Mdot[2], μ[6]*Mdot[3])

    # total element loads
    F = F + Fd
    M = M + Md

    # solve for the element resultants
    f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2, f_V, f_Ω =
        dynamic_element_equations(ΔL, Cab, CtCab, u, θ, F, M, γ, κ, V, Ω, P, H, 
        v, ω, f1, f2, m1, m2, udot, θdot, CtCabPdot, CtCabHdot)

    # insert element resultants into the residual vector
    dynamic_insert_element_residual!(resid, force_scaling, irow_e, irow_e1, irow_p1, 
        irow_e2, irow_p2, f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2, f_V, f_Ω)

    return resid
end

"""
    newmark_element_residual!(resid, x, ielem, elem, distributed_loads, point_masses, gvec,
        force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2, x0, v0, ω0, a0, α0, 
        udot_init, θdot_init, Fdot_init, Mdot_init, Vdot_init, Ωdot_init, dt)

Compute and add a beam element's contributions to the residual vector for a Newmark scheme 
time-marching analysis.

# Arguments
 - `resid`: system residual vector
 - `x`: current state vector
 - `ielem`: beam element index
 - `elem`: beam element
 - `distributed_loads`: dictionary with all distributed loads
 - `point_masses`: dictionary with all point masses
 - `gvec`: gravity vector
 - `force_scaling`: scaling parameter for forces
 - `icol`: starting index for the beam's state variables
 - `irow_e1`: row index of the first residual equation for the start of the beam element
    (a value <= 0 indicates the equations have been eliminated from the system of equations)
 - `irow_p1`: row index of the first residual equation for the point at the start of the 
    beam element
 - `irow_e2`: row index of the first residual equation for the end of the beam element
    (a value <= 0 indicates the equations have been eliminated from the system of equations)
 - `irow_p2`: row index of the first residual equation for the point at the end of the 
    beam element
 - `x0`: body frame origin (for the current time step)
 - `v0`: body frame linear velocity (for the current time step)
 - `ω0`: body frame angular velocity (for the current time step)
 - `a0`: body frame linear acceleration (for the current time step)
 - `α0`: body frame angular acceleration (for the current time step)
 - `udot_init`: `2/dt*u + udot` for the beam element from the previous time step
 - `θdot_init`: `2/dt*θ + θdot` for the beam element from the previous time step
 - `Fdot_init`: `2/dt*F + Fdot` for the beam element from the previous time step
 - `Mdot_init`: `2/dt*M + Mdot` for the beam element from the previous time step
 - `Vdot_init`: `2/dt*V + Vdot` for the beam element from the previous time step
 - `Ωdot_init`: `2/dt*Ω + Ωdot` for the beam element from the previous time step
 - `dt`: time step size
"""
@inline function newmark_element_residual!(resid, x, ielem, elem, distributed_loads, point_masses, gvec,
    force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2, x0, v0, ω0, a0, α0, 
    udot_init, θdot_init, Fdot_init, Mdot_init, Vdot_init, Ωdot_init, dt)

    # element properties
    ΔL = elem.L
    Cab = elem.Cab
    compliance = elem.compliance
    mass = elem.mass
    μ = elem.mu

    # scale mass matrix by the element length (to allow zero length elements)
    compliance *= ΔL
    mass *= ΔL

    # modify mass matrix to account for point masses, if present
    if haskey(point_masses, ielem)
        mass += transform_properties(point_masses[ielem].mass, Cab)
    end

    # element state variables
    u, θ, F, M, V, Ω = dynamic_element_state_variables(x, icol, force_scaling)
    
    # element state variable rates
    udot = 2/dt*u - udot_init
    θdot = 2/dt*θ - θdot_init
    Fdot = 2/dt*F - Fdot_init
    Mdot = 2/dt*M - Mdot_init
    Vdot = 2/dt*V - Vdot_init
    Ωdot = 2/dt*Ω - Ωdot_init

    # compliance submatrices
    S11 = compliance[SVector{3}(1:3), SVector{3}(1:3)]
    S12 = compliance[SVector{3}(1:3), SVector{3}(4:6)]
    S21 = compliance[SVector{3}(4:6), SVector{3}(1:3)]
    S22 = compliance[SVector{3}(4:6), SVector{3}(4:6)]

    # mass submatrices
    mass11 = mass[SVector{3}(1:3), SVector{3}(1:3)]
    mass12 = mass[SVector{3}(1:3), SVector{3}(4:6)]
    mass21 = mass[SVector{3}(4:6), SVector{3}(1:3)]
    mass22 = mass[SVector{3}(4:6), SVector{3}(4:6)]

    # rotation matrices
    C = get_C(θ)
    Ct = C'
    CtCab = Ct*Cab

    # rotation matrix time derivatives
    Cdot = get_C_t(C, θ, θdot)
    Ctdot = Cdot'
    CtCabdot = Ctdot*Cab

    # element strain and curvature
    γ = S11*F + S12*M
    κ = S21*F + S22*M

    # element linear and angular momentum
    P = mass11*V + mass12*Ω
    H = mass21*V + mass22*Ω

    # element linear and angular momentum rates
    Pdot = mass11*Vdot + mass12*Ωdot
    Hdot = mass21*Vdot + mass22*Ωdot

    # linear and angular momentum (in the body frame)
    CtCabPdot = CtCabdot*P + CtCab*Pdot
    CtCabHdot = CtCabdot*H + CtCab*Hdot

    # body frame velocities
    v = v0 + cross(ω0, elem.x - x0)
    ω = ω0

    # body frame accelerations
    a = a0 + cross(α0, elem.x - x0)
    α = α0 

    # effective linear and angular acceleration
    ae = a - gvec + cross(α, u)
    αe = α

    # element acceleration loads (including gravitational loads)
    f1, f2, m1, m2 = acceleration_loads(mass11, mass12, mass21, mass22, CtCab, ae, αe)

    # distributed loads
    if haskey(distributed_loads, ielem)
        dload = distributed_loads[ielem]
        f1 += dload.f1 + Ct*dload.f1_follower
        f2 += dload.f2 + Ct*dload.f2_follower
        m1 += dload.m1 + Ct*dload.m1_follower
        m2 += dload.m2 + Ct*dload.m2_follower
    end

    # damping loads
    Fd = SVector(μ[1]*Fdot[1], μ[2]*Fdot[2], μ[3]*Fdot[3])
    Md = SVector(μ[4]*Mdot[1], μ[5]*Mdot[2], μ[6]*Mdot[3])

    # total element loads
    F = F + Fd
    M = M + Md

    # solve for element resultants
    f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2, f_V, f_Ω =
        dynamic_element_equations(ΔL, Cab, CtCab, u, θ, F, M, γ, κ, V, Ω, P, H, 
        v, ω, f1, f2, m1, m2, udot, θdot, CtCabPdot, CtCabHdot)

    # insert element resultants into the residual vector
    dynamic_insert_element_residual!(resid, force_scaling, irow_e, irow_e1, irow_p1, 
        irow_e2, irow_p2, f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2, f_V, f_Ω)

    return resid
end

"""
    dynamic_element_residual!(resid, x, ielem, elem, distributed_loads, point_masses, gvec,
        force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2, x0, v0, ω0, a0, α0,
        udot, θdot, Fdot, Mdot, Vdot, Ωdot)

Compute and add a beam element's contributions to the residual vector for a general dynamic 
analysis.

# Arguments
 - `resid`: system residual vector
 - `x`: current state vector
 - `ielem`: beam element index
 - `elem`: beam element
 - `distributed_loads`: dictionary with all distributed loads
 - `point_masses`: dictionary with all point masses
 - `gvec`: gravity vector
 - `force_scaling`: scaling parameter for forces
 - `icol`: starting index for the beam's state variables
 - `irow_e1`: row index of the first residual equation for the start of the beam element
    (a value <= 0 indicates the equations have been eliminated from the system of equations)
 - `irow_p1`: row index of the first residual equation for the point at the start of the 
    beam element
 - `irow_e2`: row index of the first residual equation for the end of the beam element
    (a value <= 0 indicates the equations have been eliminated from the system of equations)
 - `irow_p2`: row index of the first residual equation for the point at the end of the 
    beam element
 - `x0`: body frame origin (for the current time step)
 - `v0`: body frame linear velocity (for the current time step)
 - `ω0`: body frame angular velocity (for the current time step)
 - `a0`: body frame linear acceleration (for the current time step)
 - `α0`: body frame angular acceleration (for the current time step)
 - `udot`: linear deflection rates for the beam element
 - `θdot`: angular deflection rates for the beam element
 - `Fdot`: resultant force rates for the beam element
 - `Mdot`: resultant moment rates for the beam element
 - `Vdot`: linear velocity rates for the beam element
 - `Ωdot`: angular velocity rates for the beam element
"""
@inline function dynamic_element_residual!(resid, x, ielem, elem, distributed_loads, 
    point_masses, gvec, force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2, 
    x0, v0, ω0, a0, α0, udot, θdot, Fdot, Mdot, Vdot, Ωdot)

    # element properties
    ΔL = elem.L
    Cab = elem.Cab
    compliance = elem.compliance
    mass = elem.mass
    μ = elem.mu

    # scale mass matrix by the element length (to allow zero length elements)
    compliance *= ΔL
    mass *= ΔL

    # modify mass matrix to account for point masses, if present
    if haskey(point_masses, ielem)
        mass += transform_properties(point_masses[ielem].mass, Cab)
    end

    # element state variables
    u, θ, F, M, V, Ω = dynamic_element_state_variables(x, icol, force_scaling)

    # compliance submatrices
    S11 = compliance[SVector{3}(1:3), SVector{3}(1:3)]
    S12 = compliance[SVector{3}(1:3), SVector{3}(4:6)]
    S21 = compliance[SVector{3}(4:6), SVector{3}(1:3)]
    S22 = compliance[SVector{3}(4:6), SVector{3}(4:6)]

    # mass submatrices
    mass11 = mass[SVector{3}(1:3), SVector{3}(1:3)]
    mass12 = mass[SVector{3}(1:3), SVector{3}(4:6)]
    mass21 = mass[SVector{3}(4:6), SVector{3}(1:3)]
    mass22 = mass[SVector{3}(4:6), SVector{3}(4:6)]

    # rotation matrices
    C = get_C(θ)
    Ct = C'
    CtCab = Ct*Cab

    # rotation matrix time derivatives
    Cdot = get_C_t(C, θ, θdot)
    Ctdot = Cdot'
    CtCabdot = Ctdot*Cab

    # element strain and curvature
    γ = S11*F + S12*M
    κ = S21*F + S22*M

    # element linear and angular momentum
    P = mass11*V + mass12*Ω
    H = mass21*V + mass22*Ω

    # element linear and angular momentum rates
    Pdot = mass11*Vdot + mass12*Ωdot
    Hdot = mass21*Vdot + mass22*Ωdot

    # linear and angular momentum (in the body frame)
    CtCabPdot = CtCabdot*P + CtCab*Pdot
    CtCabHdot = CtCabdot*H + CtCab*Hdot

    # body frame velocities
    v = v0 + cross(ω0, elem.x - x0)
    ω = ω0

    # body frame accelerations
    a = a0 + cross(α0, elem.x - x0)
    α = α0 

    # effective linear and angular acceleration
    ae = a - gvec + cross(α, u)
    αe = α

    # element acceleration loads (including gravitational loads)
    f1, f2, m1, m2 = acceleration_loads(mass11, mass12, mass21, mass22, CtCab, ae, αe)

    # distributed loads
    if haskey(distributed_loads, ielem)
        dload = distributed_loads[ielem]
        f1 += dload.f1 + Ct*dload.f1_follower
        f2 += dload.f2 + Ct*dload.f2_follower
        m1 += dload.m1 + Ct*dload.m1_follower
        m2 += dload.m2 + Ct*dload.m2_follower
    end

    # damping loads
    Fd = SVector(μ[1]*Fdot[1], μ[2]*Fdot[2], μ[3]*Fdot[3])
    Md = SVector(μ[4]*Mdot[1], μ[5]*Mdot[2], μ[6]*Mdot[3])    

    # total element loads
    F = F + Fd
    M = M + Md

    # solve for element resultants
    f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2, f_V, f_Ω =
        dynamic_element_equations(ΔL, Cab, CtCab, u, θ, F, M, γ, κ, V, Ω, P, H, 
        v, ω, f1, f2, m1, m2, udot, θdot, CtCabPdot, CtCabHdot)

    # insert element resultants into the residual vector
    dynamic_insert_element_residual!(resid, force_scaling, irow_e, irow_e1, irow_p1, 
        irow_e2, irow_p2, f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2, f_V, f_Ω)

    return resid
end

"""
    static_element_jacobian_equations(ΔL, S11, S12, S21, S22, Cab, CtCab, θ, F, M, γ, κ, 
        f1_θ, f2_θ, m1_θ, m2_θ, Ct_θ1, Ct_θ2, Ct_θ3)

Calculate the derivatives of the element resultants with respect to the state variables for 
a static analysis.

# Arguments:
 - `ΔL`: beam element length
 - `S11, S12, S21, S22`: beam element compliance matrix, divided into submatrices
 - `Cab`: transformation matrix from the undeformed beam element frame to the body frame
 - `CtCab`: transformation matrix from the deformed beam element frame to the body frame
 - `θ`: rotation variables for the element
 - `F`: Force variables for the elemen
 - `M`: Moment variables for the element [M1, M2, M3]
 - `γ`: Engineering strains in the element [γ11, 2γ12, 2γ13]
 - `κ`: Curvatures in the element [κ1, κ2, κ3]
 - `f1_θ`: Derivative of `f1` w.r.t `θ` 
 - `f2_θ`: Derivative of `f2` w.r.t `θ` 
 - `m1_θ`: Derivative of `m1` w.r.t `θ` 
 - `m2_θ`: Derivative of `m2` w.r.t `θ` 
 - `Ct_θ1`: Derivative of `C'` w.r.t. `θ[1]`
 - `Ct_θ2`: Derivative of `C'` w.r.t. `θ[2]`
 - `Ct_θ3`: Derivative of `C'` w.r.t. `θ[3]`
"""
@inline function static_element_jacobian_equations(ΔL, S11, S12, S21, S22, Cab, CtCab, 
    θ, F, M, γ, κ, f1_u, f2_u, m1_u, m2_u, f1_θ, f2_θ, m1_θ, m2_θ, Ct_θ1, Ct_θ2, Ct_θ3)

    # note that a ΔL term has been incorporated into the stiffness and mass matrix 

    # --- f_u1, f_u2 --- #

    # d_fu/d_u
    f_u1_u = -f1_u
    f_u2_u = -f2_u

    # d_fu/d_θ
    tmp = mul3(Ct_θ1, Ct_θ2, Ct_θ3, Cab*F)
    f_u1_θ = -tmp - f1_θ
    f_u2_θ =  tmp - f2_θ

    # d_fu/d_F
    f_u1_F = -CtCab
    f_u2_F =  CtCab

    # --- f_ψ1, f_ψ2 --- #

    # d_fψ_du
    f_ψ1_u = -m1_u
    f_ψ2_u = -m2_u

    # d_fψ/d_θ
    tmp1 = mul3(Ct_θ1, Ct_θ2, Ct_θ3, Cab*M)
    tmp2 = 1/2*mul3(Ct_θ1, Ct_θ2, Ct_θ3, Cab*cross(ΔL*e1 + γ, F))
    f_ψ1_θ = -tmp1 - m1_θ - tmp2
    f_ψ2_θ =  tmp1 - m2_θ - tmp2

    # d_fψ/d_F
    tmp = -1/2*CtCab*(tilde(ΔL*e1 + γ) - tilde(F)*S11)
    f_ψ1_F = tmp
    f_ψ2_F = tmp

    # d_fψ/d_M
    tmp = 1/2*CtCab*tilde(F)*S12
    f_ψ1_M = tmp - CtCab
    f_ψ2_M = tmp + CtCab

    # --- f_F1, f_F2 --- #

    # d_fF/d_u
    f_F1_u =  I3
    f_F2_u = -I3

    # d_fF/d_θ
    tmp = mul3(Ct_θ1, Ct_θ2, Ct_θ3, 1/2*Cab*(ΔL*e1 + γ))
    f_F1_θ = -tmp
    f_F2_θ = -tmp

    # d_fF/d_F
    tmp = 1/2*CtCab*S11
    f_F1_F = -tmp
    f_F2_F = -tmp

    # d_fF/d_M
    tmp = 1/2*CtCab*S12
    f_F1_M = -tmp
    f_F2_M = -tmp

    # --- f_M1, f_M2 --- #

    # d_fM/d_θ
    Qinv_θ1, Qinv_θ2, Qinv_θ3 = get_Qinv_θ(θ)
    tmp = mul3(Qinv_θ1, Qinv_θ2, Qinv_θ3, 1/2*Cab*κ)
    f_M1_θ =  I - tmp
    f_M2_θ = -I - tmp

    # d_fM/d_F
    Qinv = get_Qinv(θ)
    tmp1 = -1/2*Qinv*Cab
    tmp2 = tmp1*S21
    f_M1_F = tmp2
    f_M2_F = tmp2

    # d_fM/d_M
    tmp2 = tmp1*S22
    f_M1_M = tmp2
    f_M2_M = tmp2

    return f_u1_u, f_u2_u, f_u1_θ, f_u2_θ, f_u1_F, f_u2_F,
        f_ψ1_u, f_ψ2_u, f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M,
        f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M
end

"""
    steady_state_element_jacobian_equations(ΔL, S11, S12, S21, S22, 
        mass11, mass12, mass21, mass22, Cab, CtCab, θ, F, M, γ, κ, V, P, H, ω, f1_u, f2_u, 
        m1_u, m2_u, f1_θ, f2_θ, m1_θ, m2_θ, C_θ1, C_θ2, C_θ3, Ct_θ1, Ct_θ2, Ct_θ3)

Calculate the derivatives of the element resultants with respect to the state variables for 
a steady state analysis.

# Arguments:
 - `ΔL`: length
 - `S11, S12, S21, S22`: beam element compliance matrix, divided into submatrices
 - `mass11, mass12, mass21, mass22`: beam element mass matrix, divided into submatrices
 - `Cab`: transformation matrix from the undeformed beam element frame to the body frame
 - `CtCab`: transformation matrix from the deformed beam element frame to the body frame
 - `θ`: angular displacement
 - `F`: internal force
 - `M`: internal moment
 - `γ`: engineering strain
 - `κ`: curvature
 - `V`: linear velocity
 - `P`: linear momentum
 - `H`: angular momentum
 - `ω`: angular velocity caused by body frame velocity
 - `f1_u`: Derivative of `f1` w.r.t `u` 
 - `f2_u`: Derivative of `f2` w.r.t `u` 
 - `m1_u`: Derivative of `m1` w.r.t `u` 
 - `m2_u`: Derivative of `m2` w.r.t `u` 
 - `f1_θ`: Derivative of `f1` w.r.t `θ` 
 - `f2_θ`: Derivative of `f2` w.r.t `θ` 
 - `m1_θ`: Derivative of `m1` w.r.t `θ` 
 - `m2_θ`: Derivative of `m2` w.r.t `θ` 
 - `C_θ1`: Derivative of `C` w.r.t. `θ[1]`
 - `C_θ2`: Derivative of `C` w.r.t. `θ[2]`
 - `C_θ3`: Derivative of `C` w.r.t. `θ[3]`
 - `Ct_θ1`: Derivative of `C'` w.r.t. `θ[1]` (transpose of `C_θ1`)
 - `Ct_θ2`: Derivative of `C'` w.r.t. `θ[2]` (transpose of `C_θ2`)
 - `Ct_θ3`: Derivative of `C'` w.r.t. `θ[3]` (transpose of `C_θ3`)
"""
@inline function steady_state_element_jacobian_equations(ΔL, S11, S12, S21, S22, 
    mass11, mass12, mass21, mass22, Cab, CtCab, θ, F, M, γ, κ, V, P, H, ω, f1_u, f2_u, 
    m1_u, m2_u, f1_θ, f2_θ, m1_θ, m2_θ, C_θ1, C_θ2, C_θ3, Ct_θ1, Ct_θ2, Ct_θ3)

    f_u1_u, f_u2_u, f_u1_θ, f_u2_θ, f_u1_F, f_u2_F,
        f_ψ1_u, f_ψ2_u, f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M,
        f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M = static_element_jacobian_equations(
        ΔL, S11, S12, S21, S22, Cab, CtCab, θ, F, M, γ, κ, f1_u, f2_u, m1_u, m2_u, 
        f1_θ, f2_θ, m1_θ, m2_θ, Ct_θ1, Ct_θ2, Ct_θ3)

    # note that a ΔL term has been incorporated into the stiffness and mass matrix 

    # --- f_u1, f_u2 --- #

    # d_fu_dθ
    tmp = 1/2*tilde(ω)*mul3(Ct_θ1, Ct_θ2, Ct_θ3, Cab*P)
    f_u1_θ += tmp
    f_u2_θ += tmp

    # d_fu_dV
    tmp = 1/2*tilde(ω)*CtCab*mass11
    f_u1_V = tmp
    f_u2_V = tmp

    # d_fu_dΩ
    tmp = 1/2*tilde(ω)*CtCab*mass12
    f_u1_Ω = tmp
    f_u2_Ω = tmp

    # --- f_ψ1, f_ψ2 --- #

    # d_fψ_dθ
    tmp = 1/2*(tilde(ω)*mul3(Ct_θ1, Ct_θ2, Ct_θ3, Cab*H) + 
        mul3(Ct_θ1, Ct_θ2, Ct_θ3, Cab*cross(V, P)))
    f_ψ1_θ += tmp
    f_ψ2_θ += tmp

    # d_fψ_dV
    tmp = 1/2*(tilde(ω)*CtCab*mass21 + CtCab*(tilde(V)*mass11 - tilde(P)))
    f_ψ1_V = tmp
    f_ψ2_V = tmp

    # d_fψ_dΩ
    tmp = 1/2*(tilde(ω)*CtCab*mass22 + CtCab*(tilde(V)*mass12))
    f_ψ1_Ω = tmp
    f_ψ2_Ω = tmp

    # --- f_V --- #

    # d_fP_du
    f_V_u = -tilde(ω)

    # d_fP_dθ
    f_V_θ = mul3(Ct_θ1, Ct_θ2, Ct_θ3, Cab*V)

    # d_fP_dV
    f_V_V = CtCab

    # --- f_Ω --- #

    # d_fH_dθ
    f_Ω_θ = -Cab'*mul3(C_θ1, C_θ2, C_θ3, ω)

    # d_fH_dΩ
    f_Ω_Ω = I3

    return f_u1_u, f_u2_u, f_u1_θ, f_u2_θ, f_u1_F, f_u2_F, f_u1_V, f_u2_V, f_u1_Ω, f_u2_Ω,
        f_ψ1_u, f_ψ2_u, f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_V, f_ψ2_V, f_ψ1_Ω, f_ψ2_Ω,
        f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M,
        f_V_u, f_V_θ, f_V_V,
        f_Ω_θ, f_Ω_Ω
end

"""
    initial_condition_element_jacobian_equations(ΔL, S11, S12, S21, S22, 
        mass11, mass12, mass21, mass22, Cab, CtCab, CtCabdot, θ, F, γ, V, P, ω)

Calculate the derivatives of the element resultants with respect to the state variables for 
an initial condition analysis.
        
# Arguments:
 - `ΔL`: beam element length
 - `S11, S12, S21, S22`: beam element compliance matrix, divided into submatrices
 - `mass11, mass12, mass21, mass22`: beam element mass matrix, divided into submatrices
 - `Cab`: transformation matrix from the undeformed beam element frame to the body frame
 - `CtCab`: transformation matrix from the deformed beam element frame to the body frame
 - `CtCabdot`: time derivative of `CtCab`
 - `θ`: angular displacement
 - `F`: internal force
 - `γ`: engineering strain
 - `V`: linear velocity
 - `P`: linear momentum
 - `ω`: angular velocity caused by body frame velocity
"""
@inline function initial_condition_element_jacobian_equations(ΔL, S11, S12, S21, S22, 
    mass11, mass12, mass21, mass22, Cab, CtCab, CtCabdot, θ, F, γ, V, P, ω)

    # note that a ΔL term has been incorporated into the stiffness and mass matrix 

    # --- f_u1, f_u2 --- #

    # d_fu/d_Vdot
    tmp = 1/2*CtCab*mass11
    f_u1_Vdot = tmp
    f_u2_Vdot = tmp

    # d_fu/d_Ωdot
    tmp = 1/2*CtCab*mass12
    f_u1_Ωdot = tmp
    f_u2_Ωdot = tmp

    # d_fu/d_F
    tmp = CtCab
    f_u1_F = -tmp
    f_u2_F =  tmp

    # d_fu/dV
    tmp = 1/2*(tilde(ω)*CtCab*mass11 + CtCabdot*mass11)
    f_u1_V = tmp 
    f_u2_V = tmp 

    # d_fu/dΩ
    tmp = 1/2*(tilde(ω)*CtCab*mass12 + CtCabdot*mass12)
    f_u1_Ω = tmp 
    f_u2_Ω = tmp 

    # --- f_θ1, f_θ2 --- #

    # d_fψ/d_Ωdot
    tmp = 1/2*CtCab*mass21
    f_ψ1_Vdot = tmp 
    f_ψ2_Vdot = tmp

    # d_fψ/d_Ωdot
    tmp = 1/2*CtCab*mass22
    f_ψ1_Ωdot = tmp
    f_ψ2_Ωdot = tmp

    # d_fψ/d_F
    tmp = -1/2*CtCab*(tilde(ΔL*e1 + γ) - tilde(F)*S11)
    f_ψ1_F = tmp
    f_ψ2_F = tmp

    # d_fψ/d_M
    tmp = 1/2*CtCab*tilde(F)*S12
    f_ψ1_M = tmp - CtCab
    f_ψ2_M = tmp + CtCab

    # d_fψ_dV
    tmp = 1/2*(tilde(ω)*CtCab*mass21 + CtCabdot*mass21 + CtCab*(tilde(V)*mass11 - tilde(P)))
    f_ψ1_V = tmp 
    f_ψ2_V = tmp 

    # d_fψ_dΩ
    tmp = 1/2*(tilde(ω)*CtCab*mass22 + CtCabdot*mass22 + CtCab*(tilde(V)*mass12))
    f_ψ1_Ω = tmp
    f_ψ2_Ω = tmp

    # --- f_F1, f_F2 --- #

    # d_fF/d_F
    tmp = 1/2*CtCab*S11
    f_F1_F = -tmp
    f_F2_F = -tmp

    # d_fF/d_M
    tmp = 1/2*CtCab*S12
    f_F1_M = -tmp
    f_F2_M = -tmp

    # --- f_M1, f_M2 --- #

    # d_fM/d_F
    Qinv = get_Qinv(θ)
    tmp1 = -1/2*Qinv*Cab
    tmp2 = tmp1*S21
    f_M1_F = tmp2
    f_M2_F = tmp2

    # d_fM/d_M
    tmp2 = tmp1*S22
    f_M1_M = tmp2
    f_M2_M = tmp2

    # --- f_V --- #

    # d_fP_dV
    f_V_V = CtCab

    # --- f_Ω --- #

    # d_fH_dΩ
    f_Ω_Ω = I3

    return f_u1_Vdot, f_u2_Vdot, f_u1_Ωdot, f_u2_Ωdot, f_u1_F, f_u2_F, f_u1_V, f_u2_V, f_u1_Ω, f_u2_Ω,
        f_ψ1_Vdot, f_ψ2_Vdot, f_ψ1_Ωdot, f_ψ2_Ωdot, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_V, f_ψ2_V, f_ψ1_Ω, f_ψ2_Ω,
        f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_F, f_M2_F, f_M1_M, f_M2_M,
        f_V_V, f_Ω_Ω
end

"""
    newmark_element_jacobian_equations(ΔL, S11, S12, S21, S22, mass11, mass12, 
        mass21, mass22, μ11, μ22, Cab, CtCab, CtCabdot, θ, F, M, γ, κ, V, P, H, θdot, Pdot, Hdot, ω, 
        dt, f1_u, f2_u, m1_u, m2_u, C_θ1, C_θ2, C_θ3, Ct_θ1, Ct_θ2, Ct_θ3, Ctdot_θ1, 
        Ctdot_θ2, Ctdot_θ3, Ctdot_θdot1, Ctdot_θdot2, Ctdot_θdot3)

Calculate the derivatives of the element resultants with respect to the state variables for 
an Newmark scheme time-marching analysis.
        
# Arguments:
 - `ΔL`: length
 - `S11, S12, S21, S22`: beam element compliance matrix, divided into submatrices
 - `mass11, mass12, mass21, mass22`: beam element mass matrix, divided into submatrices
 - `Cab`: transformation matrix from the undeformed beam element frame to the body frame
 - `CtCab`: transformation matrix from the deformed beam element frame to the body frame
 - `CtCabdot`: time derivative of `CtCab`
 - `θ`: angular displacement
 - `F`: internal force
 - `M`: internal moment
 - `γ`: engineering strain
 - `κ`: curvature
 - `V`: linear velocity
 - `P`: linear momentum
 - `H`: angular momentum
 - `θdot`: angular displacement rate
 - `Pdot`: linear momentum rate
 - `Hdot`: angular momentum rate
 - `ω`: angular velocity caused by body frame velocity
 - `dt`: time step size
 - `C_θ1`: Derivative of `C` w.r.t. `θ[1]`
 - `C_θ2`: Derivative of `C` w.r.t. `θ[2]`
 - `C_θ3`: Derivative of `C` w.r.t. `θ[3]`
 - `Ct_θ1`: Derivative of `C'` w.r.t. `θ[1]` (transpose of `C_θ1`)
 - `Ct_θ2`: Derivative of `C'` w.r.t. `θ[2]` (transpose of `C_θ2`)
 - `Ct_θ3`: Derivative of `C'` w.r.t. `θ[3]` (transpose of `C_θ3`)
 - `Ctdot_θ1`: Derivative of `Cdot'` w.r.t. `θ[1]`
 - `Ctdot_θ2`: Derivative of `Cdot'` w.r.t. `θ[2]`
 - `Ctdot_θ3`: Derivative of `Cdot'` w.r.t. `θ[3]`
 - `Ctdot_θdot1`: Derivative of `Cdot'` w.r.t. `θdot[1]`
 - `Ctdot_θdot2`: Derivative of `Cdot'` w.r.t. `θdot[2]`
 - `Ctdot_θdot3`: Derivative of `Cdot'` w.r.t. `θdot[3]`
"""
@inline function newmark_element_jacobian_equations(ΔL, S11, S12, S21, S22, mass11, mass12, 
    mass21, mass22, μ11, μ22, Cab, CtCab, CtCabdot, θ, F, M, γ, κ, V, P, H, θdot, Pdot, Hdot, ω, dt,
    f1_u, f2_u, m1_u, m2_u, f1_θ, f2_θ, m1_θ, m2_θ, C_θ1, C_θ2, C_θ3, Ct_θ1, Ct_θ2, Ct_θ3, 
    Ctdot_θ1, Ctdot_θ2, Ctdot_θ3, Ctdot_θdot1, Ctdot_θdot2, Ctdot_θdot3)

    f_u1_u, f_u2_u, f_u1_θ, f_u2_θ, f_u1_F, f_u2_F, f_u1_V, f_u2_V, f_u1_Ω, f_u2_Ω,
        f_ψ1_u, f_ψ2_u, f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_V, f_ψ2_V, f_ψ1_Ω, f_ψ2_Ω,
        f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M,
        f_V_u, f_V_θ, f_V_V,
        f_Ω_θ, f_Ω_Ω = steady_state_element_jacobian_equations(ΔL, S11, S12, S21, S22, 
        mass11, mass12, mass21, mass22, Cab, CtCab, θ, F, M, γ, κ, V, P, H, ω, f1_u, f2_u, 
        m1_u, m2_u, f1_θ, f2_θ, m1_θ, m2_θ, C_θ1, C_θ2, C_θ3, Ct_θ1, Ct_θ2, Ct_θ3)

    # note that a ΔL term has been incorporated into the stiffness and mass matrix 

    # --- f_u1, f_u2 --- #

    # d_fu_dθ
    tmp = 1/2*(
        mul3(Ctdot_θ1, Ctdot_θ2, Ctdot_θ3, Cab*P) +
        2/dt*mul3(Ctdot_θdot1, Ctdot_θdot2, Ctdot_θdot3, Cab*P) + 
        mul3(Ct_θ1, Ct_θ2, Ct_θ3, Cab*Pdot)
    )
    f_u1_θ += tmp
    f_u2_θ += tmp

    # d_fu_dF
    tmp = 2/dt*CtCab*μ11
    f_u1_F += -tmp
    f_u2_F += tmp

    # d_fu_dV      
    tmp = 1/2*(CtCabdot*mass11 + 2/dt*CtCab*mass11)
    f_u1_V += tmp
    f_u2_V += tmp

    # d_fu_dΩ
    tmp = 1/2*(CtCabdot*mass12 + 2/dt*CtCab*mass12)
    f_u1_Ω += tmp
    f_u2_Ω += tmp

    # --- f_ψ1, f_ψ2 --- #

    # d_fψ_dθ
    tmp = 1/2*(
        mul3(Ctdot_θ1, Ctdot_θ2, Ctdot_θ3, Cab*H) +
        2/dt*mul3(Ctdot_θdot1, Ctdot_θdot2, Ctdot_θdot3, Cab*H) +
        mul3(Ct_θ1, Ct_θ2, Ct_θ3, Cab*Hdot)
    )
    f_ψ1_θ += tmp
    f_ψ2_θ += tmp

    # d_fψ_dF
    tmp = -1/2*CtCab*(tilde(ΔL*e1 + γ)*2/dt*μ11)
    f_ψ1_F += tmp
    f_ψ2_F += tmp

    # d_fψ_dM
    tmp = 2/dt*CtCab*μ22
    f_ψ1_M += -tmp
    f_ψ2_M += tmp
  
    # d_fψ_dV
    tmp = 1/2*(CtCabdot*mass21 + 2/dt*CtCab*mass21)
    f_ψ1_V += tmp
    f_ψ2_V += tmp

    # d_fψ_dΩ
    tmp = 1/2*(CtCabdot*mass22 + 2/dt*CtCab*mass22)
    f_ψ1_Ω += tmp
    f_ψ2_Ω += tmp

    # --- d_fP_du --- #
    f_V_u -= 2/dt*I

    # --- d_fH_dθ --- #
    Q = get_Q(θ)
    Q_θ1, Q_θ2, Q_θ3 = get_Q_θ(θ)
    f_Ω_θ -= Cab'*(mul3(Q_θ1, Q_θ2, Q_θ3, θdot) + Q*2/dt)

    return f_u1_u, f_u2_u, f_u1_θ, f_u2_θ, f_u1_F, f_u2_F, f_u1_V, f_u2_V, f_u1_Ω, f_u2_Ω,
        f_ψ1_u, f_ψ2_u, f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_V, f_ψ2_V, f_ψ1_Ω, f_ψ2_Ω,
        f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M,
        f_V_u, f_V_θ, f_V_V,
        f_Ω_θ, f_Ω_Ω
end

"""
    dynamic_element_jacobian_equations(ΔL, S11, S12, S21, S22, mass11, mass12, 
        mass21, mass22, Cab, CtCab, CtCabdot, θ, F, M, γ, κ, V, P, H, θdot, Pdot, Hdot, ω, 
        f1_u, f2_u, m1_u, m2_u, f1_θ, f2_θ, m1_θ, m2_θ, 
        C_θ1, C_θ2, C_θ3, Ct_θ1, Ct_θ2, Ct_θ3, Ctdot_θ1, Ctdot_θ2, Ctdot_θ3)

Calculate the derivatives of the element resultants with respect to the state variables for 
an Newmark scheme time-marching analysis.
        
# Arguments:
 - `ΔL`: length
 - `S11, S12, S21, S22`: beam element compliance matrix, divided into submatrices
 - `mass11, mass12, mass21, mass22`: beam element mass matrix, divided into submatrices
 - `Cab`: transformation matrix from the undeformed beam element frame to the body frame
 - `CtCab`: transformation matrix from the deformed beam element frame to the body frame
 - `CtCabdot`: time derivative of `CtCab`
 - `θ`: angular displacement
 - `F`: internal force
 - `M`: internal moment
 - `γ`: engineering strain
 - `κ`: curvature
 - `V`: linear velocity
 - `P`: linear momentum
 - `H`: angular momentum
 - `θdot`: angular displacement rate
 - `Pdot`: linear momentum rate
 - `Hdot`: angular momentum rate
 - `ω`: angular velocity caused by body frame velocity
 - `f1_u`: Derivative of `f1` w.r.t `u` 
 - `f2_u`: Derivative of `f2` w.r.t `u` 
 - `m1_u`: Derivative of `m1` w.r.t `u` 
 - `m2_u`: Derivative of `m2` w.r.t `u` 
 - `f1_θ`: Derivative of `f1` w.r.t `θ` 
 - `f2_θ`: Derivative of `f2` w.r.t `θ` 
 - `m1_θ`: Derivative of `m1` w.r.t `θ` 
 - `m2_θ`: Derivative of `m2` w.r.t `θ` 
 - `C_θ1`: Derivative of `C` w.r.t. `θ[1]`
 - `C_θ2`: Derivative of `C` w.r.t. `θ[2]`
 - `C_θ3`: Derivative of `C` w.r.t. `θ[3]`
 - `Ct_θ1`: Derivative of `C'` w.r.t. `θ[1]` (transpose of `C_θ1`)
 - `Ct_θ2`: Derivative of `C'` w.r.t. `θ[2]` (transpose of `C_θ2`)
 - `Ct_θ3`: Derivative of `C'` w.r.t. `θ[3]` (transpose of `C_θ3`)
 - `Ctdot_θ1`: Derivative of `Cdot'` w.r.t. `θ[1]`
 - `Ctdot_θ2`: Derivative of `Cdot'` w.r.t. `θ[2]`
 - `Ctdot_θ3`: Derivative of `Cdot'` w.r.t. `θ[3]`
"""
@inline function dynamic_element_jacobian_equations(ΔL, S11, S12, S21, S22, mass11, mass12, 
    mass21, mass22, Cab, CtCab, CtCabdot, θ, F, M, γ, κ, V, P, H, θdot, Pdot, Hdot, ω, 
    f1_u, f2_u, m1_u, m2_u, f1_θ, f2_θ, m1_θ, m2_θ, 
    C_θ1, C_θ2, C_θ3, Ct_θ1, Ct_θ2, Ct_θ3, Ctdot_θ1, Ctdot_θ2, Ctdot_θ3)

    f_u1_u, f_u2_u, f_u1_θ, f_u2_θ, f_u1_F, f_u2_F, f_u1_V, f_u2_V, f_u1_Ω, f_u2_Ω,
        f_ψ1_u, f_ψ2_u, f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_V, f_ψ2_V, f_ψ1_Ω, f_ψ2_Ω,
        f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M,
        f_V_u, f_V_θ, f_V_V,
        f_Ω_θ, f_Ω_Ω = steady_state_element_jacobian_equations(ΔL, S11, S12, S21, S22, 
        mass11, mass12, mass21, mass22, Cab, CtCab, θ, F, M, γ, κ, V, P, H, ω, f1_u, f2_u, 
        m1_u, m2_u, f1_θ, f2_θ, m1_θ, m2_θ, 
        C_θ1, C_θ2, C_θ3, Ct_θ1, Ct_θ2, Ct_θ3)

    # note that a ΔL term has been incorporated into the stiffness and mass matrix 

    # --- f_u1, f_u2 --- #

    # d_fu_dθ
    tmp = 1/2*(mul3(Ctdot_θ1, Ctdot_θ2, Ctdot_θ3, Cab*P) + mul3(Ct_θ1, Ct_θ2, Ct_θ3, Cab*Pdot))
    f_u1_θ += tmp
    f_u2_θ += tmp

    # d_fu_dV
    tmp = 1/2*CtCabdot*mass11
    f_u1_V += tmp
    f_u2_V += tmp

    # d_fu_dΩ
    tmp = 1/2*CtCabdot*mass12
    f_u1_Ω += tmp
    f_u2_Ω += tmp

    # --- f_ψ1, f_ψ2 --- #

    # d_fψ_dθ
    tmp = 1/2*(mul3(Ctdot_θ1, Ctdot_θ2, Ctdot_θ3, Cab*H) + mul3(Ct_θ1, Ct_θ2, Ct_θ3, Cab*Hdot))
    f_ψ1_θ += tmp
    f_ψ2_θ += tmp

    # d_fψ_dV
    tmp = 1/2*CtCabdot*mass21
    f_ψ1_V += tmp
    f_ψ2_V += tmp

    # d_fψ_dΩ
    tmp = 1/2*CtCabdot*mass22
    f_ψ1_Ω += tmp
    f_ψ2_Ω += tmp

    # --- d_fH_dθ --- #
    Q_θ1, Q_θ2, Q_θ3 = get_Q_θ(θ)
    f_Ω_θ -= Cab'*mul3(Q_θ1, Q_θ2, Q_θ3, θdot)

    return f_u1_u, f_u2_u, f_u1_θ, f_u2_θ, f_u1_F, f_u2_F, f_u1_V, f_u2_V, f_u1_Ω, f_u2_Ω,
        f_ψ1_u, f_ψ2_u, f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_V, f_ψ2_V, f_ψ1_Ω, f_ψ2_Ω,
        f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M,
        f_V_u, f_V_θ, f_V_V,
        f_Ω_θ, f_Ω_Ω
end

"""
    static_insert_element_jacobian!(jacob, force_scaling, irow_e1,
        irow_p1, irow_e2, irow_p2, icol,
        f_u1_u, f_u2_u, f_u1_θ, f_u2_θ, f_u1_F, f_u2_F,
        f_ψ1_u, f_ψ2_u, f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M,
        f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M)

Insert the the beam element jacobian entries into the jacobian matrix for a static analysis

# Arguments
 - `jacob`: System jacobian matrix
 - `force_scaling`: Scaling parameter for forces
 - `irow_e1`: row index of the first residual equation for the start of the beam element
   (a value <= 0 indicates the equations have been eliminated from the system of equations)
 - `irow_p1`: row index of the first residual equation for the point at the start of the 
   beam element
 - `irow_e2`: row index of the first residual equation for the end of the beam element
   (a value <= 0 indicates the equations have been eliminated from the system of equations)
 - `irow_p2`: row index of the first residual equation for the point at the end of the 
   beam element
 - `icol`: Row/Column index corresponding to the first beam state variable

All other arguments use the following naming convention:
 - `f_y_x`: Jacobian of resultant `y` with respect to state variable `x`
"""
@inline function static_insert_element_jacobian!(jacob, force_scaling, irow_e1,
    irow_p1, irow_e2, irow_p2, icol,
    f_u1_u, f_u2_u, f_u1_θ, f_u2_θ, f_u1_F, f_u2_F,
    f_ψ1_u, f_ψ2_u, f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M,
    f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
    f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M)

    # equilibrium equation jacobian entries for the start of the beam
    jacob[irow_p1:irow_p1+2, icol:icol+2] .= f_u1_u ./ force_scaling
    jacob[irow_p1:irow_p1+2, icol+3:icol+5] .= f_u1_θ ./ force_scaling
    jacob[irow_p1:irow_p1+2, icol+6:icol+8] .= f_u1_F

    jacob[irow_p1+3:irow_p1+5, icol:icol+2] .= f_ψ1_u ./ force_scaling
    jacob[irow_p1+3:irow_p1+5, icol+3:icol+5] .= f_ψ1_θ ./ force_scaling
    jacob[irow_p1+3:irow_p1+5, icol+6:icol+8] .= f_ψ1_F
    jacob[irow_p1+3:irow_p1+5, icol+9:icol+11] .= f_ψ1_M

    # compatability equation jacobian entries for the start of the beam
    # if irow_e1 == irow_p1 use row corresponding to compatability equations for this beam
    # if irow_e1 <= 0 use row corresponding to compatability equations for the other beam
    # otherwise use row corresponding to compatability equations for this beam
    irow = ifelse(irow_e1 == irow_p1 || irow_e1 <= 0, irow_p1+6, irow_e1)

    jacob[irow:irow+2, icol:icol+2] .= f_F1_u
    jacob[irow:irow+2, icol+3:icol+5] .= f_F1_θ
    jacob[irow:irow+2, icol+6:icol+8] .= f_F1_F .* force_scaling
    jacob[irow:irow+2, icol+9:icol+11] .= f_F1_M .* force_scaling

    jacob[irow+3:irow+5, icol+3:icol+5] .= f_M1_θ
    jacob[irow+3:irow+5, icol+6:icol+8] .= f_M1_F .* force_scaling
    jacob[irow+3:irow+5, icol+9:icol+11] .= f_M1_M .* force_scaling

    # equilibrium equation jacobian entries for the end of the beam
    jacob[irow_p2:irow_p2+2, icol:icol+2] .= f_u2_u ./ force_scaling
    jacob[irow_p2:irow_p2+2, icol+3:icol+5] .= f_u2_θ ./ force_scaling
    jacob[irow_p2:irow_p2+2, icol+6:icol+8] .= f_u2_F

    jacob[irow_p2+3:irow_p2+5, icol:icol+2] .= f_ψ2_u  ./ force_scaling
    jacob[irow_p2+3:irow_p2+5, icol+3:icol+5] .= f_ψ2_θ ./ force_scaling
    jacob[irow_p2+3:irow_p2+5, icol+6:icol+8] .= f_ψ2_F
    jacob[irow_p2+3:irow_p2+5, icol+9:icol+11] .= f_ψ2_M

    # compatability equation jacobian entries for the end of the beam
    # if irow_e2 == irow_p2 use row corresponding to compatability equations for this beam
    # if irow_e2 <= 0 use row corresponding to compatability equations for the other beam
    # otherwise use row corresponding to compatability equations for this beam
    irow = ifelse(irow_e2 == irow_p2 || irow_e2 <= 0, irow_p2 + 6, irow_e2)

    jacob[irow:irow+2, icol:icol+2] .= f_F2_u
    jacob[irow:irow+2, icol+3:icol+5] .= f_F2_θ
    jacob[irow:irow+2, icol+6:icol+8] .= f_F2_F .* force_scaling
    jacob[irow:irow+2, icol+9:icol+11] .= f_F2_M .* force_scaling

    jacob[irow+3:irow+5, icol+3:icol+5] .= f_M2_θ
    jacob[irow+3:irow+5, icol+6:icol+8] .= f_M2_F .* force_scaling
    jacob[irow+3:irow+5, icol+9:icol+11] .= f_M2_M .* force_scaling

    return jacob
end

"""
    initial_condition_insert_element_jacobian!(jacob, force_scaling,
        irow_e, irow_e1, irow_p1, irow_e2, irow_p2, icol,
        f_u1_Vdot, f_u2_Vdot, f_u1_Ωdot, f_u2_Ωdot, 
        f_u1_F, f_u2_F, f_u1_V, f_u2_V, f_u1_Ω, f_u2_Ω,
        f_ψ1_Vdot, f_ψ2_Vdot, f_ψ1_Ωdot, f_ψ2_Ωdot, 
        f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_V, f_ψ2_V, f_ψ1_Ω, f_ψ2_Ω,
        f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_F, f_M2_F, f_M1_M, f_M2_M,
        f_V_V, f_Ω_Ω)

Insert the the beam element jacobian entries into the jacobian matrix for an initial 
condition analysis

# Arguments
 - `jacob`: System jacobian matrix
 - `force_scaling`: Scaling parameter for forces
 - `irow_e1`: row index of the first residual equation for the start of the beam element
   (a value <= 0 indicates the equations have been eliminated from the system of equations)
 - `irow_p1`: row index of the first residual equation for the point at the start of the 
   beam element
 - `irow_e2`: row index of the first residual equation for the end of the beam element
   (a value <= 0 indicates the equations have been eliminated from the system of equations)
 - `irow_p2`: row index of the first residual equation for the point at the end of the 
   beam element
 - `icol`: Row/Column index corresponding to the first beam state variable

All other arguments use the following naming convention:
 - `f_y_x`: Jacobian of element equation `y` with respect to state variable `x`
"""
@inline function initial_condition_insert_element_jacobian!(jacob, force_scaling,
    irow_e, irow_e1, irow_p1, irow_e2, irow_p2, icol,
    f_u1_Vdot, f_u2_Vdot, f_u1_Ωdot, f_u2_Ωdot, 
    f_u1_F, f_u2_F, f_u1_V, f_u2_V, f_u1_Ω, f_u2_Ω,
    f_ψ1_Vdot, f_ψ2_Vdot, f_ψ1_Ωdot, f_ψ2_Ωdot, 
    f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_V, f_ψ2_V, f_ψ1_Ω, f_ψ2_Ω,
    f_F1_F, f_F2_F, f_F1_M, f_F2_M,
    f_M1_F, f_M2_F, f_M1_M, f_M2_M,
    f_V_V, f_Ω_Ω)

    # add equilibrium equation jacobian entries for the start of the beam
    jacob[irow_p1:irow_p1+2, icol:icol+2] .= f_u1_Vdot ./ force_scaling
    jacob[irow_p1:irow_p1+2, icol+3:icol+5] .= f_u1_Ωdot ./ force_scaling
    jacob[irow_p1:irow_p1+2, icol+6:icol+8] .= f_u1_F
    jacob[irow_p1:irow_p1+2, icol+12:icol+14] .= f_u1_V ./ force_scaling
    jacob[irow_p1:irow_p1+2, icol+15:icol+17] .= f_u1_Ω ./ force_scaling

    jacob[irow_p1+3:irow_p1+5, icol:icol+2] .= f_ψ1_Vdot ./ force_scaling
    jacob[irow_p1+3:irow_p1+5, icol+3:icol+5] .= f_ψ1_Ωdot ./ force_scaling
    jacob[irow_p1+3:irow_p1+5, icol+6:icol+8] .= f_ψ1_F
    jacob[irow_p1+3:irow_p1+5, icol+9:icol+11] .= f_ψ1_M
    jacob[irow_p1+3:irow_p1+5, icol+12:icol+14] .= f_ψ1_V ./ force_scaling
    jacob[irow_p1+3:irow_p1+5, icol+15:icol+17] .= f_ψ1_Ω ./ force_scaling

    # add compatability equation jacobian entries for the start of the beam
    # if irow_e1 == irow_p1 use row corresponding to compatability equations for this beam
    # if irow_e1 <= 0 use row corresponding to compatability equations for the other beam
    # otherwise use row corresponding to compatability equations for this beam
    irow = ifelse(irow_e1 == irow_p1 || irow_e1 <= 0, irow_p1+6, irow_e1)

    jacob[irow:irow+2, icol+6:icol+8] .= f_F1_F .* force_scaling
    jacob[irow:irow+2, icol+9:icol+11] .= f_F1_M .* force_scaling

    jacob[irow+3:irow+5, icol+6:icol+8] .= f_M1_F .* force_scaling
    jacob[irow+3:irow+5, icol+9:icol+11] .= f_M1_M .* force_scaling

    # add equilibrium equation jacobian entries for the end of the beam
    jacob[irow_p2:irow_p2+2, icol:icol+2] .= f_u2_Vdot ./ force_scaling
    jacob[irow_p2:irow_p2+2, icol+3:icol+5] .= f_u2_Ωdot ./ force_scaling
    jacob[irow_p2:irow_p2+2, icol+6:icol+8] .= f_u2_F
    jacob[irow_p2:irow_p2+2, icol+12:icol+14] .= f_u2_V ./ force_scaling
    jacob[irow_p2:irow_p2+2, icol+15:icol+17] .= f_u2_Ω ./ force_scaling

    jacob[irow_p2+3:irow_p2+5, icol:icol+2] .= f_ψ2_Vdot ./ force_scaling
    jacob[irow_p2+3:irow_p2+5, icol+3:icol+5] .= f_ψ2_Ωdot ./ force_scaling
    jacob[irow_p2+3:irow_p2+5, icol+6:icol+8] .= f_ψ2_F
    jacob[irow_p2+3:irow_p2+5, icol+9:icol+11] .= f_ψ2_M
    jacob[irow_p2+3:irow_p2+5, icol+12:icol+14] .= f_ψ2_V ./ force_scaling
    jacob[irow_p2+3:irow_p2+5, icol+15:icol+17] .= f_ψ2_Ω ./ force_scaling

    # add compatability equation jacobian entries for the end of the beam
    # if irow_e2 == irow_p2 use row corresponding to compatability equations for this beam
    # if irow_e2 <= 0 use row corresponding to compatability equations for the other beam
    # otherwise use row corresponding to compatability equations for this beam
    irow = ifelse(irow_e2 == irow_p2 || irow_e2 <= 0, irow_p2 + 6, irow_e2)

    jacob[irow:irow+2, icol+6:icol+8] .= f_F2_F .* force_scaling
    jacob[irow:irow+2, icol+9:icol+11] .= f_F2_M .* force_scaling

    jacob[irow+3:irow+5, icol+6:icol+8] .= f_M2_F .* force_scaling
    jacob[irow+3:irow+5, icol+9:icol+11] .= f_M2_M .* force_scaling

    # add beam residual equation jacobian entries
    jacob[irow_e:irow_e+2, icol+12:icol+14] .= f_V_V

    jacob[irow_e+3:irow_e+5, icol+15:icol+17] .= f_Ω_Ω

    return jacob
end

"""
    dynamic_insert_element_jacobian!(jacob, force_scaling, irow_e,
        irow_e1, irow_p1, irow_e2, irow_p2, icol,
        f_u1_u, f_u2_u, f_u1_θ, f_u2_θ, f_u1_F, f_u2_F, f_u1_V, f_u2_V, f_u1_Ω, f_u2_Ω,
        f_ψ1_u, f_ψ2_u, f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, 
        f_ψ1_V, f_ψ2_V, f_ψ1_Ω, f_ψ2_Ω,
        f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M,
        f_V_u, f_V_θ, f_V_V,
        f_Ω_θ, f_Ω_Ω)

Insert the the beam element jacobian entries into the jacobian matrix for a dynamic analysis
     
# Arguments
 - `jacob`: system jacobian matrix
 - `force_scaling`: scaling parameter for forces
 - `irow_e`: row index of the first linear/angular velocity residual for this element
 - `irow_e1`: row index of the first residual equation for the start of the beam element
    (a value <= 0 indicates the equations have been eliminated from the system of equations)
 - `irow_p1`: row index of the first residual equation for the point at the start of the 
    beam element
 - `irow_e2`: row index of the first residual equation for the end of the beam element
    (a value <= 0 indicates the equations have been eliminated from the system of equations)
 - `irow_p2`: row index of the first residual equation for the point at the end of the 
    beam element
 - `icol`: row/column index corresponding to the first beam state variable

All other arguments use the following naming convention:
 - `f_y_x`: jacobian of element equation `y` with respect to state variable `x`
"""
@inline function dynamic_insert_element_jacobian!(jacob, force_scaling, irow_e,
    irow_e1, irow_p1, irow_e2, irow_p2, icol,
    f_u1_u, f_u2_u, f_u1_θ, f_u2_θ, f_u1_F, f_u2_F, f_u1_V, f_u2_V, f_u1_Ω, f_u2_Ω,
    f_ψ1_u, f_ψ2_u, f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, 
    f_ψ1_V, f_ψ2_V, f_ψ1_Ω, f_ψ2_Ω,
    f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
    f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M,
    f_V_u, f_V_θ, f_V_V,
    f_Ω_θ, f_Ω_Ω)

    static_insert_element_jacobian!(jacob, force_scaling, irow_e1, irow_p1,
        irow_e2, irow_p2, icol,
        f_u1_u, f_u2_u, f_u1_θ, f_u2_θ, f_u1_F, f_u2_F,
        f_ψ1_u, f_ψ2_u, f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M,
        f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M)

    # add equilibrium equation jacobian entries for the start of the beam
    jacob[irow_p1:irow_p1+2, icol+12:icol+14] .= f_u1_V  ./ force_scaling
    jacob[irow_p1:irow_p1+2, icol+15:icol+17] .= f_u1_Ω  ./ force_scaling
    jacob[irow_p1+3:irow_p1+5, icol+12:icol+14] .= f_ψ1_V  ./ force_scaling
    jacob[irow_p1+3:irow_p1+5, icol+15:icol+17] .= f_ψ1_Ω  ./ force_scaling

    # add equilibrium equation jacobian entries for the end of the beam
    jacob[irow_p2:irow_p2+2, icol+12:icol+14] .= f_u2_V  ./ force_scaling
    jacob[irow_p2:irow_p2+2, icol+15:icol+17] .= f_u2_Ω  ./ force_scaling
    jacob[irow_p2+3:irow_p2+5, icol+12:icol+14] .= f_ψ2_V  ./ force_scaling
    jacob[irow_p2+3:irow_p2+5, icol+15:icol+17] .= f_ψ2_Ω  ./ force_scaling

    # add beam residual equation jacobian entries
    jacob[irow_e:irow_e+2, icol:icol+2] .= f_V_u
    jacob[irow_e:irow_e+2, icol+3:icol+5] .= f_V_θ
    jacob[irow_e:irow_e+2, icol+12:icol+14] .= f_V_V

    jacob[irow_e+3:irow_e+5, icol+3:icol+5] .= f_Ω_θ
    jacob[irow_e+3:irow_e+5, icol+15:icol+17] .= f_Ω_Ω

    return jacob
end

"""
    static_element_jacobian!(jacob, x, ielem, elem, distributed_loads, 
        point_masses, gvec, force_scaling, icol, irow_e1, irow_p1, irow_e2, irow_p2)

Adds a beam element's contributions to the system jacobian matrix for a static analysis.

# Arguments
 - `jacob`: System jacobian matrix
 - `x`: current state vector
 - `ielem`: beam element index
 - `elem`: beam element
 - `distributed_loads`: dictionary with all distributed loads
 - `point_masses`: dictionary with all point masses
 - `gvec`: gravity vector
 - `force_scaling`: scaling parameter for forces/moments
 - `icol`: starting index for the beam's state variables
 - `irow_e1`: row index of the first residual equation for the start of the beam element
    (a value <= 0 indicates the equations have been eliminated from the system of equations)
 - `irow_p1`: row index of the first residual equation for the point at the start of the 
    beam element
 - `irow_e2`: row index of the first residual equation for the end of the beam element
    (a value <= 0 indicates the equations have been eliminated from the system of equations)
 - `irow_p2`: row index of the first residual equation for the point at the end of the 
    beam element
"""
@inline function static_element_jacobian!(jacob, x, ielem, elem, distributed_loads, 
    point_masses, gvec, force_scaling, icol, irow_e1, irow_p1, irow_e2, irow_p2)

    # element properties
    ΔL = elem.L
    Cab = elem.Cab
    compliance = elem.compliance
    mass = elem.mass

    # scale mass matrix by the element length (to allow zero length elements)
    compliance *= ΔL
    mass *= ΔL

    # modify mass matrix to account for point masses, if present
    if haskey(point_masses, ielem)
        mass += transform_properties(point_masses[ielem].mass, Cab)
    end

    # element state variables
    u, θ, F, M = static_element_state_variables(x, icol, force_scaling)

    # compliance submatrices
    S11 = compliance[SVector{3}(1:3), SVector{3}(1:3)]
    S12 = compliance[SVector{3}(1:3), SVector{3}(4:6)]
    S21 = compliance[SVector{3}(4:6), SVector{3}(1:3)]
    S22 = compliance[SVector{3}(4:6), SVector{3}(4:6)]

    # mass submatrices
    mass11 = mass[SVector{3}(1:3), SVector{3}(1:3)]
    mass12 = mass[SVector{3}(1:3), SVector{3}(4:6)]
    mass21 = mass[SVector{3}(4:6), SVector{3}(1:3)]
    mass22 = mass[SVector{3}(4:6), SVector{3}(4:6)]

    # rotation matrices
    C = get_C(θ)
    CtCab = C'*Cab
    C_θ1, C_θ2, C_θ3 = get_C_θ(C, θ)
    Ct_θ1, Ct_θ2, Ct_θ3 = C_θ1', C_θ2', C_θ3'

    # element strain and curvature
    γ = S11*F + S12*M
    κ = S21*F + S22*M

    # effective linear and angular acceleration
    ae = -gvec
    αe = zero(ae)

    # element acceleration loads (including gravitational loads)
    f1_u, f2_u, m1_u, m2_u, f1_θ, f2_θ, m1_θ, m2_θ = acceleration_loads_jacobian(mass11, 
        mass12, mass21, mass22, ae, αe, Cab, CtCab, C_θ1, C_θ2, C_θ3, Ct_θ1, Ct_θ2, Ct_θ3)

    # distributed loads
    if haskey(distributed_loads, ielem)
        dload = distributed_loads[ielem]
        f1_θ += mul3(Ct_θ1, Ct_θ2, Ct_θ3, dload.f1_follower)
        f2_θ += mul3(Ct_θ1, Ct_θ2, Ct_θ3, dload.f2_follower)
        m1_θ += mul3(Ct_θ1, Ct_θ2, Ct_θ3, dload.m1_follower)
        m2_θ += mul3(Ct_θ1, Ct_θ2, Ct_θ3, dload.m2_follower)
    end

    # element resultant jacobians
    f_u1_u, f_u2_u, f_u1_θ, f_u2_θ, f_u1_F, f_u2_F,
        f_ψ1_u, f_ψ2_u, f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M,
        f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M = static_element_jacobian_equations(
        ΔL, S11, S12, S21, S22, Cab, CtCab, θ, F, M, γ, κ, f1_u, f2_u, m1_u, m2_u, 
        f1_θ, f2_θ, m1_θ, m2_θ, Ct_θ1, Ct_θ2, Ct_θ3)

    # insert element resultants into the jacobian matrix
    static_insert_element_jacobian!(jacob, force_scaling, irow_e1,
        irow_p1, irow_e2, irow_p2, icol,
        f_u1_u, f_u2_u, f_u1_θ, f_u2_θ, f_u1_F, f_u2_F,
        f_ψ1_u, f_ψ2_u, f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M,
        f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M)

    return jacob
end

"""
    steady_state_element_jacobian!(jacob, x, ielem, elem, distributed_loads, point_masses, 
        gvec, force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2, x0, v0, ω0, 
        a0, α0)

Adds a beam element's contributions to the system jacobian matrix for a steady state 
analysis.

# Arguments
 - `jacob`: System jacobian matrix
 - `x`: current state vector
 - `ielem`: beam element index
 - `elem`: beam element
 - `distributed_loads`: dictionary with all distributed loads
 - `point_masses`: dictionary with all point masses
 - `gvec`: gravity vector
 - `force_scaling`: scaling parameter for forces/moments
 - `icol`: starting index for the beam's state variables
 - `irow_e1`: row index of the first residual equation for the start of the beam element
    (a value <= 0 indicates the equations have been eliminated from the system of equations)
 - `irow_p1`: row index of the first residual equation for the point at the start of the 
    beam element
 - `irow_e2`: row index of the first residual equation for the end of the beam element
    (a value <= 0 indicates the equations have been eliminated from the system of equations)
 - `irow_p2`: row index of the first residual equation for the point at the end of the 
    beam element
 - `x0`: Body frame origin (for the current time step)
 - `v0`: Body frame linear velocity (for the current time step)
 - `ω0`: Body frame angular velocity (for the current time step)
 - `a0`: Body frame linear acceleration (for the current time step)
 - `α0`: Body frame angular acceleration (for the current time step)
 """
@inline function steady_state_element_jacobian!(jacob, x, ielem, elem, distributed_loads, 
    point_masses, gvec, force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2, 
    x0, v0, ω0, a0, α0)

    # element properties
    ΔL = elem.L
    Cab = elem.Cab
    compliance = elem.compliance
    mass = elem.mass

    # scale mass matrix by the element length (to allow zero length elements)
    compliance *= ΔL
    mass *= ΔL

    # modify mass matrix to account for point masses, if present
    if haskey(point_masses, ielem)
        mass += transform_properties(point_masses[ielem].mass, Cab)
    end

    # element state variables
    u, θ, F, M, V, Ω = dynamic_element_state_variables(x, icol, force_scaling)

    # compliance submatrices
    S11 = compliance[SVector{3}(1:3), SVector{3}(1:3)]
    S12 = compliance[SVector{3}(1:3), SVector{3}(4:6)]
    S21 = compliance[SVector{3}(4:6), SVector{3}(1:3)]
    S22 = compliance[SVector{3}(4:6), SVector{3}(4:6)]

    # mass submatrices
    mass11 = mass[SVector{3}(1:3), SVector{3}(1:3)]
    mass12 = mass[SVector{3}(1:3), SVector{3}(4:6)]
    mass21 = mass[SVector{3}(4:6), SVector{3}(1:3)]
    mass22 = mass[SVector{3}(4:6), SVector{3}(4:6)]

    # rotation matrices
    C = get_C(θ)
    Ct = C'
    CtCab = Ct*Cab
    C_θ1, C_θ2, C_θ3 = get_C_θ(C, θ)
    Ct_θ1, Ct_θ2, Ct_θ3 = C_θ1', C_θ2', C_θ3'

    # element strain and curvature
    γ = S11*F + S12*M
    κ = S21*F + S22*M

    # element linear and angular momentum
    P = mass11*V + mass12*Ω
    H = mass21*V + mass22*Ω

    # body frame velocities
    v = v0 + cross(ω0, elem.x - x0)
    ω = ω0

    # body frame accelerations
    a = a0 + cross(α0, elem.x - x0)
    α = α0 

    # effective linear and angular acceleration
    ae = a - gvec + cross(α, u)
    αe = α

    # element acceleration loads (including gravitational loads)
    f1_u, f2_u, m1_u, m2_u, f1_θ, f2_θ, m1_θ, m2_θ = acceleration_loads_jacobian(mass11, 
        mass12, mass21, mass22, ae, αe, Cab, CtCab, C_θ1, C_θ2, C_θ3, Ct_θ1, Ct_θ2, Ct_θ3)

    # distributed loads
    if haskey(distributed_loads, ielem)
        dload = distributed_loads[ielem]
        f1_θ += mul3(Ct_θ1, Ct_θ2, Ct_θ3, dload.f1_follower)
        f2_θ += mul3(Ct_θ1, Ct_θ2, Ct_θ3, dload.f2_follower)
        m1_θ += mul3(Ct_θ1, Ct_θ2, Ct_θ3, dload.m1_follower)
        m2_θ += mul3(Ct_θ1, Ct_θ2, Ct_θ3, dload.m2_follower)
    end

    # solve for the element resultant jacobians
    f_u1_u, f_u2_u, f_u1_θ, f_u2_θ, f_u1_F, f_u2_F, f_u1_V, f_u2_V, f_u1_Ω, f_u2_Ω,
        f_ψ1_u, f_ψ2_u, f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, 
        f_ψ1_V, f_ψ2_V, f_ψ1_Ω, f_ψ2_Ω,
        f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M,
        f_V_u, f_V_θ, f_V_V,
        f_Ω_θ, f_Ω_Ω = steady_state_element_jacobian_equations(ΔL, S11, S12, S21, S22, 
        mass11, mass12, mass21, mass22, Cab, CtCab, θ, F, M, γ, κ, V, P, H, ω, f1_u, f2_u, 
        m1_u, m2_u, f1_θ, f2_θ, m1_θ, m2_θ, C_θ1, C_θ2, C_θ3, Ct_θ1, Ct_θ2, Ct_θ3)

    # insert element resultants into the jacobian matrix
    dynamic_insert_element_jacobian!(jacob, force_scaling, irow_e,
        irow_e1, irow_p1, irow_e2, irow_p2, icol,
        f_u1_u, f_u2_u, f_u1_θ, f_u2_θ, f_u1_F, f_u2_F, f_u1_V, f_u2_V, f_u1_Ω, f_u2_Ω,
        f_ψ1_u, f_ψ2_u, f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, 
        f_ψ1_V, f_ψ2_V, f_ψ1_Ω, f_ψ2_Ω,
        f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M,
        f_V_u, f_V_θ, f_V_V,
        f_Ω_θ, f_Ω_Ω)

    return jacob
end

"""
    initial_condition_element_jacobian!(jacob, x, ielem, elem, distributed_loads, 
        point_masses, gvec, force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, 
        irow_p2, x0, v0, ω0, a0, α0, u0, θ0, udot0, θdot0, Fdot0, Mdot0)

Adds a beam element's contributions to the system jacobian matrix for an initial conditions 
analysis.

# Arguments
 - `jacob`: System jacobian matrix
 - `x`: current state vector
 - `ielem`: beam element index
 - `elem`: beam element
 - `distributed_loads`: dictionary with all distributed loads
 - `point_masses`: dictionary with all point masses
 - `gvec`: gravity vector
 - `force_scaling`: scaling parameter for forces/moments
 - `icol`: starting index for the beam's state variables
 - `irow_e1`: row index of the first residual equation for the start of the beam element
    (a value <= 0 indicates the equations have been eliminated from the system of equations)
 - `irow_p1`: row index of the first residual equation for the point at the start of the 
    beam element
 - `irow_e2`: row index of the first residual equation for the end of the beam element
    (a value <= 0 indicates the equations have been eliminated from the system of equations)
 - `irow_p2`: row index of the first residual equation for the point at the end of the 
    beam element
 - `x0`: Body frame origin (for the current time step)
 - `v0`: Body frame linear velocity (for the current time step)
 - `ω0`: Body frame angular velocity (for the current time step)
 - `a0`: Body frame linear acceleration (for the current time step)
 - `α0`: Body frame angular acceleration (for the current time step)
 - `u`: initial deflection variables for the beam element
 - `θ`: initial rotation variables for the beam element
 - `udot`: initial time derivative of u for the beam element
 - `θdot`: initial time derivative of θ for the beam element
 - `Fdot`:
 - `Mdot`:
"""
@inline function initial_condition_element_jacobian!(jacob, x, ielem, elem, 
    distributed_loads, point_masses, gvec, force_scaling, icol, irow_e, irow_e1, irow_p1, 
    irow_e2, irow_p2, x0, v0, ω0, a0, α0, u, θ, udot, θdot, Fdot, Mdot)

    # element properties
    ΔL = elem.L
    Cab = elem.Cab
    compliance = elem.compliance
    mass = elem.mass
    μ = elem.mu

    # scale mass matrix by the element length (to allow zero length elements)
    compliance *= ΔL
    mass *= ΔL

    # modify mass matrix to account for point masses, if present
    if haskey(point_masses, ielem)
        mass += transform_properties(point_masses[ielem].mass, Cab)
    end

    # element state variables
    Vdot, Ωdot, F, M, V, Ω = dynamic_element_state_variables(x, icol, force_scaling)

    # compliance submatrices
    S11 = compliance[SVector{3}(1:3), SVector{3}(1:3)]
    S12 = compliance[SVector{3}(1:3), SVector{3}(4:6)]
    S21 = compliance[SVector{3}(4:6), SVector{3}(1:3)]
    S22 = compliance[SVector{3}(4:6), SVector{3}(4:6)]

    # mass submatrices
    mass11 = mass[SVector{3}(1:3), SVector{3}(1:3)]
    mass12 = mass[SVector{3}(1:3), SVector{3}(4:6)]
    mass21 = mass[SVector{3}(4:6), SVector{3}(1:3)]
    mass22 = mass[SVector{3}(4:6), SVector{3}(4:6)]

    # rotation matrices
    C = get_C(θ)
    Ct = C'
    CtCab = Ct*Cab

    # rotation matrix time derivatives
    Cdot = get_C_t(C, θ, θdot)
    Ctdot = Cdot'
    CtCabdot = Ctdot*Cab

    # element strain and curvature
    γ = S11*F + S12*M
    κ = S21*F + S22*M

    # element linear and angular momentum
    P = mass11*V + mass12*Ω
    H = mass21*V + mass22*Ω

    # body frame velocities
    v = v0 + cross(ω0, elem.x - x0)
    ω = ω0

    # damping loads
    Fd = SVector(μ[1]*Fdot[1], μ[2]*Fdot[2], μ[3]*Fdot[3])
    Md = SVector(μ[4]*Mdot[1], μ[5]*Mdot[2], μ[6]*Mdot[3])

    # total element loads
    F = F + Fd
    M = M + Md

    # solve for the element resultant jacobians
    f_u1_Vdot, f_u2_Vdot, f_u1_Ωdot, f_u2_Ωdot, f_u1_F, f_u2_F, f_u1_V, f_u2_V, f_u1_Ω, f_u2_Ω,
        f_ψ1_Vdot, f_ψ2_Vdot, f_ψ1_Ωdot, f_ψ2_Ωdot, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_V, f_ψ2_V, f_ψ1_Ω, f_ψ2_Ω,
        f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_F, f_M2_F, f_M1_M, f_M2_M,
        f_V_V, f_Ω_Ω = initial_condition_element_jacobian_equations(ΔL, S11, S12, S21, S22, 
        mass11, mass12, mass21, mass22, Cab, CtCab, CtCabdot, θ, F, γ, V, P, ω)

    # insert element resultants into the jacobian matrix
    initial_condition_insert_element_jacobian!(jacob, force_scaling,
        irow_e, irow_e1, irow_p1, irow_e2, irow_p2, icol,
        f_u1_Vdot, f_u2_Vdot, f_u1_Ωdot, f_u2_Ωdot, 
        f_u1_F, f_u2_F, f_u1_V, f_u2_V, f_u1_Ω, f_u2_Ω,
        f_ψ1_Vdot, f_ψ2_Vdot, f_ψ1_Ωdot, f_ψ2_Ωdot, 
        f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_V, f_ψ2_V, f_ψ1_Ω, f_ψ2_Ω,
        f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_F, f_M2_F, f_M1_M, f_M2_M,
        f_V_V, f_Ω_Ω)

    return jacob
end

"""
    newmark_element_jacobian!(jacob, x, ielem, elem, distributed_loads, point_masses, gvec,
        force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2, x0, v0, ω0, a0, α0,
        udot_init, θdot_init, Fdot_init, Mdot_init, Vdot_init, Ωdot_init, dt)

Adds a beam element's contributions to the system jacobian matrix for a Newmark scheme 
time-marching simulation.

# Arguments
 - `jacob`: System jacobian matrix
 - `x`: current state vector
 - `ielem`: beam element index
 - `elem`: beam element
 - `distributed_loads`: dictionary with all distributed loads
 - `point_masses`: dictionary with all point masses
 - `gvec`: gravity vector
 - `force_scaling`: scaling parameter for forces/moments
 - `icol`: starting index for the beam's state variables
 - `irow_e1`: row index of the first residual equation for the start of the beam element
    (a value <= 0 indicates the equations have been eliminated from the system of equations)
 - `irow_p1`: row index of the first residual equation for the point at the start of the 
    beam element
 - `irow_e2`: row index of the first residual equation for the end of the beam element
    (a value <= 0 indicates the equations have been eliminated from the system of equations)
 - `irow_p2`: row index of the first residual equation for the point at the end of the 
    beam element
 - `x0`: Body frame origin (for the current time step)
 - `v0`: Body frame linear velocity (for the current time step)
 - `ω0`: Body frame angular velocity (for the current time step)
 - `a0`: Body frame linear acceleration (for the current time step)
 - `α0`: Body frame angular acceleration (for the current time step)
 - `udot_init`: `2/dt*u + udot` for the beam element from the previous time step
 - `θdot_init`: `2/dt*θ + θdot` for the beam element from the previous time step
 - `Fdot_init`: `2/dt*F + Fdot` for the beam element from the previous time step
 - `Mdot_init`: `2/dt*M + Mdot` for the beam element from the previous time step
 - `Vdot_init`: `2/dt*V + Vdot` for the beam element from the previous time step
 - `Ωdot_init`: `2/dt*Ω + Ωdot` for the beam element from the previous time step
 - `dt`: time step size
"""
@inline function newmark_element_jacobian!(jacob, x, ielem, elem, distributed_loads, 
    point_masses, gvec, force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2, 
    x0, v0, ω0, a0, α0, udot_init, θdot_init, Fdot_init, Mdot_init, Vdot_init, Ωdot_init, dt)

    # element properties
    ΔL = elem.L
    Cab = elem.Cab
    compliance = elem.compliance
    mass = elem.mass
    μ = elem.mu

    # scale mass matrix by the element length (to allow zero length elements)
    compliance *= ΔL
    mass *= ΔL

    # modify mass matrix to account for point masses, if present
    if haskey(point_masses, ielem)
        mass += transform_properties(point_masses[ielem].mass, Cab)
    end

    # element state variables
    u, θ, F, M, V, Ω = dynamic_element_state_variables(x, icol, force_scaling)
    
    # element state variable rates
    udot = 2/dt*u - udot_init
    θdot = 2/dt*θ - θdot_init
    Fdot = 2/dt*F - Fdot_init
    Mdot = 2/dt*M - Mdot_init
    Vdot = 2/dt*V - Vdot_init
    Ωdot = 2/dt*Ω - Ωdot_init

    # compliance submatrices
    S11 = compliance[SVector{3}(1:3), SVector{3}(1:3)]
    S12 = compliance[SVector{3}(1:3), SVector{3}(4:6)]
    S21 = compliance[SVector{3}(4:6), SVector{3}(1:3)]
    S22 = compliance[SVector{3}(4:6), SVector{3}(4:6)]

    # mass submatrices
    mass11 = mass[SVector{3}(1:3), SVector{3}(1:3)]
    mass12 = mass[SVector{3}(1:3), SVector{3}(4:6)]
    mass21 = mass[SVector{3}(4:6), SVector{3}(1:3)]
    mass22 = mass[SVector{3}(4:6), SVector{3}(4:6)]

    # damping coefficient submatrices
    μ11 = @SMatrix [μ[1] 0 0; 0 μ[2] 0; 0 0 μ[3]]
    μ22 = @SMatrix [μ[4] 0 0; 0 μ[5] 0; 0 0 μ[6]]

    # rotation matrices
    C = get_C(θ)
    Ct = C'
    CtCab = Ct*Cab

    # rotation matrix time derivatives
    Cdot = get_C_t(C, θ, θdot)
    Ctdot = Cdot'
    CtCabdot = Ctdot*Cab

    # rotation matrix derivatives
    C_θ1, C_θ2, C_θ3 = get_C_θ(C, θ)
    Ct_θ1, Ct_θ2, Ct_θ3 = C_θ1', C_θ2', C_θ3'
    Cdot_θ1, Cdot_θ2, Cdot_θ3 = get_C_t_θ(θ, θdot)
    Ctdot_θ1, Ctdot_θ2, Ctdot_θ3 = Cdot_θ1', Cdot_θ2', Cdot_θ3'
    Cdot_θdot1, Cdot_θdot2, Cdot_θdot3 = get_C_t_θdot(C, θ)
    Ctdot_θdot1, Ctdot_θdot2, Ctdot_θdot3 = Cdot_θdot1', Cdot_θdot2', Cdot_θdot3'

    # element strain and curvature
    γ = S11*F + S12*M
    κ = S21*F + S22*M

    # element linear and angular momentum
    P = mass11*V + mass12*Ω
    H = mass21*V + mass22*Ω

    # element linear and angular momentum rates
    Pdot = mass11*Vdot + mass12*Ωdot
    Hdot = mass21*Vdot + mass22*Ωdot

    # linear and angular momentum (in the body frame)
    CtCabPdot = CtCabdot*P + CtCab*Pdot
    CtCabHdot = CtCabdot*H + CtCab*Hdot

    # body frame velocities
    v = v0 + cross(ω0, elem.x - x0)
    ω = ω0

    # body frame accelerations
    a = a0 + cross(α0, elem.x - x0)
    α = α0 

    # effective linear and angular acceleration
    ae = a - gvec + cross(α, u)
    αe = α

    # element acceleration loads (including gravitational loads)
    f1_u, f2_u, m1_u, m2_u, f1_θ, f2_θ, m1_θ, m2_θ = acceleration_loads_jacobian(mass11, 
        mass12, mass21, mass22, ae, αe, Cab, CtCab, C_θ1, C_θ2, C_θ3, Ct_θ1, Ct_θ2, Ct_θ3)

    # distributed loads
    if haskey(distributed_loads, ielem)
        dload = distributed_loads[ielem]
        f1_θ += mul3(Ct_θ1, Ct_θ2, Ct_θ3, dload.f1_follower)
        f2_θ += mul3(Ct_θ1, Ct_θ2, Ct_θ3, dload.f2_follower)
        m1_θ += mul3(Ct_θ1, Ct_θ2, Ct_θ3, dload.m1_follower)
        m2_θ += mul3(Ct_θ1, Ct_θ2, Ct_θ3, dload.m2_follower)
    end

    # add damping loads
    Fd = SVector(μ[1]*Fdot[1], μ[2]*Fdot[2], μ[3]*Fdot[3])
    Md = SVector(μ[4]*Mdot[1], μ[5]*Mdot[2], μ[6]*Mdot[3])    

    # total element loads
    F = F + Fd
    M = M + Md

    # element resultant jacobians
    f_u1_u, f_u2_u, f_u1_θ, f_u2_θ, f_u1_F, f_u2_F, f_u1_V, f_u2_V, f_u1_Ω, f_u2_Ω,
        f_ψ1_u, f_ψ2_u, f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_V, f_ψ2_V, f_ψ1_Ω, f_ψ2_Ω,
        f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M,
        f_V_u, f_V_θ, f_V_V,
        f_Ω_θ, f_Ω_Ω = newmark_element_jacobian_equations(ΔL, S11, S12, S21, S22, 
        mass11, mass12, mass21, mass22, μ11, μ22, Cab, CtCab, CtCabdot, θ, F, M, γ, κ, V, P, H, 
        θdot, Pdot, Hdot, ω, dt, f1_u, f2_u, m1_u, m2_u, f1_θ, f2_θ, m1_θ, m2_θ, 
        C_θ1, C_θ2, C_θ3, Ct_θ1, Ct_θ2, Ct_θ3, Ctdot_θ1, Ctdot_θ2, Ctdot_θ3, 
        Ctdot_θdot1, Ctdot_θdot2, Ctdot_θdot3)

    # insert element resultant jacobians into the jacobian matrix
    dynamic_insert_element_jacobian!(jacob, force_scaling, irow_e,
        irow_e1, irow_p1, irow_e2, irow_p2, icol,
        f_u1_u, f_u2_u, f_u1_θ, f_u2_θ, f_u1_F, f_u2_F, f_u1_V, f_u2_V, f_u1_Ω, f_u2_Ω,
        f_ψ1_u, f_ψ2_u, f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_V, f_ψ2_V, f_ψ1_Ω, f_ψ2_Ω,
        f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M,
        f_V_u, f_V_θ, f_V_V,
        f_Ω_θ, f_Ω_Ω)

    return jacob
end

"""
    dynamic_element_jacobian!(jacob, x, ielem, elem, distributed_loads, point_masses, gvec,
        force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2, x0, v0, ω0, a0, α0,
        udot, θdot, Fdot, Mdot, Vdot, Ωdot)

Adds a beam element's contributions to the system jacobian matrix for a general dynamic 
simulation.

# Arguments
 - `jacob`: System jacobian matrix
 - `x`: current state vector
 - `ielem`: beam element index
 - `elem`: beam element
 - `distributed_loads`: dictionary with all distributed loads
 - `point_masses`: dictionary with all point masses
 - `gvec`: gravity vector
 - `force_scaling`: scaling parameter for forces/moments
 - `icol`: starting index for the beam's state variables
 - `irow_e1`: row index of the first residual equation for the start of the beam element
    (a value <= 0 indicates the equations have been eliminated from the system of equations)
 - `irow_p1`: row index of the first residual equation for the point at the start of the 
    beam element
 - `irow_e2`: row index of the first residual equation for the end of the beam element
    (a value <= 0 indicates the equations have been eliminated from the system of equations)
 - `irow_p2`: row index of the first residual equation for the point at the end of the 
    beam element
 - `x0`: Body frame origin (for the current time step)
 - `v0`: Body frame linear velocity (for the current time step)
 - `ω0`: Body frame angular velocity (for the current time step)
 - `a0`: Body frame linear acceleration (for the current time step)
 - `α0`: Body frame angular acceleration (for the current time step)
 - `udot_init`: `2/dt*u + udot` for the beam element from the previous time step
 - `θdot_init`: `2/dt*θ + θdot` for the beam element from the previous time step
 - `Fdot_init`: `2/dt*F + Fdot` for the beam element from the previous time step
 - `Mdot_init`: `2/dt*M + Mdot` for the beam element from the previous time step
 - `Vdot_init`: `2/dt*V + Vdot` for the beam element from the previous time step
 - `Ωdot_init`: `2/dt*Ω + Ωdot` for the beam element from the previous time step
 - `dt`: time step size
"""
@inline function dynamic_element_jacobian!(jacob, x, ielem, elem, distributed_loads, point_masses, gvec,
    force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2, x0, v0, ω0, a0, α0, 
    udot, θdot, Fdot, Mdot, Vdot, Ωdot)

    # element properties
    ΔL = elem.L
    Cab = elem.Cab
    compliance = elem.compliance
    mass = elem.mass
    μ = elem.mu

    # scale mass matrix by the element length (to allow zero length elements)
    compliance *= ΔL
    mass *= ΔL

    # modify mass matrix to account for point masses, if present
    if haskey(point_masses, ielem)
        mass += transform_properties(point_masses[ielem].mass, Cab)
    end

    # element state variables
    u, θ, F, M, V, Ω = dynamic_element_state_variables(x, icol, force_scaling)

    # compliance submatrices
    S11 = compliance[SVector{3}(1:3), SVector{3}(1:3)]
    S12 = compliance[SVector{3}(1:3), SVector{3}(4:6)]
    S21 = compliance[SVector{3}(4:6), SVector{3}(1:3)]
    S22 = compliance[SVector{3}(4:6), SVector{3}(4:6)]

    # mass submatrices
    mass11 = mass[SVector{3}(1:3), SVector{3}(1:3)]
    mass12 = mass[SVector{3}(1:3), SVector{3}(4:6)]
    mass21 = mass[SVector{3}(4:6), SVector{3}(1:3)]
    mass22 = mass[SVector{3}(4:6), SVector{3}(4:6)]

    # rotation matrices
    C = get_C(θ)
    Ct = C'
    CtCab = Ct*Cab

    # rotation matrix time derivatives
    Cdot = get_C_t(C, θ, θdot)
    Ctdot = Cdot'
    CtCabdot = Ctdot*Cab

    # rotation matrix derivatives
    C_θ1, C_θ2, C_θ3 = get_C_θ(C, θ)
    Ct_θ1, Ct_θ2, Ct_θ3 = C_θ1', C_θ2', C_θ3'
    Cdot_θ1, Cdot_θ2, Cdot_θ3 = get_C_t_θ(θ, θdot)
    Ctdot_θ1, Ctdot_θ2, Ctdot_θ3 = Cdot_θ1', Cdot_θ2', Cdot_θ3'

    # element strain and curvature
    γ = S11*F + S12*M
    κ = S21*F + S22*M

    # element linear and angular momentum
    P = mass11*V + mass12*Ω
    H = mass21*V + mass22*Ω

    # element linear and angular momentum rates
    Pdot = mass11*Vdot + mass12*Ωdot
    Hdot = mass21*Vdot + mass22*Ωdot

    # linear and angular momentum (in the body frame)
    CtCabPdot = CtCabdot*P + CtCab*Pdot
    CtCabHdot = CtCabdot*H + CtCab*Hdot

    # body frame velocities
    v = v0 + cross(ω0, elem.x - x0)
    ω = ω0

    # body frame accelerations
    a = a0 + cross(α0, elem.x - x0)
    α = α0 

    # effective linear and angular acceleration
    ae = a - gvec + cross(α, u)
    αe = α

    # element acceleration loads (including gravitational loads)
    f1_u, f2_u, m1_u, m2_u, f1_θ, f2_θ, m1_θ, m2_θ = acceleration_loads_jacobian(mass11, 
        mass12, mass21, mass22, ae, αe, Cab, CtCab, C_θ1, C_θ2, C_θ3, Ct_θ1, Ct_θ2, Ct_θ3)

    # distributed loads
    if haskey(distributed_loads, ielem)
        dload = distributed_loads[ielem]
        f1_θ += mul3(Ct_θ1, Ct_θ2, Ct_θ3, dload.f1_follower)
        f2_θ += mul3(Ct_θ1, Ct_θ2, Ct_θ3, dload.f2_follower)
        m1_θ += mul3(Ct_θ1, Ct_θ2, Ct_θ3, dload.m1_follower)
        m2_θ += mul3(Ct_θ1, Ct_θ2, Ct_θ3, dload.m2_follower)
    end

    # damping loads
    Fd = SVector(μ[1]*Fdot[1], μ[2]*Fdot[2], μ[3]*Fdot[3])
    Md = SVector(μ[4]*Mdot[1], μ[5]*Mdot[2], μ[6]*Mdot[3])        

    # total element loads
    F = F + Fd
    M = M + Md

    # solve for the element resultants
    f_u1_u, f_u2_u, f_u1_θ, f_u2_θ, f_u1_F, f_u2_F, f_u1_V, f_u2_V, f_u1_Ω, f_u2_Ω,
        f_ψ1_u, f_ψ2_u, f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_V, f_ψ2_V, f_ψ1_Ω, f_ψ2_Ω,
        f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M,
        f_V_u, f_V_θ, f_V_V,
        f_Ω_θ, f_Ω_Ω = dynamic_element_jacobian_equations(ΔL, S11, S12, S21, S22, mass11, 
        mass12, mass21, mass22, Cab, CtCab, CtCabdot, θ, F, M, γ, κ, V, P, H, θdot, Pdot, 
        Hdot, ω, f1_u, f2_u, m1_u, m2_u, f1_θ, f2_θ, m1_θ, m2_θ, C_θ1, C_θ2, C_θ3, Ct_θ1, Ct_θ2, Ct_θ3, 
        Ctdot_θ1, Ctdot_θ2, Ctdot_θ3)

    # insert element resultants into the jacobian matrix
    dynamic_insert_element_jacobian!(jacob, force_scaling, irow_e,
        irow_e1, irow_p1, irow_e2, irow_p2, icol,
        f_u1_u, f_u2_u, f_u1_θ, f_u2_θ, f_u1_F, f_u2_F, f_u1_V, f_u2_V, f_u1_Ω, f_u2_Ω,
        f_ψ1_u, f_ψ2_u, f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_V, f_ψ2_V, f_ψ1_Ω, f_ψ2_Ω,
        f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M,
        f_V_u, f_V_θ, f_V_V,
        f_Ω_θ, f_Ω_Ω)

    return jacob
end

"""
    element_mass_matrix_equations(ΔL, mass11, mass12, mass21, mass22, μ11, μ22, 
    Cab, CtCab, θ, P, H, γ, Ctdot_θdot1, Ctdot_θdot2, Ctdot_θdot3)

Calculate the derivatives of the element resultants with respect to the state rates.

# Arguments:
 - `ΔL`: length
 - `mass11, mass12, mass21, mass22`: beam element mass matrix, divided into submatrices
 - `Cab`: transformation matrix from the undeformed beam element frame to the body frame
 - `CtCab`: transformation matrix from the deformed beam element frame to the body frame
 - `θ`: angular displacement
 - `P`: linear momentum
 - `H`: angular momentum
 - `γ`: strain
 - `Ctdot_θdot1`: Derivative of `Cdot'` w.r.t. `θdot[1]`
 - `Ctdot_θdot2`: Derivative of `Cdot'` w.r.t. `θdot[2]`
 - `Ctdot_θdot3`: Derivative of `Cdot'` w.r.t. `θdot[3]`
"""
@inline function element_mass_matrix_equations(ΔL, mass11, mass12, mass21, mass22, μ11, μ22, 
    Cab, CtCab, θ, P, H, γ, Ctdot_θdot1, Ctdot_θdot2, Ctdot_θdot3)

    tmp = 1/2*mul3(Ctdot_θdot1, Ctdot_θdot2, Ctdot_θdot3, Cab*P)
    f_u1_θdot = tmp
    f_u2_θdot = tmp

    tmp = CtCab*μ11
    f_u1_Fdot = -tmp
    f_u2_Fdot = tmp

    tmp = 1/2*CtCab*mass11
    f_u1_Vdot = tmp
    f_u2_Vdot = tmp

    tmp = 1/2*CtCab*mass12
    f_u1_Ωdot = tmp
    f_u2_Ωdot = tmp

    tmp = 1/2*mul3(Ctdot_θdot1, Ctdot_θdot2, Ctdot_θdot3, Cab*H)
    f_ψ1_θdot = tmp
    f_ψ2_θdot = tmp

    tmp = 1/2*CtCab*tilde(ΔL*e1 + γ)*μ11
    f_ψ1_Fdot = -tmp
    f_ψ2_Fdot = -tmp

    tmp = CtCab*μ22
    f_ψ1_Mdot = -tmp
    f_ψ2_Mdot = tmp

    tmp = 1/2*CtCab*mass21
    f_ψ1_Vdot = tmp
    f_ψ2_Vdot = tmp

    tmp = 1/2*CtCab*mass22
    f_ψ1_Ωdot = tmp
    f_ψ2_Ωdot = tmp

    f_V_udot = -I3

    Q = get_Q(θ)
    f_Ω_θdot = -Cab'*Q

    return f_u1_θdot, f_u2_θdot, f_u1_Fdot, f_u2_Fdot, f_u1_Vdot, f_u2_Vdot, 
        f_u1_Ωdot, f_u2_Ωdot, f_ψ1_θdot, f_ψ2_θdot, f_ψ1_Fdot, f_ψ2_Fdot, 
        f_ψ1_Mdot, f_ψ2_Mdot, f_ψ1_Vdot, f_ψ2_Vdot, f_ψ1_Ωdot, f_ψ2_Ωdot,
        f_V_udot, f_Ω_θdot
end

"""
    insert_element_mass_matrix!(jacob, force_scaling, irow_e, irow_p1, 
        irow_p2, icol, 
        f_u1_θdot, f_u2_θdot, f_u1_Fdot, f_u2_Fdot, f_u1_Vdot, f_u2_Vdot, 
        f_u1_Ωdot, f_u2_Ωdot, f_ψ1_θdot, f_ψ2_θdot, f_ψ1_Fdot, f_ψ2_Fdot, 
        f_ψ1_Mdot, f_ψ2_Mdot, f_ψ1_Vdot, f_ψ2_Vdot, f_ψ1_Ωdot, f_ψ2_Ωdot,
        f_V_udot, f_Ω_θdot)

Insert the beam element's contributions into the "mass matrix": the jacobian of the
residual equations with respect to the state variable rates

# Arguments
 - `jacob`: System mass matrix
 - `force_scaling`: Scaling parameter for forces/moments
 - `irow_e`: row index of the first linear/angular velocity residual for this element
 - `irow_p1`: row index of the first residual equation for the point at the start of the 
    beam element
 - `irow_p2`: row index of the first residual equation for the point at the end of the 
    beam element
 - `icol`: Row/Column index corresponding to the first beam state variable

All other arguments use the following naming convention:
 - `f_y_x`: Jacobian of element equation "y" with respect to state variable "x"
"""
@inline function insert_element_mass_matrix!(jacob, force_scaling, irow_e, irow_p1, 
    irow_p2, icol, 
    f_u1_θdot, f_u2_θdot, f_u1_Fdot, f_u2_Fdot, f_u1_Vdot, f_u2_Vdot, 
    f_u1_Ωdot, f_u2_Ωdot, f_ψ1_θdot, f_ψ2_θdot, f_ψ1_Fdot, f_ψ2_Fdot, 
    f_ψ1_Mdot, f_ψ2_Mdot, f_ψ1_Vdot, f_ψ2_Vdot, f_ψ1_Ωdot, f_ψ2_Ωdot,
    f_V_udot, f_Ω_θdot)

    # create jacobian entries for the start of the beam
    jacob[irow_p1:irow_p1+2, icol+3:icol+5] .= f_u1_θdot ./ force_scaling
    jacob[irow_p1:irow_p1+2, icol+6:icol+8] .= f_u1_Fdot
    jacob[irow_p1:irow_p1+2, icol+12:icol+14] .= f_u1_Vdot  ./ force_scaling
    jacob[irow_p1:irow_p1+2, icol+15:icol+17] .= f_u1_Ωdot  ./ force_scaling
    jacob[irow_p1+3:irow_p1+5, icol+3:icol+5] .= f_ψ1_θdot ./ force_scaling
    jacob[irow_p1+3:irow_p1+5, icol+6:icol+8] .= f_ψ1_Fdot
    jacob[irow_p1+3:irow_p1+5, icol+9:icol+11] .= f_ψ1_Mdot
    jacob[irow_p1+3:irow_p1+5, icol+12:icol+14] .= f_ψ1_Vdot ./ force_scaling
    jacob[irow_p1+3:irow_p1+5, icol+15:icol+17] .= f_ψ1_Ωdot ./ force_scaling

    # create jacobian entries for the end of the beam
    jacob[irow_p2:irow_p2+2, icol+3:icol+5] .= f_u2_θdot ./ force_scaling
    jacob[irow_p2:irow_p2+2, icol+6:icol+8] .= f_u2_Fdot
    jacob[irow_p2:irow_p2+2, icol+12:icol+14] .= f_u2_Vdot  ./ force_scaling
    jacob[irow_p2:irow_p2+2, icol+15:icol+17] .= f_u2_Ωdot  ./ force_scaling
    jacob[irow_p2+3:irow_p2+5, icol+3:icol+5] .= f_ψ2_θdot ./ force_scaling
    jacob[irow_p2+3:irow_p2+5, icol+6:icol+8] .= f_ψ2_Fdot
    jacob[irow_p2+3:irow_p2+5, icol+9:icol+11] .= f_ψ2_Mdot
    jacob[irow_p2+3:irow_p2+5, icol+12:icol+14] .= f_ψ2_Vdot ./ force_scaling
    jacob[irow_p2+3:irow_p2+5, icol+15:icol+17] .= f_ψ2_Ωdot ./ force_scaling

    # create jacobian entries for beam residual equations
    jacob[irow_e:irow_e+2, icol:icol+2] .= f_V_udot
    jacob[irow_e+3:irow_e+5, icol+3:icol+5] .= f_Ω_θdot

    return jacob
end

"""
    insert_element_mass_matrix!(jacob, gamma, force_scaling, irow_e,
        irow_p1, irow_p2, icol, 
        f_u1_θdot, f_u2_θdot, f_u1_Fdot, f_u2_Fdot, f_u1_Vdot, f_u2_Vdot, 
        f_u1_Ωdot, f_u2_Ωdot, f_ψ1_θdot, f_ψ2_θdot, f_ψ1_Fdot, f_ψ2_Fdot, 
        f_ψ1_Mdot, f_ψ2_Mdot, f_ψ1_Vdot, f_ψ2_Vdot, f_ψ1_Ωdot, f_ψ2_Ωdot,
        f_V_udot, f_Ω_θdot)

Add the beam element's mass matrix to the system jacobian matrix `jacob`, scaled
by the scaling parameter `gamma`.
"""
@inline function insert_element_mass_matrix!(jacob, gamma, force_scaling, irow_e,
    irow_p1, irow_p2, icol, 
    f_u1_θdot, f_u2_θdot, f_u1_Fdot, f_u2_Fdot, f_u1_Vdot, f_u2_Vdot, 
    f_u1_Ωdot, f_u2_Ωdot, f_ψ1_θdot, f_ψ2_θdot, f_ψ1_Fdot, f_ψ2_Fdot, 
    f_ψ1_Mdot, f_ψ2_Mdot, f_ψ1_Vdot, f_ψ2_Vdot, f_ψ1_Ωdot, f_ψ2_Ωdot,
    f_V_udot, f_Ω_θdot)

    # create jacobian entries for the start of the beam
    jacob[irow_p1:irow_p1+2, icol+3:icol+5] .= f_u1_θdot .* (gamma/force_scaling)
    jacob[irow_p1:irow_p1+2, icol+6:icol+8] .= f_u1_Fdot .* gamma
    jacob[irow_p1:irow_p1+2, icol+12:icol+14] .= f_u1_Vdot  .* (gamma/force_scaling)
    jacob[irow_p1:irow_p1+2, icol+15:icol+17] .= f_u1_Ωdot  .* (gamma/force_scaling)
    jacob[irow_p1+3:irow_p1+5, icol+3:icol+5] .= f_ψ1_θdot .* (gamma/force_scaling)
    jacob[irow_p1+3:irow_p1+5, icol+6:icol+8] .= f_ψ1_Fdot .* gamma
    jacob[irow_p1+3:irow_p1+5, icol+9:icol+11] .= f_ψ1_Mdot .* gamma
    jacob[irow_p1+3:irow_p1+5, icol+12:icol+14] .= f_ψ1_Vdot .* (gamma/force_scaling)
    jacob[irow_p1+3:irow_p1+5, icol+15:icol+17] .= f_ψ1_Ωdot .* (gamma/force_scaling)

    # create jacobian entries for the end of the beam
    jacob[irow_p2:irow_p2+2, icol+3:icol+5] .= f_u2_θdot .* (gamma/force_scaling)
    jacob[irow_p2:irow_p2+2, icol+6:icol+8] .= f_u2_Fdot .* gamma
    jacob[irow_p2:irow_p2+2, icol+12:icol+14] .= f_u2_Vdot  .* (gamma/force_scaling)
    jacob[irow_p2:irow_p2+2, icol+15:icol+17] .= f_u2_Ωdot  .* (gamma/force_scaling)
    jacob[irow_p2+3:irow_p2+5, icol+3:icol+5] .= f_ψ2_θdot .* (gamma/force_scaling)
    jacob[irow_p2+3:irow_p2+5, icol+6:icol+8] .= f_ψ2_Fdot .* gamma
    jacob[irow_p2+3:irow_p2+5, icol+9:icol+11] .= f_ψ2_Mdot .* gamma
    jacob[irow_p2+3:irow_p2+5, icol+12:icol+14] .= f_ψ2_Vdot .* (gamma/force_scaling)
    jacob[irow_p2+3:irow_p2+5, icol+15:icol+17] .= f_ψ2_Ωdot .* (gamma/force_scaling)

    # create jacobian entries for beam residual equations
    jacob[irow_e:irow_e+2, icol:icol+2] .= f_V_udot .* gamma
    jacob[irow_e+3:irow_e+5, icol+3:icol+5] .= f_Ω_θdot .* gamma

    return jacob
end

"""
    element_mass_matrix!(jacob, x, elem, point_masses, force_scaling, icol,
        irow_e, irow_e1, irow_p1, irow_e2, irow_p2)

Insert the beam element's contributions to the "mass matrix": the jacobian of the
residual equations with respect to the time derivatives of the state variables

# Arguments
 - `jacob`: system jacobian matrix
 - `x`: current state vector
 - `ielem`: beam element index
 - `elem`: beam element
 - `point_masses`: dictionary with all point masses
 - `force_scaling`: scaling parameter for forces/moments
 - `icol`: starting index for the beam's state variables
 - `irow_e1`: row index of the first residual equation for the start of the beam element
    (a value <= 0 indicates the equations have been eliminated from the system of equations)
 - `irow_p1`: row index of the first residual equation for the point at the start of the 
    beam element
 - `irow_e2`: row index of the first residual equation for the end of the beam element
    (a value <= 0 indicates the equations have been eliminated from the system of equations)
 - `irow_p2`: row index of the first residual equation for the point at the end of the 
    beam element
 - `gamma`: Scaling parameter for scaling mass matrix contribution to `jacob`
"""
@inline function element_mass_matrix!(jacob, x, ielem, elem, point_masses, force_scaling,
    icol, irow_e, irow_p1, irow_p2)

    # element properties
    ΔL = elem.L
    Cab = elem.Cab
    compliance = elem.compliance
    mass = elem.mass
    μ = elem.mu

    # scale compliance and mass matrix by the element length (to allow zero length elements)
    compliance *= ΔL
    mass *= ΔL
    
    # modify mass matrix to account for point masses, if present
    if haskey(point_masses, ielem)
        mass += transform_properties(point_masses[ielem].mass, Cab)
    end

    # element state variables
    u, θ, F, M, V, Ω = dynamic_element_state_variables(x, icol, force_scaling)

    # compliance submatrices
    S11 = compliance[SVector{3}(1:3), SVector{3}(1:3)]
    S12 = compliance[SVector{3}(1:3), SVector{3}(4:6)]

    # mass submatrices
    mass11 = mass[SVector{3}(1:3), SVector{3}(1:3)]
    mass12 = mass[SVector{3}(1:3), SVector{3}(4:6)]
    mass21 = mass[SVector{3}(4:6), SVector{3}(1:3)]
    mass22 = mass[SVector{3}(4:6), SVector{3}(4:6)]

    # damping coefficient submatrices
    μ11 = @SMatrix [μ[1] 0 0; 0 μ[2] 0; 0 0 μ[3]]
    μ22 = @SMatrix [μ[4] 0 0; 0 μ[5] 0; 0 0 μ[6]]

    # rotation matrices
    C = get_C(θ)
    Ct = C'
    CtCab = Ct*Cab

    # rotation matrix derivatives
    Cdot_θdot1, Cdot_θdot2, Cdot_θdot3 = get_C_t_θdot(C, θ)
    Ctdot_θdot1, Ctdot_θdot2, Ctdot_θdot3 = Cdot_θdot1', Cdot_θdot2', Cdot_θdot3'

    # element strain
    γ = S11*F + S12*M

    # element linear and angular momentum
    P = mass11*V + mass12*Ω
    H = mass21*V + mass22*Ω

    # get jacobians of beam element equations
    f_u1_θdot, f_u2_θdot, f_u1_Fdot, f_u2_Fdot, f_u1_Vdot, f_u2_Vdot, 
        f_u1_Ωdot, f_u2_Ωdot, f_ψ1_θdot, f_ψ2_θdot, f_ψ1_Fdot, f_ψ2_Fdot, 
        f_ψ1_Mdot, f_ψ2_Mdot, f_ψ1_Vdot, f_ψ2_Vdot, f_ψ1_Ωdot, f_ψ2_Ωdot,
        f_V_udot, f_Ω_θdot = element_mass_matrix_equations(ΔL, mass11, mass12, mass21, 
        mass22, μ11, μ22, Cab, CtCab, θ, P, H, γ, Ctdot_θdot1, Ctdot_θdot2, Ctdot_θdot3)

    # initialize/insert into jacobian matrix for the system
    insert_element_mass_matrix!(jacob, force_scaling, irow_e,
        irow_p1, irow_p2, icol, 
        f_u1_θdot, f_u2_θdot, f_u1_Fdot, f_u2_Fdot, f_u1_Vdot, f_u2_Vdot, 
        f_u1_Ωdot, f_u2_Ωdot, f_ψ1_θdot, f_ψ2_θdot, f_ψ1_Fdot, f_ψ2_Fdot, 
        f_ψ1_Mdot, f_ψ2_Mdot, f_ψ1_Vdot, f_ψ2_Vdot, f_ψ1_Ωdot, f_ψ2_Ωdot,
        f_V_udot, f_Ω_θdot)

    return jacob
end

"""
    element_mass_matrix!(jacob, gamma, x, ielem, elem, force_scaling,
        icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2)

Add the beam element's mass matrix to the system jacobian matrix `jacob`, scaled
by the scaling parameter `gamma`.
"""
@inline function element_mass_matrix!(jacob, gamma, x, ielem, elem, point_masses, 
    force_scaling, icol, irow_e, irow_p1, irow_p2)

    # element properties
    ΔL = elem.L
    Cab = elem.Cab
    compliance = elem.compliance
    mass = elem.mass
    μ = elem.mu

    # scale compliance and mass matrix by the element length (to allow zero length elements)
    compliance *= ΔL
    mass *= ΔL
    
    # modify mass matrix to account for point masses, if present
    if haskey(point_masses, ielem)
        mass += transform_properties(point_masses[ielem].mass, Cab)
    end

    # element state variables
    u, θ, F, M, V, Ω = dynamic_element_state_variables(x, icol, force_scaling)

    # compliance submatrices
    S11 = compliance[SVector{3}(1:3), SVector{3}(1:3)]
    S12 = compliance[SVector{3}(1:3), SVector{3}(4:6)]

    # mass submatrices
    mass11 = mass[SVector{3}(1:3), SVector{3}(1:3)]
    mass12 = mass[SVector{3}(1:3), SVector{3}(4:6)]
    mass21 = mass[SVector{3}(4:6), SVector{3}(1:3)]
    mass22 = mass[SVector{3}(4:6), SVector{3}(4:6)]

    # damping coefficient submatrices
    μ11 = @SMatrix [μ[1] 0 0; 0 μ[2] 0; 0 0 μ[3]]
    μ22 = @SMatrix [μ[4] 0 0; 0 μ[5] 0; 0 0 μ[6]]

    # rotation matrices
    C = get_C(θ)
    Ct = C'
    CtCab = Ct*Cab

    # rotation matrix derivatives
    Cdot_θdot1, Cdot_θdot2, Cdot_θdot3 = get_C_t_θdot(C, θ)
    Ctdot_θdot1, Ctdot_θdot2, Ctdot_θdot3 = Cdot_θdot1', Cdot_θdot2', Cdot_θdot3'

    # element strain
    γ = S11*F + S12*M

    # element linear and angular momentum
    P = mass11*V + mass12*Ω
    H = mass21*V + mass22*Ω

    # get jacobians of beam element equations
    f_u1_θdot, f_u2_θdot, f_u1_Fdot, f_u2_Fdot, f_u1_Vdot, f_u2_Vdot, 
        f_u1_Ωdot, f_u2_Ωdot, f_ψ1_θdot, f_ψ2_θdot, f_ψ1_Fdot, f_ψ2_Fdot, 
        f_ψ1_Mdot, f_ψ2_Mdot, f_ψ1_Vdot, f_ψ2_Vdot, f_ψ1_Ωdot, f_ψ2_Ωdot,
        f_V_udot, f_Ω_θdot = element_mass_matrix_equations(ΔL, mass11, mass12, mass21, 
        mass22, μ11, μ22, Cab, CtCab, θ, P, H, γ, Ctdot_θdot1, Ctdot_θdot2, Ctdot_θdot3)

    # initialize/insert into jacobian matrix for the system
    insert_element_mass_matrix!(jacob, gamma, force_scaling, irow_e,
        irow_p1, irow_p2, icol, 
        f_u1_θdot, f_u2_θdot, f_u1_Fdot, f_u2_Fdot, f_u1_Vdot, f_u2_Vdot, 
        f_u1_Ωdot, f_u2_Ωdot, f_ψ1_θdot, f_ψ2_θdot, f_ψ1_Fdot, f_ψ2_Fdot, 
        f_ψ1_Mdot, f_ψ2_Mdot, f_ψ1_Vdot, f_ψ2_Vdot, f_ψ1_Ωdot, f_ψ2_Ωdot,
        f_V_udot, f_Ω_θdot)

    return jacob
end