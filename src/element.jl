"""
    Element{TF}

Composite type that defines a beam element's properties

# Fields
 - `L`: Length of the beam element
 - `x`: Location of the beam element (the center of the beam element)
 - `compliance`: Beam element compliance matrix
 - `mass`: Beam element mass matrix
 - `Cab`: Transformation matrix from the undeformed beam element frame to the body frame
"""
struct Element{TF}
    L::TF
    x::SVector{3,TF}
    compliance::SMatrix{6,6,TF,36}
    mass::SMatrix{6,6,TF,36}
    Cab::SMatrix{3,3,TF,9}
end

"""
    Element(L, x, compliance, mass, Cab)

Construct a beam element

# Arguments
- `L`: Length of the beam element
- `x`: Location of the beam element (the center of the beam element)
- `compliance`: Beam element compliance matrix
- `mass`: Beam element mass matrix
- `Cab`: Transformation matrix from the undeformed beam element frame to the body frame
"""
function Element(L, x, compliance, mass, Cab)
    TF = promote_type(typeof(L), eltype(x), eltype(compliance), eltype(mass), eltype(Cab))
    return Element{TF}(L, x, compliance, mass, Cab)
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
    static_element_properties(x, icol, elem, force_scaling)

Extract/calculate the properties of a beam element for a static analysis.

# Arguments
 - `x`: current state vector
 - `icol`: starting index for the beam's state variables
 - `elem`: beam element
 - `force_scaling`: scaling parameter for forces/moments
"""
@inline function static_element_properties(x, icol, elem, force_scaling)

    # state variables
    u = SVector(x[icol   ], x[icol+1 ], x[icol+2 ])
    θ = SVector(x[icol+3 ], x[icol+4 ], x[icol+5 ])
    F = SVector(x[icol+6 ], x[icol+7 ], x[icol+8 ]) .* force_scaling
    M = SVector(x[icol+9 ], x[icol+10], x[icol+11]) .* force_scaling

    # fixed element properties
    ΔL = elem.L
    compliance = elem.compliance
    mass = elem.mass

    # rotation matrices
    C = get_C(θ)
    Cab = elem.Cab
    CtCab = C'*Cab
    
    # element strain and curvature
    γ = element_strain(elem, F, M)
    κ = element_curvature(elem, F, M)

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

    return ΔL, S11, S12, S21, S22, mass11, mass12, mass21, mass22, C, Cab, CtCab, u, θ, 
        F, M, γ, κ
end

"""
    steady_state_element_properties(x, icol, elem, force_scaling, x0, v0, ω0, a0, α0)

Extract/calculate the properties of a beam element for a steady state analysis.

# Arguments
 - `x`: current state vector
 - `icol`: starting index for the beam's state variables
 - `elem`: beam element
 - `force_scaling`: scaling parameter for forces/moments
 - `x0`: body frame origin (for the current time step)
 - `v0`: body frame linear velocity (for the current time step)
 - `ω0`: body frame angular velocity (for the current time step)
 - `a0`: body frame linear acceleration (for the current time step)
 - `α0`: body frame angular acceleration (for the current time step)
"""
@inline function steady_state_element_properties(x, icol, elem, force_scaling, x0, v0, ω0, 
    a0, α0)

    ΔL, S11, S12, S21, S22, mass11, mass12, mass21, mass22, C, Cab, CtCab, u, θ, 
        F, M, γ, κ = static_element_properties(x, icol, elem, force_scaling)

    # undeformed element linear and angular velocity from body frame motion
    v = v0 + cross(ω0, elem.x - x0)
    ω = ω0

    # undeformed element linear and angular acceleration from body frame motion
    a = a0 + cross(α0, elem.x - x0)
    α = α0 

    # element linear and angular velocity
    V = SVector(x[icol+12], x[icol+13], x[icol+14])
    Ω = SVector(x[icol+15], x[icol+16], x[icol+17])

    # element linear and angular momentum
    P = element_linear_momentum(elem, V, Ω)
    H = element_angular_momentum(elem, V, Ω)

    return ΔL, S11, S12, S21, S22, mass11, mass12, mass21, mass22, C, Cab, CtCab, u, θ, 
        F, M, γ, κ, V, Ω, P, H, v, ω, a, α
end

"""
    initial_condition_element_properties(x, icol, elem, force_scaling, x0, v0, ω0, a0, α0, 
        u0, θ0, udot0, θdot0)

Extract/calculate the properties of a beam element for an initial condition analysis.

# Arguments
 - `x`: current state vector
 - `icol`: starting index for the beam's state variables
 - `elem`: beam element
 - `force_scaling`: scaling parameter for forces/moments
 - `x0`: body frame origin (for the current time step)
 - `v0`: body frame linear velocity (for the current time step)
 - `ω0`: body frame angular velocity (for the current time step)
 - `a0`: body frame linear acceleration (for the current time step)
 - `α0`: body frame angular acceleration (for the current time step)
 - `u0`: initial linear deflections for the beam element
 - `θ0`: initial angular deflections for the beam element
 - `udot0`: initial linear deflection rates for the beam element
 - `θdot0`: initial angular deflection rates for the beam element
"""
@inline function initial_condition_element_properties(x, icol, elem, force_scaling, x0, v0, 
    ω0, a0, α0, u0, θ0, udot0, θdot0)

    # state variables
    u = u0
    θ = θ0
    F = SVector(x[icol+6 ], x[icol+7 ], x[icol+8 ]) .* force_scaling
    M = SVector(x[icol+9 ], x[icol+10], x[icol+11]) .* force_scaling

    # fixed element properties
    ΔL = elem.L
    compliance = elem.compliance
    mass = elem.mass

    # rotation matrices
    C = get_C(θ)
    Cab = elem.Cab
    CtCab = C'*Cab

    # element strain and curvature
    γ = element_strain(elem, F, M)
    κ = element_curvature(elem, F, M)

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

    # element linear and angular velocity from body frame motion
    v = v0 + cross(ω0, elem.x - x0)
    ω = ω0

    # element linear and angular acceleration from body frame motion
    a = a0 + cross(α0, elem.x - x0) + cross(α0, u)
    α = α0

    # element linear and angular velocity
    V = SVector(x[icol+12], x[icol+13], x[icol+14])
    Ω = SVector(x[icol+15], x[icol+16], x[icol+17])

    # element linear and angular momentum
    P = element_linear_momentum(elem, V, Ω)
    H = element_angular_momentum(elem, V, Ω)

    # state rates
    udot = udot0
    θdot = θdot0
    Vdot = SVector(x[icol   ], x[icol+1 ], x[icol+2 ])
    Ωdot = SVector(x[icol+3 ], x[icol+4 ], x[icol+5 ])

    # momentum from local motion
    Pdot = element_linear_momentum(elem, Vdot, Ωdot)
    Hdot = element_angular_momentum(elem, Vdot, Ωdot)

    # rotation matrix time derivatives
    Cdot = get_C_t(C, θ, θdot)
    CtCabdot = Cdot'*Cab

    # linear and angular momentum in body reference frame
    CtCabPdot = CtCabdot*P + CtCab*Pdot
    CtCabHdot = CtCabdot*H + CtCab*Hdot

    return ΔL, S11, S12, S21, S22, mass11, mass12, mass21, mass22, C, Cab, CtCab, u, θ, 
        F, M, γ, κ, V, Ω, P, H, v, ω, a, α, udot, θdot, Vdot, Ωdot, Pdot, Hdot,
        Cdot, CtCabdot, CtCabPdot, CtCabHdot
end

"""
    newmark_element_properties(x, icol, elem, force_scaling, x0, v0, ω0, a0, α0, 
        udot_init, θdot_init, Vdot_init, Ωdot_init, dt)

Extract/calculate the properties of a beam element for a Newmark scheme time-marching 
analysis.

# Arguments
 - `x`: current state vector
 - `icol`: starting index for the beam's state variables
 - `elem`: beam element
 - `force_scaling`: scaling parameter for forces/moments
 - `x0`: body frame origin (for the current time step)
 - `v0`: body frame linear velocity (for the current time step)
 - `ω0`: body frame angular velocity (for the current time step)
 - `a0`: body frame linear acceleration (for the current time step)
 - `α0`: body frame angular acceleration (for the current time step)
 - `udot_init`: `2/dt*u + udot` for the beam element from the previous time step
 - `θdot_init`: `2/dt*θ + θdot` for the beam element from the previous time step
 - `Vdot_init`: `2/dt*C'*Cab*P + C'*Cab*Pdot` for the beam element from the 
    previous time step
 - `Ωdot_init`: `2/dt*C'*Cab*H + C'*Cab*Hdot` for the beam element from the 
    previous time step
 - `dt`: time step size
"""
@inline function newmark_element_properties(x, icol, elem, force_scaling, x0, v0, ω0, a0, α0, 
    udot_init, θdot_init, Vdot_init, Ωdot_init, dt)

    ΔL, S11, S12, S21, S22, mass11, mass12, mass21, mass22, C, Cab, CtCab, u, θ, F, M, 
        γ, κ, V, Ω, P, H, v, ω, a, α = steady_state_element_properties(x, icol, elem, 
        force_scaling, x0, v0, ω0, a0, α0)

    # state rates
    udot = 2/dt*u - udot_init
    θdot = 2/dt*θ - θdot_init
    Vdot = 2/dt*V - Vdot_init
    Ωdot = 2/dt*Ω - Ωdot_init

    # momentum from local motion
    Pdot = element_linear_momentum(elem, Vdot, Ωdot)
    Hdot = element_angular_momentum(elem, Vdot, Ωdot)

    # rotation matrix time derivatives
    Cdot = get_C_t(C, θ, θdot)
    CtCabdot = Cdot'*Cab

    # linear and angular momentum in body reference frame
    CtCabPdot = CtCabdot*P + CtCab*Pdot
    CtCabHdot = CtCabdot*H + CtCab*Hdot

    return ΔL, S11, S12, S21, S22, mass11, mass12, mass21, mass22, C, Cab, CtCab, u, θ, 
        F, M, γ, κ, V, Ω, P, H, v, ω, a, α, udot, θdot, Vdot, Ωdot, Pdot, Hdot, 
        Cdot, CtCabdot, CtCabPdot, CtCabHdot
end

"""
    dynamic_element_properties(x, icol, elem, force_scaling, x0, v0, ω0, a0, α0, 
        udot, θdot, Vdot, Ωdot)

Extract/calculate the properties of a beam element for a general dynamic analysis.

# Arguments
 - `x`: current state vector
 - `icol`: starting index for the beam's state variables
 - `elem`: beam element
 - `force_scaling`: scaling parameter for forces/moments
 - `x0`: body frame origin (for the current time step)
 - `v0`: body frame linear velocity (for the current time step)
 - `ω0`: body frame angular velocity (for the current time step)
 - `a0`: body frame linear acceleration (for the current time step)
 - `α0`: body frame angular acceleration (for the current time step)
 - `udot`: linear deflection rates for the beam element
 - `θdot`: angular deflection rates for the beam element
 - `Vdot`: linear velocity rates for the beam element
 - `Ωdot`: angular velocity rates for the beam element
"""
@inline function dynamic_element_properties(x, icol, elem, force_scaling, x0, v0, ω0, a0, α0, 
    udot, θdot, Vdot, Ωdot)

    ΔL, S11, S12, S21, S22, mass11, mass12, mass21, mass22, C, Cab, CtCab, u, θ, F, M, 
        γ, κ, V, Ω, P, H, v, ω, a, α = steady_state_element_properties(x, icol, elem, 
        force_scaling, x0, v0, ω0, a0, α0)

    # momentum from local motion
    Pdot = element_linear_momentum(elem, Vdot, Ωdot)
    Hdot = element_angular_momentum(elem, Vdot, Ωdot)

    # rotation matrix time derivatives
    Cdot = get_C_t(C, θ, θdot)
    CtCabdot = Cdot'*Cab

    # linear and angular momentum in body reference frame
    CtCabPdot = CtCabdot*P + CtCab*Pdot
    CtCabHdot = CtCabdot*H + CtCab*Hdot

    return ΔL, S11, S12, S21, S22, mass11, mass12, mass21, mass22, C, Cab, CtCab, u, θ, 
        F, M, γ, κ, V, Ω, P, H, v, ω, a, α, udot, θdot, Pdot, Hdot, Cdot, CtCabdot,
        CtCabPdot, CtCabHdot
        
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

    tmp = CtCab*F
    f_u1 = -tmp - f1
    f_u2 =  tmp - f2

    tmp1 = CtCab*M
    tmp2 = ΔL/2*CtCab*cross(e1 + γ, F)
    f_ψ1 = -tmp1 - m1 - tmp2
    f_ψ2 =  tmp1 - m2 - tmp2

    tmp = ΔL/2*(CtCab*(e1 + γ) - Cab*e1)
    f_F1 =  u - tmp
    f_F2 = -u - tmp

    Qinv = get_Qinv(θ)
    tmp = ΔL/2*Qinv*Cab*κ
    f_M1 =  θ - tmp
    f_M2 = -θ - tmp

    return f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2
end

"""
    steady_state_element_equations(ΔL, Cab, CtCab, u, θ, F, M, γ, κ, V, Ω, P, H, v, ω, f1, 
        f2, m1, m2)

Calculate an element's resultants for a steady state analysis.

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
    
    tmp = ΔL/2*cross(ω, CtCab*P)
    f_u1 += tmp
    f_u2 += tmp

    tmp = ΔL/2*(cross(ω, CtCab*H) + CtCab*cross(V, P))
    f_ψ1 += tmp
    f_ψ2 += tmp

    f_V = CtCab*V - v - cross(ω, u)
    f_Ω = Ω - CtCab'*ω

    return f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2, f_V, f_Ω
end

"""
    dynamic_element_equations(ΔL, Cab, CtCab, u, θ, F, M, γ, κ, V, Ω, P, H, 
        v, ω, f1, f2, m1, m2, udot, θdot, CtCabPdot, CtCabHdot)

Calculate an element's resultants for a general dynamic analysis.

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

    tmp = ΔL/2*CtCabPdot
    f_u1 += tmp
    f_u2 += tmp

    tmp = ΔL/2*CtCabHdot
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
 - `irow_e1`: row index of the first residual equation for the end of the beam element
    (a value <= 0 indicates the equations have been eliminated from the system of equations)
 - `irow_p2`: row index of the first residual equation for the point at the end of the 
    beam element
 - `f_u1`, `f_u2`: resultant forces corresponding to the start and end of the beam element, 
    respectively
 - `f_ψ1`, `f_ψ2`: resultant moments corresponding to the start and end of the beam element, 
    respectively
 - `f_F1`, `f_F2`: resultant linear displacements corresponding to the start and end of the
    beam element, respectively.
 - `f_M1`, `f_M2`: resultant angular displacements corresponding to the start and end of the 
    beam element, respectively.
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
 - `f_u1`, `f_u2`: resultant forces corresponding to the start and end of the beam element, 
   respectively
 - `f_ψ1`, `f_ψ2`: resultant moments corresponding to the start and end of the beam element, 
   respectively
 - `f_F1`, `f_F2`: resultant linear displacements corresponding to the start and end of the
   beam element, respectively.
 - `f_M1`, `f_M2`: resultant angular displacements corresponding to the start and end of the 
    beam element, respectively.
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
    ΔL, S11, S12, S21, S22, mass11, mass12, mass21, mass22, C, Cab, CtCab, u, θ, 
        F, M, γ, κ = static_element_properties(x, icol, elem, force_scaling)

    # additional precomputed quantities
    Ct = C'
    
    # gravitational loads
    f1, f2, m1, m2 = element_gravitational_loads(ΔL, mass11, mass12, CtCab, gvec)

    # distributed loads
    if haskey(distributed_loads, ielem)
        dload = distributed_loads[ielem]
        f1 += dload.f1 + Ct*dload.f1_follower
        f2 += dload.f2 + Ct*dload.f2_follower
        m1 += dload.m1 + Ct*dload.m1_follower
        m2 += dload.m2 + Ct*dload.m2_follower
    end

    # point mass loads
    if haskey(point_masses, ielem)
        # point mass properties
        massp11, massp12, massp21, massp22 = static_point_mass_properties(point_masses[ielem])
        # point mass loads
        Fp, Mp = static_point_mass_loads(massp11, massp12, C, Ct, gvec)
        # split the load between the two beam endpoints
        f1 += Fp/2
        f2 += Fp/2
        m1 += Mp/2
        m2 += Mp/2
    end

    # solve for the element resultants
    f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2 = static_element_equations(ΔL, Cab, 
        CtCab, u, θ, F, M, γ, κ, f1, f2, m1, m2)

    # insert element resultants into the residual vector
    static_insert_element_residual!(resid, force_scaling, irow_e1, irow_p1, irow_e2, 
        irow_p2, f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2)

    return resid
end

"""
    steady_state_element_residual!(resid, x, ielem, elem, distributed_loads, point_masses, gvec,
        force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2, x0, v0, ω0, a0, α0)

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
@inline function steady_state_element_residual!(resid, x, ielem, elem, distributed_loads, point_masses, gvec,
    force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2, x0, v0, ω0, a0, α0)

    # compute element properties
    ΔL, S11, S12, S21, S22, mass11, mass12, mass21, mass22, C, Cab, CtCab, u, θ, F, M, 
        γ, κ, V, Ω, P, H, v, ω, a, α = steady_state_element_properties(x, icol, elem, 
        force_scaling, x0, v0, ω0, a0, α0)
   
    # additional precomputed quantities
    Ct = C'
    
    # element acceleration loads (including gravitational loads)
    f1, f2, m1, m2 = element_acceleration_loads(ΔL, mass11, mass12, mass21, mass22, CtCab, 
        u, a, α, gvec)

    # distributed loads
    if haskey(distributed_loads, ielem)
        dload = distributed_loads[ielem]
        f1 += dload.f1 + Ct*dload.f1_follower
        f2 += dload.f2 + Ct*dload.f2_follower
        m1 += dload.m1 + Ct*dload.m1_follower
        m2 += dload.m2 + Ct*dload.m2_follower
    end

    # point mass loads
    if haskey(point_masses, ielem)
        # point mass properties
        massp11, massp12, massp21, massp22, Vp, Ωp, Pp, Hp = 
            steady_state_point_mass_properties(point_masses[ielem], Cab, V, Ω)
        # point mass loads
        Fp, Mp = steady_state_point_mass_loads(massp11, massp12, massp21, massp22, C, Ct, 
            u, Vp, Pp, Hp, ω, a, α, gvec)
        # split the load between the two beam endpoints
        f1 += Fp/2
        f2 += Fp/2
        m1 += Mp/2
        m2 += Mp/2
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
        irow_p2, x0, v0, ω0, a0, α0, u0, θ0, udot0, θdot0)

Compute and add a beam element's contributions to the residual vector for an initial 
condition analysis

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
"""
@inline function initial_condition_element_residual!(resid, x, ielem, elem, 
    distributed_loads, point_masses, gvec, force_scaling, icol, irow_e, irow_e1, irow_p1, 
    irow_e2, irow_p2, x0, v0, ω0, a0, α0, u0, θ0, udot0, θdot0)

    # compute element properties
    ΔL, S11, S12, S21, S22, mass11, mass12, mass21, mass22, C, Cab, CtCab, u, θ, 
        F, M, γ, κ, V, Ω, P, H, v, ω, a, α, udot, θdot, Vdot, Ωdot, Pdot, Hdot,
        Cdot, CtCabdot, CtCabPdot, CtCabHdot = initial_condition_element_properties(x, 
        icol, elem, force_scaling, x0, v0, ω0, a0, α0, u0, θ0, udot0, θdot0)

    # additional precomputed quantities
    Ct = C'
    Ctdot = Cdot'

    # element acceleration loads (including gravitational loads)
    f1, f2, m1, m2 = element_acceleration_loads(ΔL, mass11, mass12, mass21, mass22, CtCab, 
        u, a, α, gvec)

    # distributed loads
    if haskey(distributed_loads, ielem)
        dload = distributed_loads[ielem]
        f1 += dload.f1 + Ct*dload.f1_follower
        f2 += dload.f2 + Ct*dload.f2_follower
        m1 += dload.m1 + Ct*dload.m1_follower
        m2 += dload.m2 + Ct*dload.m2_follower
    end

    # point mass loads
    if haskey(point_masses, ielem)
        # point mass properties
        massp11, massp12, massp21, massp22, Vp, Ωp, Pp, Hp, Vpdot, Ωpdot, Ppdot, Hpdot, 
            CtPpdot, CtHpdot = dynamic_point_mass_properties(point_masses[ielem], Ct, Ctdot, 
            Cab, V, Ω, Vdot, Ωdot)
        # point mass loads
        Fp, Mp = dynamic_point_mass_loads(massp11, massp12, massp21, massp22, 
            C, Ct, Ctdot, u, Vp, Pp, Hp, Ppdot, Hpdot, ω, a, α, gvec)
        # split the load between the two beam endpoints
        f1 += Fp/2
        f2 += Fp/2
        m1 += Mp/2
        m2 += Mp/2
    end

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
        udot_init, θdot_init, Vdot_init, Ωdot_init, dt)

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
 - `Vdot_init`: `2/dt*V + Vdot` for the beam element from the previous time step
 - `Ωdot_init`: `2/dt*Ω + Ωdot` for the beam element from the previous time step
 - `dt`: time step size
"""
@inline function newmark_element_residual!(resid, x, ielem, elem, distributed_loads, point_masses, gvec,
    force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2, x0, v0, ω0, a0, α0, udot_init, θdot_init,
    Vdot_init, Ωdot_init, dt)

    # compute element properties
    ΔL, S11, S12, S21, S22, mass11, mass12, mass21, mass22, C, Cab, CtCab, u, θ, 
        F, M, γ, κ, V, Ω, P, H, v, ω, a, α, udot, θdot, Vdot, Ωdot, Pdot, Hdot, 
        Cdot, CtCabdot, CtCabPdot, CtCabHdot = newmark_element_properties(x, icol, elem, 
        force_scaling, x0, v0, ω0, a0, α0, udot_init, θdot_init, Vdot_init, Ωdot_init, dt)
    
    # additional precomputed quantities
    Ct = C'
    Ctdot = Cdot'

    # element acceleration loads (including gravitational loads)
    f1, f2, m1, m2 = element_acceleration_loads(ΔL, mass11, mass12, mass21, mass22, CtCab, 
        u, a, α, gvec)

    # distributed loads
    if haskey(distributed_loads, ielem)
        dload = distributed_loads[ielem]
        f1 += dload.f1 + Ct*dload.f1_follower
        f2 += dload.f2 + Ct*dload.f2_follower
        m1 += dload.m1 + Ct*dload.m1_follower
        m2 += dload.m2 + Ct*dload.m2_follower
    end

    # point mass loads
    if haskey(point_masses, ielem)
        # get point mass properties
        massp11, massp12, massp21, massp22, Vp, Ωp, Pp, Hp, Vpdot, Ωpdot, Ppdot, Hpdot, 
            CtPpdot, CtHpdot = dynamic_point_mass_properties(point_masses[ielem], Ct, Ctdot, 
            Cab, V, Ω, Vdot, Ωdot)
        # get point mass loads
        Fp, Mp = dynamic_point_mass_loads(massp11, massp12, massp21, massp22, 
            C, Ct, Ctdot, u, Vp, Pp, Hp, Ppdot, Hpdot, ω, a, α, gvec)
        # split the load between the two beam endpoints
        f1 += Fp/2
        f2 += Fp/2
        m1 += Mp/2
        m2 += Mp/2
    end

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
        udot, θdot, Vdot, Ωdot)

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
 - `Vdot`: linear velocity rates for the beam element
 - `Ωdot`: angular velocity rates for the beam element
"""
@inline function dynamic_element_residual!(resid, x, ielem, elem, distributed_loads, point_masses, gvec,
    force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2, x0, v0, ω0, a0, α0,
    udot, θdot, Vdot, Ωdot)

    # compute element properties
    ΔL, S11, S12, S21, S22, mass11, mass12, mass21, mass22, C, Cab, CtCab, u, θ, 
        F, M, γ, κ, V, Ω, P, H, v, ω, a, α, udot, θdot, Pdot, Hdot, Cdot, CtCabdot,
        CtCabPdot, CtCabHdot = dynamic_element_properties(x, icol, elem, force_scaling, 
        x0, v0, ω0, a0, α0, udot, θdot, Vdot, Ωdot)

    # additional precomputed quantities
    Ct = C'
    Ctdot = Cdot'

    # element acceleration loads (including gravitational loads)
    f1, f2, m1, m2 = element_acceleration_loads(ΔL, mass11, mass12, mass21, mass22, CtCab, 
        u, a, α, gvec)
        
    # distributed loads
    if haskey(distributed_loads, ielem)
        dload = distributed_loads[ielem]
        f1 += dload.f1 + Ct*dload.f1_follower
        f2 += dload.f2 + Ct*dload.f2_follower
        m1 += dload.m1 + Ct*dload.m1_follower
        m2 += dload.m2 + Ct*dload.m2_follower
    end

    # point mass loads
    if haskey(point_masses, ielem)
        # point mass properties
        massp11, massp12, massp21, massp22, Vp, Ωp, Pp, Hp, Vpdot, Ωpdot, Ppdot, Hpdot, 
            CtPpdot, CtHpdot = dynamic_point_mass_properties(point_masses[ielem], Ct, Ctdot, 
            Cab, V, Ω, Vdot, Ωdot)
        # point mass loads
        Fp, Mp = dynamic_point_mass_loads(massp11, massp12, massp21, massp22, 
            C, Ct, Ctdot, u, Vp, Pp, Hp, Ppdot, Hpdot, ω, a, α, gvec)     
        # split the load between the two beam endpoints
        f1 += Fp/2
        f2 += Fp/2
        m1 += Mp/2
        m2 += Mp/2
    end

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
    θ, F, M, γ, κ, f1_θ, f2_θ, m1_θ, m2_θ, Ct_θ1, Ct_θ2, Ct_θ3)

    # --- f_u1, f_u2 --- #

    # d_fu/d_θ
    tmp = mul3(Ct_θ1, Ct_θ2, Ct_θ3, Cab*F)
    f_u1_θ = -tmp - f1_θ
    f_u2_θ =  tmp - f2_θ

    # d_fu/d_F
    f_u1_F = -CtCab
    f_u2_F =  CtCab

    # --- f_ψ1, f_ψ2 --- #

    # d_fψ/d_θ
    tmp1 = mul3(Ct_θ1, Ct_θ2, Ct_θ3, Cab*M)
    tmp2 = ΔL/2*mul3(Ct_θ1, Ct_θ2, Ct_θ3, Cab*cross(e1 + γ, F))
    f_ψ1_θ = -tmp1 - m1_θ - tmp2
    f_ψ2_θ =  tmp1 - m2_θ - tmp2

    # d_fψ/d_F
    tmp = -ΔL/2*CtCab*(tilde(e1 + γ) - tilde(F)*S11)
    f_ψ1_F = tmp
    f_ψ2_F = tmp

    # d_fψ/d_M
    tmp = ΔL/2*CtCab*tilde(F)*S12
    f_ψ1_M = tmp - CtCab
    f_ψ2_M = tmp + CtCab

    # --- f_F1, f_F2 --- #

    # d_fF/d_u
    f_F1_u =  I3
    f_F2_u = -I3

    # d_fF/d_θ
    tmp = mul3(Ct_θ1, Ct_θ2, Ct_θ3, ΔL/2*Cab*(e1 + γ))
    f_F1_θ = -tmp
    f_F2_θ = -tmp

    # d_fF/d_F
    tmp = ΔL/2*CtCab*S11
    f_F1_F = -tmp
    f_F2_F = -tmp

    # d_fF/d_M
    tmp = ΔL/2*CtCab*S12
    f_F1_M = -tmp
    f_F2_M = -tmp

    # --- f_M1, f_M2 --- #

    # d_fM/d_θ
    Qinv_θ1, Qinv_θ2, Qinv_θ3 = get_Qinv_θ(θ)
    tmp = mul3(Qinv_θ1, Qinv_θ2, Qinv_θ3, ΔL/2*Cab*κ)
    f_M1_θ =  I - tmp
    f_M2_θ = -I - tmp

    # d_fM/d_F
    Qinv = get_Qinv(θ)
    tmp1 = -ΔL/2*Qinv*Cab
    tmp2 = tmp1*S21
    f_M1_F = tmp2
    f_M2_F = tmp2

    # d_fM/d_M
    tmp2 = tmp1*S22
    f_M1_M = tmp2
    f_M2_M = tmp2

    return f_u1_θ, f_u2_θ, f_u1_F, f_u2_F,
        f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M,
        f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M
end

"""
    steady_state_element_jacobian_equations(ΔL, S11, S12, S21, S22, 
        mass11, mass12, mass21, mass22, Cab, CtCab, θ, F, M, γ, κ, V, P, H, ω, f1_u, f2_u, 
        m1_u, m2_u, f1_θ, f2_θ, m1_θ, m2_θ, f1_V, f2_V, m1_V, m2_V, f1_Ω, f2_Ω, m1_Ω, m2_Ω, 
        C_θ1, C_θ2, C_θ3, Ct_θ1, Ct_θ2, Ct_θ3)

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
 - `f1_V`: Derivative of `f1` w.r.t `V` 
 - `f2_V`: Derivative of `f2` w.r.t `V` 
 - `m1_V`: Derivative of `m1` w.r.t `V` 
 - `m2_V`: Derivative of `m2` w.r.t `V` 
 - `f1_Ω`: Derivative of `f1` w.r.t `Ω` 
 - `f2_Ω`: Derivative of `f2` w.r.t `Ω` 
 - `m1_Ω`: Derivative of `m1` w.r.t `Ω` 
 - `m2_Ω`: Derivative of `m2` w.r.t `Ω` 
 - `C_θ1`: Derivative of `C` w.r.t. `θ[1]`
 - `C_θ2`: Derivative of `C` w.r.t. `θ[2]`
 - `C_θ3`: Derivative of `C` w.r.t. `θ[3]`
 - `Ct_θ1`: Derivative of `C'` w.r.t. `θ[1]` (transpose of `C_θ1`)
 - `Ct_θ2`: Derivative of `C'` w.r.t. `θ[2]` (transpose of `C_θ2`)
 - `Ct_θ3`: Derivative of `C'` w.r.t. `θ[3]` (transpose of `C_θ3`)
"""
@inline function steady_state_element_jacobian_equations(ΔL, S11, S12, S21, S22, 
    mass11, mass12, mass21, mass22, Cab, CtCab, θ, F, M, γ, κ, V, P, H, ω, f1_u, f2_u, 
    m1_u, m2_u, f1_θ, f2_θ, m1_θ, m2_θ, f1_V, f2_V, m1_V, m2_V, f1_Ω, f2_Ω, m1_Ω, m2_Ω, 
    C_θ1, C_θ2, C_θ3, Ct_θ1, Ct_θ2, Ct_θ3)

    f_u1_θ, f_u2_θ, f_u1_F, f_u2_F,
        f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M,
        f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M = static_element_jacobian_equations(
        ΔL, S11, S12, S21, S22, Cab, CtCab, θ, F, M, γ, κ, f1_θ, f2_θ, m1_θ, m2_θ, Ct_θ1, 
        Ct_θ2, Ct_θ3)

    # --- f_u1, f_u2 --- #

    # d_fu/d_u
    f_u1_u = -f1_u
    f_u2_u = -f2_u

    # d_fu_dθ
    tmp = ΔL/2*tilde(ω)*mul3(Ct_θ1, Ct_θ2, Ct_θ3, Cab*P)
    f_u1_θ += tmp
    f_u2_θ += tmp

    # d_fu_dV
    tmp = ΔL/2*tilde(ω)*CtCab*mass11
    f_u1_V = -f1_V + tmp
    f_u2_V = -f2_V + tmp

    # d_fu_dΩ
    tmp = ΔL/2*tilde(ω)*CtCab*mass12
    f_u1_Ω = -f1_Ω + tmp
    f_u2_Ω = -f2_Ω + tmp

    # --- f_ψ1, f_ψ2 --- #

    # d_fψ_du
    f_ψ1_u = -m1_u
    f_ψ2_u = -m2_u

    # d_fψ_dθ
    tmp = ΔL/2*(tilde(ω)*mul3(Ct_θ1, Ct_θ2, Ct_θ3, Cab*H) + 
        mul3(Ct_θ1, Ct_θ2, Ct_θ3, Cab*cross(V, P)))
    f_ψ1_θ += tmp
    f_ψ2_θ += tmp

    # d_fψ_dV
    tmp = ΔL/2*tilde(ω)*CtCab*mass21 + ΔL/2*CtCab*(tilde(V)*mass11 - tilde(P))
    f_ψ1_V = -m1_V + tmp
    f_ψ2_V = -m2_V + tmp

    # d_fψ_dΩ
    tmp = ΔL/2*tilde(ω)*CtCab*mass22 + ΔL/2*CtCab*(tilde(V)*mass12)
    f_ψ1_Ω = -m1_Ω + tmp
    f_ψ2_Ω = -m2_Ω + tmp

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
        mass11, mass12, mass21, mass22, Cab, CtCab, CtCabdot, θ, F, γ, V, P, ω, f1_Vdot,
        f2_Vdot, m1_Vdot, m2_Vdot, f1_Ωdot, f2_Ωdot, m1_Ωdot, m2_Ωdot, f1_V, f2_V, m1_V, 
        m2_V, f1_Ω, f2_Ω, m1_Ω, m2_Ω)

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
 - `f1_Vdot`: Derivative of `f1` w.r.t `Vdot` 
 - `f2_Vdot`: Derivative of `f2` w.r.t `Vdot` 
 - `m1_Vdot`: Derivative of `m1` w.r.t `Vdot` 
 - `m2_Vdot`: Derivative of `m2` w.r.t `Vdot` 
 - `f1_Ωdot`: Derivative of `f1` w.r.t `Ωdot` 
 - `f2_Ωdot`: Derivative of `f2` w.r.t `Ωdot` 
 - `m1_Ωdot`: Derivative of `m1` w.r.t `Ωdot` 
 - `m2_Ωdot`: Derivative of `m2` w.r.t `Ωdot` 
 - `f1_V`: Derivative of `f1` w.r.t `V` 
 - `f2_V`: Derivative of `f2` w.r.t `V` 
 - `m1_V`: Derivative of `m1` w.r.t `V` 
 - `m2_V`: Derivative of `m2` w.r.t `V` 
 - `f1_Ω`: Derivative of `f1` w.r.t `Ω` 
 - `f2_Ω`: Derivative of `f2` w.r.t `Ω` 
 - `m1_Ω`: Derivative of `m1` w.r.t `Ω` 
 - `m2_Ω`: Derivative of `m2` w.r.t `Ω` 
"""
@inline function initial_condition_element_jacobian_equations(ΔL, S11, S12, S21, S22, 
    mass11, mass12, mass21, mass22, Cab, CtCab, CtCabdot, θ, F, γ, V, P, ω, f1_Vdot,
    f2_Vdot, m1_Vdot, m2_Vdot, f1_Ωdot, f2_Ωdot, m1_Ωdot, m2_Ωdot, f1_V, f2_V, m1_V, 
    m2_V, f1_Ω, f2_Ω, m1_Ω, m2_Ω)

    # --- f_u1, f_u2 --- #

    # d_fu/d_Vdot
    tmp = ΔL/2*CtCab*mass11
    f_u1_Vdot = -f1_Vdot + tmp
    f_u2_Vdot = -f2_Vdot + tmp

    # d_fu/d_Ωdot
    tmp = ΔL/2*CtCab*mass12
    f_u1_Ωdot = -f1_Ωdot + tmp
    f_u2_Ωdot = -f2_Ωdot + tmp

    # d_fu/d_F
    tmp = CtCab
    f_u1_F = -tmp
    f_u2_F =  tmp

    # d_fu/dV
    tmp = ΔL/2*(tilde(ω)*CtCab*mass11 + CtCabdot*mass11)
    f_u1_V = -f1_V + tmp 
    f_u2_V = -f2_V + tmp 

    # d_fu/dΩ
    tmp = ΔL/2*(tilde(ω)*CtCab*mass12 + CtCabdot*mass12)
    f_u1_Ω = -f1_Ω + tmp 
    f_u2_Ω = -f2_Ω + tmp 

    # --- f_θ1, f_θ2 --- #

    # d_fψ/d_Ωdot
    tmp = ΔL/2*CtCab*mass21
    f_ψ1_Vdot = -m1_Vdot + tmp 
    f_ψ2_Vdot = -m2_Vdot + tmp

    # d_fψ/d_Ωdot
    tmp = ΔL/2*CtCab*mass22
    f_ψ1_Ωdot = -m1_Ωdot + tmp
    f_ψ2_Ωdot = -m2_Ωdot + tmp

    # d_fψ/d_F
    tmp = -ΔL/2*CtCab*(tilde(e1 + γ) - tilde(F)*S11)
    f_ψ1_F = tmp
    f_ψ2_F = tmp

    # d_fψ/d_M
    tmp = ΔL/2*CtCab*tilde(F)*S12
    f_ψ1_M = tmp - CtCab
    f_ψ2_M = tmp + CtCab

    # d_fψ_dV
    tmp = ΔL/2*tilde(ω)*CtCab*mass21 + ΔL/2*CtCabdot*mass21 + ΔL/2*CtCab*(tilde(V)*mass11 - tilde(P))
    f_ψ1_V = -m1_V + tmp 
    f_ψ2_V = -m2_V + tmp 

    # d_fψ_dΩ
    tmp = ΔL/2*tilde(ω)*CtCab*mass22 + ΔL/2*CtCabdot*mass22 + ΔL/2*CtCab*(tilde(V)*mass12)
    f_ψ1_Ω = -m1_Ω + tmp
    f_ψ2_Ω = -m2_Ω + tmp

    # --- f_F1, f_F2 --- #

    # d_fF/d_F
    tmp = ΔL/2*CtCab*S11
    f_F1_F = -tmp
    f_F2_F = -tmp

    # d_fF/d_M
    tmp = ΔL/2*CtCab*S12
    f_F1_M = -tmp
    f_F2_M = -tmp

    # --- f_M1, f_M2 --- #

    # d_fM/d_F
    Qinv = get_Qinv(θ)
    tmp1 = -ΔL/2*Qinv*Cab
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
        mass21, mass22, Cab, CtCab, CtCabdot, θ, F, M, γ, κ, V, P, H, θdot, Pdot, Hdot, ω, 
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
 - `f1_u`: Derivative of `f1` w.r.t `u` 
 - `f2_u`: Derivative of `f2` w.r.t `u` 
 - `m1_u`: Derivative of `m1` w.r.t `u` 
 - `m2_u`: Derivative of `m2` w.r.t `u` 
 - `f1_θ`: Derivative of `f1` w.r.t `θ` 
 - `f2_θ`: Derivative of `f2` w.r.t `θ` 
 - `m1_θ`: Derivative of `m1` w.r.t `θ` 
 - `m2_θ`: Derivative of `m2` w.r.t `θ` 
 - `f1_V`: Derivative of `f1` w.r.t `V` 
 - `f2_V`: Derivative of `f2` w.r.t `V` 
 - `m1_V`: Derivative of `m1` w.r.t `V` 
 - `m2_V`: Derivative of `m2` w.r.t `V` 
 - `f1_Ω`: Derivative of `f1` w.r.t `Ω` 
 - `f2_Ω`: Derivative of `f2` w.r.t `Ω` 
 - `m1_Ω`: Derivative of `m1` w.r.t `Ω` 
 - `m2_Ω`: Derivative of `m2` w.r.t `Ω` 
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
    mass21, mass22, Cab, CtCab, CtCabdot, θ, F, M, γ, κ, V, P, H, θdot, Pdot, Hdot, ω, dt,
    f1_u, f2_u, m1_u, m2_u, f1_θ, f2_θ, m1_θ, m2_θ, f1_V, f2_V, m1_V, m2_V, 
    f1_Ω, f2_Ω, m1_Ω, m2_Ω, C_θ1, C_θ2, C_θ3, Ct_θ1, Ct_θ2, Ct_θ3, 
    Ctdot_θ1, Ctdot_θ2, Ctdot_θ3, Ctdot_θdot1, Ctdot_θdot2, Ctdot_θdot3)

    f_u1_u, f_u2_u, f_u1_θ, f_u2_θ, f_u1_F, f_u2_F, f_u1_V, f_u2_V, f_u1_Ω, f_u2_Ω,
        f_ψ1_u, f_ψ2_u, f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_V, f_ψ2_V, f_ψ1_Ω, f_ψ2_Ω,
        f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M,
        f_V_u, f_V_θ, f_V_V,
        f_Ω_θ, f_Ω_Ω = steady_state_element_jacobian_equations(ΔL, S11, S12, S21, S22, 
        mass11, mass12, mass21, mass22, Cab, CtCab, θ, F, M, γ, κ, V, P, H, ω, f1_u, f2_u, 
        m1_u, m2_u, f1_θ, f2_θ, m1_θ, m2_θ, f1_V, f2_V, m1_V, m2_V, f1_Ω, f2_Ω, m1_Ω, m2_Ω, 
        C_θ1, C_θ2, C_θ3, Ct_θ1, Ct_θ2, Ct_θ3)

    # --- f_u1, f_u2 --- #

    # d_fu_dθ
    tmp = ΔL/2*(
        mul3(Ctdot_θ1, Ctdot_θ2, Ctdot_θ3, Cab*P) +
        2/dt*mul3(Ctdot_θdot1, Ctdot_θdot2, Ctdot_θdot3, Cab*P) + 
        mul3(Ct_θ1, Ct_θ2, Ct_θ3, Cab*Pdot)
    )
    f_u1_θ += tmp
    f_u2_θ += tmp

    # d_fu_dV      
    tmp = ΔL/2*(CtCabdot*mass11 + 2/dt*CtCab*mass11)
    f_u1_V += tmp
    f_u2_V += tmp

    # d_fu_dΩ
    tmp = ΔL/2*(CtCabdot*mass12 + 2/dt*CtCab*mass12)
    f_u1_Ω += tmp
    f_u2_Ω += tmp

    # --- f_ψ1, f_ψ2 --- #

    # d_fψ_dθ
    tmp = ΔL/2*(
        mul3(Ctdot_θ1, Ctdot_θ2, Ctdot_θ3, Cab*H) +
        2/dt*mul3(Ctdot_θdot1, Ctdot_θdot2, Ctdot_θdot3, Cab*H) +
        mul3(Ct_θ1, Ct_θ2, Ct_θ3, Cab*Hdot)
    )
    f_ψ1_θ += tmp
    f_ψ2_θ += tmp

    # d_fψ_dV
    tmp = ΔL/2*(CtCabdot*mass21 + 2/dt*CtCab*mass21)
    f_ψ1_V += tmp
    f_ψ2_V += tmp

    # d_fψ_dΩ
    tmp = ΔL/2*(CtCabdot*mass22 + 2/dt*CtCab*mass22)
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
        f1_u, f2_u, m1_u, m2_u, f1_θ, f2_θ, m1_θ, m2_θ, f1_V, f2_V, m1_V, m2_V, f1_Ω, f2_Ω, 
        m1_Ω, m2_Ω, C_θ1, C_θ2, C_θ3, Ct_θ1, Ct_θ2, Ct_θ3, Ctdot_θ1, Ctdot_θ2, Ctdot_θ3)

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
 - `f1_V`: Derivative of `f1` w.r.t `V` 
 - `f2_V`: Derivative of `f2` w.r.t `V` 
 - `m1_V`: Derivative of `m1` w.r.t `V` 
 - `m2_V`: Derivative of `m2` w.r.t `V` 
 - `f1_Ω`: Derivative of `f1` w.r.t `Ω` 
 - `f2_Ω`: Derivative of `f2` w.r.t `Ω` 
 - `m1_Ω`: Derivative of `m1` w.r.t `Ω` 
 - `m2_Ω`: Derivative of `m2` w.r.t `Ω` 
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
    f1_u, f2_u, m1_u, m2_u, f1_θ, f2_θ, m1_θ, m2_θ, f1_V, f2_V, m1_V, m2_V, f1_Ω, f2_Ω, 
    m1_Ω, m2_Ω, C_θ1, C_θ2, C_θ3, Ct_θ1, Ct_θ2, Ct_θ3, Ctdot_θ1, Ctdot_θ2, Ctdot_θ3)

    f_u1_u, f_u2_u, f_u1_θ, f_u2_θ, f_u1_F, f_u2_F, f_u1_V, f_u2_V, f_u1_Ω, f_u2_Ω,
        f_ψ1_u, f_ψ2_u, f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_V, f_ψ2_V, f_ψ1_Ω, f_ψ2_Ω,
        f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M,
        f_V_u, f_V_θ, f_V_V,
        f_Ω_θ, f_Ω_Ω = steady_state_element_jacobian_equations(ΔL, S11, S12, S21, S22, 
        mass11, mass12, mass21, mass22, Cab, CtCab, θ, F, M, γ, κ, V, P, H, ω, f1_u, f2_u, 
        m1_u, m2_u, f1_θ, f2_θ, m1_θ, m2_θ, f1_V, f2_V, m1_V, m2_V, f1_Ω, f2_Ω, m1_Ω, m2_Ω, 
        C_θ1, C_θ2, C_θ3, Ct_θ1, Ct_θ2, Ct_θ3)

    # --- f_u1, f_u2 --- #

    # d_fu_dθ
    tmp = ΔL/2*(mul3(Ctdot_θ1, Ctdot_θ2, Ctdot_θ3, Cab*P) + mul3(Ct_θ1, Ct_θ2, Ct_θ3, Cab*Pdot))
    f_u1_θ += tmp
    f_u2_θ += tmp

    # d_fu_dV
    tmp = ΔL/2*CtCabdot*mass11
    f_u1_V += tmp
    f_u2_V += tmp

    # d_fu_dΩ
    tmp = ΔL/2*CtCabdot*mass12
    f_u1_Ω += tmp
    f_u2_Ω += tmp

    # --- f_ψ1, f_ψ2 --- #

    # d_fψ_dθ
    tmp = ΔL/2*(mul3(Ctdot_θ1, Ctdot_θ2, Ctdot_θ3, Cab*H) + mul3(Ct_θ1, Ct_θ2, Ct_θ3, Cab*Hdot))
    f_ψ1_θ += tmp
    f_ψ2_θ += tmp

    # d_fψ_dV
    tmp = ΔL/2*CtCabdot*mass21
    f_ψ1_V += tmp
    f_ψ2_V += tmp

    # d_fψ_dΩ
    tmp = ΔL/2*CtCabdot*mass22
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
        f_u1_θ, f_u2_θ, f_u1_F, f_u2_F,
        f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M,
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
    f_u1_θ, f_u2_θ, f_u1_F, f_u2_F,
    f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M,
    f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
    f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M)

    # equilibrium equation jacobian entries for the start of the beam
    jacob[irow_p1:irow_p1+2, icol+3:icol+5] .= f_u1_θ ./ force_scaling
    jacob[irow_p1:irow_p1+2, icol+6:icol+8] .= f_u1_F

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
    jacob[irow_p2:irow_p2+2, icol+3:icol+5] .= f_u2_θ ./ force_scaling
    jacob[irow_p2:irow_p2+2, icol+6:icol+8] .= f_u2_F

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
        f_ψ1_u, f_ψ2_u, f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_V, f_ψ2_V, f_ψ1_Ω, f_ψ2_Ω,
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
    f_ψ1_u, f_ψ2_u, f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_V, f_ψ2_V, f_ψ1_Ω, f_ψ2_Ω,
    f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
    f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M,
    f_V_u, f_V_θ, f_V_V,
    f_Ω_θ, f_Ω_Ω)

    jacob = static_insert_element_jacobian!(jacob, force_scaling, irow_e1, irow_p1,
        irow_e2, irow_p2, icol,
        f_u1_θ, f_u2_θ, f_u1_F, f_u2_F,
        f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M,
        f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M)

    # add equilibrium equation jacobian entries for the start of the beam
    jacob[irow_p1:irow_p1+2, icol:icol+2] .= f_u1_u ./ force_scaling
    jacob[irow_p1:irow_p1+2, icol+12:icol+14] .= f_u1_V  ./ force_scaling
    jacob[irow_p1:irow_p1+2, icol+15:icol+17] .= f_u1_Ω  ./ force_scaling
    jacob[irow_p1+3:irow_p1+5, icol:icol+2] .= f_ψ1_u ./ force_scaling
    jacob[irow_p1+3:irow_p1+5, icol+12:icol+14] .= f_ψ1_V  ./ force_scaling
    jacob[irow_p1+3:irow_p1+5, icol+15:icol+17] .= f_ψ1_Ω  ./ force_scaling

    # add equilibrium equation jacobian entries for the end of the beam
    jacob[irow_p2:irow_p2+2, icol:icol+2] .= f_u2_u ./ force_scaling
    jacob[irow_p2:irow_p2+2, icol+12:icol+14] .= f_u2_V  ./ force_scaling
    jacob[irow_p2:irow_p2+2, icol+15:icol+17] .= f_u2_Ω  ./ force_scaling
    jacob[irow_p2+3:irow_p2+5, icol:icol+2] .= f_ψ2_u  ./ force_scaling
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

    # compute element properties
    ΔL, S11, S12, S21, S22, mass11, mass12, mass21, mass22, C, Cab, CtCab, u, θ, 
        F, M, γ, κ = static_element_properties(x, icol, elem, force_scaling)

    # additional precomputed quantities
    Ct = C'
    C_θ1, C_θ2, C_θ3 = get_C_θ(C, θ)
    Ct_θ1, Ct_θ2, Ct_θ3 = C_θ1', C_θ2', C_θ3'

    # gravitational load jacobians
    f1_θ, f2_θ, m1_θ, m2_θ = element_gravitational_loads_jacobian(ΔL, mass11, mass12, Cab, 
        CtCab, gvec, C_θ1, C_θ2, C_θ3, Ct_θ1, Ct_θ2, Ct_θ3)

    # follower load jacobians
    if haskey(distributed_loads, ielem)
        dload = distributed_loads[ielem]
        f1_θ += mul3(Ct_θ1, Ct_θ2, Ct_θ3, dload.f1_follower)
        f2_θ += mul3(Ct_θ1, Ct_θ2, Ct_θ3, dload.f2_follower)
        m1_θ += mul3(Ct_θ1, Ct_θ2, Ct_θ3, dload.m1_follower)
        m2_θ += mul3(Ct_θ1, Ct_θ2, Ct_θ3, dload.m2_follower)
    end

    # point mass load jacobians
    if haskey(point_masses, ielem)
        # get point mass properties
        massp11, massp12, massp21, massp22 = static_point_mass_properties(point_masses[ielem])
        # get point mass load jacobians
        Fp_θ, Mp_θ = static_point_mass_jacobian(massp11, massp12, C, Ct, gvec, 
            C_θ1, C_θ2, C_θ3, Ct_θ1, Ct_θ2, Ct_θ3)
        # split the load between the two beam endpoints
        f1_θ += Fp_θ/2
        f2_θ += Fp_θ/2
        m1_θ += Mp_θ/2
        m2_θ += Mp_θ/2
    end

    # element resultant jacobians
    f_u1_θ, f_u2_θ, f_u1_F, f_u2_F,
        f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M,
        f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M = static_element_jacobian_equations(
        ΔL, S11, S12, S21, S22, Cab, CtCab, θ, F, M, γ, κ, f1_θ, f2_θ, m1_θ, m2_θ, 
        Ct_θ1, Ct_θ2, Ct_θ3)

    # insert element resultants into the jacobian matrix
    static_insert_element_jacobian!(jacob, force_scaling, irow_e1,
        irow_p1, irow_e2, irow_p2, icol,
        f_u1_θ, f_u2_θ, f_u1_F, f_u2_F,
        f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M,
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
@inline function steady_state_element_jacobian!(jacob, x, ielem, elem, distributed_loads, point_masses, gvec,
    force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2, x0, v0, ω0, a0, α0)

    # compute element properties
    ΔL, S11, S12, S21, S22, mass11, mass12, mass21, mass22, C, Cab, CtCab, u, θ, F, M, 
        γ, κ, V, Ω, P, H, v, ω, a, α = steady_state_element_properties(x, icol, elem, 
        force_scaling, x0, v0, ω0, a0, α0)

    # additional precomputed quantities
    Ct = C'
    C_θ1, C_θ2, C_θ3 = get_C_θ(C, θ)
    Ct_θ1, Ct_θ2, Ct_θ3 = C_θ1', C_θ2', C_θ3'

    # element acceleration loads (including gravitational loads)
    f1_u, f2_u, m1_u, m2_u, f1_θ, f2_θ, m1_θ, m2_θ = element_acceleration_loads_jacobian(ΔL, 
        mass11, mass12, mass21, mass22, Cab, CtCab, u, a, α, gvec, C_θ1, C_θ2, C_θ3, 
        Ct_θ1, Ct_θ2, Ct_θ3)

    # initialize uninitialized load jacobians
    f1_V = f2_V = m1_V = m2_V = @SMatrix zeros(3,3)
    f1_Ω = f2_Ω = m1_Ω = m2_Ω = @SMatrix zeros(3,3)

    # distributed loads
    if haskey(distributed_loads, ielem)
        dload = distributed_loads[ielem]
        f1_θ += mul3(Ct_θ1, Ct_θ2, Ct_θ3, dload.f1_follower)
        f2_θ += mul3(Ct_θ1, Ct_θ2, Ct_θ3, dload.f2_follower)
        m1_θ += mul3(Ct_θ1, Ct_θ2, Ct_θ3, dload.m1_follower)
        m2_θ += mul3(Ct_θ1, Ct_θ2, Ct_θ3, dload.m2_follower)
    end

    # point mass load jacobians
    if haskey(point_masses, ielem)
        # point mass properties
        massp11, massp12, massp21, massp22, Vp, Ωp, Pp, Hp = 
            steady_state_point_mass_properties(point_masses[ielem], Cab, V, Ω)
        # point mass load jacobians
        Fp_u, Mp_u, Fp_θ, Mp_θ, Fp_V, Mp_V, Fp_Ω, Mp_Ω = steady_state_point_mass_jacobian(massp11, 
            massp12, massp21, massp22, C, Ct, Cab, u, Vp, Pp, Hp, ω, a, α, gvec, C_θ1, C_θ2, 
            C_θ3, Ct_θ1, Ct_θ2, Ct_θ3)
        # split loads between the two beam element endpoints
        f1_u += Fp_u/2
        f2_u += Fp_u/2
        m1_u += Mp_u/2
        m2_u += Mp_u/2
        f1_θ += Fp_θ/2
        f2_θ += Fp_θ/2
        m1_θ += Mp_θ/2
        m2_θ += Mp_θ/2
        f1_V += Fp_V/2
        f2_V += Fp_V/2
        m1_V += Mp_V/2
        m2_V += Mp_V/2
        f1_Ω += Fp_Ω/2
        f2_Ω += Fp_Ω/2
        m1_Ω += Mp_Ω/2
        m2_Ω += Mp_Ω/2
    end

    # solve for the element resultant jacobians
    f_u1_u, f_u2_u, f_u1_θ, f_u2_θ, f_u1_F, f_u2_F, f_u1_V, f_u2_V, f_u1_Ω, f_u2_Ω,
        f_ψ1_u, f_ψ2_u, f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_V, f_ψ2_V, f_ψ1_Ω, f_ψ2_Ω,
        f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M,
        f_V_u, f_V_θ, f_V_V,
        f_Ω_θ, f_Ω_Ω = steady_state_element_jacobian_equations(ΔL, S11, S12, S21, S22, 
        mass11, mass12, mass21, mass22, Cab, CtCab, θ, F, M, γ, κ, V, P, H, ω, f1_u, f2_u, 
        m1_u, m2_u, f1_θ, f2_θ, m1_θ, m2_θ, f1_V, f2_V, m1_V, m2_V, f1_Ω, f2_Ω, m1_Ω, m2_Ω, 
        C_θ1, C_θ2, C_θ3, Ct_θ1, Ct_θ2, Ct_θ3)

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
    initial_condition_element_jacobian!(jacob, x, ielem, elem, distributed_loads, 
        point_masses, gvec, force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, 
        irow_p2, x0, v0, ω0, a0, α0, u0, θ0, udot0, θdot0)

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
 - `u0`: initial deflection variables for the beam element
 - `θ0`: initial rotation variables for the beam element
 - `udot0`: initial time derivative of u for the beam element
 - `θdot0`: initial time derivative of θ for the beam element
"""
@inline function initial_condition_element_jacobian!(jacob, x, ielem, elem, distributed_loads, point_masses, gvec,
    force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2, x0, v0, ω0, a0, α0, u0, θ0,
    udot0, θdot0)

    # compute element properties
    ΔL, S11, S12, S21, S22, mass11, mass12, mass21, mass22, C, Cab, CtCab, u, θ, 
        F, M, γ, κ, V, Ω, P, H, v, ω, a, α, udot, θdot, Vdot, Ωdot, Pdot, Hdot,
        Cdot, CtCabdot, CtCabPdot, CtCabHdot = initial_condition_element_properties(x, icol, 
        elem, force_scaling, x0, v0, ω0, a0, α0, u0, θ0, udot0, θdot0)

    # additional precomputed quantities
    Ct = C'
    Ctdot = Cdot'

    # initialize load jacobians
    f1_Vdot = f2_Vdot = m1_Vdot = m2_Vdot = @SMatrix zeros(3,3)
    f1_Ωdot = f2_Ωdot = m1_Ωdot = m2_Ωdot = @SMatrix zeros(3,3)
    f1_V = f2_V = m1_V = m2_V = @SMatrix zeros(3,3)
    f1_Ω = f2_Ω = m1_Ω = m2_Ω = @SMatrix zeros(3,3)
    
    # point mass load jacobians
    if haskey(point_masses, ielem)
        # get point mass properties
        massp11, massp12, massp21, massp22, Vp, Ωp, Pp, Hp, Vpdot, Ωpdot, Ppdot, Hpdot, 
            CtPpdot, CtHpdot = dynamic_point_mass_properties(point_masses[ielem], Ct, Ctdot, 
            Cab, V, Ω, Vdot, Ωdot)
        # get point mass load jacobians
        Fp_Vdot, Fp_Ωdot, Fp_V, Fp_Ω, Mp_Vdot, Mp_Ωdot, Mp_V, Mp_Ω = 
            initial_condition_point_mass_jacobian(massp11, massp12, massp21, 
            massp22, Ct, Ctdot, Cab, Vp, Pp, ω)
        # split loads between the two beam element endpoints
        f1_Vdot += Fp_Vdot/2
        f2_Vdot += Fp_Vdot/2
        m1_Vdot += Mp_Vdot/2
        m2_Vdot += Mp_Vdot/2
        f1_Ωdot += Fp_Ωdot/2
        f2_Ωdot += Fp_Ωdot/2
        m1_Ωdot += Mp_Ωdot/2
        m2_Ωdot += Mp_Ωdot/2
        f1_V += Fp_V/2
        f2_V += Fp_V/2
        m1_V += Mp_V/2
        m2_V += Mp_V/2
        f1_Ω += Fp_Ω/2
        f2_Ω += Fp_Ω/2
        m1_Ω += Mp_Ω/2
        m2_Ω += Mp_Ω/2
    end

    # solve for the element resultant jacobians
    f_u1_Vdot, f_u2_Vdot, f_u1_Ωdot, f_u2_Ωdot, f_u1_F, f_u2_F, f_u1_V, f_u2_V, f_u1_Ω, f_u2_Ω,
        f_ψ1_Vdot, f_ψ2_Vdot, f_ψ1_Ωdot, f_ψ2_Ωdot, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_V, f_ψ2_V, f_ψ1_Ω, f_ψ2_Ω,
        f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_F, f_M2_F, f_M1_M, f_M2_M,
        f_V_V, f_Ω_Ω = initial_condition_element_jacobian_equations(ΔL, S11, S12, S21, S22, 
        mass11, mass12, mass21, mass22, Cab, CtCab, CtCabdot, θ, F, γ, V, P, ω, f1_Vdot,
        f2_Vdot, m1_Vdot, m2_Vdot, f1_Ωdot, f2_Ωdot, m1_Ωdot, m2_Ωdot, f1_V, f2_V, m1_V, 
        m2_V, f1_Ω, f2_Ω, m1_Ω, m2_Ω)

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
        udot_init, θdot_init, Vdot_init, Ωdot_init, dt)

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
 - `Vdot_init`: `2/dt*V + Vdot` for the beam element from the previous time step
 - `Ωdot_init`: `2/dt*Ω + Ωdot` for the beam element from the previous time step
 - `dt`: time step size
"""
@inline function newmark_element_jacobian!(jacob, x, ielem, elem, distributed_loads, point_masses, gvec,
    force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2, x0, v0, ω0, a0, α0,
    udot_init, θdot_init, Vdot_init, Ωdot_init, dt)

    # compute element properties
    ΔL, S11, S12, S21, S22, mass11, mass12, mass21, mass22, C, Cab, CtCab, u, θ, 
        F, M, γ, κ, V, Ω, P, H, v, ω, a, α, udot, θdot, Vdot, Ωdot, Pdot, Hdot, 
        Cdot, CtCabdot, CtCabPdot, CtCabHdot = newmark_element_properties(x, icol, elem, 
        force_scaling, x0, v0, ω0, a0, α0, udot_init, θdot_init, Vdot_init, Ωdot_init, dt)

    # additional precomputed quantities
    Ct = C'
    Ctdot = Cdot'

    # jacobian of rotation matrix wrt θ
    C_θ1, C_θ2, C_θ3 = get_C_θ(C, θ)
    Ct_θ1, Ct_θ2, Ct_θ3 = C_θ1', C_θ2', C_θ3'

    # jacobian of time derivative of rotation matrix wrt θ
    Cdot_θ1, Cdot_θ2, Cdot_θ3 = get_C_t_θ(θ, θdot)
    Ctdot_θ1, Ctdot_θ2, Ctdot_θ3 = Cdot_θ1', Cdot_θ2', Cdot_θ3'

    # jacobian of time derivative of rotation matrix wrt θdot
    Cdot_θdot1, Cdot_θdot2, Cdot_θdot3 = get_C_t_θdot(C, θ)
    Ctdot_θdot1, Ctdot_θdot2, Ctdot_θdot3 = Cdot_θdot1', Cdot_θdot2', Cdot_θdot3'

    # element acceleration loads (including gravitational loads)
    f1_u, f2_u, m1_u, m2_u, f1_θ, f2_θ, m1_θ, m2_θ = element_acceleration_loads_jacobian(ΔL, 
        mass11, mass12, mass21, mass22, Cab, CtCab, u, a, α, gvec, C_θ1, C_θ2, C_θ3, 
        Ct_θ1, Ct_θ2, Ct_θ3)

    # initialize remaining load jacobians
    f1_V = f2_V = m1_V = m2_V = @SMatrix zeros(3,3)
    f1_Ω = f2_Ω = m1_Ω = m2_Ω = @SMatrix zeros(3,3)

    # distributed loads
    if haskey(distributed_loads, ielem)
        dload = distributed_loads[ielem]
        f1_θ += mul3(Ct_θ1, Ct_θ2, Ct_θ3, dload.f1_follower)
        f2_θ += mul3(Ct_θ1, Ct_θ2, Ct_θ3, dload.f2_follower)
        m1_θ += mul3(Ct_θ1, Ct_θ2, Ct_θ3, dload.m1_follower)
        m2_θ += mul3(Ct_θ1, Ct_θ2, Ct_θ3, dload.m2_follower)
    end

    # point mass load jacobians
    if haskey(point_masses, ielem)
        # point mass properties
        massp11, massp12, massp21, massp22, Vp, Ωp, Pp, Hp, Vpdot, Ωpdot, Ppdot, Hpdot, 
            CtPpdot, CtHpdot = dynamic_point_mass_properties(point_masses[ielem], Ct, Ctdot, 
            Cab, V, Ω, Vdot, Ωdot)
        # point mass load jacobians
        Fp_u, Mp_u, Fp_θ, Mp_θ, Fp_V, Mp_V, Fp_Ω, Mp_Ω = newmark_point_mass_jacobian(massp11, massp12, 
            massp21, massp22, u, Vp, Pp, Hp, Ppdot, Hpdot, C, Ct, Ctdot, Cab, ω, a, α, gvec, 
            dt, C_θ1, C_θ2, C_θ3, Ct_θ1, Ct_θ2, Ct_θ3, Ctdot_θ1, Ctdot_θ2, Ctdot_θ3, 
            Ctdot_θdot1, Ctdot_θdot2, Ctdot_θdot3)
        # split loads between the two beam element endpoints
        f1_u += Fp_u/2
        f2_u += Fp_u/2
        m1_u += Mp_u/2
        m2_u += Mp_u/2
        f1_θ += Fp_θ/2
        f2_θ += Fp_θ/2
        m1_θ += Mp_θ/2
        m2_θ += Mp_θ/2
        f1_V += Fp_V/2
        f2_V += Fp_V/2
        m1_V += Mp_V/2
        m2_V += Mp_V/2
        f1_Ω += Fp_Ω/2
        f2_Ω += Fp_Ω/2
        m1_Ω += Mp_Ω/2
        m2_Ω += Mp_Ω/2
    end

    # element resultant jacobians
    f_u1_u, f_u2_u, f_u1_θ, f_u2_θ, f_u1_F, f_u2_F, f_u1_V, f_u2_V, f_u1_Ω, f_u2_Ω,
        f_ψ1_u, f_ψ2_u, f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_V, f_ψ2_V, f_ψ1_Ω, f_ψ2_Ω,
        f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M,
        f_V_u, f_V_θ, f_V_V,
        f_Ω_θ, f_Ω_Ω = newmark_element_jacobian_equations(ΔL, S11, S12, S21, S22, 
        mass11, mass12, mass21, mass22, Cab, CtCab, CtCabdot, θ, F, M, γ, κ, V, P, H, 
        θdot, Pdot, Hdot, ω, dt, f1_u, f2_u, m1_u, m2_u, f1_θ, f2_θ, m1_θ, m2_θ, f1_V, f2_V, 
        m1_V, m2_V, f1_Ω, f2_Ω, m1_Ω, m2_Ω, C_θ1, C_θ2, C_θ3, Ct_θ1, Ct_θ2, Ct_θ3, 
        Ctdot_θ1, Ctdot_θ2, Ctdot_θ3, Ctdot_θdot1, Ctdot_θdot2, Ctdot_θdot3)

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
        udot, θdot, Vdot, Ωdot)

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
 - `Vdot_init`: `2/dt*V + Vdot` for the beam element from the previous time step
 - `Ωdot_init`: `2/dt*Ω + Ωdot` for the beam element from the previous time step
 - `dt`: time step size
"""
@inline function dynamic_element_jacobian!(jacob, x, ielem, elem, distributed_loads, point_masses, gvec,
    force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2, x0, v0, ω0, a0, α0, udot, θdot,
    Vdot, Ωdot)

    # compute element properties
    ΔL, S11, S12, S21, S22, mass11, mass12, mass21, mass22, C, Cab, CtCab, u, θ, 
        F, M, γ, κ, V, Ω, P, H, v, ω, a, α, udot, θdot, Pdot, Hdot, Cdot, CtCabdot,
        CtCabPdot, CtCabHdot = dynamic_element_properties(x, icol, elem, force_scaling, 
        x0, v0, ω0, a0, α0, udot, θdot, Vdot, Ωdot)

    # additional precomputed quantities
    Ct = C'
    Ctdot = Cdot'    
    C_θ1, C_θ2, C_θ3 = get_C_θ(C, θ) # jacobian of rotation matrix wrt θ
    Ct_θ1, Ct_θ2, Ct_θ3 = C_θ1', C_θ2', C_θ3'
    Cdot_θ1, Cdot_θ2, Cdot_θ3 = get_C_t_θ(θ, θdot) # jacobian of time derivative of rotation matrix wrt θ
    Ctdot_θ1, Ctdot_θ2, Ctdot_θ3 = Cdot_θ1', Cdot_θ2', Cdot_θ3'

    # element acceleration loads (including gravitational loads)
    f1_u, f2_u, m1_u, m2_u, f1_θ, f2_θ, m1_θ, m2_θ = element_acceleration_loads_jacobian(ΔL, 
        mass11, mass12, mass21, mass22, Cab, CtCab, u, a, α, gvec, C_θ1, C_θ2, C_θ3, 
        Ct_θ1, Ct_θ2, Ct_θ3)

    # initialize remaining load jacobians
    f1_V = f2_V = m1_V = m2_V = @SMatrix zeros(3,3)
    f1_Ω = f2_Ω = m1_Ω = m2_Ω = @SMatrix zeros(3,3)

    # distributed loads
    if haskey(distributed_loads, ielem)
        dload = distributed_loads[ielem]
        f1_θ += mul3(Ct_θ1, Ct_θ2, Ct_θ3, dload.f1_follower)
        f2_θ += mul3(Ct_θ1, Ct_θ2, Ct_θ3, dload.f2_follower)
        m1_θ += mul3(Ct_θ1, Ct_θ2, Ct_θ3, dload.m1_follower)
        m2_θ += mul3(Ct_θ1, Ct_θ2, Ct_θ3, dload.m2_follower)
    end

    # point mass load jacobians
    if haskey(point_masses, ielem)
        # point mass properties
        massp11, massp12, massp21, massp22, Vp, Ωp, Pp, Hp, Vpdot, Ωpdot, Ppdot, Hpdot, 
            CtPpdot, CtHpdot = dynamic_point_mass_properties(point_masses[ielem], Ct, Ctdot, 
            Cab, V, Ω, Vdot, Ωdot)
        # point mass load jacobians
        Fp_u, Mp_u, Fp_θ, Mp_θ, Fp_V, Mp_V, Fp_Ω, Mp_Ω = dynamic_point_mass_jacobian(massp11, massp12, 
            massp21, massp22, C, Ct, Ctdot, Cab, u, Vp, Pp, Hp, Ppdot, Hpdot, ω, a, α, gvec, 
            C_θ1, C_θ2, C_θ3, Ct_θ1, Ct_θ2, Ct_θ3, Ctdot_θ1, Ctdot_θ2, Ctdot_θ3)
        # split loads between the two beam element endpoints
        f1_u += Fp_u/2
        f2_u += Fp_u/2
        m1_u += Mp_u/2
        m2_u += Mp_u/2
        f1_θ += Fp_θ/2
        f2_θ += Fp_θ/2
        m1_θ += Mp_θ/2
        m2_θ += Mp_θ/2
        f1_V += Fp_V/2
        f2_V += Fp_V/2
        m1_V += Mp_V/2
        m2_V += Mp_V/2
        f1_Ω += Fp_Ω/2
        f2_Ω += Fp_Ω/2
        m1_Ω += Mp_Ω/2
        m2_Ω += Mp_Ω/2
    end

    # solve for the element resultants
    f_u1_u, f_u2_u, f_u1_θ, f_u2_θ, f_u1_F, f_u2_F, f_u1_V, f_u2_V, f_u1_Ω, f_u2_Ω,
        f_ψ1_u, f_ψ2_u, f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_V, f_ψ2_V, f_ψ1_Ω, f_ψ2_Ω,
        f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M,
        f_V_u, f_V_θ, f_V_V,
        f_Ω_θ, f_Ω_Ω = dynamic_element_jacobian_equations(ΔL, S11, S12, S21, S22, mass11, 
        mass12, mass21, mass22, Cab, CtCab, CtCabdot, θ, F, M, γ, κ, V, P, H, θdot, Pdot, 
        Hdot, ω, f1_u, f2_u, m1_u, m2_u, f1_θ, f2_θ, m1_θ, m2_θ, f1_V, f2_V, 
        m1_V, m2_V, f1_Ω, f2_Ω, m1_Ω, m2_Ω, C_θ1, C_θ2, C_θ3, Ct_θ1, Ct_θ2, Ct_θ3, 
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
    element_mass_matrix_properties(x, icol, elem)

Extract/Compute the properties needed for mass matrix construction: `ΔL`, `C`,
`Cab`, `CtCab`, `θ`, `P`, `H`, `Ctdot_θdot1`, `Ctdot_θdot2`, and Ctdot_θdot3
"""
@inline function element_mass_matrix_properties(x, icol, elem)

    ΔL = elem.L
    mass = elem.mass
    θ = SVector(x[icol+3 ], x[icol+4 ], x[icol+5 ])
    V = SVector(x[icol+12], x[icol+13], x[icol+14])
    Ω = SVector(x[icol+15], x[icol+16], x[icol+17])
    P = element_linear_momentum(elem, V, Ω)
    H = element_angular_momentum(elem, V, Ω)
    C = get_C(θ)
    Cab = elem.Cab
    CtCab = C'*Cab
    Cdot_cdot1, Cdot_cdot2, Cdot_cdot3 = get_C_t_θdot(C, θ)
    Ctdot_θdot1, Ctdot_θdot2, Ctdot_θdot3 = Cdot_cdot1', Cdot_cdot2', Cdot_cdot3'

    mass11 = mass[SVector{3}(1:3), SVector{3}(1:3)]
    mass12 = mass[SVector{3}(1:3), SVector{3}(4:6)]
    mass21 = mass[SVector{3}(4:6), SVector{3}(1:3)]
    mass22 = mass[SVector{3}(4:6), SVector{3}(4:6)]

    return ΔL, mass11, mass12, mass21, mass22, C, Cab, CtCab, θ, V, Ω, P, H, 
        Ctdot_θdot1, Ctdot_θdot2, Ctdot_θdot3
end

"""
    element_mass_matrix_equations(ΔL, mass11, mass12, mass21, mass22, Cab, 
        CtCab, θ, P, H, Ctdot_θdot1, Ctdot_θdot2, Ctdot_θdot3)

Calculate the derivatives of the element resultants with respect to the state rates.

# Arguments:
 - `ΔL`: beam element length
 - `mass11, mass12, mass21, mass22`: beam element mass matrix, divided into submatrices
 - `Cab`: transformation matrix from the undeformed beam element frame to the body frame
 - `CtCab`: transformation matrix from the deformed beam element frame to the body frame
 - `θ`: angular displacement
 - `P`: linear momentum
 - `H`: angular momentum
 - `Ctdot_θdot1`: Derivative of `Cdot'` w.r.t. `θdot[1]`
 - `Ctdot_θdot2`: Derivative of `Cdot'` w.r.t. `θdot[2]`
 - `Ctdot_θdot3`: Derivative of `Cdot'` w.r.t. `θdot[3]`
"""
@inline function element_mass_matrix_equations(ΔL, mass11, mass12, mass21, mass22, Cab, 
    CtCab, θ, P, H, Ctdot_θdot1, Ctdot_θdot2, Ctdot_θdot3)

    tmp = ΔL/2*mul3(Ctdot_θdot1, Ctdot_θdot2, Ctdot_θdot3, Cab*P)
    f_u1_θdot = tmp
    f_u2_θdot = tmp

    tmp = ΔL/2*CtCab*mass11
    f_u1_Vdot = tmp
    f_u2_Vdot = tmp

    tmp = ΔL/2*CtCab*mass12
    f_u1_Ωdot = tmp
    f_u2_Ωdot = tmp

    tmp = ΔL/2*mul3(Ctdot_θdot1, Ctdot_θdot2, Ctdot_θdot3, Cab*H)
    f_ψ1_θdot = tmp
    f_ψ2_θdot = tmp

    tmp = ΔL/2*CtCab*mass21
    f_ψ1_Vdot = tmp
    f_ψ2_Vdot = tmp

    tmp = ΔL/2*CtCab*mass22
    f_ψ1_Ωdot = tmp
    f_ψ2_Ωdot = tmp

    f_V_udot = -I3

    Q = get_Q(θ)
    f_Ω_θdot = -Cab'*Q

    return f_u1_θdot, f_u2_θdot, f_u1_Vdot, f_u2_Vdot, f_u1_Ωdot, f_u2_Ωdot,
        f_ψ1_θdot, f_ψ2_θdot, f_ψ1_Vdot, f_ψ2_Vdot, f_ψ1_Ωdot, f_ψ2_Ωdot,
        f_V_udot, f_Ω_θdot
end

"""
    insert_element_mass_matrix!(jacob, force_scaling, irow_e, irow_p1, 
        irow_p2, icol, 
        f_u1_θdot, f_u2_θdot, f_u1_Vdot, f_u2_Vdot, f_u1_Ωdot, f_u2_Ωdot,
        f_ψ1_θdot, f_ψ2_θdot, f_ψ1_Vdot, f_ψ2_Vdot, f_ψ1_Ωdot, f_ψ2_Ωdot,
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
    f_u1_θdot, f_u2_θdot, f_u1_Vdot, f_u2_Vdot, f_u1_Ωdot, f_u2_Ωdot,
    f_ψ1_θdot, f_ψ2_θdot, f_ψ1_Vdot, f_ψ2_Vdot, f_ψ1_Ωdot, f_ψ2_Ωdot,
    f_V_udot, f_Ω_θdot)

    # create jacobian entries for the start of the beam
    jacob[irow_p1:irow_p1+2, icol+3:icol+5] .= f_u1_θdot ./ force_scaling
    jacob[irow_p1+3:irow_p1+5, icol+3:icol+5] .= f_ψ1_θdot ./ force_scaling
    jacob[irow_p1:irow_p1+2, icol+12:icol+14] .= f_u1_Vdot  ./ force_scaling
    jacob[irow_p1:irow_p1+2, icol+15:icol+17] .= f_u1_Ωdot  ./ force_scaling
    jacob[irow_p1+3:irow_p1+5, icol+12:icol+14] .= f_ψ1_Vdot ./ force_scaling
    jacob[irow_p1+3:irow_p1+5, icol+15:icol+17] .= f_ψ1_Ωdot ./ force_scaling

    # create jacobian entries for the end of the beam
    jacob[irow_p2:irow_p2+2, icol+3:icol+5] .= f_u2_θdot ./ force_scaling
    jacob[irow_p2+3:irow_p2+5, icol+3:icol+5] .= f_ψ2_θdot ./ force_scaling
    jacob[irow_p2:irow_p2+2, icol+12:icol+14] .= f_u2_Vdot  ./ force_scaling
    jacob[irow_p2:irow_p2+2, icol+15:icol+17] .= f_u2_Ωdot  ./ force_scaling
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
        f_u1_θdot, f_u2_θdot, f_u1_Vdot, f_u2_Vdot, f_u1_Ωdot, f_u2_Ωdot,
        f_ψ1_θdot, f_ψ2_θdot, f_ψ1_Vdot, f_ψ2_Vdot, f_ψ1_Ωdot, f_ψ2_Ωdot,
        f_V_udot, f_Ω_θdot)

Add the beam element's mass matrix to the system jacobian matrix `jacob`, scaled
by the scaling parameter `gamma`.
"""
@inline function insert_element_mass_matrix!(jacob, gamma, force_scaling, irow_e,
    irow_p1, irow_p2, icol, 
    f_u1_θdot, f_u2_θdot, f_u1_Vdot, f_u2_Vdot, f_u1_Ωdot, f_u2_Ωdot,
    f_ψ1_θdot, f_ψ2_θdot, f_ψ1_Vdot, f_ψ2_Vdot, f_ψ1_Ωdot, f_ψ2_Ωdot,
    f_V_udot, f_Ω_θdot)

    # create jacobian entries for the start of the beam
    jacob[irow_p1:irow_p1+2, icol+3:icol+5] .= f_u1_θdot .* (gamma/force_scaling)
    jacob[irow_p1+3:irow_p1+5, icol+3:icol+5] .= f_ψ1_θdot .* (gamma/force_scaling)
    jacob[irow_p1:irow_p1+2, icol+12:icol+14] .= f_u1_Vdot  .* (gamma/force_scaling)
    jacob[irow_p1:irow_p1+2, icol+15:icol+17] .= f_u1_Ωdot  .* (gamma/force_scaling)
    jacob[irow_p1+3:irow_p1+5, icol+12:icol+14] .= f_ψ1_Vdot .* (gamma/force_scaling)
    jacob[irow_p1+3:irow_p1+5, icol+15:icol+17] .= f_ψ1_Ωdot .* (gamma/force_scaling)

    # create jacobian entries for the end of the beam
    jacob[irow_p2:irow_p2+2, icol+3:icol+5] .= f_u2_θdot .* (gamma/force_scaling)
    jacob[irow_p2+3:irow_p2+5, icol+3:icol+5] .= f_ψ2_θdot .* (gamma/force_scaling)
    jacob[irow_p2:irow_p2+2, icol+12:icol+14] .= f_u2_Vdot  .* (gamma/force_scaling)
    jacob[irow_p2:irow_p2+2, icol+15:icol+17] .= f_u2_Ωdot  .* (gamma/force_scaling)
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

    # get beam element properties
    ΔL, mass11, mass12, mass21, mass22, C, Cab, CtCab, θ, V, Ω, P, H, Ctdot_θdot1, Ctdot_θdot2, 
        Ctdot_θdot3 = element_mass_matrix_properties(x, icol, elem)

    # additional precomputed quantities
    Ct = C'

    # get jacobians of beam element equations
    f_u1_θdot, f_u2_θdot, f_u1_Vdot, f_u2_Vdot, f_u1_Ωdot, f_u2_Ωdot,
        f_ψ1_θdot, f_ψ2_θdot, f_ψ1_Vdot, f_ψ2_Vdot, f_ψ1_Ωdot, f_ψ2_Ωdot,
        f_V_udot, f_Ω_θdot = element_mass_matrix_equations(ΔL, mass11, mass12, mass21, mass22, Cab, 
        CtCab, θ, P, H, Ctdot_θdot1, Ctdot_θdot2, Ctdot_θdot3)

    # add jacobians for point mass loads
    if haskey(point_masses, ielem)
        # get point mass properties
        massp11, massp12, massp21, massp22, Pp, Hp = point_mass_rate_jacobian_properties(
            point_masses[ielem], Cab, V, Ω)
        # get point mass load jacobians
        Fp_θdot, Fp_Vdot, Fp_Ωdot, Mp_θdot, Mp_Vdot, Mp_Ωdot = point_mass_rate_jacobian(
            massp11, massp12, massp21, massp22, Ct, Cab, Pp, Hp, Ctdot_θdot1, Ctdot_θdot2, 
            Ctdot_θdot3)
        # add to element resultant jacobians
        f_u1_θdot -= Fp_θdot/2
        f_u2_θdot -= Fp_θdot/2
        f_u1_Vdot -= Fp_Vdot/2
        f_u2_Vdot -= Fp_Vdot/2
        f_u1_Ωdot -= Fp_Ωdot/2
        f_u2_Ωdot -= Fp_Ωdot/2
        f_ψ1_θdot -= Mp_θdot/2
        f_ψ2_θdot -= Mp_θdot/2
        f_ψ1_Vdot -= Mp_Vdot/2
        f_ψ2_Vdot -= Mp_Vdot/2
        f_ψ1_Ωdot -= Mp_Ωdot/2
        f_ψ2_Ωdot -= Mp_Ωdot/2
    end

    # initialize/insert into jacobian matrix for the system
    insert_element_mass_matrix!(jacob, force_scaling, irow_e,
        irow_p1, irow_p2, icol, 
        f_u1_θdot, f_u2_θdot, f_u1_Vdot, f_u2_Vdot, f_u1_Ωdot, f_u2_Ωdot,
        f_ψ1_θdot, f_ψ2_θdot, f_ψ1_Vdot, f_ψ2_Vdot, f_ψ1_Ωdot, f_ψ2_Ωdot,
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

    # get beam element properties
    ΔL, mass11, mass12, mass21, mass22, C, Cab, CtCab, θ, V, Ω, P, H, Ctdot_θdot1, Ctdot_θdot2, 
        Ctdot_θdot3 = element_mass_matrix_properties(x, icol, elem)

    # additional precomputed quantities
    Ct = C'

    # get jacobians of beam element equations
    f_u1_θdot, f_u2_θdot, f_u1_Vdot, f_u2_Vdot, f_u1_Ωdot, f_u2_Ωdot,
        f_ψ1_θdot, f_ψ2_θdot, f_ψ1_Vdot, f_ψ2_Vdot, f_ψ1_Ωdot, f_ψ2_Ωdot,
        f_V_udot, f_Ω_θdot = element_mass_matrix_equations(ΔL, mass11, mass12, mass21, mass22, Cab, 
        CtCab, θ, P, H, Ctdot_θdot1, Ctdot_θdot2, Ctdot_θdot3)

    # add jacobians for point mass loads
    if haskey(point_masses, ielem)
        # get point mass properties
        massp11, massp12, massp21, massp22, Pp, Hp = point_mass_rate_jacobian_properties(
            point_masses[ielem], Cab, V, Ω)
        # get point mass load jacobians
        Fp_θdot, Fp_Vdot, Fp_Ωdot, Mp_θdot, Mp_Vdot, Mp_Ωdot = point_mass_rate_jacobian(
            massp11, massp12, massp21, massp22, Ct, Cab, Pp, Hp, Ctdot_θdot1, Ctdot_θdot2, 
            Ctdot_θdot3)
        # add to element resultant jacobians
        f_u1_θdot -= Fp_θdot/2
        f_u2_θdot -= Fp_θdot/2
        f_u1_Vdot -= Fp_Vdot/2
        f_u2_Vdot -= Fp_Vdot/2
        f_u1_Ωdot -= Fp_Ωdot/2
        f_u2_Ωdot -= Fp_Ωdot/2
        f_ψ1_θdot -= Mp_θdot/2
        f_ψ2_θdot -= Mp_θdot/2
        f_ψ1_Vdot -= Mp_Vdot/2
        f_ψ2_Vdot -= Mp_Vdot/2
        f_ψ1_Ωdot -= Mp_Ωdot/2
        f_ψ2_Ωdot -= Mp_Ωdot/2
    end

    # initialize/insert into jacobian matrix for the system
    insert_element_mass_matrix!(jacob, gamma, force_scaling, irow_e,
        irow_p1, irow_p2, icol, 
        f_u1_θdot, f_u2_θdot, f_u1_Vdot, f_u2_Vdot, f_u1_Ωdot, f_u2_Ωdot,
        f_ψ1_θdot, f_ψ2_θdot, f_ψ1_Vdot, f_ψ2_Vdot, f_ψ1_Ωdot, f_ψ2_Ωdot,
        f_V_udot, f_Ω_θdot)

    return jacob
end
