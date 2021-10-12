"""
    Element{TF}

Composite type that defines a beam element's properties

# Fields
 - `L`: Length of the beam element
 - `x`: Location of the beam element (the midpoint of the beam element)
 - `compliance`: Beam element compliance matrix
 - `mass`: Beam element mass matrix
 - `Cab`: Transformation matrix to the body frame from the beam element frame
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
Construct a beam element.

# Arguments
- `L`: Length of the beam element
- `x`: Location of the beam element (the midpoint of the beam element)
- `compliance`: Beam element compliance matrix
- `mass`: Beam element mass matrix
- `Cab`: Transformation matrix to the body frame from the beam element frame
"""
function Element(L, x, compliance, mass, Cab)
    TF = promote_type(typeof(L), eltype(x), eltype(compliance), eltype(mass), eltype(Cab))
    return Element{TF}(L, x, compliance, mass, Cab)
end

"""
    element_strain(element, F, M)

Calculate the strain of a beam element given the resultant force and moments
(in the deformed local frame)
"""
@inline function element_strain(element, F, M)
    C = element.compliance
    C11 = C[SVector{3}(1:3), SVector{3}(1:3)]
    C12 = C[SVector{3}(1:3), SVector{3}(4:6)]
    return C11*F + C12*M
end

@inline element_strain_F(element) = element.compliance[SVector{3}(1:3), SVector{3}(1:3)]

@inline element_strain_M(element) = element.compliance[SVector{3}(1:3), SVector{3}(4:6)]

"""
    element_curvature(element, F, M)

Calculate the curvature of a beam element given the resultant force and moments
(in the deformed local frame)
"""
@inline function element_curvature(element, F, M)
    C = element.compliance
    C21 = C[SVector{3}(4:6), SVector{3}(1:3)]
    C22 = C[SVector{3}(4:6), SVector{3}(4:6)]
    return C21*F + C22*M
end

@inline element_curvature_F(element) = element.compliance[SVector{3}(4:6), SVector{3}(1:3)]

@inline element_curvature_M(element) = element.compliance[SVector{3}(4:6), SVector{3}(4:6)]

"""
    element_linear_momentum(element, V, Ω)

Calculate the linear momentum (in the deformed local beam frame) of a beam
element given the element's linear and angular velocity
"""
@inline function element_linear_momentum(element, V, Ω)
    M = element.mass
    M11 = M[SVector{3}(1:3), SVector{3}(1:3)]
    M12 = M[SVector{3}(1:3), SVector{3}(4:6)]
    return M11*V + M12*Ω
end

@inline element_linear_momentum_V(element) = element.mass[SVector{3}(1:3), SVector{3}(1:3)]

@inline element_linear_momentum_Ω(element) = element.mass[SVector{3}(1:3), SVector{3}(4:6)]

"""
    element_angular_momentum(element, V, Ω)

Calculate the angular momentum (in the deformed local beam frame) of a beam
element given the element's linear and angular velocity
"""
@inline function element_angular_momentum(element, V, Ω)
    M = element.mass
    M21 = M[SVector{3}(4:6), SVector{3}(1:3)]
    M22 = M[SVector{3}(4:6), SVector{3}(4:6)]
    return M21*V + M22*Ω
end

@inline element_angular_momentum_V(element) = element.mass[SVector{3}(4:6), SVector{3}(1:3)]

@inline element_angular_momentum_Ω(element) = element.mass[SVector{3}(4:6), SVector{3}(4:6)]

"""
    element_properties(x, icol, elem, force_scaling)
    element_properties(x, icol, elem, force_scaling, x0, v0, ω0)
    element_properties(x, icol, elem, force_scaling, x0, v0, ω0, u0, θ0, udot0, θdot0)
    element_properties(x, icol, elem, force_scaling, x0, v0, ω0, udot_init, 
        θdot_init, Vdot_init, Ωdot_init, dt)

Extract/calculate the properties of a specific beam element.

There are four implementations corresponding to the following analysis types:
 - Static
 - Dynamic - Steady State
 - Dynamic - Initial Step (for initializing time domain simulations)
 - Dynamic - Time Marching

See "GEBT: A general-purpose nonlinear analysis tool for composite beams" by
Wenbin Yu and "Geometrically nonlinear analysis of composite beams using
Wiener-Milenković parameters" by Qi Wang and Wenbin Yu.

# Arguments
 - `x`: current state vector
 - `icol`: starting index for the beam's state variables
 - `ielem`: beam element index
 - `elem`: beam element
 - `force_scaling`: scaling parameter for forces/moments

# Additional Arguments for Dynamic Analyses
 - `x0`: Global frame origin (for the current time step)
 - `v0`: Global frame linear velocity (for the current time step)
 - `ω0`: Global frame angular velocity (for the current time step)

# Additional Arguments for Initial Step Analyses
 - `u0`: deflection variables for each beam element
 - `θ0`: rotation variables for each beam element
 - `udot0`: time derivative of u for each beam element
 - `θdot0`: time derivative of θ for each beam element

# Additional Arguments for Time Marching Analyses
 - `udot_init`: `2/dt*u + udot` for each beam element from the previous time step
 - `θdot_init`: `2/dt*θ + θdot` for each beam element from the previous time step
 - `Vdot_init`: `2/dt*C'*Cab*P + C'*Cab*Pdot` for each beam element from the 
    previous time step
 - `Ωdot_init`: `2/dt*C'*Cab*H + C'*Cab*Hdot` for each beam element from the 
    previous time step
 - `dt`: time step size
"""
element_properties

# static
@inline element_properties(x, icol, elem, force_scaling) =
    static_element_properties(x, icol, elem, force_scaling)

# dynamic - steady state
@inline function element_properties(x, icol, elem, force_scaling, x0, v0, ω0)

    return steady_state_element_properties(x, icol, elem, force_scaling, x0, v0, ω0)
end

# dynamic - initial step
@inline function element_properties(x, icol, elem, force_scaling, x0, v0, ω0,
    u0, θ0, udot0, θdot0)

    return initial_step_element_properties(x, icol, elem, force_scaling, x0, v0, ω0,
        u0, θ0, udot0, θdot0)
end

# dynamic - newmark scheme time-marching
@inline function element_properties(x, icol, elem, force_scaling, x0, v0, ω0, udot_init,
    θdot_init, Vdot_init, Ωdot_init, dt)

    return newmark_element_properties(x, icol, elem, force_scaling, x0, v0, ω0, udot_init,
        θdot_init, Vdot_init, Ωdot_init, dt)
end

# static
@inline function static_element_properties(x, icol, elem, force_scaling)

    u = SVector(x[icol   ], x[icol+1 ], x[icol+2 ])
    θ = SVector(x[icol+3 ], x[icol+4 ], x[icol+5 ])
    F = SVector(x[icol+6 ], x[icol+7 ], x[icol+8 ]) .* force_scaling
    M = SVector(x[icol+9 ], x[icol+10], x[icol+11]) .* force_scaling

    ΔL = elem.L
    Ct = get_C(θ)'
    Cab = elem.Cab
    CtCab = Ct*Cab
    γ = element_strain(elem, F, M)
    κ = element_curvature(elem, F, M)

    return ΔL, Ct, Cab, CtCab, u, θ, F, M, γ, κ
end

# dynamic - steady state
@inline function steady_state_element_properties(x, icol, elem, force_scaling, x0, v0, ω0)

    ΔL, Ct, Cab, CtCab, u, θ, F, M, γ, κ = static_element_properties(x, icol, elem, force_scaling)

    v, ω, P, H, V, Ω = element_dynamic_properties(x, icol, elem, x0, v0, ω0)

    return ΔL, Ct, Cab, CtCab, u, θ, F, M, γ, κ, v, ω, P, H, V, Ω
end

# dynamic - initial step
@inline function initial_step_element_properties(x, icol, elem, force_scaling, x0, v0, ω0,
    u0, θ0, udot0, θdot0)

    # note that CtCabPdot and CtCabHdot are state variables instead of u and θ
    u = u0
    θ = θ0
    F = SVector(x[icol+6 ], x[icol+7 ], x[icol+8 ]) .* force_scaling
    M = SVector(x[icol+9 ], x[icol+10], x[icol+11]) .* force_scaling

    ΔL = elem.L
    Ct = get_C(θ)'
    Cab = elem.Cab
    CtCab = Ct*Cab
    γ = element_strain(elem, F, M)
    κ = element_curvature(elem, F, M)

    v, ω, P, H, V, Ω = element_dynamic_properties(x, icol, elem, x0, v0, ω0)

    udot = udot0
    θdot = θdot0

    Vdot = SVector(x[icol   ], x[icol+1 ], x[icol+2 ])
    Ωdot = SVector(x[icol+3 ], x[icol+4 ], x[icol+5 ])

    Pdot = element_linear_momentum(elem, Vdot, Ωdot)
    Hdot = element_angular_momentum(elem, Vdot, Ωdot)

    CtCabdot = get_C_t(Ct', θ, θdot)'*Cab

    CtCabPdot = CtCabdot*P + CtCab*Pdot
    CtCabHdot = CtCabdot*H + CtCab*Hdot

    return ΔL, Ct, Cab, CtCab, u, θ, F, M, γ, κ, v, ω, P, H, V, Ω, udot, θdot, 
        CtCabPdot, CtCabHdot, CtCabdot
end

# dynamic - newmark scheme time-marching
@inline function newmark_element_properties(x, icol, elem, force_scaling, x0, v0, ω0, 
    udot_init, θdot_init, Vdot_init, Ωdot_init, dt)

    ΔL, Ct, Cab, CtCab, u, θ, F, M, γ, κ, v, ω, P, H, V, Ω = steady_state_element_properties(
        x, icol, elem, force_scaling, x0, v0, ω0)

    udot = 2/dt*u - udot_init
    θdot = 2/dt*θ - θdot_init

    Vdot = 2/dt*V - Vdot_init
    Ωdot = 2/dt*Ω - Ωdot_init

    Pdot = element_linear_momentum(elem, Vdot, Ωdot)
    Hdot = element_angular_momentum(elem, Vdot, Ωdot)

    CtCabdot = get_C_t(Ct', θ, θdot)'*Cab

    CtCabPdot = CtCabdot*P + CtCab*Pdot
    CtCabHdot = CtCabdot*H + CtCab*Hdot

    return ΔL, Ct, Cab, CtCab, u, θ, F, M, γ, κ, v, ω, P, H, V, Ω, udot, θdot, Vdot, Ωdot, 
        Pdot, Hdot, CtCabPdot, CtCabHdot, CtCabdot
end

# dynamic - general
@inline function dynamic_element_properties(x, icol, elem, force_scaling, x0, v0, ω0, 
    udot, θdot, Vdot, Ωdot)

    ΔL, Ct, Cab, CtCab, u, θ, F, M, γ, κ, v, ω, P, H, V, Ω = steady_state_element_properties(
        x, icol, elem, force_scaling, x0, v0, ω0)

    CtCabdot = get_C_t(Ct', θ, θdot)'*Cab
    
    Pdot = element_linear_momentum(elem, Vdot, Ωdot)
    Hdot = element_angular_momentum(elem, Vdot, Ωdot)

    CtCabPdot = CtCabdot*P + CtCab*Pdot
    CtCabHdot = CtCabdot*H + CtCab*Hdot

    return ΔL, Ct, Cab, CtCab, u, θ, F, M, γ, κ, v, ω, P, H, V, Ω, udot, θdot,
        CtCabPdot, CtCabHdot, CtCabdot
end

"""
    element_dynamic_properties(x, icol, elem, x0, v0, ω0)

Extract/Compute `v`, `ω`, `P`, `H`, `V`, and `Ω`.
"""
@inline function element_dynamic_properties(x, icol, elem, x0, v0, ω0)

    v = v0 + cross(ω0, elem.x - x0)
    ω = ω0

    V = SVector(x[icol+12], x[icol+13], x[icol+14])
    Ω = SVector(x[icol+15], x[icol+16], x[icol+17])

    P = element_linear_momentum(elem, V, Ω)
    H = element_angular_momentum(elem, V, Ω)

    return v, ω, P, H, V, Ω
end

"""
    element_equations(ΔL, Cab, CtCab, u, θ, F, M, γ, κ)
    element_equations(ΔL, Cab, CtCab, u, θ, F, M, γ, κ, v, ω, P, H, V, Ω)
    element_equations(ΔL, Cab, CtCab, u, θ, F, M, γ, κ, v, ω, P, H, V, Ω,
        udot, θdot, CtCabPdot, CtCabHdot)

Evaluate the nonlinear equations for a beam element.

There are three implementations corresponding to the following analysis types:
 - Static
 - Dynamic - Steady State
 - Dynamic - Initial Step or Time Marching

See "GEBT: A general-purpose nonlinear analysis tool for composite beams" by
Wenbin Yu and "Geometrically nonlinear analysis of composite beams using
Wiener-Milenković parameters" by Qi Wang and Wenbin Yu.

# Arguments:
 - `ΔL`: Length of the beam element
 - `Ct`: Rotation tensor of the beam deformation in the "a" frame, transposed
 - `Cab`: Direction cosine matrix from "a" to "b" frame for the element
 - `CtCab`: `C'*Cab`, precomputed for efficiency
 - `u`: Displacement variables for the element [u1, u2, u3]
 - `θ`: Rotation variables for the element [θ1, θ2, θ3]
 - `F`: Force variables for the element [F1, F2, F3]
 - `M`: Moment variables for the element [M1, M2, M3]
 - `γ`: Engineering strains in the element [γ11, 2γ12, 2γ13]
 - `κ`: Curvatures in the element [κ1, κ2, κ3]

# Additional Arguments for Dynamic Analyses
 - `v`: Linear velocity of element in global frame "a" [v1, v2, v3]
 - `ω`: Angular velocity of element in global frame "a" [ω1, ω2, ω3]
 - `P`: Linear momenta for the element [P1, P2, P3]
 - `H`: Angular momenta for the element [H1, H2, H3]
 - `V`: Velocity of the element
 - `Ω`: Rotational velocity of the element

# Additional Arguments for Initial Step Analysis
 - `udot`: user-specified time derivative of u
 - `θdot`: user-specified time derivative of θ
 - `CtCabPdot`: C'*Cab*Pdot state variable
 - `CtCabHdot`: C'*Cab*Hdot state variable

# Additional Arguments for Time Marching Analysis
 - `udot`: Time derivative of `u`
 - `θdot`: Time derivative of `θ`
 - `CtCabPdot`: Time derivative of `CtCabPdot`
 - `CtCabHdot`: Time derivative of `CtCabHdot`
"""
element_equations

# static
@inline function element_equations(ΔL, Cab, CtCab, u, θ, F, M, γ, κ)
    return static_element_equations(ΔL, Cab, CtCab, u, θ, F, M, γ, κ)
end

# dynamic - steady state
@inline function element_equations(ΔL, Cab, CtCab, u, θ, F, M, γ, κ, v, ω,
    P, H, V, Ω)
    return steady_state_element_equations(ΔL, Cab, CtCab, u, θ, F, M, γ, κ, v, ω,
        P, H, V, Ω)
end

# dynamic - general
@inline function element_equations(ΔL, Cab, CtCab, u, θ, F, M, γ, κ, v, ω, P, H,
    V, Ω, udot, θdot, CtCabPdot, CtCabHdot)
    return dynamic_element_equations(ΔL, Cab, CtCab, u, θ, F, M, γ, κ, v, ω, P, H,
        V, Ω, udot, θdot, CtCabPdot, CtCabHdot)
end

# static
@inline function static_element_equations(ΔL, Cab, CtCab, u, θ, F, M, γ, κ)

    tmp = CtCab*F
    f_u1 = -tmp
    f_u2 =  tmp

    tmp1 = CtCab*M
    tmp2 = ΔL/2*CtCab*cross(e1 + γ, F)
    f_ψ1 = -tmp1 - tmp2
    f_ψ2 =  tmp1 - tmp2

    tmp = ΔL/2*(CtCab*(e1 + γ) - Cab*e1)
    f_F1 =  u - tmp
    f_F2 = -u - tmp

    Qinv = get_Qinv(θ)
    tmp = ΔL/2*Qinv*Cab*κ
    f_M1 =  θ - tmp
    f_M2 = -θ - tmp

    return f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2
end

# dynamic - steady state
@inline function steady_state_element_equations(ΔL, Cab, CtCab, u, θ, F, M, γ, κ, v, ω,
    P, H, V, Ω)

    f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2 = static_element_equations(ΔL, Cab,
        CtCab, u, θ, F, M, γ, κ)

    tmp = cross(ω, ΔL/2*CtCab*P)
    f_u1 += tmp
    f_u2 += tmp

    tmp = cross(ω, ΔL/2*CtCab*H) + ΔL/2*CtCab*cross(V, P)
    f_ψ1 += tmp
    f_ψ2 += tmp

    f_V = CtCab*V - v - cross(ω, u)
    f_Ω = Ω - CtCab'*ω

    return f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2, f_V, f_Ω
end

# dynamic - general
@inline function dynamic_element_equations(ΔL, Cab, CtCab, u, θ, F, M, γ, κ, v, ω, P, H,
    V, Ω, udot, θdot, CtCabPdot, CtCabHdot)

    f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2, f_V, f_Ω = steady_state_element_equations(
        ΔL, Cab, CtCab, u, θ, F, M, γ, κ, v, ω, P, H, V, Ω)

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
    insert_element_residual!(resid, force_scaling, irow_e, irow_e1,
        irow_p1, irow_e2, irow_p2, f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2)
    insert_element_residual!(resid, force_scaling, irow_e, irow_e1,
        irow_p1, irow_e2, irow_p2, f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2,
        f_V, f_Ω)

Insert beam element resultants into the residual equation.  Initialize equilibrium
and constitutive equations if they are not yet initialized.

If `irow_e1 != irow_p1` and/or `irow_e2 != irow_p2`, assume the equilibrium equations
for the left and/or right side are already initialized

There are two implementations corresponding to the following analysis types:
 - Static
 - Dynamic - Steady State, Initial Step, or Time Marching

# Arguments
 - `resid`: System residual vector
 - `force_scaling`: Scaling parameter for forces
 - `irow_e1`: Row index of the first equation for the left side of the beam element
    (a value <= 0 indicates the equations have been eliminated from the system of equations)
 - `irow_p1`: Row index of the first equation for the point on the left side of the beam element
 - `irow_e1`: Row index of the first equation for the right side of the beam element
     (a value <= 0 indicates the equations have been eliminated from the system of equations)
 - `irow_p2`: Row index of the first equation for the point on the right side of the beam element
 - `f_u1`, `f_u2`: Resultant displacements for the left and right side of the beam element, respectively
 - `f_ψ1`, `f_ψ2`: Resultant rotations for the left and right side of the beam element, respectively
 - `f_F1`, `f_F2`: Resultant forces for the left and right side of the beam element, respectively
 - `f_M1`, `f_M2`: Resultant moments for the left and right side of the beam element, respectively

# Additional Arguments for Dynamic Analyses
 - `f_V`: Resultant linear velocity of the beam element
 - `f_Ω`: Resultant angular velocity of the beam element
"""
insert_element_residual!

# static
@inline function insert_element_residual!(resid, force_scaling, irow_e, irow_e1, irow_p1,
    irow_e2, irow_p2, f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2)

    return insert_static_element_residual!(resid, force_scaling, irow_e, irow_e1, irow_p1,
        irow_e2, irow_p2, f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2)
end

# dynamic
@inline function insert_element_residual!(resid, force_scaling, irow_e, irow_e1, irow_p1,
    irow_e2, irow_p2, f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2, f_V, f_Ω)

    return insert_dynamic_element_residual!(resid, force_scaling, irow_e, irow_e1, irow_p1,
        irow_e2, irow_p2, f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2, f_V, f_Ω)
end

# static
@inline function insert_static_element_residual!(resid, force_scaling, irow_e, irow_e1, irow_p1,
    irow_e2, irow_p2, f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2)

    # create/add to residual equations for left endpoint
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

    # create/add to residual equations for right endpoint
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

# dynamic
@inline function insert_dynamic_element_residual!(resid, force_scaling, irow_e, irow_e1, irow_p1,
    irow_e2, irow_p2, f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2, f_V, f_Ω)

    resid = insert_static_element_residual!(resid, force_scaling, irow_e, irow_e1, irow_p1,
        irow_e2, irow_p2, f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2)

    # residual equations for element
    resid[irow_e:irow_e+2] .= f_V
    resid[irow_e+3:irow_e+5] .= f_Ω

    return resid
end

"""
    element_residual!(resid, x, ielem, elem, distributed_loads, point_masses, gvec, force_scaling,
        icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2)
    element_residual!(resid, x, ielem, elem, distributed_loads, point_masses, gvec, force_scaling,
        icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2, x0, v0,
        ω0)
    element_residual!(resid, x, ielem, elem, distributed_loads, point_masses, gvec, force_scaling,
        icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2, x0, v0,
        ω0, u0, θ0, udot0, θdot0)
    element_residual!(resid, x, ielem, elem, distributed_loads, point_masses, gvec, force_scaling,
        icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2, x0, v0,
        ω0, udot_init, θdot_init, Vdot_init, Ωdot_init, dt)

Compute and add a beam element's contributions to the residual vector

There are four implementations corresponding to the following analysis types:
 - Static
 - Dynamic - Steady State
 - Dynamic - Initial Step (for initializing time domain simulations)
 - Dynamic - Time Marching

See "GEBT: A general-purpose nonlinear analysis tool for composite beams" by
Wenbin Yu and "Geometrically nonlinear analysis of composite beams using
Wiener-Milenković parameters" by Qi Wang and Wenbin Yu.

# Arguments
 - `resid`: System residual vector
 - `x`: current state vector
 - `ielem`: beam element index
 - `elem`: beam element
 - `distributed_loads`: dictionary with all distributed loads
 - `point_masses`: dictionary with all point masses
 - `gvec`: gravity vector
 - `force_scaling`: scaling parameter for forces
 - `icol`: starting index for the beam's state variables
 - `irow_e1`: Row index of the first equation for the left side of the beam element
    (a value <= 0 indicates the equations have been eliminated from the system of equations)
 - `irow_p1`: Row index of the first equation for the point on the left side of the beam element
 - `irow_e2`: Row index of the first equation for the right side of the beam element
    (a value <= 0 indicates the equations have been eliminated from the system of equations)
 - `irow_p2`: Row index of the first equation for the point on the right side of the beam element

# Additional Arguments for Dynamic Analyses
 - `x0`: Global frame origin (for the current time step)
 - `v0`: Global frame linear velocity (for the current time step)
 - `ω0`: Global frame angular velocity (for the current time step)

# Additional Arguments for Initial Step Analyses
 - `u0`: initial deflection variables for each beam element
 - `θ0`: initial rotation variables for each beam element
 - `udot0`: initial time derivative of u for each beam element
 - `θdot0`: initial time derivative of θ for each beam element

# Additional Arguments for Time Marching Analyses
 - `udot_init`: `2/dt*u + udot` for each beam element from the previous time step
 - `θdot_init`: `2/dt*θ + θdot` for each beam element from the previous time step
 - `Vdot_init`: `2/dt*V + Vdot` for each beam element from the previous time step
 - `Ωdot_init`: `2/dt*Ω + Ωdot` for each beam element from the previous time step
 - `dt`: time step size
"""
element_residual!

# static
@inline function element_residual!(resid, x, ielem, elem, distributed_loads, point_masses, gvec,
    force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2)

    return static_element_residual!(resid, x, ielem, elem, distributed_loads, point_masses, gvec,
        force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2)
end

# dynamic - steady state
@inline function element_residual!(resid, x, ielem, elem, distributed_loads, point_masses, gvec,
    force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2, x0, v0, ω0)

    return steady_state_element_residual!(resid, x, ielem, elem, distributed_loads, point_masses, gvec,
        force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2, x0, v0, ω0)
end

# dynamic - initial step
@inline function element_residual!(resid, x, ielem, elem, distributed_loads, point_masses, gvec,
    force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2,
    x0, v0, ω0, u0, θ0, udot0, θdot0)

    return initial_step_element_residual!(resid, x, ielem, elem, distributed_loads, point_masses, gvec,
        force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2,
        x0, v0, ω0, u0, θ0, udot0, θdot0)
end

# time marching - Newmark scheme
@inline function element_residual!(resid, x, ielem, elem, distributed_loads, point_masses, gvec,
    force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2, x0, v0, ω0, udot_init, θdot_init,
    Vdot_init, Ωdot_init, dt)

    # time marching - Newmark scheme
    return newmark_element_residual!(resid, x, ielem, elem, distributed_loads, point_masses, gvec,
        force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2, x0, v0, ω0, udot_init, θdot_init,
        Vdot_init, Ωdot_init, dt)
end

# static
@inline function static_element_residual!(resid, x, ielem, elem, distributed_loads, point_masses, gvec,
    force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2,)

    # compute element properties
    ΔL, Ct, Cab, CtCab, u, θ, F, M, γ, κ = static_element_properties(x, icol, elem, force_scaling)

    # solve for the element resultants
    f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2 = static_element_equations(ΔL, Cab,
        CtCab, u, θ, F, M, γ, κ)

    # add distributed loads to the element equations (if applicable)
    if haskey(distributed_loads, ielem)
        f_u1 -= distributed_loads[ielem].f1 + Ct*distributed_loads[ielem].f1_follower
        f_u2 -= distributed_loads[ielem].f2 + Ct*distributed_loads[ielem].f2_follower
        f_ψ1 -= distributed_loads[ielem].m1 + Ct*distributed_loads[ielem].m1_follower
        f_ψ2 -= distributed_loads[ielem].m2 + Ct*distributed_loads[ielem].m2_follower
    end

    # add element gravitational loads to the element equations
    fg1, fg2, mg1, mg2 = element_gravitational_loads(ΔL, CtCab, elem.mass, gvec)

    f_u1 -= fg1
    f_u2 -= fg2
    f_ψ1 -= mg1
    f_ψ2 -= mg2

    if haskey(point_masses, ielem)
        # add point mass gravitational loads to the element equations
        fp1, fp2, mp1, mp2 = point_mass_gravitational_loads(Ct, point_masses[ielem].mass, gvec)

        f_u1 -= fp1
        f_u2 -= fp2
        f_ψ1 -= mp1
        f_ψ2 -= mp2
    end

    # insert element resultants into residual vector
    resid = insert_static_element_residual!(resid, force_scaling, irow_e, irow_e1, irow_p1,
        irow_e2, irow_p2, f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2)

    return resid
end

# dynamic - steady state
@inline function steady_state_element_residual!(resid, x, ielem, elem, distributed_loads, point_masses, gvec,
    force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2, x0, v0, ω0)

    # compute element properties
    ΔL, Ct, Cab, CtCab, u, θ, F, M, γ, κ, v, ω, P, H, V, Ω = steady_state_element_properties(
        x, icol, elem, force_scaling, x0, v0, ω0)

    # solve for element resultants
    f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2, f_V, f_Ω = steady_state_element_equations(
        ΔL, Cab, CtCab, u, θ, F, M, γ, κ, v, ω, P, H, V, Ω)

    # add distributed loads to the element equations (if applicable)
    if haskey(distributed_loads, ielem)
        f_u1 -= distributed_loads[ielem].f1 + Ct*distributed_loads[ielem].f1_follower
        f_u2 -= distributed_loads[ielem].f2 + Ct*distributed_loads[ielem].f2_follower
        f_ψ1 -= distributed_loads[ielem].m1 + Ct*distributed_loads[ielem].m1_follower
        f_ψ2 -= distributed_loads[ielem].m2 + Ct*distributed_loads[ielem].m2_follower
    end

    # add gravitational loads to the element equations
    fg1, fg2, mg1, mg2 = element_gravitational_loads(ΔL, CtCab, elem.mass, gvec)

    f_u1 -= fg1
    f_u2 -= fg2
    f_ψ1 -= mg1
    f_ψ2 -= mg2

    if haskey(point_masses, ielem)
        # add point mass gravitational loads to the element equations
        fp1, fp2, mp1, mp2 = point_mass_gravitational_loads(Ct, point_masses[ielem].mass, gvec)

        f_u1 -= fp1
        f_u2 -= fp2
        f_ψ1 -= mp1
        f_ψ2 -= mp2

        # add point mass acceleration loads to the element equations 
        fp1, fp2, mp1, mp2 = steady_state_point_loads(Ct, point_masses[ielem].mass, V, Ω)

        f_u1 -= fp1
        f_u2 -= fp2
        f_ψ1 -= mp1
        f_ψ2 -= mp2
    end

    # insert element resultants into the residual vector
    resid = insert_dynamic_element_residual!(resid, force_scaling, irow_e, irow_e1, irow_p1,
        irow_e2, irow_p2, f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2, f_V, f_Ω)

    return resid
end

# dynamic - initial step
@inline function initial_step_element_residual!(resid, x, ielem, elem, distributed_loads, point_masses, gvec,
    force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2,
    x0, v0, ω0, u0, θ0, udot0, θdot0)

    # compute element properties
    ΔL, Ct, Cab, CtCab, u, θ, F, M, γ, κ, v, ω, P, H, V, Ω, udot, θdot,
        CtCabPdot, CtCabHdot, CtCabdot  = initial_step_element_properties(x, icol, 
        elem, force_scaling, x0, v0, ω0, u0, θ0, udot0, θdot0)

    # solve for the element resultants
    f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2, f_V, f_Ω =
        dynamic_element_equations(ΔL, Cab, CtCab, u, θ, F, M, γ, κ, v, ω, P, H, V, Ω,
        udot, θdot, CtCabPdot, CtCabHdot)

    # add distributed loads to the element equations (if applicable)
    if haskey(distributed_loads, ielem)
        f_u1 -= distributed_loads[ielem].f1 + Ct*distributed_loads[ielem].f1_follower
        f_u2 -= distributed_loads[ielem].f2 + Ct*distributed_loads[ielem].f2_follower
        f_ψ1 -= distributed_loads[ielem].m1 + Ct*distributed_loads[ielem].m1_follower
        f_ψ2 -= distributed_loads[ielem].m2 + Ct*distributed_loads[ielem].m2_follower
    end

    # add gravitational loads to the element equations
    fg1, fg2, mg1, mg2 = element_gravitational_loads(ΔL, CtCab, elem.mass, gvec)

    f_u1 -= fg1
    f_u2 -= fg2
    f_ψ1 -= mg1
    f_ψ2 -= mg2

    # add point mass gravitational loads to the element equations
    if haskey(point_masses, ielem)
        fp1, fp2, mp1, mp2 = point_mass_gravitational_loads(Ct, point_masses[ielem].mass, gvec)

        f_u1 -= fp1
        f_u2 -= fp2
        f_ψ1 -= mp1
        f_ψ2 -= mp2

        # add point mass acceleration loads to the element equations 
        fp1, fp2, mp1, mp2 = dynamic_point_loads(Ct, point_masses[ielem].mass, V, Ω)

        f_u1 -= fp1
        f_u2 -= fp2
        f_ψ1 -= mp1
        f_ψ2 -= mp2
    end

    # insert element resultants into the residual vector
    resid = insert_dynamic_element_residual!(resid, force_scaling, irow_e, irow_e1, irow_p1,
        irow_e2, irow_p2, f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2, f_V, f_Ω)

    return resid
end

# time marching - Newmark scheme
@inline function newmark_element_residual!(resid, x, ielem, elem, distributed_loads, point_masses, gvec,
    force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2, x0, v0, ω0, udot_init, θdot_init,
    Vdot_init, Ωdot_init, dt)

    # compute element properties
    ΔL, Ct, Cab, CtCab, u, θ, F, M, γ, κ, v, ω, P, H, V, Ω, udot, θdot, Vdot, Ωdot, Pdot, Hdot, 
        CtCabPdot, CtCabHdot, CtCabdot = newmark_element_properties(x, icol,
        elem, force_scaling, x0, v0, ω0, udot_init, θdot_init, Vdot_init,
        Ωdot_init, dt)

    # solve for element resultants
    f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2, f_V, f_Ω =
        dynamic_element_equations(ΔL, Cab, CtCab, u, θ, F, M, γ, κ, v, ω, P, H, V, Ω,
        udot, θdot, CtCabPdot, CtCabHdot)

    # add distributed loads to the element equations (if applicable)
    if haskey(distributed_loads, ielem)
        f_u1 -= distributed_loads[ielem].f1 + Ct*distributed_loads[ielem].f1_follower
        f_u2 -= distributed_loads[ielem].f2 + Ct*distributed_loads[ielem].f2_follower
        f_ψ1 -= distributed_loads[ielem].m1 + Ct*distributed_loads[ielem].m1_follower
        f_ψ2 -= distributed_loads[ielem].m2 + Ct*distributed_loads[ielem].m2_follower
    end

    # add gravitational loads to the element equations
    fg1, fg2, mg1, mg2 = element_gravitational_loads(ΔL, CtCab, elem.mass, gvec)

    f_u1 -= fg1
    f_u2 -= fg2
    f_ψ1 -= mg1
    f_ψ2 -= mg2

    if haskey(point_masses, ielem)
        # add point mass gravitational loads to the element equations
        fp1, fp2, mp1, mp2 = point_mass_gravitational_loads(Ct, point_masses[ielem].mass, gvec)

        f_u1 -= fp1
        f_u2 -= fp2
        f_ψ1 -= mp1
        f_ψ2 -= mp2

        # add point mass acceleration loads to the element equations 
        fp1, fp2, mp1, mp2 = dynamic_point_loads(Ct, point_masses[ielem].mass, V, Ω)

        f_u1 -= fp1
        f_u2 -= fp2
        f_ψ1 -= mp1
        f_ψ2 -= mp2
    end

    # insert element resultants into the residual vector
    resid = insert_dynamic_element_residual!(resid, force_scaling, irow_e,
        irow_e1, irow_p1, irow_e2, irow_p2, f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2,
        f_M1, f_M2, f_V, f_Ω)

    return resid
end

# dynamic - general
@inline function dynamic_element_residual!(resid, x, ielem, elem, distributed_loads, point_masses, gvec,
    force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2, x0, v0, ω0,
    udot, θdot, Vdot, Ωdot)

    # compute element properties
    ΔL, Ct, Cab, CtCab, u, θ, F, M, γ, κ, v, ω, P, H, V, Ω, udot, θdot,
        CtCabPdot, CtCabHdot, CtCabdot = dynamic_element_properties(x, icol,
        elem, force_scaling, x0, v0, ω0, udot, θdot, Vdot, Ωdot)

    # solve for element resultants
    f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2, f_V, f_Ω =
        dynamic_element_equations(ΔL, Cab, CtCab, u, θ, F, M, γ, κ, v, ω, P, H, V, Ω,
        udot, θdot, CtCabPdot, CtCabHdot)

    # add distributed loads to the element equations (if applicable)
    if haskey(distributed_loads, ielem)
        f_u1 -= distributed_loads[ielem].f1 + Ct*distributed_loads[ielem].f1_follower
        f_u2 -= distributed_loads[ielem].f2 + Ct*distributed_loads[ielem].f2_follower
        f_ψ1 -= distributed_loads[ielem].m1 + Ct*distributed_loads[ielem].m1_follower
        f_ψ2 -= distributed_loads[ielem].m2 + Ct*distributed_loads[ielem].m2_follower
    end

    # add gravitational loads to the element equations
    fg1, fg2, mg1, mg2 = element_gravitational_loads(ΔL, CtCab, elem.mass, gvec)

    f_u1 -= fg1
    f_u2 -= fg2
    f_ψ1 -= mg1
    f_ψ2 -= mg2

    if haskey(point_masses, ielem)
        # add point mass gravitational loads to the element equations
        fp1, fp2, mp1, mp2 = point_mass_gravitational_loads(Ct, point_masses[ielem].mass, gvec)

        f_u1 -= fp1
        f_u2 -= fp2
        f_ψ1 -= mp1
        f_ψ2 -= mp2

        # add point mass acceleration loads to the element equations 
        fp1, fp2, mp1, mp2 = dynamic_point_loads(Ct, point_masses[ielem].mass, V, Ω)

        f_u1 -= fp1
        f_u2 -= fp2
        f_ψ1 -= mp1
        f_ψ2 -= mp2
    end

    # insert element resultants into the residual vector
    resid = insert_dynamic_element_residual!(resid, force_scaling, irow_e,
        irow_e1, irow_p1, irow_e2, irow_p2, f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2,
        f_M1, f_M2, f_V, f_Ω)

    return resid
end

"""
    element_jacobian_equations(elem, ΔL, Ct, Cab, CtCab, θ, F, M, γ, κ, Ct_θ1,
        Ct_θ2, Ct_θ3)
    element_jacobian_equations(elem, ΔL, Cab, CtCab, θ, F, M, γ, κ, ω, P, H, V,
        Ct_θ1, Ct_θ2, Ct_θ3)
    element_jacobian_equations(elem, ΔL, Cab, CtCab, θ, F, γ, ω, P, V, CtCabdot)
    element_jacobian_equations(elem, ΔL, Cab, CtCab, θ, F, M, γ, κ, ω, P, H, V,
        θdot, dt, Ct_θ1, Ct_θ2, Ct_θ3)

Find the jacobians of the nonlinear equations for a beam element with
respect to the state variables given the distributed loads on the beam element
and the beam element's properties.

There are four implementations corresponding to the following analysis types:
 - Static
 - Dynamic - Steady State
 - Dynamic - Initial Step (for initializing time domain simulations)
 - Dynamic - Time Marching

See "GEBT: A general-purpose nonlinear analysis tool for composite beams" by
Wenbin Yu and "Geometrically nonlinear analysis of composite beams using
Wiener-Milenković parameters" by Qi Wang and Wenbin Yu.

# Arguments:
 - `elem`: Beam element
 - `ΔL`: Length of the beam element
 - `Ct`: Rotation tensor of the beam deformation in the "a" frame, transposed
 - `Cab`: Direction cosine matrix from "a" to "b" frame for the element
 - `CtCab`: `C'*Cab`, precomputed for efficiency
 - `u`: Displacement variables for the element [u1, u2, u3]
 - `θ`: Rotation variables for the element [θ1, θ2, θ3]
 - `F`: Force variables for the element [F1, F2, F3]
 - `M`: Moment variables for the element [M1, M2, M3]
 - `γ`: Engineering strains in the element [γ11, 2γ12, 2γ13]
 - `κ`: Curvatures in the element [κ1, κ2, κ3]
 - `Ct_θ1`: Gradient of `Ct` w.r.t. `θ[1]`
 - `Ct_θ2`: Gradient of `Ct` w.r.t. `θ[2]`
 - `Ct_θ3`: Gradient of `Ct` w.r.t. `θ[3]`

# Additional Arguments for Dynamic Analyses
 - `v`: Linear velocity of element in global frame "a" [v1, v2, v3]
 - `ω`: Angular velocity of element in global frame "a" [ω1, ω2, ω3]
 - `P`: Linear momenta for the element [P1, P2, P3]
 - `H`: Angular momenta for the element [H1, H2, H3]
 - `V`: Velocity of the element
 - `Ω`: Rotational velocity of the element

# Additional Arguments for Initial Step Analyses
- `udot`: user-specified time derivative of u for this beam element
- `θdot`: user-specified time derivative of θ for this beam element
 - `CtCabPdot`: `C'*Cab*Pdot` (which is a state variable for the initial step analysis)
 - `CtCabHdot`: `C'*Cab*Hdot` (which is a state variable for the initial step analysis)

# Additional Arguments for Time Marching Analyses
 - `udot_init`: `2/dt*u + udot` for each beam element from the previous time step
 - `θdot_init`: `2/dt*θ + θdot` for each beam element from the previous time step
 - `CtCabPdot`: `2/dt*C'*Cab*P + C'*Cab*Pdot` for each beam element from the previous time step
 - `CtCabHdot`: `2/dt*C'*Cab*H + C'*Cab*Hdot` for each beam element from the previous time step
 - `dt`: time step size
"""
element_jacobian_equations

# static
@inline function element_jacobian_equations(elem, ΔL, Cab, CtCab, θ, F, M, γ, κ,
    Ct_θ1, Ct_θ2, Ct_θ3)

    return static_element_jacobian_equations(elem, ΔL, Cab, CtCab, θ, F, M, γ, κ,
        Ct_θ1, Ct_θ2, Ct_θ3)
end

# dynamic - steady state
@inline function element_jacobian_equations(elem, ΔL, Cab, CtCab, θ, F, M,
    γ, κ, ω, P, H, V, Ct_θ1, Ct_θ2, Ct_θ3)

    return steady_state_element_jacobian_equations(elem, ΔL, Cab, CtCab, θ, F, M,
        γ, κ, ω, P, H, V, Ct_θ1, Ct_θ2, Ct_θ3)
end

# dynamic - initial step
@inline function element_jacobian_equations(elem, ΔL, Cab, CtCab, θ, F, γ, ω, P, V, CtCabdot)

    # dynamic - initial step
    return initial_step_element_jacobian_equations(elem, ΔL, Cab, CtCab, θ, F, γ, ω, P, V, CtCabdot)
end

# dynamic - newmark scheme time-marching
@inline function element_jacobian_equations(elem, ΔL, Cab, CtCab, θ, F, M, γ, κ,
    ω, P, H, V, θdot, Vdot, Ωdot, dt, Ct_θ1, Ct_θ2, Ct_θ3, CtCabdot)

    return newmark_element_jacobian_equations(elem, ΔL, Cab, CtCab, θ, F, M, γ, κ,
        ω, P, H, V, Ω, θdot, Vdot, Ωdot, Pdot, Hdot, dt, Ct_θ1, Ct_θ2, Ct_θ3, CtCabdot, Ctdot_θ1, 
        Ctdot_θ2, Ctdot_θ3, Ctdot_θdot1, Ctdot_θdot2, Ctdot_θdot3)
end

# static
@inline function static_element_jacobian_equations(elem, ΔL, Cab, CtCab, θ, F, M, γ, κ,
    Ct_θ1, Ct_θ2, Ct_θ3)

    C11 = element_strain_F(elem)
    C12 = element_strain_M(elem)
    C21 = element_curvature_F(elem)
    C22 = element_curvature_M(elem)

    # --- f_u1, f_u2 --- #

    # d_fu/d_θ
    tmp = mul3(Ct_θ1, Ct_θ2, Ct_θ3, Cab*F)
    f_u1_θ = -tmp
    f_u2_θ =  tmp

    # d_fu/d_F
    f_u1_F = -CtCab
    f_u2_F =  CtCab

    # --- f_ψ1, f_ψ2 --- #

    # d_fψ/d_θ
    tmp1 = mul3(Ct_θ1, Ct_θ2, Ct_θ3, Cab*M)
    tmp2 = ΔL/2*mul3(Ct_θ1, Ct_θ2, Ct_θ3, Cab*cross(e1 + γ, F))
    f_ψ1_θ = -tmp1 - tmp2
    f_ψ2_θ =  tmp1 - tmp2

    # d_fψ/d_F
    tmp = -ΔL/2*CtCab*(tilde(e1 + γ) - tilde(F)*C11)
    f_ψ1_F = tmp
    f_ψ2_F = tmp

    # d_fψ/d_M
    tmp = ΔL/2*CtCab*tilde(F)*C12
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
    tmp = ΔL/2*CtCab*C11
    f_F1_F = -tmp
    f_F2_F = -tmp

    # d_fF/d_M
    tmp = ΔL/2*CtCab*C12
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
    tmp2 = tmp1*C21
    f_M1_F = tmp2
    f_M2_F = tmp2

    # d_fM/d_M
    tmp2 = tmp1*C22
    f_M1_M = tmp2
    f_M2_M = tmp2

    return f_u1_θ, f_u2_θ, f_u1_F, f_u2_F,
        f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M,
        f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M
end

# dynamic - steady state
@inline function steady_state_element_jacobian_equations(elem, ΔL, Cab, CtCab, θ, F, M,
    γ, κ, ω, P, H, V, Ct_θ1, Ct_θ2, Ct_θ3)

    f_u1_θ, f_u2_θ, f_u1_F, f_u2_F,
        f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M,
        f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M = static_element_jacobian_equations(
        elem, ΔL, Cab, CtCab, θ, F, M, γ, κ, Ct_θ1, Ct_θ2, Ct_θ3)


    M11 = element_linear_momentum_V(elem)
    M12 = element_linear_momentum_Ω(elem)
    M21 = element_angular_momentum_V(elem)
    M22 = element_angular_momentum_Ω(elem)

    # --- f_u1, f_u2 --- #

    # d_fu_dθ
    tmp = tilde(ω)*mul3(Ct_θ1, Ct_θ2, Ct_θ3, ΔL/2*Cab*P)
    f_u1_θ += tmp
    f_u2_θ += tmp

    # d_fu_dV
    tmp = ΔL/2*tilde(ω)*CtCab*M11
    f_u1_V = tmp
    f_u2_V = tmp

    # d_fu_dΩ
    tmp = ΔL/2*tilde(ω)*CtCab*M12
    f_u1_Ω = tmp
    f_u2_Ω = tmp

    # --- f_ψ1, f_ψ2 --- #

    # d_fψ_dθ
    tmp1 = tilde(ω)*mul3(Ct_θ1, Ct_θ2, Ct_θ3, ΔL/2*Cab*H)
    tmp2 = mul3(Ct_θ1, Ct_θ2, Ct_θ3, ΔL/2*Cab*cross(V, P))
    tmp3 = tmp1 + tmp2
    f_ψ1_θ += tmp3
    f_ψ2_θ += tmp3

    # d_fψ_dV
    tmp = ΔL/2*tilde(ω)*CtCab*M21 + ΔL/2*CtCab*(tilde(V)*M11 - tilde(P))
    f_ψ1_V = tmp
    f_ψ2_V = tmp

    # d_fψ_dΩ
    tmp = ΔL/2*tilde(ω)*CtCab*M22 + ΔL/2*CtCab*(tilde(V)*M12)
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
    f_Ω_θ = -Cab'*mul3(Ct_θ1', Ct_θ2', Ct_θ3', ω)

    # d_fH_dΩ
    f_Ω_Ω = I3

    return f_u1_θ, f_u2_θ, f_u1_F, f_u2_F, f_u1_V, f_u2_V, f_u1_Ω, f_u2_Ω,
        f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_V, f_ψ2_V, f_ψ1_Ω, f_ψ2_Ω,
        f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M,
        f_V_u, f_V_θ, f_V_V,
        f_Ω_θ, f_Ω_Ω
end

# dynamic - initial step
@inline function initial_step_element_jacobian_equations(elem, ΔL, Cab, CtCab, θ, F, γ, ω, P, V, CtCabdot)

    C11 = element_strain_F(elem)
    C12 = element_strain_M(elem)
    C21 = element_curvature_F(elem)
    C22 = element_curvature_M(elem)

    M11 = element_linear_momentum_V(elem)
    M12 = element_linear_momentum_Ω(elem)
    M21 = element_angular_momentum_V(elem)
    M22 = element_angular_momentum_Ω(elem)

    # --- f_u1, f_u2 --- #

    # d_fu/d_Vdot
    tmp = ΔL/2*CtCab*M11
    f_u1_Vdot = tmp
    f_u2_Vdot = tmp

    # d_fu/d_Ωdot
    tmp = ΔL/2*CtCab*M12
    f_u1_Ωdot = tmp
    f_u2_Ωdot = tmp

    # d_fu/d_F
    tmp = CtCab
    f_u1_F = -tmp
    f_u2_F =  tmp

    # d_fu/dV
    tmp = ΔL/2*(tilde(ω)*CtCab*M11 + CtCabdot*M11)
    f_u1_V = tmp
    f_u2_V = tmp

    # d_fu/dΩ
    tmp = ΔL/2*(tilde(ω)*CtCab*M12 + CtCabdot*M12)
    f_u1_Ω = tmp
    f_u2_Ω = tmp

    # --- f_θ1, f_θ2 --- #

    # d_fψ/d_Ωdot
    tmp = ΔL/2*CtCab*M21
    f_ψ1_Vdot = tmp
    f_ψ2_Vdot = tmp

    # d_fψ/d_Ωdot
    tmp = ΔL/2*CtCab*M22
    f_ψ1_Ωdot = tmp
    f_ψ2_Ωdot = tmp

    # d_fψ/d_F
    tmp = -ΔL/2*CtCab*(tilde(e1 + γ) - tilde(F)*C11)
    f_ψ1_F = tmp
    f_ψ2_F = tmp

    # d_fψ/d_M
    tmp = ΔL/2*CtCab*tilde(F)*C12
    f_ψ1_M = tmp - CtCab
    f_ψ2_M = tmp + CtCab

    # d_fψ_dV
    tmp = ΔL/2*tilde(ω)*CtCab*M21 + ΔL/2*CtCabdot*M21 + ΔL/2*CtCab*(tilde(V)*M11 - tilde(P))
    f_ψ1_V = tmp
    f_ψ2_V = tmp

    # d_fψ_dΩ
    tmp = ΔL/2*tilde(ω)*CtCab*M22 + ΔL/2*CtCabdot*M22 + ΔL/2*CtCab*(tilde(V)*M12)
    f_ψ1_Ω = tmp
    f_ψ2_Ω = tmp

    # --- f_F1, f_F2 --- #

    # d_fF/d_F
    tmp = ΔL/2*CtCab*C11
    f_F1_F = -tmp
    f_F2_F = -tmp

    # d_fF/d_M
    tmp = ΔL/2*CtCab*C12
    f_F1_M = -tmp
    f_F2_M = -tmp

    # --- f_M1, f_M2 --- #

    # d_fM/d_F
    Qinv = get_Qinv(θ)
    tmp1 = -ΔL/2*Qinv*Cab
    tmp2 = tmp1*C21
    f_M1_F = tmp2
    f_M2_F = tmp2

    # d_fM/d_M
    tmp2 = tmp1*C22
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

# dynamic - newmark scheme time-marching
@inline function newmark_element_jacobian_equations(elem, ΔL, Cab, CtCab, θ, F, M, γ, κ,
    ω, P, H, V, Ω, θdot, Vdot, Ωdot, Pdot, Hdot, dt, Ct_θ1, Ct_θ2, Ct_θ3, CtCabdot, 
    Ctdot_θ1, Ctdot_θ2, Ctdot_θ3, Ctdot_θdot1, Ctdot_θdot2, Ctdot_θdot3)

    f_u1_θ, f_u2_θ, f_u1_F, f_u2_F, f_u1_V, f_u2_V, f_u1_Ω, f_u2_Ω,
        f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_V, f_ψ2_V, f_ψ1_Ω, f_ψ2_Ω,
        f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M,
        f_V_u, f_V_θ, f_V_V,
        f_Ω_θ, f_Ω_Ω = steady_state_element_jacobian_equations(elem, ΔL, Cab, CtCab,
        θ, F, M, γ, κ, ω, P, H, V, Ct_θ1, Ct_θ2, Ct_θ3)

    M11 = element_linear_momentum_V(elem)
    M12 = element_linear_momentum_Ω(elem)
    M21 = element_angular_momentum_V(elem)
    M22 = element_angular_momentum_Ω(elem)

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
    tmp = ΔL/2*(CtCabdot*M11 + 2/dt*CtCab*M11)
    f_u1_V += tmp
    f_u2_V += tmp

    # d_fu_dΩ
    tmp = ΔL/2*(CtCabdot*M12 + 2/dt*CtCab*M12)
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
    tmp = ΔL/2*(CtCabdot*M21 + 2/dt*CtCab*M21)
    f_ψ1_V += tmp
    f_ψ2_V += tmp

    # d_fψ_dΩ
    tmp = ΔL/2*(CtCabdot*M22 + 2/dt*CtCab*M22)
    f_ψ1_Ω += tmp
    f_ψ2_Ω += tmp

    # --- d_fP_du --- #
    f_V_u -= 2/dt*I

    # --- d_fH_dθ --- #
    Q = get_Q(θ)
    Q_θ1, Q_θ2, Q_θ3 = get_Q_θ(θ)
    f_Ω_θ -= Cab'*(mul3(Q_θ1, Q_θ2, Q_θ3, θdot) + Q*2/dt)

    return f_u1_θ, f_u2_θ, f_u1_F, f_u2_F, f_u1_V, f_u2_V, f_u1_Ω, f_u2_Ω,
        f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_V, f_ψ2_V, f_ψ1_Ω, f_ψ2_Ω,
        f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M,
        f_V_u, f_V_θ, f_V_V,
        f_Ω_θ, f_Ω_Ω
end

# dynamic - general
@inline function dynamic_element_jacobian_equations(elem, ΔL, Cab, CtCab, θ, F, M, γ, κ,
    ω, P, H, V, θdot, Pdot, Hdot, Ct_θ1, Ct_θ2, Ct_θ3, CtCabdot, Ctdot_θ1, Ctdot_θ2, Ctdot_θ3)

    f_u1_θ, f_u2_θ, f_u1_F, f_u2_F, f_u1_V, f_u2_V, f_u1_Ω, f_u2_Ω,
        f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_V, f_ψ2_V, f_ψ1_Ω, f_ψ2_Ω,
        f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M,
        f_V_u, f_V_θ, f_V_V,
        f_Ω_θ, f_Ω_Ω = steady_state_element_jacobian_equations(elem, ΔL, Cab, CtCab,
        θ, F, M, γ, κ, ω, P, H, V, Ct_θ1, Ct_θ2, Ct_θ3)

    M11 = element_linear_momentum_V(elem)
    M12 = element_linear_momentum_Ω(elem)
    M21 = element_angular_momentum_V(elem)
    M22 = element_angular_momentum_Ω(elem)

    # --- f_u1, f_u2 --- #

    # ΔL/2*Ctdot_θ*Cab*(M11*V+M21*Ω) + Ct_θ*Cab*(M11*Vdot+M21*Ωdot)

    # d_fu_dθ
    tmp = ΔL/2*(mul3(Ctdot_θ1, Ctdot_θ2, Ctdot_θ3, Cab*P) + mul3(Ct_θ1, Ct_θ2, Ct_θ3, Cab*Pdot))
    f_u1_θ += tmp
    f_u2_θ += tmp

    # d_fu_dV
    tmp = ΔL/2*CtCabdot*M11
    f_u1_V += tmp
    f_u2_V += tmp

    # d_fu_dΩ
    tmp = ΔL/2*CtCabdot*M12
    f_u1_Ω += tmp
    f_u2_Ω += tmp

    # --- f_ψ1, f_ψ2 --- #

    # d_fψ_dθ
    tmp = ΔL/2*(mul3(Ctdot_θ1, Ctdot_θ2, Ctdot_θ3, Cab*H) + mul3(Ct_θ1, Ct_θ2, Ct_θ3, Cab*Hdot))
    f_ψ1_θ += tmp
    f_ψ2_θ += tmp

    # d_fψ_dV
    tmp = ΔL/2*CtCabdot*M21
    f_ψ1_V += tmp
    f_ψ2_V += tmp

    # d_fψ_dΩ
    tmp = ΔL/2*CtCabdot*M22
    f_ψ1_Ω += tmp
    f_ψ2_Ω += tmp

    # --- d_fH_dθ --- #
    Q_θ1, Q_θ2, Q_θ3 = get_Q_θ(θ)
    f_Ω_θ -= Cab'*mul3(Q_θ1, Q_θ2, Q_θ3, θdot)

    return f_u1_θ, f_u2_θ, f_u1_F, f_u2_F, f_u1_V, f_u2_V, f_u1_Ω, f_u2_Ω,
        f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_V, f_ψ2_V, f_ψ1_Ω, f_ψ2_Ω,
        f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M,
        f_V_u, f_V_θ, f_V_V,
        f_Ω_θ, f_Ω_Ω
end

"""
    insert_element_jacobian!(jacob, force_scaling, irow_e1,
        irow_p1, irow_e2, irow_p2, icol,
        f_u1_θ, f_u2_θ, f_u1_F, f_u2_F,
        f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M,
        f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M)
    insert_element_jacobian!(jacob, force_scaling, irow_e1,
        irow_p1, irow_e2, irow_p2, icol,
        f_u1_θ, f_u2_θ, f_u1_F, f_u2_F, f_u1_V, f_u2_V, f_u1_Ω, f_u2_Ω,
        f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_V, f_ψ2_V, f_ψ1_Ω, f_ψ2_Ω,
        f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M,
        f_V_u, f_V_θ, f_V_V,
        f_Ω_θ, f_Ω_Ω)
    insert_element_jacobian!(jacob, force_scaling, irow_e1,
        irow_p1, irow_e2, irow_p2, icol,
        f_u1_CtCabPdot, f_u2_CtCabPdot, f_u1_F, f_u2_F, f_u1_V, f_u2_V, f_u1_Ω, f_u2_Ω,
        f_ψ1_CtCabHdot, f_ψ2_CtCabHdot, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_V, f_ψ2_V, f_ψ1_Ω, f_ψ2_Ω,
        f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_F, f_M2_F, f_M1_M, f_M2_M,
        f_V_V, f_Ω_Ω)

Insert the the beam element jacobian entries into the jacobian matrix

There are three implementations corresponding to the following analysis types:
 - Static
 - Dynamic - Steady State or Time Marching
 - Dynamic - Initial Step (for initializing time domain simulations)

# Arguments
 - `jacob`: System jacobian matrix
 - `force_scaling`: Scaling parameter for forces
 - `icol_p1`: Row/column index of the first unknown for the left endpoint
   (a value <= 0 indicates the unknowns have been eliminated from the system of equations)
 - `irow_e1`: Row index of the first equation for the left side of the beam
 - `irow_p1`: Row index of the first equation for the point on the left side of the beam
 - `icol_p2`: Row/column index of the first unknown for the right endpoint
     (a value <= 0 indicates the unknowns have been eliminated from the system of equations)
 - `irow_e2`: Row index of the first equation for the right side of the beam
 - `irow_p2`: Row index of the first equation for the point on the right side of the beam
 - `icol`: Row/Column index corresponding to the first beam state variable

All other arguments use the following naming convention:
 - `f_y_x`: Jacobian of element equation "y" with respect to state variable "x"
"""
insert_element_jacobian!

# static
@inline function insert_element_jacobian!(jacob, force_scaling, irow_e1, irow_p1,
    irow_e2, irow_p2, icol,
    f_u1_θ, f_u2_θ, f_u1_F, f_u2_F,
    f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M,
    f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
    f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M)

    return insert_static_element_jacobian!(jacob, force_scaling, irow_e1,
        irow_p1, irow_e2, irow_p2, icol,
        f_u1_θ, f_u2_θ, f_u1_F, f_u2_F,
        f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M,
        f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M)
end

# dynamic - general
@inline function insert_element_jacobian!(jacob, force_scaling, irow_e, irow_e1,
    irow_p1, irow_e2, irow_p2, icol,
    f_u1_θ, f_u2_θ, f_u1_F, f_u2_F, f_u1_V, f_u2_V, f_u1_Ω, f_u2_Ω,
    f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_V, f_ψ2_V, f_ψ1_Ω, f_ψ2_Ω,
    f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
    f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M,
    f_V_u, f_V_θ, f_V_V,
    f_Ω_θ, f_Ω_Ω)

    return insert_dynamic_element_jacobian!(jacob, force_scaling, irow_e, irow_e1,
        irow_p1, irow_e2, irow_p2, icol,
        f_u1_θ, f_u2_θ, f_u1_F, f_u2_F, f_u1_V, f_u2_V, f_u1_Ω, f_u2_Ω,
        f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_V, f_ψ2_V, f_ψ1_Ω, f_ψ2_Ω,
        f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M,
        f_V_u, f_V_θ, f_V_V,
        f_Ω_θ, f_Ω_Ω)
end

# dynamic - initial step
@inline function insert_element_jacobian!(jacob, force_scaling, irow_e, irow_e1,
    irow_p1, irow_e2, irow_p2, icol,
    f_u1_Vdot, f_u2_Vdot, f_u1_Ωdot, f_u2_Ωdot, f_u1_F, f_u2_F, f_u1_V, f_u2_V, f_u1_Ω, f_u2_Ω,
    f_ψ1_Vdot, f_ψ2_Vdot, f_ψ1_Ωdot, f_ψ2_Ωdot, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_V, f_ψ2_V, f_ψ1_Ω, f_ψ2_Ω,
    f_F1_F, f_F2_F, f_F1_M, f_F2_M,
    f_M1_F, f_M2_F, f_M1_M, f_M2_M,
    f_V_V, f_Ω_Ω)

    return insert_initial_step_element_jacobian!(jacob, force_scaling, irow_e,
        irow_e1, irow_p1, irow_e2, irow_p2, icol,
        f_u1_Vdot, f_u2_Vdot, f_u1_Ωdot, f_u2_Ωdot, f_u1_F, f_u2_F, f_u1_V, f_u2_V, f_u1_Ω, f_u2_Ω,
        f_ψ1_Vdot, f_ψ2_Vdot, f_ψ1_Ωdot, f_ψ2_Ωdot, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_V, f_ψ2_V, f_ψ1_Ω, f_ψ2_Ω,
        f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_F, f_M2_F, f_M1_M, f_M2_M,
        f_V_V, f_Ω_Ω)
end

# static
@inline function insert_static_element_jacobian!(jacob, force_scaling, irow_e1,
    irow_p1, irow_e2, irow_p2, icol,
    f_u1_θ, f_u2_θ, f_u1_F, f_u2_F,
    f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M,
    f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
    f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M)

    # equilibrium equation jacobian entries for left endpoint
    jacob[irow_p1:irow_p1+2, icol+3:icol+5] .= f_u1_θ ./ force_scaling
    jacob[irow_p1:irow_p1+2, icol+6:icol+8] .= f_u1_F

    jacob[irow_p1+3:irow_p1+5, icol+3:icol+5] .= f_ψ1_θ ./ force_scaling
    jacob[irow_p1+3:irow_p1+5, icol+6:icol+8] .= f_ψ1_F
    jacob[irow_p1+3:irow_p1+5, icol+9:icol+11] .= f_ψ1_M

    # compatability equation jacobian entries for left endpoint
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

    # equilibrium equation jacobian entries for right endpoint
    jacob[irow_p2:irow_p2+2, icol+3:icol+5] .= f_u2_θ ./ force_scaling
    jacob[irow_p2:irow_p2+2, icol+6:icol+8] .= f_u2_F

    jacob[irow_p2+3:irow_p2+5, icol+3:icol+5] .= f_ψ2_θ ./ force_scaling
    jacob[irow_p2+3:irow_p2+5, icol+6:icol+8] .= f_ψ2_F
    jacob[irow_p2+3:irow_p2+5, icol+9:icol+11] .= f_ψ2_M

    # compatability equation jacobian entries for right endpoint
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

# dynamic - general
@inline function insert_dynamic_element_jacobian!(jacob, force_scaling, irow_e,
    irow_e1, irow_p1, irow_e2, irow_p2, icol,
    f_u1_θ, f_u2_θ, f_u1_F, f_u2_F, f_u1_V, f_u2_V, f_u1_Ω, f_u2_Ω,
    f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_V, f_ψ2_V, f_ψ1_Ω, f_ψ2_Ω,
    f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
    f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M,
    f_V_u, f_V_θ, f_V_V,
    f_Ω_θ, f_Ω_Ω)

    jacob = insert_element_jacobian!(jacob, force_scaling, irow_e1, irow_p1,
        irow_e2, irow_p2, icol,
        f_u1_θ, f_u2_θ, f_u1_F, f_u2_F,
        f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M,
        f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M)

    # add equilibrium equation jacobian entries for left endpoint
    jacob[irow_p1:irow_p1+2, icol+12:icol+14] .= f_u1_V  ./ force_scaling
    jacob[irow_p1:irow_p1+2, icol+15:icol+17] .= f_u1_Ω  ./ force_scaling
    jacob[irow_p1+3:irow_p1+5, icol+12:icol+14] .= f_ψ1_V  ./ force_scaling
    jacob[irow_p1+3:irow_p1+5, icol+15:icol+17] .= f_ψ1_Ω  ./ force_scaling

    # add equilibrium equation jacobian entries for right endpoint
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

# dynamic - initial step
@inline function insert_initial_step_element_jacobian!(jacob, force_scaling,
    irow_e, irow_e1, irow_p1, irow_e2, irow_p2, icol,
    f_u1_Vdot, f_u2_Vdot, f_u1_Ωdot, f_u2_Ωdot, f_u1_F, f_u2_F, f_u1_V, f_u2_V, f_u1_Ω, f_u2_Ω,
    f_ψ1_Vdot, f_ψ2_Vdot, f_ψ1_Ωdot, f_ψ2_Ωdot, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_V, f_ψ2_V, f_ψ1_Ω, f_ψ2_Ω,
    f_F1_F, f_F2_F, f_F1_M, f_F2_M,
    f_M1_F, f_M2_F, f_M1_M, f_M2_M,
    f_V_V, f_Ω_Ω)

    # add equilibrium equation jacobian entries for left endpoint
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

    # add compatability equation jacobian entries for left endpoint
    # if irow_e1 == irow_p1 use row corresponding to compatability equations for this beam
    # if irow_e1 <= 0 use row corresponding to compatability equations for the other beam
    # otherwise use row corresponding to compatability equations for this beam
    irow = ifelse(irow_e1 == irow_p1 || irow_e1 <= 0, irow_p1+6, irow_e1)

    jacob[irow:irow+2, icol+6:icol+8] .= f_F1_F .* force_scaling
    jacob[irow:irow+2, icol+9:icol+11] .= f_F1_M .* force_scaling

    jacob[irow+3:irow+5, icol+6:icol+8] .= f_M1_F .* force_scaling
    jacob[irow+3:irow+5, icol+9:icol+11] .= f_M1_M .* force_scaling

    # add equilibrium equation jacobian entries for right endpoint
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

    # add compatability equation jacobian entries for right endpoint
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
    element_jacobian!(jacob, x, ielem, elem, distributed_loads, point_masses, gvec, force_scaling,
        icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2)
    element_jacobian!(jacob, x, ielem, elem, distributed_loads, point_masses, gvec, force_scaling,
        icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2, x0, v0,
        ω0)
    element_jacobian!(jacob, x, ielem, elem, distributed_loads, point_masses, gvec, force_scaling,
        icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2, x0, v0,
        ω0, u0, θ0, udot0, θdot0)
    element_jacobian!(jacob, x, ielem, elem, distributed_loads, point_masses, gvec, force_scaling,
        icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2, x0, v0,
        ω0, udot_init, θdot_init, Vdot_init, Ωdot_init, dt)

Adds a beam element's contributions to the jacobian matrix

There are four implementations corresponding to the following analysis types:
 - Static
 - Dynamic - Steady State
 - Dynamic - Initial Step (for initializing time domain simulations)
 - Dynamic - Time Marching

See "GEBT: A general-purpose nonlinear analysis tool for composite beams" by
Wenbin Yu and "Geometrically nonlinear analysis of composite beams using
Wiener-Milenković parameters" by Qi Wang and Wenbin Yu.

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
 - `irow_e1`: Row index of the first equation for the left side of the beam element
    (a value <= 0 indicates the equations have been eliminated from the system of equations)
 - `irow_p1`: Row index of the first equation for the point on the left side of the beam element
 - `irow_e2`: Row index of the first equation for the right side of the beam element
    (a value <= 0 indicates the equations have been eliminated from the system of equations)
 - `irow_p2`: Row index of the first equation for the point on the right side of the beam element

# Additional Arguments for Dynamic Analyses
 - `x0`: Global frame origin (for the current time step)
 - `v0`: Global frame linear velocity (for the current time step)
 - `ω0`: Global frame angular velocity (for the current time step)

# Additional Arguments for Initial Step Analyses
 - `u0`: initial deflection variables for each beam element
 - `θ0`: initial rotation variables for each beam element
 - `udot0`: initial time derivative of u for each beam element
 - `θdot0`: initial time derivative of θ for each beam element

# Additional Arguments for Time Marching Analyses
 - `udot_init`: `2/dt*u + udot` for each beam element from the previous time step
 - `θdot_init`: `2/dt*θ + θdot` for each beam element from the previous time step
 - `Vdot_init`: `2/dt*V + Vdot` for each beam element from the previous time step
 - `Ωdot_init`: `2/dt*Ω + Ωdot` for each beam element from the previous time step
 - `dt`: time step size
"""
element_jacobian!

# static
@inline function element_jacobian!(jacob, x, ielem, elem, distributed_loads, point_masses, gvec,
    force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2)

    # static
    return static_element_jacobian!(jacob, x, ielem, elem, distributed_loads, point_masses, gvec,
        force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2)
end

# dynamic - steady state
@inline function element_jacobian!(jacob, x, ielem, elem, distributed_loads, point_masses, gvec,
    force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2, x0, v0, ω0)

    return steady_state_element_jacobian!(jacob, x, ielem, elem, distributed_loads, point_masses, gvec,
        force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2, x0, v0, ω0)
end

# dynamic - initial step
@inline function element_jacobian!(jacob, x, ielem, elem, distributed_loads, point_masses, gvec,
    force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2, x0, v0, ω0, u0, θ0,
    udot0, θdot0)

    return initial_step_element_jacobian!(jacob, x, ielem, elem, distributed_loads, point_masses, gvec,
        force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2, x0, v0, ω0, u0, θ0,
        udot0, θdot0)
end

# dynamic - time marching
@inline function element_jacobian!(jacob, x, ielem, elem, distributed_loads, point_masses, gvec,
    force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2, x0, v0, ω0,
    udot_init, θdot_init, Vdot_init, Ωdot_init, dt)

    return newmark_element_jacobian!(jacob, x, ielem, elem, distributed_loads, point_masses, gvec,
        force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2, x0, v0, ω0,
        udot_init, θdot_init, Vdot_init, Ωdot_init, dt)
end

# static
@inline function static_element_jacobian!(jacob, x, ielem, elem, distributed_loads, point_masses, gvec,
    force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2)

    # compute element properties
    ΔL, Ct, Cab, CtCab, u, θ, F, M, γ, κ = static_element_properties(x, icol, elem, force_scaling)

    # pre-calculate jacobian of rotation matrix wrt θ
    C_θ1, C_θ2, C_θ3 = get_C_θ(Ct', θ)
    Ct_θ1, Ct_θ2, Ct_θ3 = C_θ1', C_θ2', C_θ3'

    # solve for the element resultants
    f_u1_θ, f_u2_θ, f_u1_F, f_u2_F,
        f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M,
        f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M = static_element_jacobian_equations(
        elem, ΔL, Cab, CtCab, θ, F, M, γ, κ, Ct_θ1, Ct_θ2, Ct_θ3)

    # add jacobians for follower loads (if applicable)
    if haskey(distributed_loads, ielem)
        f_u1_θ -= mul3(Ct_θ1, Ct_θ2, Ct_θ3, distributed_loads[ielem].f1_follower)
        f_u2_θ -= mul3(Ct_θ1, Ct_θ2, Ct_θ3, distributed_loads[ielem].f2_follower)
        f_ψ1_θ -= mul3(Ct_θ1, Ct_θ2, Ct_θ3, distributed_loads[ielem].m1_follower)
        f_ψ2_θ -= mul3(Ct_θ1, Ct_θ2, Ct_θ3, distributed_loads[ielem].m2_follower)
    end

    # add jacobians for distributed gravitational loads
    fg1_θ, fg2_θ, mg1_θ, mg2_θ = element_gravitational_loads_jacobian(ΔL, Cab, CtCab, Ct_θ1, Ct_θ2, Ct_θ3, elem, gvec)
    f_u1_θ -= fg1_θ
    f_u2_θ -= fg2_θ
    f_ψ1_θ -= mg1_θ
    f_ψ2_θ -= mg2_θ

    # add jacobians for point mass gravitational loads
    if haskey(point_masses, ielem)
        fp1_θ, fp2_θ, mp1_θ, mp2_θ = point_mass_gravitational_loads_jacobian(Ct, Ct_θ1, Ct_θ2, Ct_θ3, point_masses[ielem].mass, gvec)

        f_u1 -= fp1_θ
        f_u2 -= fp2_θ
        f_ψ1 -= mp1_θ
        f_ψ2 -= mp2_θ
    end

    # insert element resultants into the jacobian matrix
    jacob = insert_static_element_jacobian!(jacob, force_scaling, irow_e1,
        irow_p1, irow_e2, irow_p2, icol,
        f_u1_θ, f_u2_θ, f_u1_F, f_u2_F,
        f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M,
        f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M)

    return jacob
end

# dynamic - steady state
@inline function steady_state_element_jacobian!(jacob, x, ielem, elem, distributed_loads, point_masses, gvec,
    force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2, x0, v0, ω0)

    # compute element properties
    ΔL, Ct, Cab, CtCab, u, θ, F, M, γ, κ, v, ω, P, H, V, Q = steady_state_element_properties(
        x, icol, elem, force_scaling, x0, v0, ω0)

    # pre-calculate jacobian of rotation matrix wrt θ
    C_θ1, C_θ2, C_θ3 = get_C_θ(Ct', θ)
    Ct_θ1, Ct_θ2, Ct_θ3 = C_θ1', C_θ2', C_θ3'

    # solve for the element resultants
    f_u1_θ, f_u2_θ, f_u1_F, f_u2_F, f_u1_V, f_u2_V, f_u1_Ω, f_u2_Ω,
        f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_V, f_ψ2_V, f_ψ1_Ω, f_ψ2_Ω,
        f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M,
        f_V_u, f_V_θ, f_V_V,
        f_Ω_θ, f_Ω_Ω = steady_state_element_jacobian_equations(elem, ΔL, Cab, CtCab,
        θ, F, M, γ, κ, ω, P, H, V, Ct_θ1, Ct_θ2, Ct_θ3)

    # add jacobians for follower loads (if applicable)
    if haskey(distributed_loads, ielem)
        f_u1_θ -= mul3(Ct_θ1, Ct_θ2, Ct_θ3, distributed_loads[ielem].f1_follower)
        f_u2_θ -= mul3(Ct_θ1, Ct_θ2, Ct_θ3, distributed_loads[ielem].f2_follower)
        f_ψ1_θ -= mul3(Ct_θ1, Ct_θ2, Ct_θ3, distributed_loads[ielem].m1_follower)
        f_ψ2_θ -= mul3(Ct_θ1, Ct_θ2, Ct_θ3, distributed_loads[ielem].m2_follower)
    end

    # add jacobians for distributed gravitational loads
    fg1_θ, fg2_θ, mg1_θ, mg2_θ = element_gravitational_loads_jacobian(ΔL, Cab, CtCab, Ct_θ1, Ct_θ2, Ct_θ3, elem, gvec)
    f_u1_θ -= fg1_θ
    f_u2_θ -= fg2_θ
    f_ψ1_θ -= mg1_θ
    f_ψ2_θ -= mg2_θ

    # add jacobians for point mass gravitational loads
    if haskey(point_masses, ielem)
        fp1_θ, fp2_θ, mp1_θ, mp2_θ = point_mass_gravitational_loads_jacobian(Ct, Ct_θ1, Ct_θ2, Ct_θ3, point_masses[ielem].mass, gvec)

        f_u1 -= fp1_θ
        f_u2 -= fp2_θ
        f_ψ1 -= mp1_θ
        f_ψ2 -= mp2_θ
    end

    # insert element resultants into the jacobian matrix
    jacob = insert_dynamic_element_jacobian!(jacob, force_scaling, irow_e,
        irow_e1, irow_p1, irow_e2, irow_p2, icol,
        f_u1_θ, f_u2_θ, f_u1_F, f_u2_F, f_u1_V, f_u2_V, f_u1_Ω, f_u2_Ω,
        f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_V, f_ψ2_V, f_ψ1_Ω, f_ψ2_Ω,
        f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M,
        f_V_u, f_V_θ, f_V_V,
        f_Ω_θ, f_Ω_Ω)

    return jacob
end

# dynamic - initial step
@inline function initial_step_element_jacobian!(jacob, x, ielem, elem, distributed_loads, point_masses, gvec,
    force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2, x0, v0, ω0, u0, θ0,
    udot0, θdot0)

    # compute element properties
    ΔL, Ct, Cab, CtCab, u, θ, F, M, γ, κ, v, ω, P, H, V, Ω, udot, θdot,
        CtCabPdot, CtCabHdot, CtCabdot = initial_step_element_properties(x, icol, elem,
        force_scaling, x0, v0, ω0, u0, θ0, udot0, θdot0)

    # solve for the element resultants
    f_u1_Vdot, f_u2_Vdot, f_u1_Ωdot, f_u2_Ωdot, f_u1_F, f_u2_F, f_u1_V, f_u2_V, f_u1_Ω, f_u2_Ω,
        f_ψ1_Vdot, f_ψ2_Vdot, f_ψ1_Ωdot, f_ψ2_Ωdot, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_V, f_ψ2_V, f_ψ1_Ω, f_ψ2_Ω,
        f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_F, f_M2_F, f_M1_M, f_M2_M,
        f_V_V, f_Ω_Ω = initial_step_element_jacobian_equations(elem, ΔL, Cab, CtCab, θ, F, γ,
        ω, P, V, CtCabdot)

    # insert element resultants into the jacobian matrix
    jacob = insert_initial_step_element_jacobian!(jacob, force_scaling, irow_e,
        irow_e1, irow_p1, irow_e2, irow_p2, icol,
        f_u1_Vdot, f_u2_Vdot, f_u1_Ωdot, f_u2_Ωdot, f_u1_F, f_u2_F, f_u1_V, f_u2_V, f_u1_Ω, f_u2_Ω,
        f_ψ1_Vdot, f_ψ2_Vdot, f_ψ1_Ωdot, f_ψ2_Ωdot, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_V, f_ψ2_V, f_ψ1_Ω, f_ψ2_Ω,
        f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_F, f_M2_F, f_M1_M, f_M2_M,
        f_V_V, f_Ω_Ω)

    return jacob
end

# dynamic - newmark scheme time marching
@inline function newmark_element_jacobian!(jacob, x, ielem, elem, distributed_loads, point_masses, gvec,
    force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2, x0, v0, ω0,
    udot_init, θdot_init, Vdot_init, Ωdot_init, dt)

    # compute element properties
    ΔL, Ct, Cab, CtCab, u, θ, F, M, γ, κ, v, ω, P, H, V, Ω, udot, θdot, Vdot, Ωdot, Pdot, Hdot,
        CtCabPdot, CtCabHdot, CtCabdot = newmark_element_properties(x, icol, elem,
        force_scaling, x0, v0, ω0, udot_init, θdot_init, Vdot_init,
        Ωdot_init, dt)

    # pre-calculate jacobian of rotation matrix wrt θ
    C_θ1, C_θ2, C_θ3 = get_C_θ(Ct', θ)
    Ct_θ1, Ct_θ2, Ct_θ3 = C_θ1', C_θ2', C_θ3'

    Cdot_θ1, Cdot_θ2, Cdot_θ3 = get_C_t_θ(θ, θdot)
    Ctdot_θ1, Ctdot_θ2, Ctdot_θ3 = Cdot_θ1', Cdot_θ2', Cdot_θ3'

    Cdot_θdot1, Cdot_θdot2, Cdot_θdot3 = get_C_θdot(Ct', θ)
    Ctdot_θdot1, Ctdot_θdot2, Ctdot_θdot3 = Cdot_θdot1', Cdot_θdot2', Cdot_θdot3'

    # solve for the element resultants
    f_u1_θ, f_u2_θ, f_u1_F, f_u2_F, f_u1_V, f_u2_V, f_u1_Ω, f_u2_Ω,
        f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_V, f_ψ2_V, f_ψ1_Ω, f_ψ2_Ω,
        f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M,
        f_V_u, f_V_θ, f_V_V,
        f_Ω_θ, f_Ω_Ω = newmark_element_jacobian_equations(elem, ΔL, Cab, CtCab,
        θ, F, M, γ, κ, ω, P, H, V, Ω, θdot, Vdot, Ωdot, Pdot, Hdot, dt, Ct_θ1, Ct_θ2, Ct_θ3, CtCabdot,
        Ctdot_θ1, Ctdot_θ2, Ctdot_θ3, Ctdot_θdot1, Ctdot_θdot2, Ctdot_θdot3)

    # add jacobians for follower loads (if applicable)
    if haskey(distributed_loads, ielem)
        f_u1_θ -= mul3(Ct_θ1, Ct_θ2, Ct_θ3, distributed_loads[ielem].f1_follower)
        f_u2_θ -= mul3(Ct_θ1, Ct_θ2, Ct_θ3, distributed_loads[ielem].f2_follower)
        f_ψ1_θ -= mul3(Ct_θ1, Ct_θ2, Ct_θ3, distributed_loads[ielem].m1_follower)
        f_ψ2_θ -= mul3(Ct_θ1, Ct_θ2, Ct_θ3, distributed_loads[ielem].m2_follower)
    end

    # add jacobians for distributed gravitational loads
    fg1_θ, fg2_θ, mg1_θ, mg2_θ = element_gravitational_loads_jacobian(ΔL, Cab, CtCab, Ct_θ1, Ct_θ2, Ct_θ3, elem, gvec)
    f_u1_θ -= fg1_θ
    f_u2_θ -= fg2_θ
    f_ψ1_θ -= mg1_θ
    f_ψ2_θ -= mg2_θ

    # add jacobians for point mass gravitational loads
    if haskey(point_masses, ielem)
        fp1_θ, fp2_θ, mp1_θ, mp2_θ = point_mass_gravitational_loads_jacobian(Ct, Ct_θ1, Ct_θ2, Ct_θ3, point_masses[ielem].mass, gvec)

        f_u1 -= fp1_θ
        f_u2 -= fp2_θ
        f_ψ1 -= mp1_θ
        f_ψ2 -= mp2_θ
    end

    # insert element resultants into the jacobian matrix
    jacob = insert_dynamic_element_jacobian!(jacob, force_scaling, irow_e,
        irow_e1, irow_p1, irow_e2, irow_p2, icol,
        f_u1_θ, f_u2_θ, f_u1_F, f_u2_F, f_u1_V, f_u2_V, f_u1_Ω, f_u2_Ω,
        f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_V, f_ψ2_V, f_ψ1_Ω, f_ψ2_Ω,
        f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M,
        f_V_u, f_V_θ, f_V_V,
        f_Ω_θ, f_Ω_Ω)

    return jacob
end

# dynamic - general
@inline function dynamic_element_jacobian!(jacob, x, ielem, elem, distributed_loads, point_masses, gvec,
    force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2, x0, v0, ω0, udot, θdot,
    Vdot, Ωdot)

    # compute element properties
    ΔL, Ct, Cab, CtCab, u, θ, F, M, γ, κ, v, ω, P, H, V, Ω, udot, θdot,
        CtCabPdot, CtCabHdot, CtCabdot = dynamic_element_properties(x, icol, elem,
        force_scaling, x0, v0, ω0, udot, θdot, Vdot, Ωdot)

    Pdot = element_linear_momentum(elem, Vdot, Ωdot)
    Hdot = element_angular_momentum(elem, Vdot, Ωdot)

    # pre-calculate jacobian of rotation matrix wrt θ
    C_θ1, C_θ2, C_θ3 = get_C_θ(Ct', θ)
    Ct_θ1, Ct_θ2, Ct_θ3 = C_θ1', C_θ2', C_θ3'

    Cdot_θ1, Cdot_θ2, Cdot_θ3 = get_C_t_θ(θ, θdot)
    Ctdot_θ1, Ctdot_θ2, Ctdot_θ3 = Cdot_θ1', Cdot_θ2', Cdot_θ3'

    # solve for the element resultants
    f_u1_θ, f_u2_θ, f_u1_F, f_u2_F, f_u1_V, f_u2_V, f_u1_Ω, f_u2_Ω,
        f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_V, f_ψ2_V, f_ψ1_Ω, f_ψ2_Ω,
        f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M,
        f_V_u, f_V_θ, f_V_V,
        f_Ω_θ, f_Ω_Ω = dynamic_element_jacobian_equations(elem, ΔL, Cab, CtCab,
        θ, F, M, γ, κ, ω, P, H, V, θdot, Pdot, Hdot, Ct_θ1, Ct_θ2, Ct_θ3, CtCabdot,
        Ctdot_θ1, Ctdot_θ2, Ctdot_θ3)

    # add jacobians for follower loads (if applicable)
    if haskey(distributed_loads, ielem)
        f_u1_θ -= mul3(Ct_θ1, Ct_θ2, Ct_θ3, distributed_loads[ielem].f1_follower)
        f_u2_θ -= mul3(Ct_θ1, Ct_θ2, Ct_θ3, distributed_loads[ielem].f2_follower)
        f_ψ1_θ -= mul3(Ct_θ1, Ct_θ2, Ct_θ3, distributed_loads[ielem].m1_follower)
        f_ψ2_θ -= mul3(Ct_θ1, Ct_θ2, Ct_θ3, distributed_loads[ielem].m2_follower)
    end

    # add jacobians for distributed gravitational loads
    fg1_θ, fg2_θ, mg1_θ, mg2_θ = element_gravitational_loads_jacobian(ΔL, Cab, CtCab, Ct_θ1, Ct_θ2, Ct_θ3, elem, gvec)
    f_u1_θ -= fg1_θ
    f_u2_θ -= fg2_θ
    f_ψ1_θ -= mg1_θ
    f_ψ2_θ -= mg2_θ

    # add jacobians for point mass gravitational loads
    if haskey(point_masses, ielem)
        fp1_θ, fp2_θ, mp1_θ, mp2_θ = point_mass_gravitational_loads_jacobian(Ct, Ct_θ1, Ct_θ2, Ct_θ3, point_masses[ielem].mass, gvec)

        f_u1 -= fp1_θ
        f_u2 -= fp2_θ
        f_ψ1 -= mp1_θ
        f_ψ2 -= mp2_θ
    end

    # insert element resultants into the jacobian matrix
    jacob = insert_dynamic_element_jacobian!(jacob, force_scaling, irow_e,
        irow_e1, irow_p1, irow_e2, irow_p2, icol,
        f_u1_θ, f_u2_θ, f_u1_F, f_u2_F, f_u1_V, f_u2_V, f_u1_Ω, f_u2_Ω,
        f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_V, f_ψ2_V, f_ψ1_Ω, f_ψ2_Ω,
        f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
        f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M,
        f_V_u, f_V_θ, f_V_V,
        f_Ω_θ, f_Ω_Ω)

    return jacob
end

"""
    element_mass_matrix_properties(x, icol, elem)

Extract/Compute the properties needed for mass matrix construction: `ΔL`, `Ct`,
`Cab`, `CtCab`, `θ`, `P`, `H`, `Ctdot_cdot1`, `Ctdot_cdot2`, and Ctdot_cdot3
"""
@inline function element_mass_matrix_properties(x, icol, elem)

    ΔL = elem.L
    θ = SVector(x[icol+3 ], x[icol+4 ], x[icol+5 ])
    V = SVector(x[icol+12], x[icol+13], x[icol+14])
    Ω = SVector(x[icol+15], x[icol+16], x[icol+17])
    P = element_linear_momentum(elem, V, Ω)
    H = element_angular_momentum(elem, V, Ω)
    C = get_C(θ)
    Ct = C'
    Cab = elem.Cab
    CtCab = Ct*Cab
    Cdot_cdot1, Cdot_cdot2, Cdot_cdot3 = get_C_θdot(C, θ)
    Ctdot_cdot1, Ctdot_cdot2, Ctdot_cdot3 = Cdot_cdot1', Cdot_cdot2', Cdot_cdot3'

    return ΔL, Ct, Cab, CtCab, θ, P, H, Ctdot_cdot1, Ctdot_cdot2, Ctdot_cdot3
end

"""
    element_mass_matrix_equations(elem, ΔL, Ct, Cab, CtCab, θ, P, H)

Calculates the jacobians of the nonlinear equations for a beam element with
respect to the time derivatives of the state variables.

See "GEBT: A general-purpose nonlinear analysis tool for composite beams" by
Wenbin Yu and "Geometrically nonlinear analysis of composite beams using
Wiener-Milenković parameters" by Qi Wang and Wenbin Yu.

# Arguments:
 - `elem`: Beam element
 - `ΔL`: Length of the beam element
 - `Ct`: Rotation tensor of the beam deformation in the "a" frame, transposed
 - `Cab`: Direction cosine matrix from "a" to "b" frame for the element
 - `θ`: Rotation variables for the element [θ1, θ2, θ3]
 - `P`: Linear momenta for the element [P1, P2, P3]
 - `H`: Angular momenta for the element [H1, H2, H3]
"""
@inline function element_mass_matrix_equations(elem, ΔL, Ct, Cab, CtCab, θ, P, H,
    Ctdot_cdot1, Ctdot_cdot2, Ctdot_cdot3)

    M11 = element_linear_momentum_V(elem)
    M12 = element_linear_momentum_Ω(elem)
    M21 = element_angular_momentum_V(elem)
    M22 = element_angular_momentum_Ω(elem)

    tmp = ΔL/2*mul3(Ctdot_cdot1, Ctdot_cdot2, Ctdot_cdot3, Cab*P)
    f_u1_θdot = tmp
    f_u2_θdot = tmp

    tmp = ΔL/2*CtCab*M11
    f_u1_Vdot = tmp
    f_u2_Vdot = tmp

    tmp = ΔL/2*CtCab*M12
    f_u1_Ωdot = tmp
    f_u2_Ωdot = tmp

    tmp = ΔL/2*mul3(Ctdot_cdot1, Ctdot_cdot2, Ctdot_cdot3, Cab*H)
    f_ψ1_θdot = tmp
    f_ψ2_θdot = tmp

    tmp = ΔL/2*CtCab*M21
    f_ψ1_Vdot = tmp
    f_ψ2_Vdot = tmp

    tmp = ΔL/2*CtCab*M22
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
    insert_element_mass_matrix!(jacob, force_scaling, irow_e1,
        irow_p1, irow_e2, irow_p2, icol, 
        f_u1_θdot, f_u2_θdot, f_u1_Vdot, f_u2_Vdot, f_u1_Ωdot, f_u2_Ωdot,
        f_ψ1_θdot, f_ψ2_θdot, f_ψ1_Vdot, f_ψ2_Vdot, f_ψ1_Ωdot, f_ψ2_Ωdot,
        f_V_udot, f_Ω_θdot)

Insert the beam element's contributions into the "mass matrix": the jacobian of the
residual equations with respect to the time derivatives of the state variables

# Arguments
 - `jacob`: System mass matrix
 - `force_scaling`: Scaling parameter for forces/moments
 - `irow_e1`: Row index of the first equation for the left side of the beam
 - `irow_p1`: Row index of the first equation for the point on the left side of the beam
 - `irow_e2`: Row index of the first equation for the right side of the beam
 - `irow_p2`: Row index of the first equation for the point on the right side of the beam
 - `icol`: Row/Column index corresponding to the first beam state variable

All other arguments use the following naming convention:
 - `f_y_x`: Jacobian of element equation "y" with respect to state variable "x"
"""
@inline function insert_element_mass_matrix!(jacob, force_scaling, irow_e, irow_e1, irow_p1,
    irow_e2, irow_p2, icol, 
    f_u1_θdot, f_u2_θdot, f_u1_Vdot, f_u2_Vdot, f_u1_Ωdot, f_u2_Ωdot,
    f_ψ1_θdot, f_ψ2_θdot, f_ψ1_Vdot, f_ψ2_Vdot, f_ψ1_Ωdot, f_ψ2_Ωdot,
    f_V_udot, f_Ω_θdot)

    # create jacobian entries for left endpoint
    jacob[irow_p1:irow_p1+2, icol+3:icol+5] .= f_u1_θdot ./ force_scaling
    jacob[irow_p1+3:irow_p1+5, icol+3:icol+5] .= f_ψ1_θdot ./ force_scaling
    jacob[irow_p1:irow_p1+2, icol+12:icol+14] .= f_u1_Vdot  ./ force_scaling
    jacob[irow_p1:irow_p1+2, icol+15:icol+17] .= f_u1_Ωdot  ./ force_scaling
    jacob[irow_p1+3:irow_p1+5, icol+12:icol+14] .= f_ψ1_Vdot ./ force_scaling
    jacob[irow_p1+3:irow_p1+5, icol+15:icol+17] .= f_ψ1_Ωdot ./ force_scaling

    # create jacobian entries for right endpoint
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
    insert_element_mass_matrix!(jacob, gamma, force_scaling,
        irow_e1, irow_p1, irow_e2, irow_p2, icol, 
        f_u1_θdot, f_u2_θdot, f_u1_Vdot, f_u2_Vdot, f_u1_Ωdot, f_u2_Ωdot,
        f_ψ1_θdot, f_ψ2_θdot, f_ψ1_Vdot, f_ψ2_Vdot, f_ψ1_Ωdot, f_ψ2_Ωdot,
        f_V_udot, f_Ω_θdot

Add the beam element's mass matrix to the system jacobian matrix `jacob`, scaled
by the scaling parameter `gamma`.
"""
@inline function insert_element_mass_matrix!(jacob, gamma, force_scaling, irow_e,
    irow_e1, irow_p1, irow_e2, irow_p2, icol, 
    f_u1_θdot, f_u2_θdot, f_u1_Vdot, f_u2_Vdot, f_u1_Ωdot, f_u2_Ωdot,
    f_ψ1_θdot, f_ψ2_θdot, f_ψ1_Vdot, f_ψ2_Vdot, f_ψ1_Ωdot, f_ψ2_Ωdot,
    f_V_udot, f_Ω_θdot)

    # create jacobian entries for left endpoint
    jacob[irow_p1:irow_p1+2, icol+3:icol+5] .= f_u1_θdot .* (gamma/force_scaling)
    jacob[irow_p1+3:irow_p1+5, icol+3:icol+5] .= f_ψ1_θdot .* (gamma/force_scaling)
    jacob[irow_p1:irow_p1+2, icol+12:icol+14] .= f_u1_Vdot  .* (gamma/force_scaling)
    jacob[irow_p1:irow_p1+2, icol+15:icol+17] .= f_u1_Ωdot  .* (gamma/force_scaling)
    jacob[irow_p1+3:irow_p1+5, icol+12:icol+14] .= f_ψ1_Vdot .* (gamma/force_scaling)
    jacob[irow_p1+3:irow_p1+5, icol+15:icol+17] .= f_ψ1_Ωdot .* (gamma/force_scaling)

    # create jacobian entries for right endpoint
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
    element_mass_matrix!(jacob, x, elem, force_scaling, icol,
        irow_e, irow_e1, irow_p1, irow_e2, irow_p2)

Insert the beam element's contributions to the "mass matrix": the jacobian of the
residual equations with respect to the time derivatives of the state variables

See "GEBT: A general-purpose nonlinear analysis tool for composite beams" by
Wenbin Yu and "Geometrically nonlinear analysis of composite beams using
Wiener-Milenković parameters" by Qi Wang and Wenbin Yu.

# Arguments
 - `jacob`: System jacobian matrix to add mass matrix jacobian to
 - `x`: current state vector
 - `beam`: beam element
 - `force_scaling`: scaling parameter for forces/moments
 - `icol`: starting index for the beam's state variables
 - `irow_e1`: Row index of the first equation for the left side of the beam element
    (a value <= 0 indicates the equations have been eliminated from the system of equations)
 - `irow_p1`: Row index of the first equation for the point on the left side of the beam element
 - `irow_e2`: Row index of the first equation for the right side of the beam element
    (a value <= 0 indicates the equations have been eliminated from the system of equations)
 - `irow_p2`: Row index of the first equation for the point on the right side of the beam element
 - `gamma`: Scaling parameter for scaling mass matrix contribution to `jacob`
"""
@inline function element_mass_matrix!(jacob, x, elem, force_scaling,
    icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2)

    # get beam element properties
    ΔL, Ct, Cab, CtCab, θ, P, H, Ctdot_cdot1, Ctdot_cdot2, Ctdot_cdot3 =
        element_mass_matrix_properties(x, icol, elem)

    # get jacobians of beam element equations
    f_u1_θdot, f_u2_θdot, f_u1_Vdot, f_u2_Vdot, f_u1_Ωdot, f_u2_Ωdot,
        f_ψ1_θdot, f_ψ2_θdot, f_ψ1_Vdot, f_ψ2_Vdot, f_ψ1_Ωdot, f_ψ2_Ωdot,
        f_V_udot, f_Ω_θdot = element_mass_matrix_equations(elem, ΔL, Ct, Cab, CtCab,
        θ, P, H, Ctdot_cdot1, Ctdot_cdot2, Ctdot_cdot3)

    # initialize/insert into jacobian matrix for the system
    jacob = insert_element_mass_matrix!(jacob, force_scaling, irow_e, irow_e1, irow_p1,
        irow_e2, irow_p2, icol, f_u1_θdot, f_u2_θdot, f_u1_Vdot, f_u2_Vdot, f_u1_Ωdot, f_u2_Ωdot,
        f_ψ1_θdot, f_ψ2_θdot, f_ψ1_Vdot, f_ψ2_Vdot, f_ψ1_Ωdot, f_ψ2_Ωdot,
        f_V_udot, f_Ω_θdot)

    return jacob
end

"""
    element_mass_matrix!(jacob, gamma, x, elem, force_scaling,
        icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2)

Add the beam element's mass matrix to the system jacobian matrix `jacob`, scaled
by the scaling parameter `gamma`.
"""
@inline function element_mass_matrix!(jacob, gamma, x, elem, force_scaling,
    icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2)

    # get beam element properties
    ΔL, Ct, Cab, CtCab, θ, P, H, Ctdot_cdot1, Ctdot_cdot2, Ctdot_cdot3 =
        element_mass_matrix_properties(x, icol, elem)

    # get jacobians of beam element equations
    f_u1_θdot, f_u2_θdot, f_u1_Vdot, f_u2_Vdot, f_u1_Ωdot, f_u2_Ωdot,
        f_ψ1_θdot, f_ψ2_θdot, f_ψ1_Vdot, f_ψ2_Vdot, f_ψ1_Ωdot, f_ψ2_Ωdot,
        f_V_udot, f_Ω_θdot = element_mass_matrix_equations(elem, ΔL, Ct, Cab, CtCab,
        θ, P, H, Ctdot_cdot1, Ctdot_cdot2, Ctdot_cdot3)

    # initialize/insert into jacobian matrix for the system
    jacob = insert_element_mass_matrix!(jacob, gamma, force_scaling,
        irow_e, irow_e1, irow_p1, irow_e2, irow_p2, icol, 
        f_u1_θdot, f_u2_θdot, f_u1_Vdot, f_u2_Vdot, f_u1_Ωdot, f_u2_Ωdot,
        f_ψ1_θdot, f_ψ2_θdot, f_ψ1_Vdot, f_ψ2_Vdot, f_ψ1_Ωdot, f_ψ2_Ωdot,
        f_V_udot, f_Ω_θdot)

    return jacob
end
