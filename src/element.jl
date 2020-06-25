struct Element{TF}
	ΔL::TF
	x::SVector{3,TF}
	pt1::Int
	pt2::Int
	C11::SMatrix{3,3,TF,9}
	C12::SMatrix{3,3,TF,9}
	C22::SMatrix{3,3,TF,9}
	m11::SMatrix{3,3,TF,9}
	m12::SMatrix{3,3,TF,9}
	m22::SMatrix{3,3,TF,9}
	Cab::SMatrix{3,3,TF,9}
end

Element(ΔL, x, pt1, pt2, C, m, Cab) = Element(ΔL, x, pt1, pt2, C[1:3,1:3],
	C[1:3,4:6], C[4:6,4:6], m[1:3,1:3], m[1:3,4:6], m[4:6,4:6], Cab)

strain(element, F, M) = element.C11*F + element.C12*M
curvature(element, F, M) = element.C12'*F + element.C22*M
linear_velocity(element, P, H) = element.m11*P + element.m12*H
angular_velocity(element, P, H) = element.m12'*P + element.m22*H

function element_state_variables(x, icol_beam, ibeam, ::Steady)
	icol = icol_beam[ibeam]
	u = SVector(x[icol  ], x[icol+1 ], x[icol+2 ])
	ψ = SVector(x[icol+3], x[icol+4 ], x[icol+5 ])
	F = SVector(x[icol+6], x[icol+7 ], x[icol+8 ])
	M = SVector(x[icol+9], x[icol+10], x[icol+11])
	return u, ψ, F, M
end

function element_state_variables(x, icol_beam, ibeam, ::Unsteady)

	icol = icol_beam[ibeam]
	u = SVector(x[icol   ], x[icol+1 ], x[icol+2 ])
	ψ = SVector(x[icol+3 ], x[icol+4 ], x[icol+5 ])
	F = SVector(x[icol+6 ], x[icol+7 ], x[icol+8 ])
	M = SVector(x[icol+9 ], x[icol+10], x[icol+11])
	P = SVector(x[icol+12], x[icol+13], x[icol+14])
	H = SVector(x[icol+15], x[icol+16], x[icol+17])

	return u, ψ, F, M, P, H
end

"""
	element_properties(element, x0, v0, ω0, u, θ, F, M)
	element_properties(element, x0, v0, ω0, u, θ, F, M, P, H)
	element_properties(element, x0, v0, ω0, u, θ, F, M, P, H, udot_p, θdot_p, CtCabPdot_p, CtCabHdot_p)
	element_properties(element, x0, v0, ω0, u, θ, F, M, P, H, u_p, θ_p, P_p, H_p, udot_p, θdot_p, CtCabPdot_p, CtCabHdot_p, dt)

Calculates element properties.

# Arguments
 - element: Element from which properties are extracted
 - x0: Global frame origin
 - v0: Global frame linear velocity
 - ω0: Global frame angular velocity
 - u: Displacement variables for the element [u1, u2, u3]
 - θ: Rotation variables for the element [θ1, θ2, θ3]
 - F: Force variables for the element [F1, F2, F3]
 - M: Moment variables for the element [M1, M2, M3]
 - P: Linear momenta for the element [P1, P2, P3] (Only for unsteady simulations)
 - H: Angular momenta for the element [H1, H2, H3] (Only for unsteady simulations)
 - u_p: `u` from previous time step (Only for time-marching simulations)
 - θ_p: `θ` from previous time step (Only for time-marching simulations)
 - P_p: `P` from previous time step (Only for time-marching simulations)
 - H_p: `H` from previous time step (Only for time-marching simulations)
 - udot_p: `udot` from previous time step (Only for time-marching simulations)
 - θdot_p: `θdot` from previous time step (Only for time-marching simulations)
 - CtCabPdot_p: `C'*Cab*Pdot` from previous time step (Only for time-marching simulations)
 - CtCabHdot_p: `C'*Cab*Hdot` from previous time step (Only for time-marching simulations)
 - dt: time step size
"""
element_properties

# for steady simulations
function element_properties(element, x0, v0, ω0, u, θ, F, M)

	ΔL = element.ΔL
	Ct = wiener_milenkovic(θ)
	Cab = element.Cab
	CtCab = Ct*Cab
	v = v0 + cross(ω0, element.x - x0)
	ω = ω0
	γ = strain(element, F, M)
	κ = curvature(element, F, M)

	return ΔL, Ct, Cab, CtCab, u, θ, F, M, v, ω, γ, κ
end

# for unsteady simulations
function element_properties(element, x0, v0, ω0, u, θ, F, M, P, H)

	ΔL, Ct, Cab, CtCab, u, θ, F, M, v, ω, γ, κ = steady_element_properties(element, x0, v0, ω0, u, θ, F, M)

	V = linear_velocity(element, P)
	Ω = angular_velocity(element, H)

	return ΔL, Ct, Cab, CtCab, u, θ, F, M, P, H, v, ω, γ, κ, V, Ω
end

# for time-marching simulations
function element_properties(element, x0, v0, ω0, u, θ, F, M, P, H,
	u_p, θ_p, P_p, H_p, udot_p, θdot_p, CtCabPdot_p, CtCabHdot_p, dt)

	ΔL, Ct, Cab, CtCab, u, θ, F, M, P, H, v, ω, γ, κ, V, Ω = element_properties(element, x0, v0, ω0, u, θ, F, M, P, H)

	udot = -2/dt*u_p - udot_p
	θdot = -2/dt*(θ - θ_p) + θdot_p
	CtCabPdot = -2/dt*CtCabP_p - CtCabPdot_p
	CtCabHdot = -2/dt*CtCabH_p - CtCabHdot_p

	return ΔL, Ct, Cab, CtCab, u, θ, F, M, v, ω, γ, κ, P, H, V, Ω, udot, θdot, CtCabPdot, CtCabHdot
end

"""
	element_equations(ΔL, Ct, Cab, CtCab, u, θ, F, M, v, ω, γ, κ)
	element_equations(ΔL, Ct, Cab, CtCab, u, θ, F, M, v, ω, γ, κ, P, H, V, Ω)
    element_equations(ΔL, Ct, Cab, CtCab, u, θ, F, M, v, ω, γ, κ, P, H, V, Ω,
		udot, θdot, CtCabPdot, CtCabHdot)
	element_equations(ΔL, Ct, Cab, CtCab, u, θ, F, M, v, ω, γ, κ, P, H, V, Ω,
		udot, θdot, CtCabPdot, CtCabHdot, dt)

Evaluates nonlinear equations for an element using shape functions of the lowest
possible order in order to avoid numerical integration.

See "GEBT: A general-purpose nonlinear analysis tool for composite beams" by
Wenbin Yu

The Rodriguez parameters from the original paper have been swapped out for
Wiener-Milenkovic parameters. See "Geometrically nonlinear analysis of composite
beams using Wiener-Milenković parameters" by Qi Wang and Wenbin Yu.

# Arguments:
 - ΔL: Length of the element
 - Ct: Rotation tensor of the beam deformation in the "a" frame, transposed
 - Cab: Direction cosine matrix from "a" to "b" frame for the element
 - CtCab: Ct*Cab, precomputed for efficiency
 - u: Displacement variables for the element [u1, u2, u3]
 - θ: Rotation variables for the element [θ1, θ2, θ3]
 - F: Force variables for the element [F1, F2, F3]
 - M: Moment variables for the element [M1, M2, M3]
 - v: Linear velocity of element in global frame "a" [v1, v2, v3]
 - ω: Angular velocity of element in global frame "a" [ω1, ω2, ω3]
 - γ: Engineering strains in the element [γ11, 2γ12, 2γ13]
 - κ: Curvatures in the element [κ1, κ2, κ3]
 - P: Linear momenta for the element [P1, P2, P3] (unsteady simulations only)
 - H: Angular momenta for the element [H1, H2, H3] (unsteady simulations only)
 - V: Velocity of the element (unsteady simulations only)
 - Ω: Rotational velocity of the element (unsteady simulations only)
 - udot: time derivative of u, evaluated at the previous iteration (time-marching simulations only)
 - θdot: time derivative of θ, evaluated at the previous iteration (time-marching simulations only)
 - CtCabPdot: C'*Cab*Pdot evaluated at the previous iteration (time-marching simulations only)
 - CtCabHdot: C'*Cab*Hdot evaluated at the previous iteration (time-marching simulations only)
 - dt: time step size (for time-marching simulations, during time marching)
"""
element_equations

# for steady simulations
function element_equations(ΔL, Ct, Cab, CtCab, u, θ, F, M, v, ω, γ, κ)

	e1 = SVector(1, 0, 0)
	θtilde = tilde(θ)

	#TODO: is ∓ or ± correct? currently using ∓
    f_u1 = -CtCab*F
    f_u2 =  CtCab*F
    f_ψ1 = -CtCab*M + ΔL/2*CtCab*cross(e1 + γ, F)
    f_ψ2 =  CtCab*M + ΔL/2*CtCab*cross(e1 + γ, F)
    f_F1 =  u - ΔL/2*(CtCab*(e1 + γ) - Cab*e1)
    f_F2 = -u - ΔL/2*(CtCab*(e1 + γ) - Cab*e1)
    f_M1 =  θ - ΔL/2*(θ'*θ*I/16 - θtilde/2 - θ*θ'/8)*Cab*κ
    f_M2 = -θ - ΔL/2*(θ'*θ*I/16 - θtilde/2 - θ*θ'/8)*Cab*κ

    if loads
        f1, f2, m1, m2 = element_loads(ΔL, Ct, load, follower)
        f_u1 -= f1
        f_u2 -= f2
        f_ψ1 -= m1
        f_ψ2 -= m2
    end

	return f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2
end

# for unsteady simulations
function element_equations(ΔL, Ct, Cab, CtCab, u, θ, F, M, v, ω, γ, κ, P, H, V, Ω)

    f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2 = steady_element_equations(
		ΔL, Ct, Cab, CtCab, u, θ, F, M, v, ω, γ, κ)

    f_u1 += cross(ω, ΔL/2*CtCab*P)
    f_u2 += cross(ω, ΔL/2*CtCab*P)

    f_ψ1 += cross(ω, ΔL/2*CtCab*H) + ΔL/2*CtCab*cross(V, P)
    f_ψ2 += cross(ω, ΔL/2*CtCab*H) + ΔL/2*CtCab*cross(V, P)

    f_P = CtCab*V - v - cross(ω, u)
	f_H = Ω - CtCab'*ω

	return f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2, f_P, f_H
end

# for initial time step of time stepping simulation
function element_equations(ΔL, Ct, Cab, CtCab, u, θ, F, M, v, ω, γ, κ, P, H, V, Ω,
	udot, θdot, CtCabPdot, CtCabHdot)

    f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2, f_P, f_H = unsteady_element_equations(
		ΔL, Ct, Cab, CtCab, u, θ, F, M, v, ω, γ, κ, P, H, V, Ω)

    f_u1 += ΔL/2*CtCabPdot
    f_u2 += ΔL/2*CtCabPdot
    f_ψ1 += ΔL/2*CtCabHdot
    f_ψ2 += ΔL/2*CtCabHdot

    f_P -= udot
    f_H -= Cab'*(((4 - θ'*θ/4 - 2*θtilde + θ*θ'/2)*θdot)/(4-(2-θ'*θ))^2)

	return f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2, f_P, f_H
end

# for time-marching simulations
function element_equations(ΔL, Ct, Cab, CtCab, u, θ, F, M, P, H,
	v, ω, γ, κ, V, Ω, udot, θdot, CtCabPdot, CtCabHdot, dt)

    f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2, f_P, f_H = element_equations(
		ΔL, Ct, Cab, CtCab, u, θ, F, M, v, ω, γ, κ, P, H, V, Ω, udot, θdot, CtCabPdot, CtCabHdot)

    f_u1 += 2/dt*ΔL/2*CtCab*P
    f_u2 += 2/dt*ΔL/2*CtCab*P
    m_ψ1 += 2/dt*ΔL/2*CtCab*H
    m_ψ2 += 2/dt*ΔL/2*CtCab*H

	f_P -= 2/dt*u

	return f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2, f_P, f_H
end
