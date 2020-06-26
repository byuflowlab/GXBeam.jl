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

"""
	element_properties(x, icol, ibeam, beam)
	element_properties(x, icol, ibeam, beam, x0, v0, ω0)
	element_properties(x, icol, ibeam, beam, x0, v0, ω0, udot, θdot, CtCabPdot, CtCabHdot)
	element_properties(x, icol, ibeam, beam, x0, v0, ω0, udot_p, θdot_p, CtCabPdot_p, CtCabHdot_p, u_p, θ_p, P_p, H_p, dt)

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
 - x: current state vector
 - icol: starting index for the beam's state variables
 - ibeam: beam element number
 - beam: beam element

# Additional Arguments for Dynamic Analyses
 - x0: Global frame origin
 - v0: Global frame linear velocity
 - ω0: Global frame angular velocity

# Additional Arguments for Initial Step Analyses
 - udot: time derivative of u for each beam element (used only for initial step of time-marching simulations)
 - θdot: time derivative of θ for each beam element (used only for initial step of time-marching simulations)
 - CtCabPdot: C'*Cab*Pdot for each beam element (used only for initial step of time-marching simulations)
 - CtCabHdot: C'*Cab*Hdot for each beam element (used only for initial step of time-marching simulations)

# Additional Arguments for Time Marching Analyses
 - u_p: `u` for each beam element from the previous time step (used for time-marching)
 - θ_p: `θ` for each beam element from the previous time step (used for time-marching)
 - CtCabP_p: `C'*Cab*P` for each beam element from the previous time step (used for time-marching)
 - CtCabH_p: `C'*Cab*H` for each beam element from the previous time step (used for time-marching)
 - udot_p: `udot` for each beam element from the previous time step (used for time-marching)
 - θdot_p: `θdot` for each beam element from the previous time step (used for time-marching)
 - CtCabPdot_p: `C'*Cab*Pdot` for each beam element from the previous time step (used for time-marching)
 - CtCabHdot_p: `C'*Cab*Hdot` for each beam element from the previous time step (used for time-marching)
 - dt: time step size (used for time-marching)
"""
element_properties

# static
function element_properties(x, icol, ibeam, beam)

	u = SVector(x[icol   ], x[icol+1 ], x[icol+2 ])
	θ = SVector(x[icol+3 ], x[icol+4 ], x[icol+5 ])
	F = SVector(x[icol+6 ], x[icol+7 ], x[icol+8 ])
	M = SVector(x[icol+9 ], x[icol+10], x[icol+11])

	ΔL = beam.ΔL
	Ct = wiener_milenkovic(θ)
	Cab = beam.Cab
	CtCab = Ct*Cab
	γ = strain(beam, F, M)
	κ = curvature(beam, F, M)

	return ΔL, Ct, Cab, CtCab, u, θ, F, M, γ, κ
end

# dynamic
function element_properties(x, icol, ibeam, beam, x0, v0, ω0)

	ΔL, Ct, Cab, CtCab, u, θ, F, M, γ, κ = element_properties(x, icol, ibeam, beam)

	v = v0 + cross(ω0, element.x - x0)
	ω = ω0

	P = SVector(x[icol+12], x[icol+13], x[icol+14])
	H = SVector(x[icol+15], x[icol+16], x[icol+17])

	V = linear_velocity(beam, P)
	Ω = angular_velocity(beam, H)

	return ΔL, Ct, Cab, CtCab, u, θ, F, M, γ, κ, v, ω, P, H, V, Ω
end

# initial step
function element_properties(x, icol, ibeam, beam, x0, v0, ω0, udot, θdot, CtCabPdot, CtCabHdot)

	ΔL, Ct, Cab, CtCab, u, θ, F, M, γ, κ, v, ω, P, H, V, Ω = element_properties(x, icol, ibeam, beam, x0, v0, ω0)

	return ΔL, Ct, Cab, CtCab, u, θ, F, M, γ, κ, v, ω, P, H, V, Ω,
		udot[ibeam], θdot[ibeam], CtCabPdot[ibeam], CtCabHdot[ibeam]
end

# time-marching
function element_properties(x, icol, ibeam, beam, x0, v0, ω0, udot_p, θdot_p,
	CtCabPdot_p, CtCabHdot_p, u_p, θ_p, P_p, H_p, dt)

	ΔL, Ct, Cab, CtCab, u, θ, F, M, P, H, v, ω, γ, κ, V, Ω = element_properties(x, icol, ibeam, beam, x0, v0, ω0)

	udot = -2/dt*u_p[ibeam] - udot_p[ibeam]
	θdot = -2/dt*(θ - θ_p[ibeam]) + θdot_p[ibeam]
	CtCabPdot = -2/dt*CtCabP_p[ibeam] - CtCabPdot_p[ibeam]
	CtCabHdot = -2/dt*CtCabH_p[ibeam] - CtCabHdot_p[ibeam]

	return ΔL, Ct, Cab, CtCab, u, θ, F, M, γ, κ, v, ω, P, H, V, Ω,
		udot, θdot, CtCabPdot, CtCabHdot, dt
end

"""
	element_equations(distributed_loads, ΔL, Ct, Cab, CtCab, u, θ, F, M, γ, κ)
	element_equations(distributed_loads, ΔL, Ct, Cab, CtCab, u, θ, F, M, γ, κ, v, ω, P, H, V, Ω)
    element_equations(distributed_loads, ΔL, Ct, Cab, CtCab, u, θ, F, M, γ, κ, v, ω, P, H, V, Ω,
		udot, θdot, CtCabPdot, CtCabHdot)

Evaluate the nonlinear equations for a beam element given the distributed loads
	on the beam element and its properties.

There are four implementations corresponding to the following analysis types:
 - Static
 - Dynamic - Steady State
 - Dynamic - Initial Step (for initializing time domain simulations)
 - Dynamic - Time Marching

See "GEBT: A general-purpose nonlinear analysis tool for composite beams" by
Wenbin Yu and "Geometrically nonlinear analysis of composite beams using
Wiener-Milenković parameters" by Qi Wang and Wenbin Yu.

# Arguments:
 - distributed_load: Distributed load on the beam element
 - ΔL: Length of the beam element
 - Ct: Rotation tensor of the beam deformation in the "a" frame, transposed
 - Cab: Direction cosine matrix from "a" to "b" frame for the element
 - CtCab: Ct*Cab, precomputed for efficiency
 - u: Displacement variables for the element [u1, u2, u3]
 - θ: Rotation variables for the element [θ1, θ2, θ3]
 - F: Force variables for the element [F1, F2, F3]
 - M: Moment variables for the element [M1, M2, M3]
 - γ: Engineering strains in the element [γ11, 2γ12, 2γ13]
 - κ: Curvatures in the element [κ1, κ2, κ3]

# Additional Arguments for Dynamic Analyses
 - v: Linear velocity of element in global frame "a" [v1, v2, v3]
 - ω: Angular velocity of element in global frame "a" [ω1, ω2, ω3]
 - P: Linear momenta for the element [P1, P2, P3] (unsteady simulations only)
 - H: Angular momenta for the element [H1, H2, H3] (unsteady simulations only)
 - V: Velocity of the element (unsteady simulations only)
 - Ω: Rotational velocity of the element (unsteady simulations only)

# Additional Arguments for Initial Step Analyses
 - udot: time derivative of u, evaluated at the previous iteration (time-marching simulations only)
 - θdot: time derivative of θ, evaluated at the previous iteration (time-marching simulations only)
 - CtCabPdot: C'*Cab*Pdot evaluated at the previous iteration (time-marching simulations only)
 - CtCabHdot: C'*Cab*Hdot evaluated at the previous iteration (time-marching simulations only)
"""
element_equations

# static
function element_equations(distributed_load, ΔL, Ct, Cab, CtCab, u, θ, F, M, γ, κ)

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
        f1, f2, m1, m2 = element_loads(ΔL, Ct, distributed_loads)
        f_u1 -= f1
        f_u2 -= f2
        f_ψ1 -= m1
        f_ψ2 -= m2
    end

	return f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2
end

# dynamic
function element_equations(distributed_loads, ΔL, Ct, Cab, CtCab, u, θ, F, M,
	v, ω, γ, κ, P, H, V, Ω)

    f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2 = element_equations(
		distributed_load, ΔL, Ct, Cab, CtCab, u, θ, F, M, γ, κ)

    f_u1 += cross(ω, ΔL/2*CtCab*P)
    f_u2 += cross(ω, ΔL/2*CtCab*P)

    f_ψ1 += cross(ω, ΔL/2*CtCab*H) + ΔL/2*CtCab*cross(V, P)
    f_ψ2 += cross(ω, ΔL/2*CtCab*H) + ΔL/2*CtCab*cross(V, P)

    f_P = CtCab*V - v - cross(ω, u)
	f_H = Ω - CtCab'*ω

	return f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2, f_P, f_H
end

# initial step
function element_equations(distributed_loads, ΔL, Ct, Cab, CtCab, u, θ, F, M,
	v, ω, γ, κ, P, H, V, Ω, udot, θdot, CtCabPdot, CtCabHdot)

    f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2, f_P, f_H = element_equations(
		distributed_loads, ΔL, Ct, Cab, CtCab, u, θ, F, M, v, ω, γ, κ, P, H, V, Ω)

    f_u1 += ΔL/2*CtCabPdot
    f_u2 += ΔL/2*CtCabPdot
    f_ψ1 += ΔL/2*CtCabHdot
    f_ψ2 += ΔL/2*CtCabHdot

    f_P -= udot
    f_H -= Cab'*(((4 - θ'*θ/4 - 2*θtilde + θ*θ'/2)*θdot)/(4-(2-θ'*θ))^2)

	return f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2, f_P, f_H
end

# time-marching
function element_equations(distributed_loads, ΔL, Ct, Cab, CtCab, u, θ, F, M,
	v, ω, γ, κ, P, H, V, Ω, udot, θdot, CtCabPdot, CtCabHdot, dt)

    f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2, f_P, f_H = element_equations(
		distributed_loads, ΔL, Ct, Cab, CtCab, u, θ, F, M, v, ω, γ, κ, P, H, V, Ω,
		udot, θdot, CtCabPdot, CtCabHdot)

    f_u1 += 2/dt*ΔL/2*CtCab*P
    f_u2 += 2/dt*ΔL/2*CtCab*P
    m_ψ1 += 2/dt*ΔL/2*CtCab*H
    m_ψ2 += 2/dt*ΔL/2*CtCab*H

	f_P -= 2/dt*u

	return f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2, f_P, f_H
end
