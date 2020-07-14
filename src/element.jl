"""
    Element{TF}

Composite type that defines a beam element's properties

# Fields
 - `ΔL`: Length of the beam element
 - `x`: Location of the beam element (the midpoint of the beam element)
 - `C11`: Upper left portion of the beam element's compliance matrix
 - `C12`: Upper right portion of the beam element's compliance matrix
 - `C22`: Lower right portion of the beam element's compliance matrix
 - `minv11`: Upper left portion of the inverse of the beam element's mass matrix
 - `minv12`: Upper right portion of the inverse of the beam element's mass matrix
 - `minv22`: Lower right portion of the inverse of the beam element's mass matrix
 - `Cab`: Rotation matrix from the global frame to beam element's frame
"""
struct Element{TF}
	ΔL::TF
	x::SVector{3,TF}
	C11::SMatrix{3,3,TF,9}
	C12::SMatrix{3,3,TF,9}
	C22::SMatrix{3,3,TF,9}
	minv11::SMatrix{3,3,TF,9}
	minv12::SMatrix{3,3,TF,9}
	minv22::SMatrix{3,3,TF,9}
	Cab::SMatrix{3,3,TF,9}
end

"""
    Element(ΔL, x, C, m, Cab)
Construct a beam element.

# Arguments
- `ΔL`: Length of the beam element
- `x`: Location of the beam element (the midpoint of the beam element)
- `C`: 6 x 6 compliance matrix of the beam element
- `m`: 6 x 6 mass matrix of the beam element
- `Cab`: Rotation matrix from the global frame to beam element's frame
"""
function Element(ΔL, x, C, m, Cab)
	TF = promote_type(typeof(ΔL), eltype(x), eltype(C), eltype(m), eltype(Cab))
	return Element{TF}(ΔL, x, C, m, Cab)
end

function Element{TF}(ΔL, x, C, m, Cab) where TF

	minv = m^-1

	return Element{TF}(ΔL, SVector{3,TF}(x),
	SMatrix{3,3,TF}(C[1:3,1:3]), SMatrix{3,3,TF}(C[1:3,4:6]),
	SMatrix{3,3,TF}(C[4:6,4:6]), SMatrix{3,3,TF}(minv[1:3,1:3]),
	SMatrix{3,3,TF}(minv[1:3,4:6]), SMatrix{3,3,TF}(minv[4:6,4:6]),
	SMatrix{3,3,TF}(Cab))
end

"""
    element_strain(element, F, M)

Calculate the strain of a beam element given the resultant force and moments
"""
@inline element_strain(element, F, M) = element.C11*F + element.C12*M

"""
    element_curvature(element, F, M)

Calculate the curvature of a beam element given the resultant force and moments
"""
@inline element_curvature(element, F, M) = element.C12'*F + element.C22*M

"""
    element_linear_velocity(element, P, H)

Calculate the linear velocity of a beam element given the linear and angular momenta
"""
@inline element_linear_velocity(element, P, H) = element.minv11*P + element.minv12*H

"""
    element_angular_velocity(element, P, H)

Calculate the angular velocity of a beam element given the linear and angular momenta
"""
@inline element_angular_velocity(element, P, H) = element.minv12'*P + element.minv22*H

"""
	element_properties(x, icol, ibeam, beam)
	element_properties(x, icol, ibeam, beam, x0, v0, ω0)
	element_properties(x, icol, ibeam, beam, x0, v0, ω0, u, θ, udot, θdot)
	element_properties(x, icol, ibeam, beam, x0, v0, ω0, udot, θdot_init, CtCabPdot, CtCabHdot, dt)

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
 - `ibeam`: beam element index
 - `beam`: beam element

# Additional Arguments for Dynamic Analyses
 - `x0`: Global frame origin (for the current time step)
 - `v0`: Global frame linear velocity (for the current time step)
 - `ω0`: Global frame angular velocity (for the current time step)

# Additional Arguments for Initial Step Analyses
 - `u`: deflection variables for each beam element
 - `θ`: rotation variables for each beam element
 - `udot`: time derivative of u for each beam element
 - `θdot`: time derivative of θ for each beam element

# Additional Arguments for Time Marching Analyses
 - `udot`: `-2/dt*u - udot` for each beam element from the previous time step
 - `θdot_init`: `-2/dt*θ - θdot` for each beam element from the previous time step
 - `CtCabPdot`: `-2/dt*C'*Cab*P - C'*Cab*Pdot` for each beam element from the previous time step
 - `CtCabHdot`: `-2/dt*C'*Cab*H - C'*Cab*Hdot` for each beam element from the previous time step
 - `dt`: time step size
"""
element_properties

# static
@inline function element_properties(x, icol, ibeam, beam)

	u = SVector(x[icol   ], x[icol+1 ], x[icol+2 ])
	θ = SVector(x[icol+3 ], x[icol+4 ], x[icol+5 ])
	F = SVector(x[icol+6 ], x[icol+7 ], x[icol+8 ])
	M = SVector(x[icol+9 ], x[icol+10], x[icol+11])

	ΔL = beam.ΔL
	Ct = wiener_milenkovic(θ)'
	Cab = beam.Cab
	CtCab = Ct*Cab
	γ = element_strain(beam, F, M)
	κ = element_curvature(beam, F, M)

	return ΔL, Ct, Cab, CtCab, u, θ, F, M, γ, κ
end

# dynamic
@inline function element_properties(x, icol, ibeam, beam, x0, v0, ω0)

	ΔL, Ct, Cab, CtCab, u, θ, F, M, γ, κ = element_properties(x, icol, ibeam, beam)

	v = v0 + cross(ω0, beam.x - x0)
	ω = ω0

	P = SVector(x[icol+12], x[icol+13], x[icol+14])
	H = SVector(x[icol+15], x[icol+16], x[icol+17])

	V = element_linear_velocity(beam, P, H)
	Ω = element_angular_velocity(beam, P, H)

	return ΔL, Ct, Cab, CtCab, u, θ, F, M, γ, κ, v, ω, P, H, V, Ω
end

# initial step
@inline function element_properties(x, icol, ibeam, beam, x0, v0, ω0, u, θ, udot, θdot)

	# note that CtCabPdot and CtCabHdot are state variables instead of u and θ

	# static variables
	CtCabPdot = SVector(x[icol   ], x[icol+1 ], x[icol+2 ])
	CtCabHdot = SVector(x[icol+3 ], x[icol+4 ], x[icol+5 ])
	F = SVector(x[icol+6 ], x[icol+7 ], x[icol+8 ])
	M = SVector(x[icol+9 ], x[icol+10], x[icol+11])

	ΔL = beam.ΔL
	Ct = wiener_milenkovic(θ[ibeam])'
	Cab = beam.Cab
	CtCab = Ct*Cab
	γ = element_strain(beam, F, M)
	κ = element_curvature(beam, F, M)

	# dynamic variables
	v = v0 + cross(ω0, beam.x - x0)
	ω = ω0

	P = SVector(x[icol+12], x[icol+13], x[icol+14])
	H = SVector(x[icol+15], x[icol+16], x[icol+17])

	V = element_linear_velocity(beam, P, H)
	Ω = element_angular_velocity(beam, P, H)

	return ΔL, Ct, Cab, CtCab, u[ibeam], θ[ibeam], F, M, γ, κ, v, ω, P, H, V, Ω,
		udot[ibeam], θdot[ibeam], CtCabPdot, CtCabHdot
end

# time-marching
@inline function element_properties(x, icol, ibeam, beam, x0, v0, ω0, udot_init, θdot_init,
	CtCabPdot_init, CtCabHdot_init, dt)

	ΔL, Ct, Cab, CtCab, u, θ, F, M, γ, κ, v, ω, P, H, V, Ω = element_properties(x, icol, ibeam, beam, x0, v0, ω0)

	udot = 2/dt*u - udot_init[ibeam]
	θdot = 2/dt*θ - θdot_init[ibeam]
	CtCabPdot = 2/dt*CtCab*P - CtCabPdot_init[ibeam]
	CtCabHdot = 2/dt*CtCab*H - CtCabHdot_init[ibeam]

	return ΔL, Ct, Cab, CtCab, u, θ, F, M, γ, κ, v, ω, P, H, V, Ω,
		udot, θdot, CtCabPdot, CtCabHdot, dt
end

"""
	element_equations(ΔL, Ct, Cab, CtCab, u, θ, F, M, γ, κ, f1, f2, m1, m2)
	element_equations(ΔL, Ct, Cab, CtCab, u, θ, F, M, γ, κ, v, ω, P, H, V, Ω,
		f1, f2, m1, m2)
    element_equations(ΔL, Ct, Cab, CtCab, u, θ, F, M, γ, κ, v, ω, P, H, V, Ω,
		udot, θdot, CtCabPdot, CtCabHdot, f1, f2, m1, m2)

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
 - `CtCab`: C'*Cab, precomputed for efficiency
 - `u`: Displacement variables for the element [u1, u2, u3]
 - `θ`: Rotation variables for the element [θ1, θ2, θ3]
 - `F`: Force variables for the element [F1, F2, F3]
 - `M`: Moment variables for the element [M1, M2, M3]
 - `γ`: Engineering strains in the element [γ11, 2γ12, 2γ13]
 - `κ`: Curvatures in the element [κ1, κ2, κ3]
 - `f1, f2`: Integrated distributed forces on the left and right side of the element
 - `m1, m2`: Integrated distributed moments on the left and right side of the element

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
 - `udot`: `2/dt*(u-u_p) - udot_p` for this beam element where `_p` denotes values
 	taken from the previous time step
 - `θdot`: `2/dt*(θ-θ_p) - θdot_p` for this beam element where `_p` denotes values
 	taken from the previous time step
 - `CtCabPdot`: `2/dt*(C'*Cab*P - (C'*Cab*P)_p) - (C'*Cab*Pdot)_p` for this beam
 	element where `_p` denotes values taken from the previous time step
 - `CtCabHdot`: `2/dt*(C'*Cab*H - (C'*Cab*H)_p) - (C'*Cab*Hdot)_p` for this beam
 	element where `_p` denotes values taken from the previous time step
"""
element_equations

# static
@inline function element_equations(ΔL, Ct, Cab, CtCab, u, θ, F, M, γ, κ, f1, f2, m1, m2)

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

# dynamic
@inline function element_equations(ΔL, Ct, Cab, CtCab, u, θ, F, M, γ, κ, v, ω, P, H, V, Ω,
	f1, f2, m1, m2)

    f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2 = element_equations(
		ΔL, Ct, Cab, CtCab, u, θ, F, M, γ, κ, f1, f2, m1, m2)

	tmp = cross(ω, ΔL/2*CtCab*P)
    f_u1 += tmp
    f_u2 += tmp

	tmp = cross(ω, ΔL/2*CtCab*H) + ΔL/2*CtCab*cross(V, P)
    f_ψ1 += tmp
    f_ψ2 += tmp

    f_P = CtCab*V - v - cross(ω, u)
	f_H = Ω - CtCab'*ω

	return f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2, f_P, f_H
end

# initial step
@inline function element_equations(ΔL, Ct, Cab, CtCab, u, θ, F, M, γ, κ, v, ω, P, H, V, Ω,
	udot, θdot, CtCabPdot, CtCabHdot, f1, f2, m1, m2)

    f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2, f_P, f_H = element_equations(
		ΔL, Ct, Cab, CtCab, u, θ, F, M, γ, κ, v, ω, P, H, V, Ω, f1, f2, m1, m2)

	tmp = ΔL/2*CtCabPdot
    f_u1 += tmp
    f_u2 += tmp

	tmp = ΔL/2*CtCabHdot
    f_ψ1 += tmp
    f_ψ2 += tmp

    f_P -= udot

	Q = get_Q(θ)
	f_H -= Cab'*Q*θdot

	return f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2, f_P, f_H
end

# time-marching
@inline element_equations(ΔL, Ct, Cab, CtCab, u, θ, F, M, γ, κ, v, ω, P, H, V, Ω,
	udot, θdot, CtCabPdot, CtCabHdot, dt, f1, f2, m1, m2) = element_equations(
	ΔL, Ct, Cab, CtCab, u, θ, F, M, γ, κ, v, ω, P, H, V, Ω, udot, θdot, CtCabPdot,
	CtCabHdot, f1, f2, m1, m2)

"""
	element_residual!(resid, irow_b, irow_b1, irow_p1, irow_b2, irow_p2, f_u1, f_u2, f_ψ1,
		f_ψ2, f_F1, f_F2, f_M1, f_M2)
	element_residual!(resid, irow_b, irow_b1, irow_p1, irow_b2, irow_p2, f_u1, f_u2, f_ψ1,
		f_ψ2, f_F1, f_F2, f_M1, f_M2, f_P, f_H)

Add the beam element's contributions to the residual equations.  Initialize
equilibrium and constitutive equations if they are not yet initialized.

If `irow_b1 != irow_p1` and/or `irow_b2 != irow_p2`, assume the equilibrium equations
for the left and/or right side are already initialized

There are two implementations, one for static and one for dynamic analyses

# Arguments
 - `resid`: System residual vector
 - `irow_b1`: Row index of the first equation for the left side of the beam element
    (a value <= 0 indicates the equations have been eliminated from the system of equations)
 - `irow_p1`: Row index of the first equation for the point on the left side of the beam element
 - `irow_b1`: Row index of the first equation for the right side of the beam element
 	(a value <= 0 indicates the equations have been eliminated from the system of equations)
 - `irow_p2`: Row index of the first equation for the point on the right side of the beam element
 - `f_u1`, `f_u2`: Resultant displacements for the left and right side of the beam element, respectively
 - `f_ψ1`, `f_ψ2`: Resultant rotations for the left and right side of the beam element, respectively
 - `f_F1`, `f_F2`: Resultant forces for the left and right side of the beam element, respectively
 - `f_M1`, `f_M2`: Resultant moments for the left and right side of the beam element, respectively

# Additional Arguments for Dynamic Analyses
 - `f_P`: Resultant linear momenta of the beam element
 - `f_H`: Resultant angular momenta of the beam element
"""
element_residual

# static
@inline function element_residual!(resid, irow_b, irow_b1, irow_p1, irow_b2, irow_p2,
	f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2)

	# create/add to residual equations for left endpoint
	if irow_b1 == irow_p1
		# add equilibrium equations
		resid[irow_p1:irow_p1+2] .= f_u1./FORCE_SCALING
		resid[irow_p1+3:irow_p1+5] .= f_ψ1./FORCE_SCALING
		# add compatability equations
		resid[irow_p1+6:irow_p1+8] .= f_F1
		resid[irow_p1+9:irow_p1+11] .= f_M1
	else
		# add to existing equilibrium equations
		resid[irow_p1] += f_u1[1]/FORCE_SCALING
		resid[irow_p1+1] += f_u1[2]/FORCE_SCALING
		resid[irow_p1+2] += f_u1[3]/FORCE_SCALING
		resid[irow_p1+3] += f_ψ1[1]/FORCE_SCALING
		resid[irow_p1+4] += f_ψ1[2]/FORCE_SCALING
		resid[irow_p1+5] += f_ψ1[3]/FORCE_SCALING
		if irow_b1 <= 0
			# compatability equations have been combined with other beam
			resid[irow_p1+6] += f_F1[1]
			resid[irow_p1+7] += f_F1[2]
			resid[irow_p1+8] += f_F1[3]
			resid[irow_p1+9] += f_M1[1]
			resid[irow_p1+10] += f_M1[2]
			resid[irow_p1+11] += f_M1[3]
		else
			# create compatability equations for this beam
			resid[irow_b1:irow_b1+2] .= f_F1
			resid[irow_b1+3:irow_b1+5] .= f_M1
		end
	end

	# create/add to residual equations for right endpoint
	if irow_b2 == irow_p2
		# add equilibrium equations
		resid[irow_p2:irow_p2+2] .= f_u2./FORCE_SCALING
		resid[irow_p2+3:irow_p2+5] .= f_ψ2./FORCE_SCALING
		# add compatability equations
		resid[irow_p2+6:irow_p2+8] .= f_F2
		resid[irow_p2+9:irow_p2+11] .= f_M2
	else
		# add to existing equilibrium equations
		resid[irow_p2] += f_u2[1]/FORCE_SCALING
		resid[irow_p2+1] += f_u2[2]/FORCE_SCALING
		resid[irow_p2+2] += f_u2[3]/FORCE_SCALING
		resid[irow_p2+3] += f_ψ2[1]/FORCE_SCALING
		resid[irow_p2+4] += f_ψ2[2]/FORCE_SCALING
		resid[irow_p2+5] += f_ψ2[3]/FORCE_SCALING
		if irow_b2 <= 0
			# compatability equations have been combined with other beam
			resid[irow_p2+6] += f_F2[1]
			resid[irow_p2+7] += f_F2[2]
			resid[irow_p2+8] += f_F2[3]
			resid[irow_p2+9] += f_M2[1]
			resid[irow_p2+10] += f_M2[2]
			resid[irow_p2+11] += f_M2[3]
		else
			# create compatability equations for this beam
			resid[irow_b2:irow_b2+2] .= f_F2
			resid[irow_b2+3:irow_b2+5] .= f_M2
		end
	end

	return resid
end

# dynamic
@inline function element_residual!(resid, irow_b, irow_b1, irow_p1, irow_b2, irow_p2,
	f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2, f_P, f_H)

	resid = element_residual!(resid, irow_b, irow_b1, irow_p1, irow_b2, irow_p2, f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2)

	# residual equations for element
	resid[irow_b:irow_b+2] .= f_P
	resid[irow_b+3:irow_b+5] .= f_H

	return resid
end

"""
	element_jacobian_equations(beam, ΔL, Ct, Cab, CtCab, u, θ, F, M,
		γ, κ, Ct_θ1, Ct_θ2, Ct_θ3, f1_θ, f2_θ, m1_θ, m2_θ)
	element_jacobian_equations(beam, ΔL, Ct, Cab, CtCab, u, θ, F, M,
		γ, κ, v, ω, P, H, V, Ω, Ct_θ1, Ct_θ2, Ct_θ3, f1_θ, f2_θ, m1_θ, m2_θ)
    element_jacobian_equations(beam, ΔL, Ct, Cab, CtCab, u, θ, F, M,
		γ, κ, v, ω, P, H, V, Ω, udot, θdot, CtCabPdot, CtCabHdot, Ct_θ1, Ct_θ2,
		Ct_θ3, f1_θ, f2_θ, m1_θ, m2_θ)
	element_jacobian_equations(beam, ΔL, Ct, Cab, CtCab, u, θ, F, M,
		γ, κ, v, ω, P, H, V, Ω, udot, θdot, CtCabPdot, CtCabHdot, dt, Ct_θ1,
		Ct_θ2, Ct_θ3, f1_θ, f2_θ, m1_θ, m2_θ)

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
 - `beam`: Beam element
 - `ΔL`: Length of the beam element
 - `Ct`: Rotation tensor of the beam deformation in the "a" frame, transposed
 - `Cab`: Direction cosine matrix from "a" to "b" frame for the element
 - `CtCab`: Ct*Cab, precomputed for efficiency
 - `u`: Displacement variables for the element [u1, u2, u3]
 - `θ`: Rotation variables for the element [θ1, θ2, θ3]
 - `F`: Force variables for the element [F1, F2, F3]
 - `M`: Moment variables for the element [M1, M2, M3]
 - `γ`: Engineering strains in the element [γ11, 2γ12, 2γ13]
 - `κ`: Curvatures in the element [κ1, κ2, κ3]
 - `Ct_θ1`: Gradient of `Ct` w.r.t. `θ[1]`
 - `Ct_θ2`: Gradient of `Ct` w.r.t. `θ[2]`
 - `Ct_θ3`: Gradient of `Ct` w.r.t. `θ[3]`
 - `f1_θ`: Gradient w.r.t. θ of integrated distributed force for left side of element
 - `m1_θ`: Gradient w.r.t. θ of integrated distributed moment for left side of element
 - `f2_θ`: Gradient w.r.t. θ of integrated distributed force for right side of element
 - `m2_θ`: Gradient w.r.t. θ of integrated distributed moment for right side of element

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
 - `udot`: `-2/dt*u - udot` for each beam element from the previous time step
 - `θdot_init`: `-2/dt*θ - θdot` for each beam element from the previous time step
 - `CtCabPdot`: `-2/dt*C'*Cab*P - C'*Cab*Pdot` for each beam element from the previous time step
 - `CtCabHdot`: `-2/dt*C'*Cab*H - C'*Cab*Hdot` for each beam element from the previous time step
 - `dt`: time step size
"""
element_jacobian_equations

# static
@inline function element_jacobian_equations(beam, ΔL, Ct, Cab, CtCab, u, θ, F, M,
	γ, κ, Ct_θ1, Ct_θ2, Ct_θ3, f1_θ, f2_θ, m1_θ, m2_θ)

	# --- f_u1, f_u2 --- #

	# d_fu/d_θ
	tmp = matrix_jacobian_product(Ct_θ1, Ct_θ2, Ct_θ3, Cab*F)
	f_u1_θ = -tmp - f1_θ
	f_u2_θ =  tmp - f2_θ

	# d_fu/d_F
	f_u1_F = -CtCab
	f_u2_F =  CtCab

	# --- f_ψ1, f_ψ2 --- #

	# d_fψ/d_θ
	tmp1 = matrix_jacobian_product(Ct_θ1, Ct_θ2, Ct_θ3, Cab*M)
	tmp2 = ΔL/2*matrix_jacobian_product(Ct_θ1, Ct_θ2, Ct_θ3, Cab*cross(e1 + γ, F))
	f_ψ1_θ = -tmp1 - m1_θ - tmp2
	f_ψ2_θ =  tmp1 - m2_θ - tmp2

	# d_fψ/d_F
	tmp = -ΔL/2*CtCab*(tilde(e1 + γ) - tilde(F)*beam.C11)
	f_ψ1_F = tmp
	f_ψ2_F = tmp

	# d_fψ/d_M
	tmp = ΔL/2*CtCab*tilde(F)*beam.C12
    f_ψ1_M = tmp - CtCab
    f_ψ2_M = tmp + CtCab

	# --- f_F1, f_F2 --- #

	# d_fF/d_u
    f_F1_u =  I3
    f_F2_u = -I3

	# d_fF/d_θ
	tmp = matrix_jacobian_product(Ct_θ1, Ct_θ2, Ct_θ3, ΔL/2*Cab*(e1 + γ))
	f_F1_θ = -tmp
	f_F2_θ = -tmp

	# d_fF/d_F
	tmp = ΔL/2*CtCab*beam.C11
	f_F1_F = -tmp
	f_F2_F = -tmp

	# d_fF/d_M
	tmp = ΔL/2*CtCab*beam.C12
	f_F1_M = -tmp
	f_F2_M = -tmp

	# --- f_M1, f_M2 --- #

	# d_fM/d_θ
	Qinv_θ1, Qinv_θ2, Qinv_θ3 = get_Qinv_jacobian(θ)
	tmp = matrix_jacobian_product(Qinv_θ1, Qinv_θ2, Qinv_θ3, ΔL/2*Cab*κ)
	f_M1_θ =  I - tmp
	f_M2_θ = -I - tmp

	# d_fM/d_F
	Qinv = get_Qinv(θ)
	tmp1 = -ΔL/2*Qinv*Cab
	tmp2 = tmp1*beam.C12'
	f_M1_F = tmp2
    f_M2_F = tmp2

	# d_fM/d_M
	tmp2 = tmp1*beam.C22
	f_M1_M = tmp2
    f_M2_M = tmp2

	return f_u1_θ, f_u2_θ, f_u1_F, f_u2_F,
		f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M,
		f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
		f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M
end

# dynamic
@inline function element_jacobian_equations(beam, ΔL, Ct, Cab, CtCab, u, θ, F, M,
	γ, κ, v, ω, P, H, V, Ω, Ct_θ1, Ct_θ2, Ct_θ3, f1_θ, f2_θ, m1_θ, m2_θ)

	f_u1_θ, f_u2_θ, f_u1_F, f_u2_F,
	f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M,
    f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
	f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M =
	    element_jacobian_equations(beam, ΔL, Ct, Cab, CtCab, u, θ, F, M,
		γ, κ, Ct_θ1, Ct_θ2, Ct_θ3, f1_θ, f2_θ, m1_θ, m2_θ)

	# --- f_u1, f_u2 --- #

	# d_fu_dθ
	tmp = tilde(ω)*matrix_jacobian_product(Ct_θ1, Ct_θ2, Ct_θ3, ΔL/2*Cab*P)
	f_u1_θ += tmp
	f_u2_θ += tmp

	# d_fu_dP
	tmp = ΔL/2*tilde(ω)*CtCab
	f_u1_P = tmp
	f_u2_P = tmp

	# --- f_ψ1, f_ψ2 --- #

	# d_fψ_dθ
	tmp1 = tilde(ω)*matrix_jacobian_product(Ct_θ1, Ct_θ2, Ct_θ3, ΔL/2*Cab*H)
	tmp2 = matrix_jacobian_product(Ct_θ1, Ct_θ2, Ct_θ3, ΔL/2*Cab*cross(V, P))
	tmp3 = tmp1 + tmp2
	f_ψ1_θ += tmp3
	f_ψ2_θ += tmp3

	# d_fψ_dP
	tmp = ΔL/2*CtCab*(tilde(V) - tilde(P)*beam.minv11)
	f_ψ1_P = tmp
	f_ψ2_P = tmp

	# d_fψ_dH
	tmp = tilde(ω)*ΔL/2*CtCab - ΔL/2*CtCab*(tilde(P)*beam.minv12)
	f_ψ1_H = tmp
	f_ψ2_H = tmp

	# --- f_P --- #

	# d_fP_du
	f_P_u = -tilde(ω)

	# d_fP_dθ
	f_P_θ = matrix_jacobian_product(Ct_θ1, Ct_θ2, Ct_θ3, Cab*V)

	# d_fP_dP
	f_P_P = CtCab*beam.minv11

	# d_fP_dH
	f_P_H = CtCab*beam.minv12

	# --- f_H --- #

	# d_fH_dθ
	f_H_θ = -Cab'*matrix_jacobian_product(Ct_θ1', Ct_θ2', Ct_θ3', ω)

	# d_fH_dP
	f_H_P = beam.minv12'

	# d_fH_dH
	f_H_H = beam.minv22

	return f_u1_θ, f_u2_θ, f_u1_F, f_u2_F, f_u1_P, f_u2_P,
		f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_P, f_ψ2_P, f_ψ1_H, f_ψ2_H,
    	f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
		f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M,
		f_P_u, f_P_θ, f_P_P, f_P_H,
		f_H_θ, f_H_P, f_H_H
end

# initial step
@inline function element_jacobian_equations(beam, ΔL, Ct, Cab, CtCab, u, θ, F, M,
	γ, κ, v, ω, P, H, V, Ω, udot, θdot, CtCabPdot, CtCabHdot,  Ct_θ1, Ct_θ2,
	Ct_θ3, f1_θ, f2_θ, m1_θ, m2_θ)

	# --- f_u1, f_u2 --- #

	# d_fu/d_CtCabPdot
	tmp = ΔL/2*I3
	f_u1_CtCabPdot = tmp
	f_u2_CtCabPdot = tmp

	# d_fu/d_F
	f_u1_F = -CtCab
	f_u2_F =  CtCab

	# d_fu/d_P
	tmp = ΔL/2*tilde(ω)*CtCab
	f_u1_P = tmp
	f_u2_P = tmp

	# --- f_θ1, f_θ2 --- #

	# d_fm/d_CtCabHdot
	tmp = ΔL/2*I3
	f_ψ1_CtCabHdot = tmp
	f_ψ2_CtCabHdot = tmp

	# d_fψ/d_F
	tmp = -ΔL/2*CtCab*(tilde(e1 + γ) - tilde(F)*beam.C11)
	f_ψ1_F = tmp
	f_ψ2_F = tmp

	# d_fψ/d_M
	tmp = ΔL/2*CtCab*tilde(F)*beam.C12
    f_ψ1_M = tmp - CtCab
    f_ψ2_M = tmp + CtCab

	# d_fψ_dP
	tmp = ΔL/2*CtCab*(tilde(V) - tilde(P)*beam.minv11)
	f_ψ1_P = tmp
	f_ψ2_P = tmp

	# d_fψ_dH
	tmp = tilde(ω)*ΔL/2*CtCab - ΔL/2*CtCab*(tilde(P)*beam.minv12)
	f_ψ1_H = tmp
	f_ψ2_H = tmp

	# --- f_F1, f_F2 --- #

	# d_fF/d_F
	tmp = ΔL/2*CtCab*beam.C11
	f_F1_F = -tmp
	f_F2_F = -tmp

	# d_fF/d_M
	tmp = ΔL/2*CtCab*beam.C12
	f_F1_M = -tmp
	f_F2_M = -tmp

	# --- f_M1, f_M2 --- #

	# d_fM/d_F
	Qinv = get_Qinv(θ)
	tmp1 = -ΔL/2*Qinv*Cab
	tmp2 = tmp1*beam.C12'
	f_M1_F = tmp2
    f_M2_F = tmp2

	# d_fM/d_M
	tmp2 = tmp1*beam.C22
	f_M1_M = tmp2
    f_M2_M = tmp2

	# --- f_P --- #

	# d_fP_dP
	f_P_P = CtCab*beam.minv11

	# d_fP_dH
	f_P_H = CtCab*beam.minv12

	# --- f_H --- #

	# d_fH_dP
	f_H_P = beam.minv12'

	# d_fH_dH
	f_H_H = beam.minv22

	return f_u1_CtCabPdot, f_u2_CtCabPdot, f_u1_F, f_u2_F, f_u1_P, f_u2_P,
		f_ψ1_CtCabHdot, f_ψ2_CtCabHdot, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_P, f_ψ2_P, f_ψ1_H, f_ψ2_H,
		f_F1_F, f_F2_F, f_F1_M, f_F2_M,
		f_M1_F, f_M2_F, f_M1_M, f_M2_M,
		f_P_P, f_P_H,
		f_H_P, f_H_H
end

# time-marching
@inline function element_jacobian_equations(beam, ΔL, Ct, Cab, CtCab, u, θ, F, M, γ, κ,
	v, ω, P, H, V, Ω, udot, θdot, CtCabPdot, CtCabHdot, dt, Ct_θ1, Ct_θ2, Ct_θ3,
	f1_θ, f2_θ, m1_θ, m2_θ)

	f_u1_θ, f_u2_θ, f_u1_F, f_u2_F, f_u1_P, f_u2_P,
	f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_P, f_ψ2_P, f_ψ1_H, f_ψ2_H,
	f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
	f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M,
	f_P_u, f_P_θ, f_P_P, f_P_H,
	f_H_θ, f_H_P, f_H_H =
		element_jacobian_equations(beam, ΔL, Ct, Cab, CtCab, u, θ, F, M,
		γ, κ, v, ω, P, H, V, Ω, Ct_θ1, Ct_θ2, Ct_θ3, f1_θ, f2_θ, m1_θ, m2_θ)

	# --- f_u1, f_u2 --- #

	# d_fu_dθ
	tmp = ΔL/dt*matrix_jacobian_product(Ct_θ1, Ct_θ2, Ct_θ3, Cab*P)
	f_u1_θ += tmp
	f_u2_θ += tmp

	# d_fu_dP
	tmp = ΔL/dt*CtCab
	f_u1_P += tmp
	f_u2_P += tmp

	# --- f_ψ1, f_ψ2 --- #

	# d_fψ_dθ
	tmp = ΔL/dt*matrix_jacobian_product(Ct_θ1, Ct_θ2, Ct_θ3, Cab*H)
	f_ψ1_θ += tmp
	f_ψ2_θ += tmp

	# d_fψ_dH
	tmp = ΔL/dt*CtCab
	f_ψ1_H += tmp
	f_ψ2_H += tmp

	# --- d_fP_du --- #
	f_P_u -= 2/dt*I

	# --- d_fH_dθ --- #
	Q = get_Q(θ)
	Q_θ1, Q_θ2, Q_θ3 = get_Q_jacobian(θ)
	f_H_θ -= Cab'*(matrix_jacobian_product(Q_θ1, Q_θ2, Q_θ3, θdot) + Q*2/dt)

	return f_u1_θ, f_u2_θ, f_u1_F, f_u2_F, f_u1_P, f_u2_P,
		f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_P, f_ψ2_P, f_ψ1_H, f_ψ2_H,
		f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
		f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M,
		f_P_u, f_P_θ, f_P_P, f_P_H,
		f_H_θ, f_H_P, f_H_H
end


"""
	element_jacobian!(jacob, irow_b1, irow_p1, irow_b2, irow_p2, icol,
		f_u1_θ, f_u2_θ, f_u1_F, f_u2_F,
		f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M,
		f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
		f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M)
	element_jacobian!(jacob, irow_b1, irow_p1, irow_b2, irow_p2, icol,
		f_u1_CtCabPdot, f_u2_CtCabPdot, f_u1_F, f_u2_F, f_u1_P, f_u2_P,
		f_ψ1_CtCabHdot, f_ψ2_CtCabHdot, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_P, f_ψ2_P, f_ψ1_H, f_ψ2_H,
		f_F1_F, f_F2_F, f_F1_M, f_F2_M,
		f_M1_F, f_M2_F, f_M1_M, f_M2_M,
		f_P_P, f_P_H,
		f_H_P, f_H_H)
	element_jacobian!(jacob, irow_b1, irow_p1, irow_b2, irow_p2, icol,
		f_u1_θ, f_u2_θ, f_u1_F, f_u2_F, f_u1_P, f_u2_P,
		f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_P, f_ψ2_P, f_ψ1_H, f_ψ2_H,
		f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
		f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M,
		f_P_u, f_P_θ, f_P_P, f_P_H,
		f_H_θ, f_H_P, f_H_H)

Add the beam element's contributions to the jacobian matrix

There are three implementations:
 - Static
 - Initial Step for Time-Marching Simulations
 - Dynamic

# Arguments
 - `jacob`: System jacobian matrix
 - `icol_p1`: Row/column index of the first unknown for the left endpoint
   (a value <= 0 indicates the unknowns have been eliminated from the system of equations)
 - `irow_b1`: Row index of the first equation for the left side of the beam
 - `irow_p1`: Row index of the first equation for the point on the left side of the beam
 - `icol_p2`: Row/column index of the first unknown for the right endpoint
 	(a value <= 0 indicates the unknowns have been eliminated from the system of equations)
 - `irow_b2`: Row index of the first equation for the right side of the beam
 - `irow_p2`: Row index of the first equation for the point on the right side of the beam
 - `icol_b`: Row/Column index corresponding to the first beam state variable

All other arguments use the following naming convention:
 - `f_y_x`: Jacobian of element equation "y" with respect to state variable "x"

"""
element_jacobian

# static
@inline function element_jacobian!(jacob, irow_b, irow_b1, irow_p1, irow_b2, irow_p2, icol_b,
	f_u1_θ, f_u2_θ, f_u1_F, f_u2_F,
	f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M,
	f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
	f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M)

	# create/add to residual jacobian entries for left endpoint
	if irow_b1 == irow_p1
		# initialize equilibrium and compatability equation jacobian entries
		jacob[irow_p1:irow_p1+2, icol_b+3:icol_b+5] .= f_u1_θ./FORCE_SCALING
		jacob[irow_p1:irow_p1+2, icol_b+6:icol_b+8] .= f_u1_F./FORCE_SCALING

		jacob[irow_p1+3:irow_p1+5, icol_b+3:icol_b+5] .= f_ψ1_θ./FORCE_SCALING
		jacob[irow_p1+3:irow_p1+5, icol_b+6:icol_b+8] .= f_ψ1_F./FORCE_SCALING
		jacob[irow_p1+3:irow_p1+5, icol_b+9:icol_b+11] .= f_ψ1_M./FORCE_SCALING

		jacob[irow_p1+6:irow_p1+8, icol_b:icol_b+2] .= f_F1_u
		jacob[irow_p1+6:irow_p1+8, icol_b+3:icol_b+5] .= f_F1_θ
		jacob[irow_p1+6:irow_p1+8, icol_b+6:icol_b+8] .= f_F1_F
		jacob[irow_p1+6:irow_p1+8, icol_b+9:icol_b+11] .= f_F1_M

		jacob[irow_p1+9:irow_p1+11, icol_b+3:icol_b+5] .= f_M1_θ
		jacob[irow_p1+9:irow_p1+11, icol_b+6:icol_b+8] .= f_M1_F
		jacob[irow_p1+9:irow_p1+11, icol_b+9:icol_b+11] .= f_M1_M
	else
		# add to existing equilibrium equation jacobian entries
		jacob[irow_p1:irow_p1+2, icol_b+3:icol_b+5] .= f_u1_θ./FORCE_SCALING
		jacob[irow_p1:irow_p1+2, icol_b+6:icol_b+8] .= f_u1_F./FORCE_SCALING

		jacob[irow_p1+3:irow_p1+5, icol_b+3:icol_b+5] .= f_ψ1_θ./FORCE_SCALING
		jacob[irow_p1+3:irow_p1+5, icol_b+6:icol_b+8] .= f_ψ1_F./FORCE_SCALING
		jacob[irow_p1+3:irow_p1+5, icol_b+9:icol_b+11] .= f_ψ1_M./FORCE_SCALING

		if irow_b1 <= 0
			# compatability equations have been combined with other beam
			jacob[irow_p1+6:irow_p1+8, icol_b:icol_b+2] .= f_F1_u
			jacob[irow_p1+6:irow_p1+8, icol_b+3:icol_b+5] .= f_F1_θ
			jacob[irow_p1+6:irow_p1+8, icol_b+6:icol_b+8] .= f_F1_F
			jacob[irow_p1+6:irow_p1+8, icol_b+9:icol_b+11] .= f_F1_M

			jacob[irow_p1+9:irow_p1+11, icol_b+3:icol_b+5] .= f_M1_θ
			jacob[irow_p1+9:irow_p1+11, icol_b+6:icol_b+8] .= f_M1_F
			jacob[irow_p1+9:irow_p1+11, icol_b+9:icol_b+11] .= f_M1_M
		else
			# initialize compatability equation jacobian entries for this beam
			jacob[irow_b1:irow_b1+2, icol_b:icol_b+2] .= f_F1_u
			jacob[irow_b1:irow_b1+2, icol_b+3:icol_b+5] .= f_F1_θ
			jacob[irow_b1:irow_b1+2, icol_b+6:icol_b+8] .= f_F1_F
			jacob[irow_b1:irow_b1+2, icol_b+9:icol_b+11] .= f_F1_M

			jacob[irow_b1+3:irow_b1+5, icol_b+3:icol_b+5] .= f_M1_θ
			jacob[irow_b1+3:irow_b1+5, icol_b+6:icol_b+8] .= f_M1_F
			jacob[irow_b1+3:irow_b1+5, icol_b+9:icol_b+11] .= f_M1_M
		end
	end

	# create/add to residual jacobian entries for right endpoint
	if irow_b2 == irow_p2
		# initialize equilibrium equation jacobian entries
		jacob[irow_p2:irow_p2+2, icol_b+3:icol_b+5] .= f_u2_θ./FORCE_SCALING
		jacob[irow_p2:irow_p2+2, icol_b+6:icol_b+8] .= f_u2_F./FORCE_SCALING

		jacob[irow_p2+3:irow_p2+5, icol_b+3:icol_b+5] .= f_ψ2_θ./FORCE_SCALING
		jacob[irow_p2+3:irow_p2+5, icol_b+6:icol_b+8] .= f_ψ2_F./FORCE_SCALING
		jacob[irow_p2+3:irow_p2+5, icol_b+9:icol_b+11] .= f_ψ2_M./FORCE_SCALING

		# initialize compatability equation jacobian entries
		jacob[irow_p2+6:irow_p2+8, icol_b:icol_b+2] .= f_F2_u
		jacob[irow_p2+6:irow_p2+8, icol_b+3:icol_b+5] .= f_F2_θ
		jacob[irow_p2+6:irow_p2+8, icol_b+6:icol_b+8] .= f_F2_F
		jacob[irow_p2+6:irow_p2+8, icol_b+9:icol_b+11] .= f_F2_M

		jacob[irow_p2+9:irow_p2+11, icol_b+3:icol_b+5] .= f_M2_θ
		jacob[irow_p2+9:irow_p2+11, icol_b+6:icol_b+8] .= f_M2_F
		jacob[irow_p2+9:irow_p2+11, icol_b+9:icol_b+11] .= f_M2_M
	else
		# add to existing equilibrium equation jacobian entries
		jacob[irow_p2:irow_p2+2, icol_b+3:icol_b+5] .= f_u2_θ./FORCE_SCALING
		jacob[irow_p2:irow_p2+2, icol_b+6:icol_b+8] .= f_u2_F./FORCE_SCALING

		jacob[irow_p2+3:irow_p2+5, icol_b+3:icol_b+5] .= f_ψ2_θ./FORCE_SCALING
		jacob[irow_p2+3:irow_p2+5, icol_b+6:icol_b+8] .= f_ψ2_F./FORCE_SCALING
		jacob[irow_p2+3:irow_p2+5, icol_b+9:icol_b+11] .= f_ψ2_M./FORCE_SCALING

		if irow_b2 <= 0
			jacob[irow_p2+6:irow_p2+8, icol_b:icol_b+2] .= f_F2_u
			jacob[irow_p2+6:irow_p2+8, icol_b+3:icol_b+5] .= f_F2_θ
			jacob[irow_p2+6:irow_p2+8, icol_b+6:icol_b+8] .= f_F2_F
			jacob[irow_p2+6:irow_p2+8, icol_b+9:icol_b+11] .= f_F2_M

			jacob[irow_p2+9:irow_p2+11, icol_b+3:icol_b+5] .= f_M2_θ
			jacob[irow_p2+9:irow_p2+11, icol_b+6:icol_b+8] .= f_M2_F
			jacob[irow_p2+9:irow_p2+11, icol_b+9:icol_b+11] .= f_M2_M
		else
			# initialize compatability equation jacobian entries
			jacob[irow_b2:irow_b2+2, icol_b:icol_b+2] .= f_F2_u
			jacob[irow_b2:irow_b2+2, icol_b+3:icol_b+5] .= f_F2_θ
			jacob[irow_b2:irow_b2+2, icol_b+6:icol_b+8] .= f_F2_F
			jacob[irow_b2:irow_b2+2, icol_b+9:icol_b+11] .= f_F2_M

			jacob[irow_b2+3:irow_b2+5, icol_b+3:icol_b+5] .= f_M2_θ
			jacob[irow_b2+3:irow_b2+5, icol_b+6:icol_b+8] .= f_M2_F
			jacob[irow_b2+3:irow_b2+5, icol_b+9:icol_b+11] .= f_M2_M
		end
	end

	return jacob
end

# initial step
@inline function element_jacobian!(jacob, irow_b, irow_b1, irow_p1, irow_b2, irow_p2, icol_b,
	f_u1_CtCabPdot, f_u2_CtCabPdot, f_u1_F, f_u2_F, f_u1_P, f_u2_P,
	f_ψ1_CtCabHdot, f_ψ2_CtCabHdot, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_P, f_ψ2_P, f_ψ1_H, f_ψ2_H,
	f_F1_F, f_F2_F, f_F1_M, f_F2_M,
	f_M1_F, f_M2_F, f_M1_M, f_M2_M,
	f_P_P, f_P_H,
	f_H_P, f_H_H)

	# create/add to residual jacobian entries for left endpoint
	if irow_b1 == irow_p1
		# initialize equilibrium equation jacobian entries
		jacob[irow_p1:irow_p1+2, icol_b:icol_b+2] .= f_u1_CtCabPdot./FORCE_SCALING
		jacob[irow_p1:irow_p1+2, icol_b+6:icol_b+8] .= f_u1_F./FORCE_SCALING
		jacob[irow_p1:irow_p1+2, icol_b+12:icol_b+14] .= f_u1_P./FORCE_SCALING

		jacob[irow_p1+3:irow_p1+5, icol_b+3:icol_b+5] .= f_ψ1_CtCabHdot./FORCE_SCALING
		jacob[irow_p1+3:irow_p1+5, icol_b+6:icol_b+8] .= f_ψ1_F./FORCE_SCALING
		jacob[irow_p1+3:irow_p1+5, icol_b+9:icol_b+11] .= f_ψ1_M./FORCE_SCALING
		jacob[irow_p1+3:irow_p1+5, icol_b+12:icol_b+14] .= f_ψ1_P./FORCE_SCALING
		jacob[irow_p1+3:irow_p1+5, icol_b+15:icol_b+17] .= f_ψ1_H./FORCE_SCALING

		# initialize compatability equation jacobian entries
		jacob[irow_p1+6:irow_p1+8, icol_b+6:icol_b+8] .= f_F1_F
		jacob[irow_p1+6:irow_p1+8, icol_b+9:icol_b+11] .= f_F1_M

		jacob[irow_p1+9:irow_p1+11, icol_b+6:icol_b+8] .= f_M1_F
		jacob[irow_p1+9:irow_p1+11, icol_b+9:icol_b+11] .= f_M1_M

	else
		# add to existing equilibrium equation jacobian entries
		jacob[irow_p1:irow_p1+2, icol_b:icol_b+2] .= f_u1_CtCabPdot./FORCE_SCALING
		jacob[irow_p1:irow_p1+2, icol_b+6:icol_b+8] .= f_u1_F./FORCE_SCALING
		jacob[irow_p1:irow_p1+2, icol_b+12:icol_b+14] .= f_u1_P./FORCE_SCALING

		jacob[irow_p1+3:irow_p1+5, icol_b+3:icol_b+5] .= f_ψ1_CtCabHdot./FORCE_SCALING
		jacob[irow_p1+3:irow_p1+5, icol_b+6:icol_b+8] .= f_ψ1_F./FORCE_SCALING
		jacob[irow_p1+3:irow_p1+5, icol_b+9:icol_b+11] .= f_ψ1_M./FORCE_SCALING
		jacob[irow_p1+3:irow_p1+5, icol_b+12:icol_b+14] .= f_ψ1_P./FORCE_SCALING
		jacob[irow_p1+3:irow_p1+5, icol_b+15:icol_b+17] .= f_ψ1_H./FORCE_SCALING

		if irow_b1 <= 0
			# compatability equations have been combined with other beam
			jacob[irow_p1+6:irow_p1+8, icol_b+6:icol_b+8] .= f_F1_F
			jacob[irow_p1+6:irow_p1+8, icol_b+9:icol_b+11] .= f_F1_M

			jacob[irow_p1+9:irow_p1+11, icol_b+6:icol_b+8] .= f_M1_F
			jacob[irow_p1+9:irow_p1+11, icol_b+9:icol_b+11] .= f_M1_M
		else
			# initialize compatability equation jacobian entries
			jacob[irow_b1:irow_b1+2, icol_b+6:icol_b+8] .= f_F1_F
			jacob[irow_b1:irow_b1+2, icol_b+9:icol_b+11] .= f_F1_M

			jacob[irow_b1+3:irow_b1+5, icol_b+6:icol_b+8] .= f_M1_F
			jacob[irow_b1+3:irow_b1+5, icol_b+9:icol_b+11] .= f_M1_M
		end

	end

	# create/add to residual jacobian entries for right endpoint
	if irow_b2 == irow_p2
		# initialize equilibrium equation jacobian entries
		jacob[irow_p2:irow_p2+2, icol_b:icol_b+2] .= f_u2_CtCabPdot./FORCE_SCALING
		jacob[irow_p2:irow_p2+2, icol_b+6:icol_b+8] .= f_u2_F./FORCE_SCALING
		jacob[irow_p2:irow_p2+2, icol_b+12:icol_b+14] .= f_u2_P./FORCE_SCALING

		jacob[irow_p2+3:irow_p2+5, icol_b+3:icol_b+5] .= f_ψ2_CtCabHdot./FORCE_SCALING
		jacob[irow_p2+3:irow_p2+5, icol_b+6:icol_b+8] .= f_ψ2_F./FORCE_SCALING
		jacob[irow_p2+3:irow_p2+5, icol_b+9:icol_b+11] .= f_ψ2_M./FORCE_SCALING
		jacob[irow_p2+3:irow_p2+5, icol_b+12:icol_b+14] .= f_ψ2_P./FORCE_SCALING
		jacob[irow_p2+3:irow_p2+5, icol_b+15:icol_b+17] .= f_ψ2_H./FORCE_SCALING

		# initialize compatability equation jacobian entries
		jacob[irow_p2+6:irow_p2+8, icol_b+6:icol_b+8] .= f_F2_F
		jacob[irow_p2+6:irow_p2+8, icol_b+9:icol_b+11] .= f_F2_M

		jacob[irow_p2+9:irow_p2+11, icol_b+6:icol_b+8] .= f_M2_F
		jacob[irow_p2+9:irow_p2+11, icol_b+9:icol_b+11] .= f_M2_M
	else
		# add to existing equilibrium equation jacobian entries
		jacob[irow_p2:irow_p2+2, icol_b:icol_b+2] .= f_u2_CtCabPdot./FORCE_SCALING
		jacob[irow_p2:irow_p2+2, icol_b+6:icol_b+8] .= f_u2_F./FORCE_SCALING
		jacob[irow_p2:irow_p2+2, icol_b+12:icol_b+14] .= f_u2_P./FORCE_SCALING

		jacob[irow_p2+3:irow_p2+5, icol_b+3:icol_b+5] .= f_ψ2_CtCabHdot./FORCE_SCALING
		jacob[irow_p2+3:irow_p2+5, icol_b+6:icol_b+8] .= f_ψ2_F./FORCE_SCALING
		jacob[irow_p2+3:irow_p2+5, icol_b+9:icol_b+11] .= f_ψ2_M./FORCE_SCALING
		jacob[irow_p2+3:irow_p2+5, icol_b+12:icol_b+14] .= f_ψ2_P./FORCE_SCALING
		jacob[irow_p2+3:irow_p2+5, icol_b+15:icol_b+17] .= f_ψ2_H./FORCE_SCALING

		if irow_b2 <= 0
			# compatability equations have been combined with other beam
			jacob[irow_p2+6:irow_p2+8, icol_b+6:icol_b+8] .= f_F2_F
			jacob[irow_p2+6:irow_p2+8, icol_b+9:icol_b+11] .= f_F2_M

			jacob[irow_p2+9:irow_p2+11, icol_b+6:icol_b+8] .= f_M2_F
			jacob[irow_p2+9:irow_p2+11, icol_b+9:icol_b+11] .= f_M2_M
		else
			# initialize compatability equation jacobian entries
			jacob[irow_b2+6:irow_b2+8, icol_b+6:icol_b+8] .= f_F2_F
			jacob[irow_b2+6:irow_b2+8, icol_b+9:icol_b+11] .= f_F2_M

			jacob[irow_b2+9:irow_b2+11, icol_b+6:icol_b+8] .= f_M2_F
			jacob[irow_b2+9:irow_b2+11, icol_b+9:icol_b+11] .= f_M2_M
		end

	end

	# initialize beam residual equation jacobian entries
	jacob[irow_b:irow_b+2, icol_b+12:icol_b+14] .= f_P_P
	jacob[irow_b:irow_b+2, icol_b+15:icol_b+17] .= f_P_H

	jacob[irow_b+3:irow_b+5, icol_b+12:icol_b+14] .= f_H_P
	jacob[irow_b+3:irow_b+5, icol_b+15:icol_b+17] .= f_H_H


	return jacob
end

# dynamic
@inline function element_jacobian!(jacob, irow_b, irow_b1, irow_p1, irow_b2, irow_p2, icol_b,
	f_u1_θ, f_u2_θ, f_u1_F, f_u2_F, f_u1_P, f_u2_P,
	f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M, f_ψ1_P, f_ψ2_P, f_ψ1_H, f_ψ2_H,
	f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
	f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M,
	f_P_u, f_P_θ, f_P_P, f_P_H,
	f_H_θ, f_H_P, f_H_H)

	# add jacobian entries corresponding to static case
	jacob = element_jacobian!(jacob, irow_b, irow_b1, irow_p1, irow_b2, irow_p2, icol_b,
		f_u1_θ, f_u2_θ, f_u1_F, f_u2_F,
		f_ψ1_θ, f_ψ2_θ, f_ψ1_F, f_ψ2_F, f_ψ1_M, f_ψ2_M,
		f_F1_u, f_F2_u, f_F1_θ, f_F2_θ, f_F1_F, f_F2_F, f_F1_M, f_F2_M,
		f_M1_θ, f_M2_θ, f_M1_F, f_M2_F, f_M1_M, f_M2_M)

	# create/add to jacobian entries for left endpoint
	jacob[irow_p1:irow_p1+2, icol_b+12:icol_b+14] .= f_u1_P./FORCE_SCALING
	jacob[irow_p1+3:irow_p1+5, icol_b+12:icol_b+14] .= f_ψ1_P./FORCE_SCALING
	jacob[irow_p1+3:irow_p1+5, icol_b+15:icol_b+17] .= f_ψ1_H./FORCE_SCALING

	# create/add to residual jacobian entries for right endpoint
	jacob[irow_p2:irow_p2+2, icol_b+12:icol_b+14] .= f_u2_P./FORCE_SCALING
	jacob[irow_p2+3:irow_p2+5, icol_b+12:icol_b+14] .= f_ψ2_P./FORCE_SCALING
	jacob[irow_p2+3:irow_p2+5, icol_b+15:icol_b+17] .= f_ψ2_H./FORCE_SCALING

	# initialize beam residual equation jacobian entries
	jacob[irow_b:irow_b+2, icol_b:icol_b+2] .= f_P_u
	jacob[irow_b:irow_b+2, icol_b+3:icol_b+5] .= f_P_θ
	jacob[irow_b:irow_b+2, icol_b+12:icol_b+14] .= f_P_P
	jacob[irow_b:irow_b+2, icol_b+15:icol_b+17] .= f_P_H

	jacob[irow_b+3:irow_b+5, icol_b+3:icol_b+5] .= f_H_θ
	jacob[irow_b+3:irow_b+5, icol_b+12:icol_b+14] .= f_H_P
	jacob[irow_b+3:irow_b+5, icol_b+15:icol_b+17] .= f_H_H

	return jacob
end
