"""
	system_residual!(resid, x, beams, prescribed_conditions, distributed_loads,
	irow_pt, irow_beam1, irow_beam2, icol_pt, icol_beam)

Populate the residual vector `resid` with the results of the residual equations
for the system

# Arguments
 - resid: Residual vector
 - x: Vector containing current state variables of the system
 - prescribed_conditions: Vector of prescribed conditions for each point with prescribed conditions
 - distributed_loads: Vector of distribued loads on each beam
 - irow_pt: Row index of first equilibrium equation for each point
 - irow_beam1: Row index of first equation for the left side of each beam
 - irow_beam2: Row index of first equation for the right side of each beam
 - icol_pt: Column index of first state variable for each point
 - icol_beam: Column index of first state variable for each beam element
"""
function system_residual!(resid, x, beams, prescribed_conditions, distributed_loads,
	irow_pt, irow_beam1, irow_beam2, icol_pt, icol_beam)

	nbeam = length(beams)

	# add contributions to residual equations from the elements
	for ibeam = 1:nbeam
		irow_b1 = irow_beam1[ibeam]
		irow_p1 = irow_pt[beams[ibeam].pt1]
		irow_b2 = irow_beam2[ibeam]
		irow_p2 = irow_pt[beams[ibeam].pt2]
		state_variables = element_state_variables(solution, icol_beam, ibeam, mode)
		properties = element_properties(beam, x0, v0, ω0, state_variables..., )
		internal_resultants = element_equations(properties...)
		resid = element_residual!(resid, irow_b1, irow_p1, irow_b2, irow_p2, internal_resultants...)
	end

	# add contributions to the residual equations from the prescribed point conditions
	for prescribed_condition in prescribed_conditions
		# get variables resulting from the prescribed point condition
		ipt = prescribed_condition.ipt
		u, θ, F, M = point_variables(solution, icol_pt, prescribed_condition)
		# search for beams that are connected to the specified point
		for ibeam = 1:nbeam
			# check left side of beam
			if ipt == beams[ibeam].pt1
				# add to residual equations for the beam endpoint
				irow_b = irow_beam1[ibeam]
				irow_p = irow_pt[ipt]
				point_residual!(resid, irow_b, irow_p, u, θ, F, M)
			end
			# check right side of beam
			if ipt == beams[ibeam].pt2
				# add to residual equations for the beam endpoint
				irow_b = irow_beam2[ibeam]
				irow_p = irow_pt[ipt]
				point_residual!(resid, irow_b, irow_p, u, θ, F, M)
			end
		end
	end
	return resid
end

"""
	element_residual!(resid, irow_b1, irow_p1, irow_b2, irow_p2, f_u1, f_u2, f_ψ1,
		f_ψ2, f_F1, f_F2, f_M1, f_M2, f_P, f_H)

Initialize the equilibrium and constitutive equations at irow_p1 and irow_b1,
for the left side of the beam element and irow_p2 and irow_b2 for the right side
of the beam element, respectively, to account for the beam element internal resultants given by
f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2, f_P, f_H

If irow_b1 != irow_p1 and/or irow_b2 != irow_p2, assume the equilibrium equations
have already been initialized for the left and/or right internal resultants, respectively.
"""
element_residual

function element_residual!(resid, irow_b1, irow_p1, irow_b2, irow_p2, f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2)

	# create/add to residual equations for left endpoint
	if irow_b1 == irow_p1
		# add equilibrium equations
		resid[irow_b1:irow_b1+2] = -f_u1
		resid[irow_b1+3:irow_b1+5] = -f_ψ1
		# add compatability equations
		resid[irow_b1+6:irow_b1+8] = -f_F1
		resid[irow_b1+9:irow_b1+11] = -f_M1
	else
		# add to existing equilibrium equations
		resid[irow_p1:irow_p1+2] -= f_u1
		resid[irow_p1+3:irow_p1+5] -= f_ψ1
		# add compatability equations
		resid[irow_b1:irow_b1+2] = -f_F1
		resid[irow_b1+3:irow_b1+5] = -f_M1
	end

	# create/add to residual equations for right endpoint
	if irow_b2 == irow_p2
		# add equilibrium equations
		resid[irow_b2:irow_b2+2] = -f_u2
		resid[irow_b2+3:irow_b2+5] = -f_ψ2
		# add compatability equations
		resid[irow_b2+6:irow_b2+8] = f_F2
		resid[irow_b2+9:irow_b2+11] = f_M2
	else
		# add to existing equilibrium equations
		resid[irow_p2:irow_p2+2] -= f_u2
		resid[irow_p2+3:irow_p2+5] -= f_ψ2
		# add compatability equations
		resid[irow_b2:irow_b2+2] = f_F2
		resid[irow_b2+3:irow_b2+5] = f_M2
	end

	return resid
end

function element_residual!(resid, irow_b1, irow_p1, irow_b2, irow_p2, f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2, f_P, f_H)

	resid = element_residual!(resid, irow_b1, irow_p1, irow_b2, irow_p2, f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2)

	# residual equations for element
	if irow_b1 == irow_p1
		resid[irow_b1+12:irow_b1+14] = f_P
		resid[irow_b1+15:irow_b1+17] = f_H
	else
		resid[irow_b1+6:irow_b1+8] = f_P
		resid[irow_b1+9:irow_b1+11] = f_H
	end

	return resid
end

"""
	point_residual!(resid, irow_b, irow_p, u, θ, F, M)

Modify the equilibrium and constitutive equations specified at irow_p and irow_b
respectively to account for the point variables given by u, θ, F, M

If irow_b == irow_p, assume that the equilibrium equations have already
been modified
"""
function point_residual!(resid, irow_b, irow_p, u, θ, F, M)

	if irow_b == irow_p
		# add to equilibrium and compatability equations
		resid[irow_b:irow_b+2] += F
		resid[irow_b+3:irow_b+5] += M
		resid[irow_b+6:irow_b+8] += u
		resid[irow_b+9:irow_b+11] += θ
	else
		# add to compatability equations
		resid[irow_b:irow_b+2] += u
		resid[irow_b+3:irow_b+5] += θ
	end

	return resid
end

"""
	point_variables(x, icol_p, prescribed_condition)

Extract u, θ, F, M for the point described by the state variables at `icol_p` in
x after incorporating the prescribed conditions in `prescribed_condition`

Note that the degrees of freedom that are not specified in `prescribed_condition`
are used as state variables (e.g. prescribing F[2] would mean u[2] = x[icol+1])
"""
function point_variables(x, icol_p, prescribed_condition)

	force = prescribed_condition.force

	# get the displacement and rotations of the point
	u = SVector(ifelse(force, x[icol_p  ], prescribed_condition.val[1]),
				ifelse(force, x[icol_p+1], prescribed_condition.val[2]),
				ifelse(force, x[icol_p+2], prescribed_condition.val[3]))
	θ = SVector(ifelse(force, x[icol_p+3], prescribed_condition.val[4]),
				ifelse(force, x[icol_p+4], prescribed_condition.val[5]),
				ifelse(force, x[icol_p+5], prescribed_condition.val[6]))

	# get the rotation matrix for the point
	Ct = wiener_milenkovic(θ)

	# solve for the force applied at the point due to the prescribed loads
	Fp = zero(u)
	for i = 1:3
		if cond.force[i]
			if cond.follower[i]
				Fp += SVector(CT[1,i], Ct[2,i], Ct[3,i])*cond.val[i]
			else
				Fp += SVector(I3[1,i], I3[2,i], I3[3,i])*cond.val[i]
			end
		end
	end

	# solve for the moment applied at the point due to the prescribed loads
	Mp = zero(θ)
	for i = 4:6
		if cond.force[i]
			if cond.follower[i]
				Mp += SVector(CT[1,i], Ct[2,i], Ct[3,i])*cond.val[i]
			else
				Mp += SVector(I3[1,i], I3[2,i], I3[3,i])*cond.val[i]
			end
		end
	end

	# overwrite external forces/moments with solved for forces/moments
	F = SVector(ifelse(force, Fp[1], x[icol  ]),
				ifelse(force, Fp[2], x[icol+1]),
				ifelse(force, Fp[3], x[icol+2]))
	M = SVector(ifelse(force, Mp[1], x[icol+3]),
				ifelse(force, Mp[2], x[icol+4]),
				ifelse(force, Mp[3], x[icol+5]))

	return u, θ, F, M
end
