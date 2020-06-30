"""
	system_jacobian!(jacob, x, beams, prescribed_conditions, distributed_loads,
		irow_pt, irow_beam1, irow_beam2, icol_pt, icol_beam)
	system_jacobian!(jacob, x, beams, prescribed_conditions, distributed_loads,
		irow_pt, irow_beam1, irow_beam2, icol_pt, icol_beam, x0, v0, ω0)
	system_jacobian!(jacob, x, beams, prescribed_conditions, distributed_loads,
		irow_pt, irow_beam1, irow_beam2, icol_pt, icol_beam, x0, v0, ω0,
		udot, θdot, CtCabPdot, CtCabHdot)
	system_jacobian!(jacob, x, beams, prescribed_conditions, distributed_loads,
		irow_pt, irow_beam1, irow_beam2, icol_pt, icol_beam, x0, v0, ω0,
		u_p, θ_p, CtCabP_p, CtCabH_p, udot_p, θdot_p, CtCabPdot_p, CtCabHdot_p, dt)

Populate the jacobian matrix `jacob` with the gradient of the residual equations
with respect to the state variables.

There are four implementations corresponding to the following analysis types:
 - Static
 - Dynamic - Steady State
 - Dynamic - Initial Step (for initializing time domain simulations)
 - Dynamic - Time Marching

 See "GEBT: A general-purpose nonlinear analysis tool for composite beams" by
 Wenbin Yu and "Geometrically nonlinear analysis of composite beams using
 Wiener-Milenković parameters" by Qi Wang and Wenbin Yu.

# Arguments
 - jacob: Jacobian matrix
 - x: Vector containing current state variables of the system
 - beams: Vector containing beam elements
 - prescribed_conditions: Vector of prescribed conditions for each point
 - distributed_loads: Vector of distributed loads for each beam element
 - irow_pt: Row index of first equilibrium equation for each point
 - irow_beam1: Row index of first equation for the left side of each beam
 - irow_beam2: Row index of first equation for the right side of each beam
 - icol_pt: Column index of first state variable for each point
 - icol_beam: Column index of first state variable for each beam element

# Additional Arguments for Dynamic Analyses
 - x0: Global frame origin (only for dynamic analyses)
 - v0: Global frame linear velocity (only for dynamic analyses)
 - ω0: Global frame angular velocity (only for dynamic analyses)

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
function system_jacobian!(jacob, x, beams, prescribed_conditions, distributed_loads,
	irow_pt, irow_beam1, irow_beam2, icol_pt, icol_beam, args...)

	nbeam = length(beams)

	# add contributions to residual equations from the elements
	for ibeam = 1:nbeam
		icol = icol_beam[ibeam]
		irow_b1 = irow_beam1[ibeam]
		irow_p1 = irow_pt[beams[ibeam].pt1]
		irow_b2 = irow_beam2[ibeam]
		irow_p2 = irow_pt[beams[ibeam].pt2]
		properties = element_properties(x, icol, ibeam, beam, args...)
		internal_resultants = element_equations(distributed_loads[ibeam], properties...)
		jacob = element_jacobian!(jacob, irow_b1, irow_p1, irow_b2, irow_p2, internal_resultants...)
	end

	# add contributions to the residual equations from the prescribed point conditions
	for (point, prescribed_condition) in enumerate(prescribed_conditions)
		# get variables resulting from the prescribed point condition
		u, θ, F, M = point_variables(solution, icol_pt[point], prescribed_condition)
		# search for beams that are connected to the specified point
		for ibeam = 1:nbeam
			# check left side of beam
			if point == beams[ibeam].pt1
				# add to residual equations for the beam endpoint
				irow_b = irow_beam1[ibeam]
				irow_p = irow_pt[point]
				point_residual!(resid, irow_b, irow_p, u, θ, F, M)
			end
			# check right side of beam
			if point == beams[ibeam].pt2
				# add to residual equations for the beam endpoint
				irow_b = irow_beam2[ibeam]
				irow_p = irow_pt[point]
				point_residual!(resid, irow_b, irow_p, u, θ, F, M)
			end
		end
	end
	return resid
end
