"""
	point_connections(assembly)

Count the number of beams connected to each point
"""
function point_connections(assembly)

	npt = length(assembly.points)

	n_connections = Array{Int, 1}(undef, npt)

	for ipt = 1:npt
		n_connections[ipt] = count(x -> x == ipt, assembly.start) + count(x -> x == ipt, assembly.stop)
	end

	return n_connections
end

"""
    system_indices(assembly, static)

Solve for the row indices of the first equilibrium or compatability equations for
each point and side of each beam element.  Also solve for the row/column index of
each point and beam state variable.

If only two beams meet at a point, and no prescribed conditions are applied to
that point, the 6 unknowns associated with that point as well as 6 compatability
equations are eliminated from the system.  Points for which unknowns have been
eliminated are assigned a column index of -1.  Beams for which the compatability
equations have been eliminated are also assigned an index of -1

# Arguments:
 - `assembly`: object of type Assembly
 - `prescribed_conditions`: Dictionary of prescribed point conditions
 - `n_connections`: Number of connections to each point
 - `static`: flag indicating whether analysis is static

# Return Arguments:
 - `n`: total number of equations/unknowns in the system
 - `irow_pt`: Row index of first equilibrium equation for each point
 - `irow_beam`: Row index of first equation for just this beam element
 - `irow_beam1`: Row index of first equation for the left side of each beam
 - `irow_beam2`: Row index of first equation for the right side of each beam
 - `icol_pt`: Column index of first state variable for each point
 - `icol_beam`: Column index of first state variable for each beam element
"""
function system_indices(assembly, prescribed_conditions, n_connections, static)

    npt = length(assembly.points)
    nbeam = length(assembly.elements)

	assigned = fill(false, npt)

    irow_pt = Array{Int, 1}(undef, npt)
	irow_beam = Array{Int, 1}(undef, nbeam)
    irow_beam1 = Array{Int, 1}(undef, nbeam)
    irow_beam2 = Array{Int, 1}(undef, nbeam)

    icol_pt = Array{Int, 1}(undef, npt)
    icol_beam = Array{Int, 1}(undef, nbeam)

    irow = 1
    icol = 1
    for ibeam = 1:nbeam
        ipt = assembly.start[ibeam]
        if !assigned[ipt]
			assigned[ipt] = true
            # 6 equilibrium equations + 6 compatability equations
            irow_pt[ipt] = irow
            irow_beam1[ibeam] = irow
            irow += 12
			# add unknowns for each point
			if n_connections[ipt] == 2 && !haskey(prescribed_conditions, ipt)
				# no additional unknowns when only two beams meet at a point
				icol_pt[ipt] = -1
			else
				# 6 additional unknowns for each point
				icol_pt[ipt] = icol
            	icol += 6
			end
        else
			# add compatability equations
			if n_connections[ipt] == 2 && !haskey(prescribed_conditions, ipt)
				# no additional equations when only two beams meet at a point
				irow_beam1[ibeam] = -1
			else
				# 6 additional compatibility equations
            	irow_beam1[ibeam] = irow
            	irow += 6
			end
        end

        # 12 unknowns for each element
        icol_beam[ibeam] = icol
        icol += 12

        if static
			irow_beam[ibeam] = -1
		else
            # 6 linear and angular momentum residual equations for each element
			irow_beam[ibeam] = irow
			irow += 6
            # 6 additional unknowns for each member for unsteady analyses
            icol += 6
        end

        ipt = assembly.stop[ibeam]
        if !assigned[ipt]
			assigned[ipt] = true
            # 6 equilibrium equations + 6 compatability equations
            irow_pt[ipt] = irow
            irow_beam2[ibeam] = irow
            irow += 12
			# add unknowns for each point
			if n_connections[ipt] == 2 && !haskey(prescribed_conditions, ipt)
				# no additional unknowns when only two beams meet at a point
				icol_pt[ipt] = -1
			else
				# 6 additional unknowns for each point
				icol_pt[ipt] = icol
				icol += 6
			end
        else
			if n_connections[ipt] == 2 && !haskey(prescribed_conditions, ipt)
				# no additional compatability equations when only two beams meet at a point
				irow_beam2[ibeam] = -1
			else
				# 6 additional compatibility equations
	            irow_beam2[ibeam] = irow
	            irow += 6
			end
        end
    end

	# number of unknowns/equations
	n = irow - 1

    return n, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam
end

# function scale_state_variables!(x, prescribed_conditions, icol_pt, icol_beam)
#
# 	npoint = length(icol_pt)
# 	nbeam = length(icol_beam)
#
# 	for ipoint = 1:npoint
# 		if icol_pt[ipoint] > 0
# 			prescribed = haskey(prescribed_conditions, ipoint)
# 			if prescribed
# 				force_dof = prescribed_conditions[ipoint].force_dof
# 				for i in findall(x -> x == false, force_dof)
# 					# force is a state variable
# 					x[icol_pt[ipoint]+i-1] *= FORCE_SCALING
# 				end
# 			end
# 		end
# 	end
#
# 	for ibeam = 1:nbeam
# 		icol = icol_beam[ibeam]
# 		x[icol+6] *= FORCE_SCALING
# 		x[icol+7] *= FORCE_SCALING
# 		x[icol+8] *= FORCE_SCALING
# 		x[icol+9] *= FORCE_SCALING
# 		x[icol+10] *= FORCE_SCALING
# 		x[icol+11] *= FORCE_SCALING
# 	end
#
# 	return x
# end
#
# function descale_state_variables!(x, prescribed_conditions, icol_pt, icol_beam)
#
# 	npoint = length(icol_pt)
# 	nbeam = length(icol_beam)
#
# 	for ipoint = 1:npoint
# 		if icol_pt[ipoint] > 0
# 			prescribed = haskey(prescribed_conditions, ipoint)
# 			if prescribed
# 				force_dof = prescribed_conditions[ipoint].force_dof
# 				for i in findall(x -> x == false, force_dof)
# 					# force is a state variable
# 					x[icol_pt[ipoint]+i-1] /= FORCE_SCALING
# 				end
# 			end
# 		end
# 	end
#
# 	for ibeam = 1:nbeam
# 		icol = icol_beam[ibeam]
# 		x[icol+6] /= FORCE_SCALING
# 		x[icol+7] /= FORCE_SCALING
# 		x[icol+8] /= FORCE_SCALING
# 		x[icol+9] /= FORCE_SCALING
# 		x[icol+10] /= FORCE_SCALING
# 		x[icol+11] /= FORCE_SCALING
# 	end
#
# 	return x
# end

"""
	system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads,
		time_function_values, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam)
	system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads,
		time_function_values, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam, x0, v0, ω0)
	system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads,
		time_function_values, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam, x0, v0, ω0,
		u, θ, udot, θdot)
	system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads,
		time_function_values, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam, x0, v0, ω0,
		udot, θdot_init, CtCabPdot, CtCabHdot, dt)

Populate the residual vector `resid` with the results of the residual equations
for the system.

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
 - `x`: Current state variables of the system
 - `assembly`: Assembly of rigidly connected nonlinear beam elements
 - `prescribed_conditions`: Dictionary of prescribed conditions
 - `distributed_loads`: Dictionary of distributed loads
 - `time_function_values`: Time functions, evaluated at the current time step
 - `irow_pt`: Row index of first equilibrium equation for each point
 - `irow_beam`: Row index of first equation for just this beam element
 - `irow_beam1`: Row index of first equation for the left side of each beam
 - `irow_beam2`: Row index of first equation for the right side of each beam
 - `icol_pt`: Column index of first state variable for each point
 - `icol_beam`: Column index of first state variable for each beam element

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
function system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads,
	time_function_values, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam, args...)

	npoint = length(assembly.points)
	nbeam = length(assembly.elements)

	# add contributions to residual equations from the elements
	for ibeam = 1:nbeam

		icol_b = icol_beam[ibeam]
		irow_b = irow_beam[ibeam]
		irow_b1 = irow_beam1[ibeam]
		irow_p1 = irow_pt[assembly.start[ibeam]]
		irow_b2 = irow_beam2[ibeam]
		irow_p2 = irow_pt[assembly.stop[ibeam]]

		# get element properties
		properties = element_properties(x, icol_b, ibeam, assembly.elements[ibeam], args...)
		# integrate distributed loads if applicable
		loads = haskey(distributed_loads, ibeam)
		if loads
			ΔL, Ct, _ = properties
			f1, f2, m1, m2 = integrate_element_loads(ΔL, Ct, distributed_loads, time_function_values)
		else
			f1 = @SVector zeros(3)
			f2 = @SVector zeros(3)
			m1 = @SVector zeros(3)
			m2 = @SVector zeros(3)
		end
		# get the resultants from the element equations
		internal_resultants = element_equations(properties..., f1, f2, m1, m2)
		# assemble into the system residual vector
		resid = element_residual!(resid, irow_b, irow_b1, irow_p1, irow_b2, irow_p2, internal_resultants...)
	end

	# add contributions to the residual equations from the prescribed point conditions
	for ipoint = 1:npoint

		# skip if the unknowns have been eliminated from the system of equations
		if icol_pt[ipoint] <= 0
			continue
		end

		# incorporate prescribed condition if applicable
		prescribed = haskey(prescribed_conditions, ipoint)
		if prescribed
			u, θ, F, M = point_variables(x, icol_pt[ipoint], prescribed_conditions[ipoint], time_function_values)
		else
			u, θ, F, M = point_variables(x, icol_pt[ipoint])
		end

		# search for beams that are connected to the specified point
		for ibeam = 1:nbeam
			# check left side of beam
			if ipoint == assembly.start[ibeam]
				# add to residual equations for the beam endpoint
				side = -1
                irow_b = irow_beam1[ibeam]
				irow_p = irow_pt[ipoint]
				point_residual!(resid, irow_b, irow_p, u, θ, F, M, side)
			end
			# check right side of beam
			if ipoint == assembly.stop[ibeam]
				# add to residual equations for the beam endpoint
				side = 1
				irow_b = irow_beam2[ibeam]
				irow_p = irow_pt[ipoint]
				point_residual!(resid, irow_b, irow_p, u, θ, F, M, side)
			end
		end
	end

	return resid
end

"""
	system_jacobian!(jacob, x, assembly, prescribed_conditions, distributed_loads,
		time_function_values, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam)
	system_jacobian!(jacob, x, assembly, prescribed_conditions, distributed_loads,
		time_function_values, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam,
		x0, v0, ω0)
	system_jacobian!(jacob, x, assembly, prescribed_conditions, distributed_loads,
		time_function_values, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam,
		x0, v0, ω0, u, θ, udot, θdot)
	system_jacobian!(jacob, x, assembly, prescribed_conditions, distributed_loads,
		time_function_values, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam,
		x0, v0, ω0, udot_init, θdot_init, CtCabPdot_init, CtCabHdot_init, dt)

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
 - `jacob`: Jacobian matrix
 - `x`: Vector containing current state variables of the system
 - `assembly`: Vector containing beam elements
 - `prescribed_conditions`: Vector of prescribed conditions for each point
 - `distributed_loads`: Vector of distributed loads for each beam element
 - `time_function_values`: Time functions, evaluated at the current time step
 - `irow_pt`: Row index of first equilibrium equation for each point
 - `irow_beam`: Row index of first equation for just this beam element
 - `irow_beam1`: Row index of first equation for the left side of each beam
 - `irow_beam2`: Row index of first equation for the right side of each beam
 - `icol_pt`: Column index of first state variable for each point
 - `icol_beam`: Column index of first state variable for each beam element

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
function system_jacobian!(jacob, x, assembly, prescribed_conditions, distributed_loads,
	time_function_values, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam, args::Vararg{<:Any, N}) where N

	npoint = length(assembly.points)
	nbeam = length(assembly.elements)

	initial_step = N == 7

	# add contributions to residual equations from the elements
	for ibeam = 1:nbeam

		icol = icol_beam[ibeam]
		irow_b = irow_beam[ibeam]
		irow_b1 = irow_beam1[ibeam]
		irow_p1 = irow_pt[assembly.start[ibeam]]
		irow_b2 = irow_beam2[ibeam]
		irow_p2 = irow_pt[assembly.stop[ibeam]]

		# get beam element properties
		properties = element_properties(x, icol, ibeam, assembly.elements[ibeam], args...)

		# pre-calculate jacobian of wiener-milenkovic rotation matrix wrt θ
		ΔL, Ct, θ = properties[1], properties[2], properties[6]
		C_θ1, C_θ2, C_θ3 = wiener_milenkovic_jacobian(Ct', θ)
		Ct_θ1, Ct_θ2, Ct_θ3 = C_θ1', C_θ2', C_θ3'

		# get jacobian of integrated loads if applicable
		loads = haskey(distributed_loads, ibeam)
		if loads
			f1_θ, f2_θ, m1_θ, m2_θ = follower_load_jacobians(ΔL, Ct_θ1, Ct_θ2, Ct_θ3, distributed_loads, time_function_values)
		else
			f1_θ = @SMatrix zeros(3, 3)
			f2_θ = @SMatrix zeros(3, 3)
			m1_θ = @SMatrix zeros(3, 3)
			m2_θ = @SMatrix zeros(3, 3)
		end

		# get jacobians of beam element equations
		internal_resultant_jacobians = element_jacobian_equations(assembly.elements[ibeam], properties..., Ct_θ1, Ct_θ2, Ct_θ3, f1_θ, f2_θ, m1_θ, m2_θ)

		# initialize/insert into jacobian matrix for the system
		jacob = element_jacobian!(jacob, irow_b, irow_b1, irow_p1, irow_b2, irow_p2, icol, internal_resultant_jacobians...)
	end

	# add contributions to the system jacobian matrix from the prescribed point conditions
	for ipoint = 1:npoint

		# skip if the unknowns have been eliminated from the system of equations
		if icol_pt[ipoint] <= 0
			continue
		end

		# incorporate prescribed condition if applicable
		prescribed = haskey(prescribed_conditions, ipoint)
		if prescribed
			F_θ, M_θ = point_follower_jacobians(x, icol_pt[ipoint], prescribed_conditions[ipoint], time_function_values)
			force_dof = prescribed_conditions[ipoint].force_dof
		else
			F_θ = @SMatrix zeros(3,3)
			M_θ = @SMatrix zeros(3,3)
			force_dof = @SVector ones(Bool, 6)
		end

		# search for beams that are connected to the specified point
		for ibeam = 1:nbeam
			# check left side of beam
			if ipoint == assembly.start[ibeam]
				# add to jacobian entries for the beam endpoint
				side = -1
				irow_b = irow_beam1[ibeam]
				irow_p = irow_pt[ipoint]
				icol = icol_pt[ipoint]
				point_jacobian!(jacob, irow_b, irow_p, icol, side, force_dof, F_θ, M_θ)
			end
			# check right side of beam
			if ipoint == assembly.stop[ibeam]
				# add to jacobian entries for the beam endpoint
				side = 1
				irow_b = irow_beam2[ibeam]
				irow_p = irow_pt[ipoint]
				icol = icol_pt[ipoint]
				point_jacobian!(jacob, irow_b, irow_p, icol, side, force_dof, F_θ, M_θ)
			end
		end
	end
	return jacob
end
