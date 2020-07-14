"""
	steady_state_analysis(assembly; kwargs...)

Performs a steady-state analysis for the system of nonlinear beams contained in `assembly`

# Keyword Arguments
 - `prescribed_conditions=Dict{Int,PrescribedConditions{Float64}}()`: Dictionary
 	holding PrescribedConditions composite types for the points in `keys(prescribed_conditions)`
 - `distributed_loads=Dict{Int,DistributedLoads{Float64}}()`: Dictionary holding
 	DistributedLoads composite types for the beam elements in `keys(distributed_loads)`
 - `time_functions = TimeFunction{Float64}[]`: Array of time functions (of type TimeFunction)
 - `origin = zeros(3)`: Global frame origin (for dynamic simulations)
 - `linear_velocity = zeros(3)`: Global frame linear velocity (for dynamic simulations)
 - `angular_velocity = zeros(3)`: Global frame angular velocity (for dynamic simulations)
 - `linear_velocity_tf = zeros(Int, 3)`: Time function for global frame linear velocity
 - `angular_velocity_tf = zeros(Int, 3)`: Time function for global frame angular velocity
 - `method=:newton`: Method (as defined in NLsolve) to solve nonlinear system of equations
 - `linesearch=LineSearches.BackTracking()`: Line search used to solve nonlinear system of equations
 - `ftol=1e-12`: tolerance for solving nonlinear system of equations
 - `iterations=1000`: maximum iterations for solving the nonlinear system of equations
 - `nstep=1`: Number of time steps. May be used in conjunction with `time_functions`
   to gradually increase displacements/loads.
"""
function steady_state_analysis(assembly;
	prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}(),
	distributed_loads = Dict{Int,DistributedLoads{Float64}}(),
	time_functions = TimeFunction{Float64}[],
	origin = (@SVector zeros(3)),
	linear_velocity = (@SVector zeros(3)),
	angular_velocity = (@SVector zeros(3)),
	linear_velocity_tf = (@SVector zeros(Int, 3)),
	angular_velocity_tf = (@SVector zeros(Int, 3)),
	method = :newton,
	linesearch = BackTracking(),
	ftol = 1e-12,
	iterations = 1000,
	nstep = 1,
	)

	TF = eltype(assembly)

	static = linear_velocity == zero(linear_velocity) || angular_velocity == zero(angular_velocity)

	n_connections = point_connections(assembly)

	n, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam = system_indices(
		assembly, prescribed_conditions, n_connections, static)

	F = zeros(TF, n)
	J = spzeros(TF, n, n)
	x = rand(TF, n)
	converged = true

	# initialize evaluated time functions
	n_tf = length(time_functions)
	time_function_values = OffsetArray(ones(eltype(eltype(time_functions)), n_tf+1), 0:n_tf)

	# default time function does nothing
	time_function_values[0] = 1

	# convert inputs for global frame motion to static arrays
	x0 = SVector{3}(origin)
	v0 = SVector{3}(linear_velocity)
	ω0 = SVector{3}(angular_velocity)
	v0_tf = SVector{3}(linear_velocity_tf)
	ω0_tf = SVector{3}(angular_velocity_tf)

	for istep = 1:nstep

		# update time functions (except default time function)
		time_function_values[1:n_tf] .= getindex.(time_functions, istep)

		# update global frame motion
		v0 = linear_velocity .* time_function_values[v0_tf]
		ω0 = angular_velocity .* time_function_values[ω0_tf]

		# set extra arguments for static or dynamic case
		if static
			args = ()
		else
			args = (x0, v0, ω0)
		end

		# solve the nonlinear system of equations
		f! = (F, x) -> system_residual!(F, x, assembly, prescribed_conditions,
			distributed_loads, time_function_values, irow_pt, irow_beam,
			irow_beam1, irow_beam2, icol_pt, icol_beam, args...)

		j! = (J, x) -> system_jacobian!(J, x, assembly, prescribed_conditions,
			distributed_loads, time_function_values, irow_pt, irow_beam,
			irow_beam1, irow_beam2, icol_pt, icol_beam, args...)

		df = NLsolve.OnceDifferentiable(f!, j!, x, F, J)

		result = NLsolve.nlsolve(f!, x,
			method=method,
			linesearch=linesearch,
			ftol=ftol,
			iterations=iterations)

		# update solution
		x .= result.zero

		# stop early if unconverged
		if !result.f_converged
			converged = false
			break
		end
	end

	if static
		args = ()
	else
		args = (x0, v0, ω0)
	end

	return AssemblyState(converged, x,
		assembly, prescribed_conditions, distributed_loads, time_function_values,
		irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam, args...)
end

function eigenvalue_analysis(assembly, prescribed_conditions, distributed_loads,
	x0, ω0, v0, niter, nev)

	TF = eltype(assembly)

	nbeam = length(assembly.elements)
	npoint = length(assembly.points)

	static = false

	n, irow_pt, irow_beam1, irow_beam2, icol_pt, icol_beam = system_indices(
		assembly.points, assembly.elements, static)

	# solve steady state problem
	x = zeros(TF, n)
	F = zeros(TF, n)
	J = spzeros(TF, n, n)

	f! = (F, x) -> system_residual!(F, x, assembly, prescribed_conditions,
		distributed_loads, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam,
		x0, v0, ω0)

	j! = (J, x) -> system_jacobian!(J, x, assembly, prescribed_conditions,
		distributed_loads, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam,
		x0, v0, ω0)

	df = OnceDifferentiable(f!, j!, x, F, J)

	nlsolve!(df, initial_x)

	# solve eigenvalue problem
	K = J
	M = spzeros(TF, n, n)

	M = mass_matrix!(M, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam,
		prescribed_conditions, distributed_loads, x0, ω0, v0)

	λ, ϕ, nconv, niter, nmult, resid = Arpack.eigs(K, M; nev=nev, which=:LM)

	return λ, [AssemblyState(ϕ[:,i], prescribed_conditions, icol_pt, icol_beam, static) for i = 1:size(ϕ, 2)]
end

"""
	time_domain_analysis(assembly, dt; kwargs...)

Performs a steady-state analysis for the system of nonlinear beams contained in `assembly`

# Keyword Arguments
 - `prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}()`: Dictionary
 	holding PrescribedConditions composite types for the points in `keys(prescribed_conditions)`
 - `distributed_loads = Dict{Int,DistributedLoads{Float64}}()`: Dictionary holding
 	DistributedLoads composite types for the beam elements in `keys(distributed_loads)`
 - `time_functions = TimeFunction{Float64}[]`: Array of time functions (of type TimeFunction)
 - `origin = zeros(3)`: Global frame origin (for dynamic simulations)
 - `linear_velocity = zeros(3)`: Global frame linear velocity (for dynamic simulations)
 - `angular_velocity = zeros(3)`: Global frame angular velocity (for dynamic simulations)
 - `linear_velocity_tf = zeros(Int, 3)`: Time function for global frame linear velocity
 - `angular_velocity_tf = zeros(Int, 3)`: Time function for global frame angular velocity
 - `u0=fill(zeros(3), length(assembly.elements))`: initial displacment of each beam element,
 - `theta0=fill(zeros(3), length(assembly.elements))`: initial angular displacement of each beam element,
 - `udot0=fill(zeros(3), length(assembly.elements))`: initial time derivative with respect to `u`
 - `thetadot0=fill(zeros(3), length(assembly.elements))`: initial time derivative with respect to `theta`
 - `method=:newton`: Method (as defined in NLsolve) to solve nonlinear system of equations
 - `linesearch=LineSearches.BackTracking()`: Line search used to solve nonlinear system of equations
 - `ftol=1e-12`: tolerance for solving nonlinear system of equations
 - `iterations=1000`: maximum iterations for solving the nonlinear system of equations
 - `nstep=1`: Number of time steps.
"""
function time_domain_analysis(assembly, dt;
	prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}(),
	distributed_loads = Dict{Int,DistributedLoads{Float64}}(),
	time_functions = TimeFunction{Float64}[],
	origin = (@SVector zeros(3)),
	linear_velocity = (@SVector zeros(3)),
	angular_velocity = (@SVector zeros(3)),
	linear_velocity_tf = (@SVector zeros(Int, 3)),
	angular_velocity_tf = (@SVector zeros(Int, 3)),
	u0 = fill((@SVector zeros(3)), length(assembly.elements)),
	theta0 = fill((@SVector zeros(3)), length(assembly.elements)),
	udot0 = fill((@SVector zeros(3)), length(assembly.elements)),
	thetadot0 = fill((@SVector zeros(3)), length(assembly.elements)),
	method = :newton,
	linesearch = BackTracking(),
	ftol = 1e-12,
	iterations = 1000,
	nstep = 1,
	)

	TF = eltype(assembly)

	static = false

	nbeam = length(assembly.elements)

	n_connections = point_connections(assembly)

	n, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam = system_indices(
		assembly, prescribed_conditions, n_connections, static)

	F = zeros(TF, n)
	J = spzeros(TF, n, n)
	x = zeros(TF, n)

	# initialize evaluated time functions
	n_tf = length(time_functions)
	time_function_values = OffsetArray(ones(eltype(eltype(time_functions)), n_tf+1), 0:n_tf)

	# default time function does nothing
	time_function_values[0] = 1

	# convert inputs for global frame motion to static arrays
	x0 = SVector{3}(origin)
	v0 = SVector{3}(linear_velocity)
	ω0 = SVector{3}(angular_velocity)
	v0_tf = SVector{3}(linear_velocity_tf)
	ω0_tf = SVector{3}(angular_velocity_tf)

	# convert initial displacments/velocities to static arrays
	u = SVector{3, TF}.(u0)
	θ = SVector{3, TF}.(theta0)
	udot = SVector{3, TF}.(udot0)
	θdot = SVector{3, TF}.(thetadot0)

	# update time functions
	time_function_values[1:n_tf] .= getindex.(time_functions, 1)

	# update global frame motion
	v0 = linear_velocity .* time_function_values[v0_tf]
	ω0 = angular_velocity .* time_function_values[ω0_tf]

	# set extra arguments for initial step case
	args = (x0, v0, ω0, u, θ, udot, θdot)

	# construct residual and jacobian functions
	f! = (F, x) -> system_residual!(F, x, assembly, prescribed_conditions,
		distributed_loads, time_function_values, irow_pt, irow_beam,
		irow_beam1, irow_beam2, icol_pt, icol_beam, args...)
	j! = (J, x) -> system_jacobian!(J, x, assembly, prescribed_conditions,
		distributed_loads, time_function_values, irow_pt, irow_beam,
		irow_beam1, irow_beam2, icol_pt, icol_beam, args...)

	# solve nonlinear system of equations
	df = OnceDifferentiable(f!, j!, x, F, J)

	result = NLsolve.nlsolve(df, x,
		method=method,
		linesearch=linesearch,
		ftol=ftol,
		iterations=iterations)

	# update solution
	x .= result.zero

	# set converged flag
	converged = result.f_converged

	# now set up for the time-domain run
	udot_init = Vector{SVector{3,TF}}(undef, nbeam)
	θdot_init = Vector{SVector{3,TF}}(undef, nbeam)
	CtCabPdot_init = Vector{SVector{3, TF}}(undef, nbeam)
	CtCabHdot_init = Vector{SVector{3, TF}}(undef, nbeam)
	for ibeam = 1:nbeam
		icol = icol_beam[ibeam]
		# calculate udot
		udot_init[ibeam] = 2/dt*u[ibeam] + udot[ibeam]
		# calculate θdot + 2/dt*θ (remaining term will be added back in later)
		θdot_init[ibeam] = 2/dt*θ[ibeam] + θdot[ibeam]
		# extract rotation parameters
		Ct = wiener_milenkovic(θ[ibeam])'
		Cab = assembly.elements[ibeam].Cab
		# calculate CtCabPdot
		CtCabP = Ct*Cab*SVector{3, TF}(x[icol+12], x[icol+13], x[icol+14])
		CtCabPdot = SVector{3, TF}(x[icol], x[icol+1], x[icol+2])
		CtCabPdot_init[ibeam] = 2/dt*CtCabP + CtCabPdot
		# calculate CtCabHdot
		CtCabH = Ct*Cab*SVector{3, TF}(x[icol+15], x[icol+16], x[icol+17])
		CtCabHdot = SVector{3, TF}(x[icol+3], x[icol+4], x[icol+5])
		CtCabHdot_init[ibeam] = 2/dt*CtCabH + CtCabHdot
		# insert initial conditions for time-domain analysis
		x[icol:icol+2] .= u[ibeam]
		x[icol+3:icol+5] .= θ[ibeam]
	end

	# initialize storage for each time step
	history = Vector{AssemblyState{TF}}(undef, nstep)

	# add state for initial conditions to the solution history
	history[1] = AssemblyState(converged, x,
		assembly, prescribed_conditions, distributed_loads, time_function_values,
		irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam, args...)

	# loop for each time step
	for istep = 2:nstep
		# update time functions (except default time function)
		time_function_values[1:n_tf] .= getindex.(time_functions, istep)

		# update global frame motion
		v0 = linear_velocity .* time_function_values[v0_tf]
		ω0 = angular_velocity .* time_function_values[ω0_tf]

		# set extra arguments
		args = (x0, v0, ω0, udot_init, θdot_init, CtCabPdot_init, CtCabHdot_init, dt)

		# solve for the state variables at next time step
		f! = (F, x) -> system_residual!(F, x, assembly, prescribed_conditions,
			distributed_loads, time_function_values, irow_pt, irow_beam,
			irow_beam1, irow_beam2, icol_pt, icol_beam, args...)

		j! = (J, x) -> system_jacobian!(J, x, assembly, prescribed_conditions,
			distributed_loads, time_function_values, irow_pt, irow_beam,
			irow_beam1, irow_beam2, icol_pt, icol_beam, args...)

		df = OnceDifferentiable(f!, j!, x, F, J)

		result = NLsolve.nlsolve(df, x,
			method=method,
			linesearch=linesearch,
			ftol=ftol,
			iterations=iterations)

		# update solution
		x .= result.zero

		# add state to history
		history[istep] = AssemblyState(converged, x,
			assembly, prescribed_conditions, distributed_loads, time_function_values,
			irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam, args...)

		# stop early if unconverged
		if !converged
			break
		end

		# exit loop if done iterating
		if istep == nstep
			break
		end

		# update time derivative terms for the next time step
		for ibeam = 1:nbeam
			icol = icol_beam[ibeam]
			# calculate udot for next time step
			u = SVector{3, TF}(x[icol], x[icol+1], x[icol+2])
			udot = 2/dt*u - udot_init[ibeam]
			udot_init[ibeam] = 2/dt*u + udot
			# calculate θdot for next time step
			θ = SVector{3, TF}(x[icol+3], x[icol+4], x[icol+5])
			θdot = 2/dt*θ - θdot_init[ibeam]
			θdot_init[ibeam] = 2/dt*θ + θdot
			# extract rotation parameters
			Ct = wiener_milenkovic(θ)'
			Cab = assembly.elements[ibeam].Cab
			# calculate CtCabPdot for next time step
			CtCabP = Ct*Cab*SVector{3, TF}(x[icol+12], x[icol+13], x[icol+14])
			CtCabPdot = 2/dt*CtCabP - CtCabPdot_init[ibeam]
			CtCabPdot_init[ibeam] = 2/dt*CtCabP + CtCabPdot
			# calculate CtCabHdot for next time step
			CtCabH = Ct*Cab*SVector{3, TF}(x[icol+15], x[icol+16], x[icol+17])
			CtCabHdot = 2/dt*CtCabH - CtCabHdot_init[ibeam]
			CtCabHdot_init[ibeam] = 2/dt*CtCabH + CtCabHdot
		end
	end

	return history
end
