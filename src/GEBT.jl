module GEBT

using LinearAlgebra
using StaticArrays
using SparseArrays
using FLOWMath

include("math.jl")
include("setup.jl")
include("assembly.jl")
include("loads.jl")

function static_analysis(assembly, prescribed_conditions, distributed_loads)

	TF = eltype(assembly)

	static = true

	N, irow_pt, irow_beam1, irow_beam2, icol_pt, icol_beam = system_indices(
		assembly.points, assembly.elements, static)

	x = zeros(TF, N)
	A = spzeros(TF, N, N)
	b = zeros(TF, N)

	A = system_jacobian!(A, x, assembly.elements, prescribed_conditions, distributed_loads,
		irow_pt, irow_beam1, irow_beam2, icol_pt, icol_beam)

	b = system_residual!(b, x, assembly.elements, prescribed_conditions, distributed_loads,
		irow_pt, irow_beam1, irow_beam2, icol_pt, icol_beam)

	x = A\b

	return AssemblyState(x, prescribed_conditions, icol_pt, icol_beam, static)
end

function steady_state_analysis(assembly, prescribed_conditions, distributed_loads,
	x0, ω0, v0, niter)

	TF = eltype(assembly)

	static = false

	N, irow_pt, irow_beam1, irow_beam2, icol_pt, icol_beam = system_indices(
		assembly.points, assembly.elements, static)

	x = zeros(TF, N)
	A = spzeros(TF, N, N)
	b = zeros(TF, N)

	x, A, b = newton_solve!(x, A, b, N, irow_pt, irow_beam1, irow_beam2,
		icol_pt, icol_beam, prescribed_conditions, distributed_loads, x0, ω0, v0, niter)

	return AssemblyState(x, prescribed_conditions, icol_pt, icol_beam, static)
end

function eigenvalue_analysis(assembly, prescribed_conditions, distributed_loads, x0, ω0, v0, niter, nev)

	TF = eltype(assembly)

	x = zeros(TF, N)
	A = spzeros(TF, N, N)
	M = spzeros(TF, N, N)
	b = zeros(TF, N)

	N, irow_pt, irow_beam1, irow_beam2, icol_pt, icol_beam = system_indices(
		assembly.points, assembly.elements, static)

	x, A, b = newton_solve!(x, A, b, N, irow_pt, irow_beam1, irow_beam2,
		icol_pt, icol_beam, prescribed_conditions, distributed_loads, x0, ω0, v0, niter)

	M = mass_matrix!(M, N, irow_pt, irow_beam1, irow_beam2, icol_pt, icol_beam,
		prescribed_conditions, distributed_loads, x0, ω0, v0)

	λ, v = eigenvalue_solve!(x, A, M, nev)

	return λ, [AssemblyState(v[i], prescribed_conditions, icol_pt, icol_beam, static) for i = 1:length(v)]
end

function time_domain_analysis(assembly, prescribed_conditions, distributed_loads,
	u, θ, udot, θdot, x0, v0, ω0, dt, niter, nstep)

	TF = eltype(assembly)

	nbeam = length(icol_beam)

	static = false

	N, irow_pt, irow_beam1, irow_beam2, icol_pt, icol_beam = system_indices(
		assembly.points, assembly.elements, static)

	x = zeros(TF, N)
	A = spzeros(TF, N, N)
	M = spzeros(TF, N, N)
	b = zeros(TF, N)

	# solve for the initial conditions given u, udot, θ, θdot
	u = SVector{3, TF}.(u)
	θ = SVector{3, TF}.(θ)
	udot = SVector{3, TF}.(udot)
	θdot = SVector{3, TF}.(θdot)
	x, A, b = newton_solve!(x, A, b, N, irow_pt, irow_beam1, irow_beam2,
		icol_pt, icol_beam, prescribed_conditions, distributed_loads,
		u, θ, udot, θdot, x0, ω0, v0, niter)

	# set up for time-domain run
	udot = Vector{SVector{3,TF}}(undef, nbeam)
	θdot = Vector{SVector{3,TF}}(undef, nbeam)
	CtCabPdot = Vector{SVector{3, TF}}(undef, nbeam)
	CtCabHdot = Vector{SVector{3, TF}}(undef, nbeam)
	for i = 1:nbeam
		icol = icol_beam[i]
		# calculate udot
		u_p = SVector{3, TF}(u[ibeam])
		udot_p = SVector{3, TF}(udot[ibeam])
		udot[ibeam] = 2/dt*u_p - udot_p
		# calculate θdot
		θ_p = SVector{3, TF}(θ[ibeam])
		θdot_p = SVector{3, TF}(θdot[ibeam])
		θdot[ibeam] = 2/dt*θ_p - θdot_p
		# extract rotation parameters
		Ct = wiener_milenkovic(θ_p[i])'
		Cab = assembly.beams[ibeam].Cab
		# calculate CtCabPdot
		CtCabP_p = Ct*Cab*SVector{3, TF}(x[icol+6], x[icol+7], x[icol+8])
		CtCabPdot_p = SVector{3, TF}(x[icol], x[icol+1], x[icol+2])
		CtCabPdot[ibeam] = 2/dt*CtCabP_p - CtCabPdot_p
		# calculate CtCabHdot
		CtCabH_p = Ct*Cab*SVector{3, TF}(x[icol+9], x[icol+10], x[icol+11])
		CtCabHdot_p = SVector{3, TF}(x[icol+3], x[icol+4], x[icol+5])
		CtCabHdot[ibeam] = 2/dt*CtCabH_p - CtCabHdot_p
		# insert initial conditions for time-domain analysis
		x[icol:icol+2] .= u[i]
		x[icol+3:icol+5] .= θ[i]
	end

	# initialize storage for each time step
	history = Vector{AssemblyState{TF}}(undef, nstep)

	# loop for each time step
	for i = 1:nstep

		# solve for the state variables at next time step
		x, A, b = newton_solve!(x, A, b, N, irow_pt, irow_beam1, irow_beam2, icol_pt, icol_beam,
			prescribed_conditions, distributed_loads, x0, v0, ω0,
			udot, θdot, CtCabPdot, CtCabHdot, dt)

		# add state to history
		history[i] = AssemblyState(x, prescribed_conditions, icol_pt, icol_beam, static)

		# exit loop if done iterating
		if i == nstep
			break
		end

		# update time derivative terms for the next time step
		for ibeam = 1:nbeam
			icol = icol_beam[i]
			# calculate udot for next time step
			u_p = SVector{3, TF}(x[icol], x[icol+1], x[icol+2])
			udot_p = udot[ibeam]
			udot[ibeam] = 2/dt*u_p - udot_p
			# calculate θdot for next time step
			θ_p = SVector{3, TF}(x[icol+3], x[icol+4], x[icol+5])
			θdot_p = θdot[ibeam]
			θdot[ibeam] = 2/dt*θ_p - θdot_p
			# extract rotation parameters
			Ct = wiener_milenkovic(θ_p[i])'
			Cab = assembly.beams[ibeam].Cab
			# calculate CtCabPdot for next time step
			CtCabP_p = Ct*Cab*SVector{3, TF}(x[icol+6], x[icol+7], x[icol+8])
			CtCabPdot_p = CtCabPdot[ibeam]
			CtCabPdot[ibeam] = 2/dt*CtCabP_p - CtCabPdot_p
			# calculate CtCabHdot for next time step
			CtCabH_p = Ct*Cab*SVector{3, TF}(x[icol+9], x[icol+10], x[icol+11])
			CtCabHdot_p = CtCabPdot[ibeam]
			CtCabHdot[ibeam] = 2/dt*CtCabH_p - CtCabHdot_p
		end
	end

	return history
end

end
