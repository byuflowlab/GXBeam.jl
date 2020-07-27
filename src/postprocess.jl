struct PointState{TF}
	u::SVector{3, TF}
	theta::SVector{3, TF}
	F::SVector{3, TF}
	M::SVector{3, TF}
end

struct ElementState{TF}
	u::SVector{3, TF}
	theta::SVector{3, TF}
	F::SVector{3, TF}
	M::SVector{3, TF}
	P::SVector{3, TF}
	H::SVector{3, TF}
end

struct AssemblyState{TF, TP<:AbstractVector{PointState{TF}}, TE<:AbstractVector{ElementState{TF}}}
	points::TP
	elements::TE
end
Base.eltype(::AssemblyState{TF, TP, TE}) where {TF, TP, TE} = TF

"""
    AssemblyState(system, assembly, x = system.x;
		prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}())

Post-process the system state given the solution vector `x`.  Return an object
of type `AssemblyState` that defines the state of the assembly for the time step.

If `prescribed_conditions` is not provided, all point state variables are assumed
to be displacements/rotations, rather than their actual identities as used in the
analysis.
"""
function AssemblyState(system, assembly, x = system.x;
	prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}())

	irow_pt = system.irow_pt
	irow_beam = system.irow_beam
	irow_beam1 = system.irow_beam1
	irow_beam2 = system.irow_beam2
	icol_pt = system.icol_pt
	icol_beam = system.icol_beam


	TF = eltype(x)

	istep = system.current_step[]

	npoint = length(assembly.points)
	nbeam = length(assembly.elements)

	points = Vector{PointState{TF}}(undef, npoint)
	for ipoint = 1:npoint

		if icol_pt[ipoint] <= 0
			# point variables are not state variables, solve for point variables

			# find the first beam connected to the point, checking both sides
			ibeam = findfirst(x -> x == ipoint, assembly.start)
			if !isnothing(ibeam)
				side = -1
			else
				ibeam = findfirst(x -> x == ipoint, assembly.stop)
				side = 1
			end

			beam = assembly.elements[ibeam]

			# get beam state variables
			icol_b = icol_beam[ibeam]
			u_b = SVector(x[icol_b   ], x[icol_b+1 ], x[icol_b+2 ])
			θ_b = SVector(x[icol_b+3 ], x[icol_b+4 ], x[icol_b+5 ])
			F_b = SVector(x[icol_b+6 ], x[icol_b+7 ], x[icol_b+8 ]) .* FORCE_SCALING
			M_b = SVector(x[icol_b+9 ], x[icol_b+10], x[icol_b+11]) .* FORCE_SCALING

			# get beam element properties
			ΔL_b = beam.L
			Ct_b = get_C(θ_b)'
			Cab_b = beam.Cab
			γ_b = element_strain(beam, F_b, M_b)
			κ_b = element_curvature(beam, F_b, M_b)

			# solve for the displacements/rotations of the point
			u = u_b + side*ΔL_b/2*(Ct_b*Cab_b*(e1 + γ_b) - Cab_b*e1)
			theta = θ_b + side*ΔL_b/2*get_Qinv(θ_b)*Cab_b*κ_b
			F = @SVector zeros(3)
			M = @SVector zeros(3)
		else
			# point variables are state variables, extract point state variables
			prescribed = haskey(prescribed_conditions, ipoint)
			if prescribed
				u, theta, F, M = point_variables(x, icol_pt[ipoint],
					prescribed_conditions[ipoint], istep)
			else
				u, theta, F, M = point_variables(x, icol_pt[ipoint])
			end
		end

		points[ipoint] = PointState{TF}(u, theta, F, M)
	end

	# all element variables are state variables
	elements = Vector{ElementState{TF}}(undef, nbeam)
	for ibeam = 1:nbeam
		icol = icol_beam[ibeam]
		static = irow_beam[ibeam] <= 0
		u = SVector{3, TF}(x[icol], x[icol+1], x[icol+2])
		theta = SVector{3, TF}(x[icol+3], x[icol+4], x[icol+5])
		F = SVector{3, TF}(x[icol+6], x[icol+7], x[icol+8]) .* FORCE_SCALING
		M = SVector{3, TF}(x[icol+9], x[icol+10], x[icol+11]) .* FORCE_SCALING
		P = ifelse(static, zero(u), SVector{3, TF}(x[icol+12], x[icol+13], x[icol+14]))
		H = ifelse(static, zero(u), SVector{3, TF}(x[icol+15], x[icol+16], x[icol+17]))
		elements[ibeam] = ElementState{TF}(u, theta, F, M, P, H)
	end

	return AssemblyState(points, elements)
end

"""
    write_vtk(name, assembly::Assembly, [state::AssemblyState, scaling=1.0])

Write the deformed geometry (and associated data) to a VTK file for visualization
using ParaView.  The deflections can be scaled with `scaling`.

The `state` and `scaling` parameters may be omitted to write the original geometry
to a VTK file without any associated data.
"""
function write_vtk(name, assembly::Assembly)

	# get problem dimensions
	npoint = length(assembly.points)
	nbeam = length(assembly.elements)

	# extract point locations
	points = Matrix{eltype(assembly)}(undef, 3, npoint)
	for ipoint = 1:npoint
		for i = 1:3
			points[i,ipoint] = assembly.points[ipoint][i]
		end
	end

	# create cells
	cells = [MeshCell(PolyData.Lines(), [assembly.start[i], assembly.stop[i]]) for i = 1:nbeam]

	# write vtk file
	vtk_grid(name, points, cells) do vtkfile

		# no associated data to add

	end

	return nothing
end

function write_vtk(name, assembly::Assembly, state::AssemblyState, scaling=1.0)

	# get problem dimensions
	npoint = length(assembly.points)
	nbeam = length(assembly.elements)

	# extract point locations
	points = Matrix{eltype(assembly)}(undef, 3, npoint)
	for ipoint = 1:npoint
		for i = 1:3
			points[i,ipoint] = assembly.points[ipoint][i] + scaling*state.points[ipoint].u[i]
		end
	end

	global points

	# create cells
	cells = [MeshCell(PolyData.Lines(), [assembly.start[i], assembly.stop[i]]) for i = 1:nbeam]

	# write vtk file
	vtk_grid(name, points, cells) do vtkfile

		# add point data
		for field in fieldnames(PointState)
			data = Matrix{eltype(assembly)}(undef, 3, npoint)
			for ipoint = 1:npoint
				data[:, ipoint] .= getproperty(state.points[ipoint], field)
			end
			vtkfile[string(field), VTKPointData()] = data
		end

		# add cell data
		for field in fieldnames(ElementState)
			data = Matrix{eltype(assembly)}(undef, 3, nbeam)
			for ibeam = 1:nbeam
				data[:, ibeam] .= getproperty(state.elements[ibeam], field)
			end
			vtkfile[string(field), VTKCellData()] = data
		end

	end

	return nothing
end

"""
	left_eigenvectors(system, λ, V)
	left_eigenvectors(K, M, λ, V)

Compute the left eigenvector matrix `Us` for the `system` using inverse power
iteration given the eigenvalues `λ` and the corresponding right eigenvector
matrix `V`.

The complex conjugate of each left eigenvector is stored in each row of the
matrix `Us`

Left and right eigenvectors satisfy the following M-orthogonality condition:
 - u'*M*v = 1 if u and v correspond to the same eigenvalue
 - u'*M*v = 0 if u and v correspond to different eigenvalues

This function assumes that `system` has not been modified since the eigenvalues
and right eigenvectors were computed.
"""
left_eigenvectors(system, λ, V) = left_eigenvectors(system.K, system.M, λ, V)

function left_eigenvectors(K, M, λ, V)

	# problem type and dimensions
	TC = eltype(V)
	nx = size(V,1)
	nev = size(V,2)

	# allocate storage
	Us = rand(TC, nev, nx)
	u = Vector{TC}(undef, nx)
	tmp = Vector{TC}(undef, nx)

	# get entries in M
	iM, jM, valM = findnz(M)

	# compute eigenvectors for each eigenvalue
	for iλ = 1:nev

		# factorize (K + λ*M)'
		KmλMfact = factorize(K' - λ[iλ]'*M')

		# initialize left eigenvector
		for i = 1:nx
			u[i] = Us[iλ,i]
		end

		# perform a few iterations to converge the left eigenvector
		for ipass = 1:3
			# get updated u
			mul!(tmp, M, u)
			ldiv!(u, KmλMfact, tmp)
			# normalize u
			unorm = zero(TC)
			for k = 1:length(valM)
				unorm += conj(u[iM[k]])*valM[k]*V[jM[k],iλ]
			end
			rdiv!(u, conj(unorm))
		end

		# store conjugate of final eigenvector
		for i = 1:nx
			Us[iλ,i] = conj(u[i])
		end

	end

	return Us
end
