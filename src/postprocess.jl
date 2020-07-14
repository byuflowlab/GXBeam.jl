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
	converged::Bool
end
Base.eltype(::AssemblyState{TF, TP, TE}) where {TF, TP, TE} = TF

function AssemblyState(converged, x, assembly, prescribed_conditions, distributed_loads,
	time_function_values, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt,
	icol_beam, args::Vararg{<:Any,N}) where N

	TF = eltype(x)

	static = N == 0

	npoint = length(assembly.points)
	nbeam = length(assembly.elements)

	points = Vector{PointState{TF}}(undef, npoint)
	for ipoint = 1:npoint

		# unknowns have been eliminated from the system of equations
		if icol_pt[ipoint] <= 0
			# solve for point variables
			ibeam = findfirst(x -> x == ipoint, assembly.start)
			# get indices
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
				f1, f2, m1, m2 = integrate_element_loads(ΔL, Ct, distributed_loads)
			else
				f1 = @SVector zeros(3)
				f2 = @SVector zeros(3)
				m1 = @SVector zeros(3)
				m2 = @SVector zeros(3)
			end
			# get the resultants from the element equations
			internal_resultants = element_equations(properties..., f1, f2, m1, m2)
			# extract resultants of interest
			f_F1 = internal_resultants[5]
			f_M1 = internal_resultants[7]
			# assign displacements and rotations of point
			u = f_F1
			theta = f_M1
			F = @SVector zeros(3)
			M = @SVector zeros(3)
		else
			prescribed = haskey(prescribed_conditions, ipoint)
			if prescribed
				u, theta, F, M = point_variables(x, icol_pt[ipoint],
					prescribed_conditions[ipoint], time_function_values)
			else
				u, theta, F, M = point_variables(x, icol_pt[ipoint])
			end
		end

		points[ipoint] = PointState{TF}(u, theta, F, M)
	end

	elements = Vector{ElementState{TF}}(undef, nbeam)
	for ibeam = 1:nbeam
		icol = icol_beam[ibeam]
		u = SVector{3, TF}(x[icol], x[icol+1], x[icol+2])
		theta = SVector{3, TF}(x[icol+3], x[icol+4], x[icol+5])
		F = SVector{3, TF}(x[icol+6], x[icol+7], x[icol+8])
		M = SVector{3, TF}(x[icol+9], x[icol+10], x[icol+11])
		P = ifelse(static, zero(u), SVector{3, TF}(x[icol+12], x[icol+13], x[icol+14]))
		H = ifelse(static, zero(u), SVector{3, TF}(x[icol+15], x[icol+16], x[icol+17]))
		elements[ibeam] = ElementState{TF}(u, theta, F, M, P, H)
	end

	return AssemblyState(points, elements, converged)
end
