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
	time_function_values, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam)

	TF = eltype(x)

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
			ΔL_b = beam.ΔL
			Ct_b = get_C(θ_b)'
			Cab_b = beam.Cab
			γ_b = element_strain(beam, F_b, M_b)
			κ_b = element_curvature(beam, F_b, M_b)

			# solve for the displacements/rotations of the point
			u =  side*u_b - ΔL_b/2*(Ct_b*Cab_b*(e1 + γ_b) - Cab_b*e1)
			theta =  side*θ_b - ΔL_b/2*get_Qinv(θ_b)*Cab_b*κ_b
			F = @SVector zeros(3)
			M = @SVector zeros(3)
		else
			# point variables are state variables, extract point state variables
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

	return AssemblyState(points, elements, converged)
end
