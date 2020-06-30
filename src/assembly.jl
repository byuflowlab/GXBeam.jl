struct Assembly{TF, TP<:AbstractVector{<:AbstractVector{TF}}, TE<:AbstractVector{<:AbstractElement{TF}}}
	points::TP
	elements::TE
end
Base.eltype(::Assembly{TF, TP, TE}) = TF

"""
	Assembly(points, connectivity, compliance, mass, frames, [lengths, midpoints])

Construct an assembly of connected nonlinear beam elements for analysis.  Beam lengths
and midpoints may be manually specified in case beam elements are curved rather than
straight.

# Arguments
 - points: Array of all beam element endpoints
 - pt1: Array containing index of start of each beam element
 - pt2: Array containing index of end of each beam element
 - compliance: Array of (6 x 6) compliance matrices at the midpoint of each beam element
 - mass: Array of (6 x 6) mass matrices at the midpoint of each beam element
 - frames: Array of (3 x 3) direction cosine matrices at the midpoint of each beam
 - lengths: Array containing length of each beam, defaults to distance between beam endpoints
 - midpoints: Array containing the midpoints of each beam element, defaults to average of the beam element endpoints
"""
Assembly

function Assembly(points, pt1, pt2, compliance, mass, frames, lengths, midpoints)

	elements = Element.(lengths, midpoints, pt1, pt2, compliance, mass, frames)

	return Assembly(points, elements)
end

function Assembly(points, pt1, pt2, compliance, mass, frames)

	lengths = norm.(points[pt2] - points[pt1])
	midpoints = (points[pt2] + points[pt1])/2

	return Assembly(points, connectivity, compliance, mass, frames, lengths, midpoints)
end

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
Base.eltype(::AssemblyState{TF, TP, TE}) = TF

function AssemblyState(x, prescribed_conditions, icol_pt, icol_beam, static)

	TF = eltype(x)
	npt = length(icol_pt)
	nbeam = length(icol_beam)

	points = AbstractVector{PointState{TF}}(undef, npt)
	for ipt = 1:npt
		u, theta, F, M = point_variables(x, icol_pt[ipt], prescribed_conditions[ipt])
		points[ipt] = PointState{TF}(u, theta, F, M)
	end

	elements = AbstractVector{ElementState{TF}}(undef, nbeam)
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

	return AssemblyState(points, elements)
end
