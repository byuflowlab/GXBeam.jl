"""
    Assembly{TF, TP<:AbstractVector{<:AbstractVector{TF}}, TC<:AbstractVector{<:Integer}, TE<:AbstractVector{Element{TF}}}

Composite type that defines an assembly of connected nonlinear beam elements for
analysis.

# Fields
 - `points`: Array of all beam element endpoints
 - `start`: Array containing point index where each beam element starts
 - `stop`: Array containing point index where each beam element stops
 - `element`: Array of `Element`s
"""
struct Assembly{TF, TP<:AbstractVector{<:AbstractVector{TF}}, TC<:AbstractVector{<:Integer}, TE<:AbstractVector{Element{TF}}}
	points::TP
	start::TC
	stop::TC
	elements::TE
end
Base.eltype(::Assembly{TF, TP, TC, TE}) where {TF, TP, TC, TE} = TF

"""
	Assembly(points, start, stop; kwargs...)

Construct an assembly of connected nonlinear beam elements for analysis.  Beam lengths
and midpoints may be manually specified in case beam elements are curved rather than
straight.

# Arguments
 - `points`: Array of all beam element endpoints
 - `start`: Array containing point index where each beam element starts
 - `stop`: Array containing point index where each beam element stops

# Keyword Arguments
 - `compliance = fill(Diagonal((@SMatrix zeros(6,6))))), length(start))`: Array of (6 x 6)
 	compliance matrices for each beam element,
 - `minv = fill((Diagonal((@SVector ones(6)))), length(start))`: Array of (6 x 6) inverse
 	mass matrix inverses for each beam element
 - `frames = fill(SMatrix{3,3}(I))`: Array of (3 x 3) direction cosine matrices for each beam element
 - `lengths = norm.(points[stop] - points[start])`: Array containing the length of each beam, defaults to the distance between beam endpoints
 - `midpoints = (points[stop] + points[start])/2`: Array containing the midpoint of each beam element, defaults to the average of the beam element endpoints
"""
function Assembly(points, start, stop;
	compliance = fill((@SMatrix zeros(6,6)), length(start)),
	minv = fill(Diagonal((@SVector ones(6))), length(start)),
	frames = fill(SMatrix{3,3}(I)),
	lengths = norm.(points[stop] - points[start]),
	midpoints = (points[stop] + points[start])/2)

	TF = promote_type(
		eltype(eltype(points)),
		eltype(eltype(compliance)),
		eltype(eltype(minv)),
		eltype(eltype(frames)),
		eltype(eltype(lengths)),
		eltype(eltype(midpoints))
		)

	elements = Element{TF}.(lengths, midpoints, compliance, minv, frames)

	return Assembly(points, promote(start, stop)..., elements)
end

"""
	discretize_beam(L, r, d; Cab = Matrix(I,3,3)), k = zeros(3))

Discretize a beam according to the discretization provided in `d`, an array that
ranges from 0 to 1, with 0 representing the beginning of the beam and 1 representing
the end of the beam.

If `d` is an integer, the beam is discretized into `d` uniformly spaced elements.

Return the lengths, endpoints, midpoints, and rotation matrices of the beam elements.

# Arguments
 - `L`: Beam length
 - `r`: Beam starting point
 - `d`: Discretization vector
 - `Cab`: 3x3 beam rotation matrix at the starting point
 - `k`: curvature vector
"""
function discretize_beam(L, r, ndiv::Integer; Cab = I3, k = (@SVector zeros(3)))

	d = range(0, 1, length=ndiv+1)

	return discretize_beam(L, r, d; Cab=Cab, k=k)
end

function discretize_beam(L, r, d; Cab = I3, k = (@SVector zeros(3)))

	# discretize beam
	sp = L*d
	sm = (sp[1:end-1] .+ sp[2:end])/2
	ΔL = (sp[2:end] .- sp[1:end-1])

	# precompute some curvature quantities
	kkt = k*k'
	ktilde = tilde(k)
	kn = sqrt(k'*k)

	if k == zero(k)
		triads = fill(Cab, length(sm))
		xp = [r + s*Cab*e1 for s in sp]
		xm = [r + s*Cab*e1 for s in sm]
	else
		triads = curve_triad.(Ref(Cab), Ref(kkt), Ref(ktilde), Ref(kn), sm)
		xp = curve_coordinates.(Ref(r), Ref(Cab), Ref(kkt), Ref(ktilde), Ref(kn), sp)
		xm = curve_coordinates.(Ref(r), Ref(Cab), Ref(kkt), Ref(ktilde), Ref(kn), sm)
	end

	return ΔL, xp, xm, triads
end


"""
	curve_length(r1, r2, k)

Calculate the length of a curve given its endpoints and its curvature vector
"""
function curve_length(r1, r2, k)
	kn = sqrt(k'*k)
	r = r2 - r1
	rn = sqrt(r'*r)
	if kn == 0 || (k[2] == 0 && k[3] == 0)
		ΔL = rn
	else
		if k[1] == 0
			ΔL = 2*asin(kn*rn/2)/kn
		else
			lower_bound = rn
			upper_bound = rn*kn/k
			residual = (L) -> (2*(kn^2-k[1]^2)*(1-cos(kn*L)) + k[1]^2*(kn*L)^2)/kn^4 - rnorm
			ΔL = fzero(residual, lower_bound, upper_bound)
		end
	end
	return ΔL
end

@inline curve_triad(Cab, k, s) = curve_triad(Cab, k*k', tilde(k), sqrt(k'*k), s)
@inline curve_triad(Cab, kkt, ktilde, kn, s) = SMatrix{3,3}(Cab*((I - kkt/kn^2)*cos(kn*s) +
	ktilde/kn*sin(kn*s) + kkt/kn^2))

@inline curve_coordinates(r, Cab, k, s) = curve_coordinates(r, Cab, k*k', tilde(k), sqrt(k'*k), s)
@inline curve_coordinates(r, Cab, kkt, ktilde, kn, s) = r + SVector{3}(Cab*((I/kn - kkt/kn^2)*sin(kn*s) + ktilde/kn^2*(1-cos(kn*s)) + kkt/kn^2*s)*e1)
