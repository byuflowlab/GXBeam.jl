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
 - `compliance = fill((@SMatrix zeros(6,6)), length(start))`: Array of (6 x 6)
 	compliance matrices for each beam element,
 - `mass = fill((@SMatrix zeros(6,6)), length(start))`: Array of (6 x 6) inverse
 	mass matrices for each beam element
 - `frames = fill(SMatrix{3,3}(I))`: Array of (3 x 3) direction cosine matrices for each beam element
 - `lengths = norm.(points[stop] - points[start])`: Array containing the length of each beam, defaults to the distance between beam endpoints
 - `midpoints = (points[stop] + points[start])/2`: Array containing the midpoint of each beam element, defaults to the average of the beam element endpoints
"""
function Assembly(points, start, stop;
	compliance = fill((@SMatrix zeros(6,6)), length(start)),
	minv = fill((@SMatrix zeros(6,6)), length(start)),
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
