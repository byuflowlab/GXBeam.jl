struct Assembly{TF}
	points::Array{TF, 2}
	elements::Array{Element{TF}, 1}
end

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

# function curve_length(r1, r2, k)
# 	kn = sqrt(k'*k)
# 	r = r2 - r1
# 	rn = sqrt(r'*r)
# 	if kn == 0 || (k[2] == 0 && k[3] == 0)
# 		ΔL = rn
# 	else
# 		if k[1] == 0
# 			ΔL = 2*asin(kn*rn/2)/kn
# 		else
# 			lower_bound = rn
# 			upper_bound = rn*kn/k
# 			residual = (L) -> (2*(kn^2-k[1]^2)*(1-cos(kn*L)) + k[1]^2*(kn*L)^2)/kn^4 - rnorm
# 			ΔL = fzero(residual, lower_bound, upper_bound)
# 		end
# 	return ΔL
# end
