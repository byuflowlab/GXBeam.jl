"""
    Assembly{TF, TP<:AbstractVector{<:AbstractVector{TF}},
        TC<:AbstractVector{<:Integer}, TE<:AbstractVector{Element{TF}}}

Composite type that defines an assembly of connected nonlinear beam elements for
analysis.

# Fields
 - `points`: Array of all beam element endpoints
 - `start`: Array containing point index where each beam element starts
 - `stop`: Array containing point index where each beam element stops
 - `elements`: Array of `Element`s
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
 - `start`: Array containing point indices where each beam element starts
 - `stop`: Array containing point indices where each beam element stops

# Keyword Arguments
 - `mu`: Mass per unit length used for gravitational loads (if applicable). Will be derived 
        from the mass matrix or inverse mass matrix if not specified, or set to zero if 
        `mass` or `minv` are not provided.
 - `xm2`: Mass center location along axis 2 used for gravitational loads (if applicable). 
        Will be derived from the mass matrix or inverse mass matrix if not specified, 
        or set to zero if `mass` or `minv` are not provided.
 - `xm3`: Mass center location along axis 3 used for gravitational loads (if applicable). 
        Will be derived from the mass matrix or inverse mass matrix if not specified, 
        or set to zero if `mass` or `minv` are not provided.
 - `stiffness`: Array of (6 x 6) stiffness matrices for each beam element, 
        acts as an alternative to providing `compliance`
 - `compliance`: Array of (6 x 6) compliance matrices for each beam element, 
        defaults to `zeros(6,6)` for each beam element
 - `mass`: Array of (6 x 6) mass matrices for each beam element, 
        acts as an alternative to providing `minv`
 - `minv`: Array of (6 x 6) mass matrices for each beam element, 
        defaults to the identity matrix for each beam element
 - `frames`: Array of (3 x 3) tranformation matrices for each beam element.
        Transforms from the local undeformed beam frame to the global frame) and defaults
        to the identity matrix for each beam element
 - `lengths`: Array containing the length of each beam, defaults to the distance between 
        beam endpoints
 - `midpoints`: Array containing the midpoint of each beam element, defaults to the average 
        of the beam element endpoints
"""
function Assembly(points, start, stop;
    mu = nothing, 
    xm2 = nothing,
    xm3 = nothing,
    stiffness = nothing,
    compliance = nothing,
    mass = nothing,
    minv = nothing,
    frames = nothing,
    lengths = norm.(points[stop] - points[start]),
    midpoints = (points[stop] + points[start])/2)

    nelem = length(start)

    if isnothing(mu)
        if isnothing(mass) && isnothing(minv)
            μ = zeros(nelem)
        else
            if isnothing(mass)
                mass = inv.(minv)
            end
            μ = getindex.(mass, 1, 1)
        end
    end

    if isnothing(xm2)
        if isnothing(mass) && isnothing(minv)
            xm2 = zeros(nelem)
        else 
            if isnothing(mass)
                xm2 = -getindex.(minv,2,3)./getindex.(minv,2,4)
            else
                xm2 = -getindex.(mass,1,6)./getindex.(mass,1,1)
            end 
        end
    end

    if isnothing(xm3)
        if isnothing(mass) && isnothing(minv)
            xm3 = zeros(nelem)
        else 
            if isnothing(mass)
                xm3 = getindex.(minv,2,3)./getindex.(minv,4,3)
            else
                xm3 = getindex.(mass,1,5)./getindex.(mass,1,1)
            end 
        end
    end

    if isnothing(compliance)
        if isnothing(stiffness)
            compliance = fill((@SMatrix zeros(6,6)), nelem)
        else
            compliance = [(@MMatrix zeros(eltype(eltype(stiffness)), 6,6)) for i=1:nelem] #can't use fill because it copies the reference. Need a different value for every i.
            for i = 1:nelem
                filled_cols = findall(vec(mapslices(col -> any(row -> !isapprox(row, 0), col), stiffness[i], dims = 1)))
                compliance[i][filled_cols,filled_cols] .= inv(Matrix(stiffness[i][filled_cols, filled_cols]))
            end
            compliance = SMatrix.(compliance)
        end
    end

    if isnothing(minv)
        if isnothing(mass)
            minv = fill(Diagonal(@SVector ones(6)), nelem)
        else
            minv = [MMatrix{6,6}(Diagonal(@SVector ones(eltype(eltype(mass)), 6))) for i=1:nelem] # can't use infill because it copies the reference. Need a different value for every i.
            for i = 1:nelem
                filled_cols = findall(vec(mapslices(col -> any(row -> !isapprox(row, 0), col), mass[i], dims = 1)))
                minv[i][filled_cols,filled_cols] .= inv(Matrix(mass[i][filled_cols, filled_cols]))
            end
            minv = SMatrix{6,6}.(minv)
        end
    end

    if isnothing(frames)
        frames = fill(I3, nelem)
    end

    elements = Element.(lengths, midpoints, μ, xm2, xm3, compliance, minv, frames)

    return Assembly(points, promote(start, stop)..., elements)
end

"""
    discretize_beam(L, start, discretization; frame = Matrix(I,3,3)),
        curvature = zeros(3))

Discretize a beam according to the discretization provided in `discretization`
given the beam length (`L`), and starting point (`start`).

Return the lengths, endpoints, midpoints, and trasformation matrices of the beam elements.

# Arguments
 - `L`: Beam length
 - `start`: Beam starting point
 - `discretization`: May be either an integer, representing the number of
        elements that the beam should be discretized into, or a vector containing
        the normalized endpoints of each beam element, where 0 is the beginning
        of the beam and 1 is the end of the beam.
 - `frame`: 3x3 tranformation matrix which transforms from the local beam
        coordinate frame at the start of the beam to the global coordinate frame.
 - `curvature`: curvature vector
"""
discretize_beam(L, start, discretization::Integer; frame = I3, curvature = (@SVector zeros(3))) =
    discretize_beam(L, start, range(0, 1, length=discretization+1); frame=frame, curvature=curvature)

function discretize_beam(L, start, discretization; frame = I3, curvature = (@SVector zeros(3)))

    r1 = SVector{3}(start)
    Cab = SMatrix{3,3}(frame)
    k = SVector{3}(curvature)

    # discretize beam
    sp = L*discretization
    sm = (sp[1:end-1] .+ sp[2:end])/2
    ΔL = (sp[2:end] .- sp[1:end-1])

    # precompute some curvature quantities
    kkt = k*k'
    ktilde = tilde(k)
    kn = sqrt(k'*k)

    if curvature == zero(curvature)
        triads = fill(Cab, length(sm))
        xp = [r1 + s*Cab*e1 for s in sp]
        xm = [r1 + s*Cab*e1 for s in sm]
    else
        triads = curve_triad.(Ref(Cab), Ref(kkt), Ref(ktilde), Ref(kn), sm)
        xp = curve_coordinates.(Ref(r1), Ref(Cab), Ref(kkt), Ref(ktilde), Ref(kn), sp)
        xm = curve_coordinates.(Ref(r1), Ref(Cab), Ref(kkt), Ref(ktilde), Ref(kn), sm)
    end

    return ΔL, xp, xm, triads
end


"""
    curve_length(start, stop, curvature)

Calculate the length of a curve given its endpoints and its curvature vector
"""
function curve_length(start, stop, curvature)
    r1, r2, k = SVector{3}(start), SVector{3}(stop), SVector{3}(curvature)
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

"""
    curve_triad(Cab, k, s)
    curve_triad(Cab, kkt, ktilde, kn, s)

Return the transformation matrix at `s` along the length of the beam given
the curvature vector `k` and the initial transformation matrix `Cab`.
"""
@inline curve_triad(Cab, k, s) = curve_triad(Cab, k*k', tilde(k), sqrt(k'*k), s)
@inline curve_triad(Cab, kkt, ktilde, kn, s) = SMatrix{3,3}(Cab*((I - kkt/kn^2)*cos(kn*s) +
    ktilde/kn*sin(kn*s) + kkt/kn^2))

"""
    curve_coordiantes(r, Cab, k, s)
    curve_coordinates(r, Cab, kkt, ktilde, kn, s)

Return the coordinates at `s` along the length of the beam given the starting
point `r`, initial transformation matrix `Cab`, and curvature vector `k`.
"""
@inline curve_coordinates(r, Cab, k, s) = curve_coordinates(r, Cab, k*k', tilde(k), sqrt(k'*k), s)
@inline curve_coordinates(r, Cab, kkt, ktilde, kn, s) = r + SVector{3}(Cab*((I/kn - kkt/kn^2)*sin(kn*s) + ktilde/kn^2*(1-cos(kn*s)) + kkt/kn^2*s)*e1)
