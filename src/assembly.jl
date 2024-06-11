"""
    Element{TF}

Composite type that defines a beam element's properties

# Fields
 - `L`: Beam element length
 - `x`: Beam element location
 - `compliance`: Beam element compliance matrix
 - `mass`: Beam element mass matrix
 - `Cab`: Transformation matrix from the undeformed beam element frame to the body frame
 - `mu`: Beam element damping coefficients
"""
struct Element{TF}
    L::TF
    x::SVector{3,TF}
    compliance::SMatrix{6,6,TF,36}
    mass::SMatrix{6,6,TF,36}
    Cab::SMatrix{3,3,TF,9}
    mu::SVector{6,TF}
end

Base.eltype(::Element{TF}) where TF = TF
Base.eltype(::Type{Element{TF}}) where TF = TF

Element{TF}(e::Element) where {TF} = Element{TF}(e.L, e.x, e.compliance, e.mass, e.Cab, e.mu)
Base.convert(::Type{Element{TF}}, e::Element) where {TF} = Element{TF}(e)

"""
    Element(L, x, compliance, mass, Cab, mu)

Construct a beam element.

# Arguments
- `L`: Length of the beam element
- `x`: Location of the beam element (the center of the beam element)
- `compliance`: Beam element compliance matrix
- `mass`: Beam element mass matrix
- `Cab`: Transformation matrix from the undeformed beam element frame to the body frame
- `mu`: Beam element damping coefficients
"""
function Element(L, x, compliance, mass, Cab, mu)
    TF = promote_type(typeof(L), eltype(x), eltype(compliance), eltype(mass), eltype(Cab), eltype(mu))
    return Element{TF}(L, x, compliance, mass, Cab, mu)
end

"""
    Element(points, start, stop; kwargs...)

Construct a beam element.  Beam lengths and midpoints may be manually specified in case 
the beam element is curved rather than straight.

# Arguments
 - `points`: Array of all beam element endpoints in the assembly
 - `start`: Point index where this beam element starts
 - `stop`: Point index where this beam element stops

# Keyword Arguments
 - `frame`: 3 x 3 tranformation matrix for this beam element.  Transforms from the local 
        undeformed beam frame to the global frame. Defaults to the identity matrix.
 - `compliance`: 6 x 6 compliance matrix for this beam element. Defaults to `zeros(6,6)`. 
 - `mass`: 6 x 6 mass matrix for this beam element. Defaults to `zeros(6,6)`
 - `damping`: structural damping coefficients for this beam element.
        Defaults to `fill(0.01, 6)`.
 - `length`: Beam length, defaults to the distance between beam endpoints
 - `midpoint`: Beam midpoint, defaults to the average of the beam element endpoints
"""
function Element(points, start, stop;
    frame = nothing,
    stiffness = nothing,
    compliance = nothing,
    mass = nothing,
    damping = nothing,
    length = nothing,
    midpoint = nothing)

    if isnothing(length)
        length = norm(points[stop] - points[start])
    end

    if isnothing(midpoint)
        midpoint = (points[stop] + points[start])/2
    end

    if isnothing(compliance)
        if isnothing(stiffness)
            compliance = @SMatrix zeros(6,6)
        else
            compliance = @MMatrix zeros(eltype(eltype(stiffness)), 6, 6)
            filled_cols = findall(vec(mapslices(col -> any(row -> !isapprox(row, 0), col), stiffness, dims = 1)))
            compliance[filled_cols,filled_cols] .= inv(Matrix(stiffness[filled_cols, filled_cols]))
            compliance = SMatrix(compliance)
        end
    end

    if isnothing(mass)
        mass = @SMatrix zeros(6,6)
    end

    if isnothing(frame)
        frame = I3
    end

    if isnothing(damping)
        damping = 0.01*(@SVector ones(6))
    end

    return Element(length, midpoint, compliance, mass, frame, damping)
end

# this definition is for broadcasting
function Element(points, start, stop, frame, stiffness, compliance, mass, damping, length, midpoint)
    return Element(points, start, stop; frame, stiffness, compliance, mass, damping, length, midpoint)
end

"""
    Assembly{TF, TP<:AbstractVector{<:AbstractVector{TF}},
        TC<:AbstractVector{<:Integer}, TE<:AbstractVector{Element{TF}}}

Composite type that defines an assembly of connected nonlinear beam elements.

# Fields
 - `points`: Array of all beam element endpoints
 - `start`: Array containing point index where each beam element starts
 - `stop`: Array containing point index where each beam element stops
 - `elements`: Array containing beam element definitions (see [`Element`](@ref))
"""
struct Assembly{TF, TP<:AbstractVector{<:AbstractVector{TF}}, TC<:AbstractVector{<:Integer}, TE<:AbstractVector{Element{TF}}}
    points::TP
    start::TC
    stop::TC
    elements::TE
end

Base.eltype(::Assembly{TF, TP, TC, TE}) where {TF, TP, TC, TE} = TF
Base.eltype(::Type{Assembly{TF, TP, TC, TE}}) where {TF, TP, TC, TE} = TF

Assembly{TF, TP, TC, TE}(a::Assembly) where {TF, TP, TC, TE} = Assembly{TF, TP, TC, TE}(a.points, a.start, a.stop, a.elements)
Base.convert(::Type{Assembly{TF, TP, TC, TE}}, a::Assembly) where {TF, TP, TC, TE} = Assembly{TF, TP, TC, TE}(a)

"""
    Assembly{TF}(points, start, stop, elements)

Construct an assembly of connected nonlinear beam elements.  

# Arguments
 - `points`: Array of all beam element endpoints
 - `start`: Array containing point indices where each beam element starts
 - `stop`: Array containing point indices where each beam element stops
 - `elements`: Array containing beam element definitions (see [`Element`](@ref))
"""
function Assembly(points, start, stop, elements)
    TF = promote_type(eltype(eltype(points)), eltype(eltype(elements)))
    return Assembly(SVector{3,TF}.(points), promote(start, stop)..., Element{TF}.(elements))
end

"""
    Assembly(points, start, stop; kwargs...)

Construct an assembly of connected nonlinear beam elements.  Beam lengths and midpoints 
may be manually specified in case beam elements are curved rather than straight.

# Arguments
 - `points`: Array of all beam element endpoints
 - `start`: Array containing point indices where each beam element starts
 - `stop`: Array containing point indices where each beam element stops

# Keyword Arguments
 - `frames`: Array of (3 x 3) tranformation matrices for each beam element.
        Transforms from the local undeformed beam frame to the global frame. 
        Defaults to the identity matrix.
 - `compliance`: Array of (6 x 6) compliance matrices for each beam element.
        Defaults to `zeros(6,6)` for each beam element.
 - `mass`: Array of (6 x 6) mass matrices for each beam element. 
        Defaults to `zeros(6,6)` for each beam element.
 - `damping`: Array of (6) structural damping coefficients for each beam element. 
        Defaults to `fill(0.01, 6)` for each beam element.
 - `lengths`: Array containing the length of each beam element.
        Defaults to the distance between beam endpoints.
 - `midpoints`: Array containing the midpoint of each beam element.
        Defaults to the average of the beam element endpoints.
"""
function Assembly(points, start, stop;
    frames = nothing,
    stiffness = nothing,
    compliance = nothing,
    mass = nothing,
    damping = nothing,
    lengths = nothing,
    midpoints = nothing)

    elements = Element.(Ref(points), start, stop, frames, stiffness, compliance, mass, damping, lengths, midpoints)

    return Assembly(SVector{3}.(points), promote(start, stop)..., elements)
end

"""
    discretize_beam(L, start, discretization; frame, curvature)

Discretize a beam of length `L` located at `start` according to the discretization provided 
in `discretization`

Return the lengths, endpoints, midpoints, and reference frame of the beam elements.

# Arguments
 - `L`: Beam length
 - `start`: Beam starting point
 - `discretization`: Number of beam elements, or the normalized endpoints of each beam 
        element, with values ranging from 0 to 1.

# Keyword Arguments
 - `frame`: Reference frame at the start of the beam element, represented by a 3x3
         transformation matrix from the undeformed local frame to the body frame.
 - `curvature`: curvature vector
"""
function discretize_beam(L::Number, start::AbstractVector, discretization::Integer; 
    frame=I3, curvature=(@SVector zeros(3)))
  
    return discretize_beam(L, start, range(0, 1, length=discretization+1); 
        frame=frame, curvature=curvature)
end

function discretize_beam(L::Number, start::AbstractVector, discretization; 
    frame=I3, curvature=(@SVector zeros(3)))

    r1 = SVector{3}(start)
    Cab = SMatrix{3,3}(frame)
    k = SVector{3}(curvature)

    # discretize beam
    sp = L*discretization
    sm = (sp[1:end-1] .+ sp[2:end])/2
    ΔL = sp[2:end] .- sp[1:end-1]

    # precompute some curvature quantities
    kkt = k*k'
    ktilde = tilde(k)
    kn = sqrt(k'*k)

    if iszero(curvature)
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
    discretize_beam(start, stop, discretization; frame, curvature)

Discretize a beam from `start` to `stop` according to the discretization provided in 
`discretization`.

Return the lengths, endpoints, midpoints, and reference frame of the beam elements.

# Arguments
 - `start`: Beam starting point
 - `stop`: Beam ending point
 - `discretization`: Number of beam elements, or the normalized endpoints of each beam 
        element, with values ranging from 0 to 1.

# Keyword Arguments
 - `frame`: Reference frame at the start of the beam element, represented by a 3x3
        transformation matrix from the undeformed local frame to the body frame.
 - `curvature`: curvature vector
"""
function discretize_beam(start::AbstractVector, stop::AbstractVector, 
    discretization::Integer; frame=I3, curvature=(@SVector zeros(3)))
    
    return discretize_beam(start, stop, range(0, 1, length=discretization+1); 
        frame=frame, curvature=curvature)
end

function discretize_beam(start::AbstractVector, stop::AbstractVector, discretization; 
    frame=I3, curvature=(@SVector zeros(3)))

    # compute curve length
    L = curve_length(start, stop, curvature)

    # discretize beam
    ΔL, xp, xm, triads = discretize_beam(L, start, discretization; frame=frame, curvature=curvature)

    return ΔL, xp, xm, triads
end

"""
    curve_length(start, stop, curvature)

Calculate the length of a curve given its endpoints and its curvature vector.
"""
function curve_length(start, stop, curvature)

    r1 = SVector{3}(start)
    r2 = SVector{3}(stop)
    k = SVector{3}(curvature)

    # precomputed quantities
    kn = sqrt(k'*k)
    r = r2 - r1
    rn2 = r'*r
    rn = sqrt(rn2)

    # compute curve length
    if iszero(kn) || (iszero(k[2]) && iszero(k[3]))
        # the beam is straight or twisted, but not curved
        ΔL = rn
    elseif iszero(k[1])
        # the beam is curved, but not twisted
        ΔL = 2*asin(kn*rn/2)/kn
    else 
        # the beam is curved and twisted
        f = (L) -> (2*(kn^2-k[1]^2)*(1-cos(kn*L)) + k[1]^2*(kn*L)^2)/kn^4 - rn2
        ΔL = Roots.find_zero(f, (rn, rn*kn/k), Roots.Brent())
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
@inline curve_coordinates(r, Cab, kkt, ktilde, kn, s) = r + SVector{3}(Cab*((I/kn - kkt/kn^3)*sin(kn*s) + ktilde/kn^2*(1-cos(kn*s)) + kkt/kn^2*s)*e1)

function total_mass(elements::Vector{Element{TF}}) where TF
    mass = 0.0
    for elem in elements
        mass += elem.mass[1]*elem.L
    end
    return mass
end

function total_mass(assembly::Assembly)
    return total_mass(assembly.elements)
end