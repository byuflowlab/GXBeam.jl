"""
    PrescribedConditions{T}

Type which defines the prescribed displacements and loads at a point.

# Fields:
 - `pd`: Flag for each degree of freedom indicating whether displacements are prescribed
 - `pl`: Flag for each degree of freedom indicating whether loads are prescribed
 - `u`: Linear displacement
 - `theta`: Angular displacement
 - `F`: External forces
 - `M`: External moments
 - `Ff`: Follower forces
 - `Mf`: Follower moments
"""
struct PrescribedConditions{T}
    pd::SVector{6, Bool}
    pl::SVector{6, Bool}
    u::SVector{3, T}
    theta::SVector{3, T}
    F::SVector{3, T}
    M::SVector{3, T}
    Ff::SVector{3, T}
    Mf::SVector{3, T}
end
Base.eltype(::PrescribedConditions{T}) where T = T

function PrescribedConditions{T}(p::PrescribedConditions) where T
    PrescribedConditions{T}(p.pd, p.pl, p.u, p.theta, p.F, p.M, p.Ff, p.Mf)
end

PrescribedConditions(pd, pl, u, theta, F, M, Ff, Mf) = PrescribedConditions(
    promote(pd, pl)..., promote(u, theta, F, M, Ff, Mf)...)

"""
    PrescribedConditions(; kwargs...)

Define the prescribed conditions at a point.  Individual prescribed conditions 
may be assigned as either a scalar parameter or as a function of time.

Prescribed Wiener-Milenkovic parameters must satisfy the following inequality:
sqrt(theta_x^2 + theta_y^2 + theta_z^2) <= 4.  Note that this restriction still
allows all possible rotations to be represented.

Note that if displacements and loads corresponding to the same degree of freedom are 
prescribed at the same point, the global body-fixed acceleration corresponding to the 
same degree of freedom will be modified to attempt to satisfy both conditions.

# Keyword Arguments
 - `ux`: Prescribed x-displacement (in the body frame)
 - `uy`: Prescribed y-displacement (in the body frame)
 - `uz`: Prescribed z-displacement (in the body frame)
 - `theta_x`: Prescribed first Wiener-Milenkovic parameter
 - `theta_y`: Prescribed second Wiener-Milenkovic parameter
 - `theta_z`: Prescribed third Wiener-Milenkovic parameter
 - `Fx`: Prescribed x-direction force
 - `Fy`: Prescribed y-direction force
 - `Fz`: Prescribed z-direction force
 - `Mx`: Prescribed x-direction moment
 - `My`: Prescribed y-direction moment
 - `Mz`: Prescribed z-direction moment
 - `Fx_follower`: Prescribed x-direction follower force
 - `Fy_follower`: Prescribed y-direction follower force
 - `Fz_follower`: Prescribed z-direction follower force
 - `Mx_follower`: Prescribed x-direction follower moment
 - `My_follower`: Prescribed y-direction follower moment
 - `Mz_follower`: Prescribed z-direction follower moment
"""
function PrescribedConditions(;
    ux = nothing,
    uy = nothing,
    uz = nothing,
    theta_x = nothing,
    theta_y = nothing,
    theta_z = nothing,
    Fx = nothing,
    Fy = nothing,
    Fz = nothing,
    Mx = nothing,
    My = nothing,
    Mz = nothing,
    Fx_follower = nothing,
    Fy_follower = nothing,
    Fz_follower = nothing,
    Mx_follower = nothing,
    My_follower = nothing,
    Mz_follower = nothing)

    # first degree of freedom
    if isnothing(ux)
        # prescribed load, displacement state variable
        pd1 = false
        pl1 = true
        ux = NaN
        Fx = isnothing(Fx) ? 0.0 : Fx
        Fx_follower = isnothing(Fx_follower) ? 0.0 : Fx_follower
    elseif isnothing(Fx) && isnothing(Fx_follower)
        # prescribed displacement, load state variable
        pd1 = true
        pl1 = false
        ux = isnothing(ux) ? 0.0 : ux
        Fx = NaN
        Fx_follower = 0.0
    else
        # prescribed displacement and load, body-fixed acceleration state variable
        pd1 = true
        pl1 = true
        ux = isnothing(ux) ? 0.0 : ux
        Fx = isnothing(Fx) ? 0.0 : Fx
        Fx_follower = isnothing(Fx_follower) ? 0.0 : Fx_follower
    end

    # second degree of freedom
    if isnothing(uy)
        # prescribed load, displacement state variable
        pd2 = false
        pl2 = true
        uy = NaN
        Fy = isnothing(Fy) ? 0.0 : Fy
        Fy_follower = isnothing(Fy_follower) ? 0.0 : Fy_follower
    elseif isnothing(Fz) && isnothing(Fz_follower)
        # prescribed displacement, load state variable
        pd2 = true
        pl2 = false
        uy = isnothing(uy) ? 0.0 : uy
        Fy = NaN
        Fy_follower = 0.0
    else
        # prescribed displacement and load, body-fixed acceleration state variable
        pd2 = true
        pl2 = true
        uy = isnothing(uy) ? 0.0 : uy
        Fy = isnothing(Fy) ? 0.0 : Fy
        Fy_follower = isnothing(Fy_follower) ? 0.0 : Fy_follower
    end

    # third degree of freedom
    if isnothing(uz)
        # prescribed load, displacement state variable
        pd3 = false
        pl3 = true
        uz = NaN
        Fz = isnothing(Fz) ? 0.0 : Fz
        Fz_follower = isnothing(Fz_follower) ? 0.0 : Fz_follower
    elseif isnothing(Fz) && isnothing(Fz_follower)
        # prescribed displacement, load state variable
        pd3 = true
        pl3 = false
        uz = isnothing(uz) ? 0.0 : uz
        Fz = NaN
        Fz_follower = 0.0
    else
        # prescribed displacement and load, body-fixed acceleration state variable
        pd3 = true
        pl3 = true
        uz = isnothing(uz) ? 0.0 : uz
        Fz = isnothing(Fz) ? 0.0 : Fz
        Fz_follower = isnothing(Fz_follower) ? 0.0 : Fz_follower
    end

    # fourth degree of freedom
    if isnothing(theta_x)
        # prescribed load, displacement state variable
        pd4 = false
        pl4 = true
        theta_x = NaN
        Mx = isnothing(Mx) ? 0.0 : Mx
        Mx_follower = isnothing(Mx_follower) ? 0.0 : Mx_follower
    elseif isnothing(Mx) && isnothing(Mx_follower)
        # prescribed displacement, load state variable
        pd4 = true
        pl4 = false
        theta_x = isnothing(theta_x) ? 0.0 : theta_x
        Mx = NaN
        Mx_follower = 0.0
    else
        # prescribed displacement and load, body-fixed acceleration state variable
        pd4 = true
        pl4 = true
        theta_x = isnothing(theta_x) ? 0.0 : theta_x
        Mx = isnothing(Mx) ? 0.0 : Mx
        Mx_follower = isnothing(Mx_follower) ? 0.0 : Mx_follower
    end

    # fifth degree of freedom
    if isnothing(theta_y)
        # prescribed load, displacement state variable
        pd5 = false
        pl5 = true
        theta_y = NaN
        My = isnothing(My) ? 0.0 : My
        My_follower = isnothing(My_follower) ? 0.0 : My_follower
    elseif isnothing(My) && isnothing(My_follower)
        # prescribed displacement, load state variable
        pd5 = true
        pl5 = false
        theta_y = isnothing(theta_y) ? 0.0 : theta_y
        My = NaN
        My_follower = 0.0
    else
        # prescribed displacement and load, body-fixed acceleration state variable
        pd5 = true
        pl5 = true
        theta_y = isnothing(theta_y) ? 0.0 : theta_y
        My = isnothing(My) ? 0.0 : My
        My_follower = isnothing(My_follower) ? 0.0 : My_follower
    end

    # sixth degree of freedom
    if isnothing(theta_z)
        # prescribed load, displacement state variable
        pd6 = false
        pl6 = true
        theta_z = NaN
        Mz = isnothing(Mz) ? 0.0 : Mz
        Mz_follower = isnothing(Mz_follower) ? 0.0 : Mz_follower
    elseif isnothing(Mz) && isnothing(Mz_follower)
        # prescribed displacement, load state variable
        pd6 = true
        pl6 = false
        theta_z = isnothing(theta_z) ? 0.0 : theta_z
        Mz = NaN
        Mz_follower = 0.0
    else
        # prescribed displacement and load, body-fixed acceleration state variable
        pd6 = true
        pl6 = true
        theta_z = isnothing(theta_z) ? 0.0 : theta_z
        Mz = isnothing(Mz) ? 0.0 : Mz
        Mz_follower = isnothing(Mz_follower) ? 0.0 : Mz_follower
    end

    # define prescribed conditions
    pd = SVector(pd1, pd2, pd3, pd4, pd5, pd6)
    pl = SVector(pl1, pl2, pl3, pl4, pl5, pl6)
    u = SVector(ux, uy, uz)
    theta = SVector(theta_x, theta_y, theta_z)
    F = SVector(Fx, Fy, Fz)
    M = SVector(Mx, My, Mz)
    Ff = SVector(Fx_follower, Fy_follower, Fz_follower)
    Mf = SVector(Mx_follower, My_follower, Mz_follower)

    return PrescribedConditions(pd, pl, u, theta, F, M, Ff, Mf)
end


"""
    DistributedLoads{T}

Type which contains pre-integrated distributed forces and moments applied to a beam element.

# Fields
- f1: Integrated non-follower distributed force corresponding to the start of the beam element.
- f2: Integrated non-follower distributed force corresponding to the end of the beam element.
- m1: Integrated non-follower distributed moment corresponding to the start of the beam element.
- m2: Integrated non-follower distributed moment corresponding to the end of the beam element.
- f1_follower: Integrated follower distributed force corresponding to the start of the beam element.
- f2_follower: Integrated follower distributed force corresponding to the end of the beam element.
- m1_follower: Integrated follower distributed moment corresponding to the start of the beam element.
- m2_follower: Integrated follower distributed moment corresponding to the end of the beam element.
"""
struct DistributedLoads{T}
    f1::SVector{3, T}
    f2::SVector{3, T}
    m1::SVector{3, T}
    m2::SVector{3, T}
    f1_follower::SVector{3, T}
    f2_follower::SVector{3, T}
    m1_follower::SVector{3, T}
    m2_follower::SVector{3, T}
end
Base.eltype(::DistributedLoads{T}) where T = T

function DistributedLoads{T}(d::DistributedLoads) where T
    return DistributedLoads{T}(d.f1, d.f2, d.m1, d.m2, d.f1_follower, d.f2_follower,
        d.m1_follower, d.m2_follower)
end

DistributedLoads(f1, f2, m1, m2, f1_follower, f2_follower, m1_follower, m2_follower) =
DistributedLoads(promote(f1, f2, m1, m2, f1_follower, f2_follower, m1_follower, m2_follower)...)

"""
    DistributedLoads(assembly, ielem; kwargs...)

Pre-integrate distributed loads on a beam element for use in an analysis.

# Arguments
- `assembly`: Beam element assembly (of type [`Assembly`](@ref))
- `ielem`: Beam element index

# Keyword Arguments
- `s1 = 0.0`: Start of the beam element (used solely for integrating the distributed loads)
- `s2 = 1.0`: End of the beam element (used solely for integrating the distributed loads)
- `method = (f, s1, s2) -> gauss_quadrature(f, s1, s2)`: Method which integrates function
`f` from `s1` to `s2`. Defaults to the Gauss-Legendre quadrature with 4 points 
on each element.
- `fx = (s) -> 0.0`: Distributed x-direction force
- `fy = (s) -> 0.0`: Distributed y-direction force
- `fz = (s) -> 0.0`: Distributed z-direction force
- `mx = (s) -> 0.0`: Distributed x-direction moment
- `my = (s) -> 0.0`: Distributed y-direction moment
- `mz = (s) -> 0.0`: Distributed z-direction moment
- `fx_follower = (s) -> 0.0`: Distributed x-direction follower force
- `fy_follower = (s) -> 0.0`: Distributed y-direction follower force
- `fz_follower = (s) -> 0.0`: Distributed z-direction follower force
- `mx_follower = (s) -> 0.0`: Distributed x-direction follower moment
- `my_follower = (s) -> 0.0`: Distributed y-direction follower moment
- `mz_follower = (s) -> 0.0`: Distributed z-direction follower moment
"""
function DistributedLoads(assembly, ielem;
    s1 = 0.0,
    s2 = 1.0,
    method = gauss_quadrature,
    fx = (s) -> 0.0,
    fy = (s) -> 0.0,
    fz = (s) -> 0.0,
    mx = (s) -> 0.0,
    my = (s) -> 0.0,
    mz = (s) -> 0.0,
    fx_follower = (s) -> 0.0,
    fy_follower = (s) -> 0.0,
    fz_follower = (s) -> 0.0,
    mx_follower = (s) -> 0.0,
    my_follower = (s) -> 0.0,
    mz_follower = (s) -> 0.0,
    )

    return DistributedLoads(assembly, ielem, 0.0;
        s1 = s1,
        s2 = s2,
        method = method,
        fx = (s, t) -> fx(s),
        fy = (s, t) -> fy(s),
        fz = (s, t) -> fz(s),
        mx = (s, t) -> mx(s),
        my = (s, t) -> my(s),
        mz = (s, t) -> mz(s),
        fx_follower = (s, t) -> fx_follower(s),
        fy_follower = (s, t) -> fy_follower(s),
        fz_follower = (s, t) -> fz_follower(s),
        mx_follower = (s, t) -> mx_follower(s),
        my_follower = (s, t) -> my_follower(s),
        mz_follower = (s, t) -> mz_follower(s),
        )
end

"""
    DistrubutedLoads(assembly, ielem, t; kwargs...)

Pre-integrate distributed loads on a beam element for use in an analysis.

# Arguments
- `assembly`: Beam element assembly (of type [`Assembly`](@ref))
- `ielem`: Beam element index
- `t`: Current time

# Keyword Arguments
- `s1 = 0.0`: Start of the beam element (used solely for integrating the distributed loads)
- `s2 = 1.0`: End of the beam element (used solely for integrating the distributed loads)
- `method = (f, s1, s2) -> gauss_quadrature(f, s1, s2)`: Method which integrates function
`f` from `s1` to `s2`. Defaults to the Gauss-Legendre quadrature with 4 points on each element.
- `fx = (s, t) -> 0.0`: Distributed x-direction force
- `fy = (s, t) -> 0.0`: Distributed y-direction force
- `fz = (s, t) -> 0.0`: Distributed z-direction force
- `mx = (s, t) -> 0.0`: Distributed x-direction moment
- `my = (s, t) -> 0.0`: Distributed y-direction moment
- `mz = (s, t) -> 0.0`: Distributed z-direction moment
- `fx_follower = (s, t) -> 0.0`: Distributed x-direction follower force
- `fy_follower = (s, t) -> 0.0`: Distributed y-direction follower force
- `fz_follower = (s, t) -> 0.0`: Distributed z-direction follower force
- `mx_follower = (s, t) -> 0.0`: Distributed x-direction follower moment
- `my_follower = (s, t) -> 0.0`: Distributed y-direction follower moment
- `mz_follower = (s, t) -> 0.0`: Distributed z-direction follower moment
"""
function DistributedLoads(assembly, ielem, t;
    s1 = 0.0,
    s2 = 1.0,
    method = gauss_quadrature,
    fx = (s, t) -> 0.0,
    fy = (s, t) -> 0.0,
    fz = (s, t) -> 0.0,
    mx = (s, t) -> 0.0,
    my = (s, t) -> 0.0,
    mz = (s, t) -> 0.0,
    fx_follower = (s, t) -> 0.0,
    fy_follower = (s, t) -> 0.0,
    fz_follower = (s, t) -> 0.0,
    mx_follower = (s, t) -> 0.0,
    my_follower = (s, t) -> 0.0,
    mz_follower = (s, t) -> 0.0,
    )

    # get element length
    ΔL = assembly.elements[ielem].L

    # create function for general coordinate
    ξ = (s) -> (s-s1)/(s2-s1)

    # element length scaling
    scaling = ΔL/(s2-s1)

    # integrate to get f1 at this point in time
    f1 = SVector(
        method((s)->(1-ξ(s))*fx(s, t)*scaling, s1, s2),
        method((s)->(1-ξ(s))*fy(s, t)*scaling, s1, s2),
        method((s)->(1-ξ(s))*fz(s, t)*scaling, s1, s2)
    )

    # integrate to get f2 at this point in time
    f2 = SVector(
        method((s)->ξ(s)*fx(s, t)*scaling, s1, s2),
        method((s)->ξ(s)*fy(s, t)*scaling, s1, s2),
        method((s)->ξ(s)*fz(s, t)*scaling, s1, s2)
    )

    # integrate to get m1 at this point in time
    m1 = SVector(
        method((s)->(1-ξ(s))*mx(s, t)*scaling, s1, s2),
        method((s)->(1-ξ(s))*my(s, t)*scaling, s1, s2),
        method((s)->(1-ξ(s))*mz(s, t)*scaling, s1, s2)
    )

    # integrate to get m2 at this point in time
    m2 = SVector(
        method((s)->ξ(s)*mx(s, t)*scaling, s1, s2),
        method((s)->ξ(s)*my(s, t)*scaling, s1, s2),
        method((s)->ξ(s)*mz(s, t)*scaling, s1, s2)
    )

    # integrate to get f1_follower at this point in time
    f1_follower =SVector(
        method((s)->(1-ξ(s))*fx_follower(s, t)*scaling, s1, s2),
        method((s)->(1-ξ(s))*fy_follower(s, t)*scaling, s1, s2),
        method((s)->(1-ξ(s))*fz_follower(s, t)*scaling, s1, s2)
    )

    # integrate to get f2_follower at this point in time
    f2_follower = SVector(
        method((s)->ξ(s)*fx_follower(s, t)*scaling, s1, s2),
        method((s)->ξ(s)*fy_follower(s, t)*scaling, s1, s2),
        method((s)->ξ(s)*fz_follower(s, t)*scaling, s1, s2)
    )

    # integrate to get m1_follower at this point in time
    m1_follower = SVector(
        method((s)->(1-ξ(s))*mx_follower(s, t)*scaling, s1, s2),
        method((s)->(1-ξ(s))*my_follower(s, t)*scaling, s1, s2),
        method((s)->(1-ξ(s))*mz_follower(s, t)*scaling, s1, s2)
    )

    # integrate to get m2_follower at this point in time
    m2_follower = SVector(
        method((s)->ξ(s)*mx_follower(s, t)*scaling, s1, s2),
        method((s)->ξ(s)*my_follower(s, t)*scaling, s1, s2),
        method((s)->ξ(s)*mz_follower(s, t)*scaling, s1, s2)
    )

    return DistributedLoads(f1, f2, m1, m2, f1_follower, f2_follower, m1_follower, m2_follower)
end

"""
    combine_loads(loads)

Combine the distributed loads in the iterable collection `loads`
"""
function combine_loads(loads)
    f1 = @SVector zeros(3)
    f2 = @SVector zeros(3)
    m1 = @SVector zeros(3)
    m2 = @SVector zeros(3)
    f1_follower = @SVector zeros(3)
    f2_follower = @SVector zeros(3)
    m1_follower = @SVector zeros(3)
    m2_follower = @SVector zeros(3)
    for load in loads
        f1 += load.f1
        f2 += load.f2
        m1 += load.m1
        m2 += load.m2
        f1_follower += load.f1_follower
        f2_follower += load.f2_follower
        m1_follower += load.m1_follower
        m2_follower += load.m2_follower
    end
    return DistributedLoads(f1, f2, m1, m2, f1_follower, f2_follower, m1_follower, m2_follower)
end

"""
    PointMass{T}

Type which contains the aggregated inertial properties of one or more point masses which 
are rigidly attached to the center of an element.

# Fields:
 - `mass`: Mass matrix corresponding to the point masses.
"""
struct PointMass{T}
    mass::SMatrix{6,6,T,36}
end
Base.eltype(::PointMass{T}) where T = T

function PointMass{T}(p::PointMass) where T
    PointMass{T}(p.mass)
end

"""
    PointMass(mass)

Define a point mass given its mass matrix
"""
PointMass(mass::AbstractMatrix{T}) where T = PointMass{T}(SMatrix{6,6,T,36}(mass))

"""
    PointMass(m, p, I)

Define a point mass given its mass `m`, offset `p`, and inertia matrix `I`
"""
PointMass(m, p, I) = PointMass([m*I3 -m*tilde(p); m*tilde(p) I-m*tilde(p)*tilde(p)])

"""
    combine_masses(masses)

Combine the point masses in the iterable collection `masses`
"""
function combine_masses(masses)
    M = @SMatrix zeros(6,6)
    for mass in masses
        M += mass.M
    end
    return PointMass(M)
end