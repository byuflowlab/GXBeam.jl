
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

