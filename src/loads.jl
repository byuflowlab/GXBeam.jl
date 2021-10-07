"""
    PrescribedConditions{T}

Describes the forces, moments, displacements, and/or rotations prescribed at a
point at a specific point in time.
"""
struct PrescribedConditions{T}
    force::SVector{6, Bool}
    value::SVector{6, T}
    follower::SVector{6, T}
end
Base.eltype(::PrescribedConditions{T}) where T = T

function PrescribedConditions{T}(p::PrescribedConditions) where T
    PrescribedConditions{T}(p.force, p.value, p.follower)
end

PrescribedConditions(force, value, follower) = PrescribedConditions(force, promote(value, follower)...)

"""
    PrescribedConditions(t = 0.0; kwargs...)

Return the combined prescribed conditions at a point (as an object of type
[`PrescribedConditions`](@ref)) at time `t`.

Individual prescribed conditions may be assigned as either a scalar parameter or
as a function of time.

Prescribed Wiener-Milenkovic parameters must satisfy the following inequality:
sqrt(theta_x^2 + theta_y^2 + theta_z^2) <= 4.  Note that this restriction still
allows all possible rotations to be represented.

# Keyword Arguments
 - `ux`: Prescribed x-direction displacement of the point
 - `uy`: Prescribed y-direction displacement of the point
 - `uz`: Prescribed z-direction displacement of the point
 - `theta_x`: Prescribed first Wiener-Milenkovic parameter of the point
 - `theta_y`: Prescribed second Wiener-Milenkovic parameter of the point
 - `theta_z`: Prescribed third Wiener-Milenkovic parameter of the point
 - `Fx`: Prescribed force in x-direction applied on the point
 - `Fy`: Prescribed force in y-direction applied on the point
 - `Fz`: Prescribed force in z-direction applied on the point
 - `Mx`: Prescribed moment about x-axis applied on the point
 - `My`: Prescribed moment about y-axis applied on the point
 - `Mz`: Prescribed moment about z-axis applied on the point
 - `Fx_follower`: Prescribed follower force in x-direction applied on the point
 - `Fy_follower`: Prescribed follower force in y-direction applied on the point
 - `Fz_follower`: Prescribed follower force in z-direction applied on the point
 - `Mx_follower`: Prescribed follower moment about x-axis applied on the point
 - `My_follower`: Prescribed follower moment about y-axis applied on the point
 - `Mz_follower`: Prescribed follower moment about z-axis applied on the point
"""
function PrescribedConditions(t = 0.0;
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

    # Set function for first slot
    if isnothing(ux)
        f1 = true
        v1 = isnothing(Fx) ? 0.0 : typeof(Fx) <: Number ? Fx : Fx(t)
        v1_f = isnothing(Fx_follower) ? 0.0 : typeof(Fx_follower) <: Number ? Fx_follower : Fx_follower(t)
    elseif !isnothing(ux) && (isnothing(Fx) && isnothing(Fx_follower))
        f1 = false
        v1 = isnothing(ux) ? 0.0 : typeof(ux) <: Number ? ux : ux(t)
        v1_f = 0.0
    else
        error("Both `ux` and `Fx` or `Fx_follower` cannot be specified at the same time")
    end

    # Set function for second slot
    if isnothing(uy)
        f2 = true
        v2 = isnothing(Fy) ? 0.0 : typeof(Fy) <: Number ? Fy : Fy(t)
        v2_f = isnothing(Fy_follower) ? 0.0 : typeof(Fy_follower) <: Number ? Fy_follower : Fy_follower(t)
    elseif !isnothing(uy) && (isnothing(Fy) && isnothing(Fy_follower))
        f2 = false
        v2 = isnothing(uy) ? 0.0 : typeof(uy) <: Number ? uy : uy(t)
        v2_f = 0.0
    else
        error("Both `uy` and `Fy` or `Fy_follower` cannot be specified at the same time")
    end

    # Set function for third slot
    if isnothing(uz)
        f3 = true
        v3 = isnothing(Fz) ? 0.0 : typeof(Fz) <: Number ? Fz : Fz(t)
        v3_f = isnothing(Fz_follower) ? 0.0 : typeof(Fz_follower) <: Number ? Fz_follower : Fz_follower(t)
    elseif !isnothing(uz) && (isnothing(Fz) && isnothing(Fz_follower))
        f3 = false
        v3 = isnothing(uz) ? 0.0 : typeof(uz) <: Number ? uz : uz(t)
        v3_f = 0.0
    else
        error("Both `uz` and `Fz` or `Fz_follower` cannot be specified at the same time")
    end

    # Set function for fourth slot
    if isnothing(theta_x)
        f4 = true
        v4 = isnothing(Mx) ? 0.0 : typeof(Mx) <: Number ? Mx : Mx(t)
        v4_f = isnothing(Mx_follower) ? 0.0 : typeof(Mx_follower) <: Number ? Mx_follower : Mx_follower(t)
    elseif !isnothing(theta_x) && (isnothing(Mx) && isnothing(Mx_follower))
        f4 = false
        v4 = isnothing(theta_x) ? 0.0 : typeof(theta_x) <: Number ? theta_x : theta_x(t)
        v4_f = 0.0
    else
        error("Both `theta_x` and `Mx` or `Mx_follower` cannot be specified at the same time")
    end

    # Set function for fifth slot
    if isnothing(theta_y)
        f5 = true
        v5 = isnothing(My) ? 0.0 : typeof(My) <: Number ? My : My(t)
        v5_f = isnothing(My_follower) ? 0.0 : typeof(My_follower) <: Number ? My_follower : My_follower(t)
    elseif !isnothing(theta_y) && (isnothing(My) && isnothing(My_follower))
        f5 = false
        v5 = isnothing(theta_y) ? 0.0 : typeof(theta_y) <: Number ? theta_y : theta_y(t)
        v5_f = 0.0
    else
        error("Both `theta_y` and `My` or `My_follower` cannot be specified at the same time")
    end

    # Set function for sixth slot
    if isnothing(theta_z)
        f6 = true
        v6 = isnothing(Mz) ? 0.0 : typeof(Mz) <: Number ? Mz : Mz(t)
        v6_f = isnothing(Mz_follower) ? 0.0 : typeof(Mz_follower) <: Number ? Mz_follower : Mz_follower(t)
    elseif !isnothing(theta_z) && (isnothing(Mz) && isnothing(Mz_follower))
        f6 = false
        v6 = isnothing(theta_z) ? 0.0 : typeof(theta_z) <: Number ? theta_z : theta_z(t)
        v6_f = 0.0
    else
        error("Both `theta_z` and `Mz` or `Mz_follower` cannot be specified at the same time")
    end

    force = SVector(f1, f2, f3, f4, f5, f6)
    value = SVector(v1, v2, v3, v4, v5, v6)
    follower = SVector(v1_f, v2_f, v3_f, v4_f, v5_f, v6_f)

    return PrescribedConditions(force, value, follower)
end

"""
    DistributedLoads{T}

Contains the integrated distributed forces and moments for each beam element at
a specific time step.

# Fields
 - f1: Integrated non-follower distributed force for the beam element's left endpoint for each time step
 - f2: Integrated non-follower distributed force for the beam element's right endpoint for each time step
 - m1: Integrated non-follower distributed moment for the beam element's left endpoint for each time step
 - m2: Integrated non-follower distributed moment for the beam element's right endpoint for each time step
 - f1_follower: Integrated follower distributed force for the beam element's left endpoint for each time step
 - f2_follower: Integrated follower distributed force for the beam element's right endpoint for each time step
 - m1_follower: Integrated follower distributed moment for the beam element's left endpoint for each time step
 - m2_follower: Integrated follower distributed moment for the beam element's right endpoint for each time step
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

Return the integrated distributed loads at a point (as an object of type [`DistributedLoads`](@ref)).

# Arguments
 - `assembly`: The beam element assembly
 - `ielem`: The index of the beam element which the distributed load is assigned to

# Keyword Arguments
 - `s1 = 0.0`: Start of the beam element (used solely for integrating the distributed loads)
 - `s2 = 1.0`: End of the beam element (used solely for integrating the distributed loads)
 - `method = (f, a, b) -> gauss_quadrature(f, a, b)`: Method which integrates function
    `f` from `a` to `b`. Defaults to the Gauss-Legendre quadrature with 4 points on each element.
 - `fx = (s) -> 0.0`: Distributed non-follower force on beam element in x-direction
 - `fy = (s) -> 0.0`: Distributed non-follower force on beam element in y-direction
 - `fz = (s) -> 0.0`: Distributed non-follower force on beam element in z-direction
 - `mx = (s) -> 0.0`: Distributed non-follower moment on beam element in x-direction
 - `my = (s) -> 0.0`: Distributed non-follower moment on beam element in y-direction
 - `mz = (s) -> 0.0`: Distributed non-follower moment on beam element in z-direction
 - `fx_follower = (s) -> 0.0`: Distributed follower force on beam element in x-direction
 - `fy_follower = (s) -> 0.0`: Distributed follower force on beam element in y-direction
 - `fz_follower = (s) -> 0.0`: Distributed follower force on beam element in z-direction
 - `mx_follower = (s) -> 0.0`: Distributed follower moment on beam element in x-direction
 - `my_follower = (s) -> 0.0`: Distributed follower moment on beam element in y-direction
 - `mz_follower = (s) -> 0.0`: Distributed follower moment on beam element in z-direction
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

Return the integrated distributed loads on a beam element (as an object of type
[`DistributedLoads`](@ref)) at time `t`.

# Arguments
 - `assembly`: The beam element assembly
 - `ielem`: The index of the beam element which the distributed load is assigned to
 - `t`: time at which to evaluate the distributed loads

# Keyword Arguments
 - `s1 = 0.0`: Start of the beam element (used solely for integrating the distributed loads)
 - `s2 = 1.0`: End of the beam element (used solely for integrating the distributed loads)
 - `method = (f, a, b) -> gauss_quadrature(f, a, b)`: Method which integrates function
    `f` from `a` to `b`. Defaults to the Gauss-Legendre quadrature with 4 points on each element.
 - `fx = (s, t) -> 0.0`: Distributed non-follower force on beam element in x-direction
 - `fy = (s, t) -> 0.0`: Distributed non-follower force on beam element in y-direction
 - `fz = (s, t) -> 0.0`: Distributed non-follower force on beam element in z-direction
 - `mx = (s, t) -> 0.0`: Distributed non-follower moment on beam element in x-direction
 - `my = (s, t) -> 0.0`: Distributed non-follower moment on beam element in y-direction
 - `mz = (s, t) -> 0.0`: Distributed non-follower moment on beam element in z-direction
 - `fx_follower = (s, t) -> 0.0`: Distributed follower force on beam element in x-direction
 - `fy_follower = (s, t) -> 0.0`: Distributed follower force on beam element in y-direction
 - `fz_follower = (s, t) -> 0.0`: Distributed follower force on beam element in z-direction
 - `mx_follower = (s, t) -> 0.0`: Distributed follower moment on beam element in x-direction
 - `my_follower = (s, t) -> 0.0`: Distributed follower moment on beam element in y-direction
 - `mz_follower = (s, t) -> 0.0`: Distributed follower moment on beam element in z-direction
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
    PointMass{T}

Describes the inertial properties of a mass attached at a point.
"""
struct PointMass{T}
    m::T
    r::SVector{3,T}
    I::SMatrix{3,3,T,9}
end
Base.eltype(::PointMass{T}) where T = T

function PointMass{T}(p::PointMass) where T
    PointMass{T}(p.m, p.r, p.I)
end

"""
    element_gravitational_loads(CtCab, elem, gvec)

Calculate the integrated distributed loads on an element due to gravity.
"""
@inline function element_gravitational_loads(ΔL, CtCab, mass, gvec)
    # extract mass matrix submatrices
    M11 = mass[SVector{3}(1:3), SVector{3}(1:3)]
    M12 = mass[SVector{3}(1:3), SVector{3}(4:6)]
    # calculate force and moment per unit length due to graviational forces
    f = CtCab*M11*CtCab'*gvec
    m = -CtCab*M12*CtCab'*gvec
    # calculate integrated force and moment per unit length
    f1 = f2 = ΔL*f/2
    m1 = m2 = ΔL*m/2
    # return result
    return f1, f2, m1, m2
end

@inline function element_gravitational_loads_jacobian(ΔL, Cab, CtCab, Ct_θ1, Ct_θ2, Ct_θ3, elem, gvec)   
    C_θ1, C_θ2, C_θ3 = Ct_θ1', Ct_θ2', Ct_θ3'
    # extract mass matrix submatrices
    M11 = elem.mass[SVector{3}(1:3), SVector{3}(1:3)]
    M12 = elem.mass[SVector{3}(1:3), SVector{3}(4:6)]
    # calculate force and moment per unit length due to gravitational forces
    f_θ = mul3(Ct_θ1, Ct_θ2, Ct_θ3, Cab*M11*CtCab'*gvec) + CtCab*M11*Cab'*mul3(C_θ1, C_θ2, C_θ3, gvec)
    m_θ = -mul3(Ct_θ1, Ct_θ2, Ct_θ3, Cab*M12*CtCab'*gvec) + CtCab*M12*Cab'*mul3(C_θ1, C_θ2, C_θ3, gvec)
    # calculate integrated force and moment per unit length
    f1_θ = f2_θ = ΔL*f_θ/2
    m1_θ = m2_θ = ΔL*m_θ/2
    # return result
    return f1_θ, f2_θ, m1_θ, m2_θ
end

"""
    combine_loads(l1, l2)

Combine the distributed loads in `l1` and `l2`
"""
function combine_loads(l1, l2)
    f1 = l1.f1 + l2.f1
    f2 = l1.f2 + l2.f2
    m1 = l1.m1 + l2.m1
    m2 = l1.m2 + l2.m2
    f1_follower = l1.f1_follower + l2.f1_follower
    f2_follower = l1.f2_follower + l2.f2_follower
    m1_follower = l1.m1_follower + l2.m1_follower
    m2_follower = l1.m2_follower + l2.m2_follower
    return DistributedLoads(f1, f2, m1, m2, f1_follower, f2_follower, m1_follower, m2_follower)
end
