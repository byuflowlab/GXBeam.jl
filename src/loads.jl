"""
    PrescribedConditions{T}

Describes the forces, moments, displacements, and/or rotations prescribed at a
point for each time step
"""
struct PrescribedConditions{T}
    force::Vector{Bool}
    value::Vector{SVector{6, T}}
    follower::Vector{SVector{6, T}}
end
Base.eltype(::PrescribedConditions{T}) where T = T

"""
    PrescribedConditions(dt=0.0; kwargs...)

Construct an object of type PrescribedConditions which stores the prescribed
conditions for a point at each time step.

Prescribed conditions may be assigned as either a scalar parameter or as a
function of time.

Prescribed Wiener-Milenkovic parameters must satisfy the following inequality:
sqrt(theta_x^2 + theta_y^2 + theta_z^2) <= 4.  Note that this restriction still
allows all possible rotations to be represented.

# Arguments
 - `dt`: Time step size.
 - `nstep`: The total length of the time vector
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
function PrescribedConditions(dt=0.0;
    nstep = 1,
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

    force = Vector{Bool}(undef, 6)

    t = range(0, dt*(nstep-1), length=nstep)

    # Set function for first slot
    if isnothing(ux)
        force[1] = true
        v1 = isnothing(Fx) ? (t) -> 0.0 : typeof(Fx) <: Number ? (t) -> Fx : Fx
        v1_follower = isnothing(Fx_follower) ? (t) -> 0.0 : typeof(Fx_follower) <: Number ? (t) -> Fx_follower : Fx
    elseif !isnothing(ux) && (isnothing(Fx) && isnothing(Fx_follower))
        force[1] = false
        v1 = isnothing(ux) ? (t) -> 0.0 : typeof(ux) <: Number ? (t) -> ux : ux
        v1_follower = (t) -> 0.0
    else
        error("Both `ux` and `Fx` or `Fx_follower` cannot be specified at the same time")
    end

    # Set function for second slot
    if isnothing(uy)
        force[2] = true
        v2 = isnothing(Fy) ? (t) -> 0.0 : typeof(Fy) <: Number ? (t) -> Fy : Fy
        v2_follower = isnothing(Fy_follower) ? (t) -> 0.0 : typeof(Fy_follower) <: Number ? (t) -> Fy_follower : Fy
    elseif !isnothing(uy) && (isnothing(Fy) && isnothing(Fy_follower))
        force[2] = false
        v2 = isnothing(uy) ? (t) -> 0.0 : typeof(uy) <: Number ? (t) -> uy : uy
        v2_follower = (t) -> 0.0
    else
        error("Both `uy` and `Fy` or `Fy_follower` cannot be specified at the same time")
    end

    # Set function for third slot
    if isnothing(uz)
        force[3] = true
        v3 = isnothing(Fz) ? (t) -> 0.0 : typeof(Fz) <: Number ? (t) -> Fz : Fz
        v3_follower = isnothing(Fz_follower) ? (t) -> 0.0 : typeof(Fz_follower) <: Number ? (t) -> Fz_follower : Fz
    elseif !isnothing(uz) && (isnothing(Fz) && isnothing(Fz_follower))
        force[3] = false
        v3 = isnothing(uz) ? (t) -> 0.0 : typeof(uz) <: Number ? (t) -> uz : uz
        v3_follower = (t) -> 0.0
    else
        error("Both `uz` and `Fz` or `Fz_follower` cannot be specified at the same time")
    end

    # Set function for fourth slot
    if isnothing(theta_x)
        force[4] = true
        v4 = isnothing(Mx) ? (t) -> 0.0 : typeof(Mx) <: Number ? (t) -> Mx : Mx
        v4_follower = isnothing(Mx_follower) ? (t) -> 0.0 : typeof(Mx_follower) <: Number ? (t) -> Mx_follower : Mx
    elseif !isnothing(theta_x) && (isnothing(Mx) && isnothing(Mx_follower))
        force[4] = false
        v4 = isnothing(theta_x) ? (t) -> 0.0 : typeof(theta_x) <: Number ? (t) -> theta_x : theta_x
        v4_follower = (t) -> 0.0
    else
        error("Both `theta_x` and `Mx` or `Mx_follower` cannot be specified at the same time")
    end

    # Set function for fifth slot
    if isnothing(theta_y)
        force[5] = true
        v5 = isnothing(My) ? (t) -> 0.0 : typeof(My) <: Number ? (t) -> My : My
        v5_follower = isnothing(My_follower) ? (t) -> 0.0 : typeof(My_follower) <: Number ? (t) -> My_follower : My
    elseif !isnothing(theta_y) && (isnothing(My) && isnothing(My_follower))
        force[5] = false
        v5 = isnothing(theta_y) ? (t) -> 0.0 : typeof(theta_y) <: Number ? (t) -> theta_y : theta_y
        v5_follower = (t) -> 0.0
    else
        error("Both `theta_y` and `My` or `My_follower` cannot be specified at the same time")
    end

    # Set function for sixth slot
    if isnothing(theta_z)
        force[6] = true
        v6 = isnothing(Mz) ? (t) -> 0.0 : typeof(Mz) <: Number ? (t) -> Mz : Mz
        v6_follower = isnothing(Mz_follower) ? (t) -> 0.0 : typeof(Mz_follower) <: Number ? (t) -> Mz_follower : Mz
    elseif !isnothing(theta_z) && (isnothing(Mz) && isnothing(Mz_follower))
        force[6] = false
        v6 = isnothing(theta_z) ? (t) -> 0.0 : typeof(theta_z) <: Number ? (t) -> theta_z : theta_z
        v6_follower = (t) -> 0.0
    else
        error("Both `theta_z` and `Mz` or `Mz_follower` cannot be specified at the same time")
    end

    # now populate values for each time step
    value = [
        SVector{6}(
            v1(t[i]), v2(t[i]), v3(t[i]),
            v4(t[i]), v5(t[i]), v6(t[i])
            ) for i in eachindex(t)
    ]

    follower = [
        SVector{6}(
            v1_follower(t[i]), v2_follower(t[i]), v3_follower(t[i]),
            v4_follower(t[i]), v5_follower(t[i]), v6_follower(t[i])
        ) for i in eachindex(t)
    ]

    return PrescribedConditions(force, promote(value, follower)...)
end

"""
    DistributedLoads{T}

Contains the integrated distributed forces and moments for each beam element for each time step.

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
    f1::Vector{SVector{3, T}}
    f2::Vector{SVector{3, T}}
    m1::Vector{SVector{3, T}}
    m2::Vector{SVector{3, T}}
    f1_follower::Vector{SVector{3, T}}
    f2_follower::Vector{SVector{3, T}}
    m1_follower::Vector{SVector{3, T}}
    m2_follower::Vector{SVector{3, T}}
end
Base.eltype(::DistributedLoads{T}) where T = T

"""
    DistributedLoads(assembly, ibeam[, dt]; kwargs...)

Integrates the specified distributed loads over the element for each time step.

# Arguments
 - `assembly`: The beam element assembly
 - `ibeam`: The index of the beam element which the distributed load is assigned to
 - `dt`: Time step size.  If omitted a single time step is assumed and specified
     functions become a function of `s` only.
 - `s1 = 0.0`: Start of beam element (used for integrating the distributed loads)
 - `s2 = 1.0`: End of beam element (used for integrating the distributed loads)
 - `nstep`: The total length of the time vector
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
function DistributedLoads(assembly, ibeam;
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

    return DistributedLoads(assembly, ibeam, 0.0;
        s1 = s1,
        s2 = s2,
        nstep = 1,
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

function DistributedLoads(assembly, ibeam, dt;
    s1 = 0.0,
    s2 = 1.0,
    nstep = 1,
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
    ΔL = assembly.elements[ibeam].L

    # create function for general coordinate
    ξ = (s) -> (s-s1)/(s2-s1)

    # create time range
    t = range(0.0, step=dt, length=nstep+1)

    # integrate to get f1 for each time step
    f1 = [
        SVector(
            method((s)->(1-ξ(s))*fx(s, t[i])*ΔL/(s2-s1), s1, s2),
            method((s)->(1-ξ(s))*fy(s, t[i])*ΔL/(s2-s1), s1, s2),
            method((s)->(1-ξ(s))*fz(s, t[i])*ΔL/(s2-s1), s1, s2)
        ) for i in eachindex(t)
    ]

    # integrate to get f2 for each time step
    f2 = [
        SVector(
            method((s)->ξ(s)*fx(s, t[i])*ΔL/(s2-s1), s1, s2),
            method((s)->ξ(s)*fy(s, t[i])*ΔL/(s2-s1), s1, s2),
            method((s)->ξ(s)*fz(s, t[i])*ΔL/(s2-s1), s1, s2)
        ) for i in eachindex(t)
    ]

    # integrate to get m1 for each time step
    m1 = [
        SVector(
            method((s)->(1-ξ(s))*mx(s, t[i])*ΔL/(s2-s1), s1, s2),
            method((s)->(1-ξ(s))*my(s, t[i])*ΔL/(s2-s1), s1, s2),
            method((s)->(1-ξ(s))*mz(s, t[i])*ΔL/(s2-s1), s1, s2)
        ) for i in eachindex(t)
    ]

    # integrate to get m2 for each time step
    m2 = [
        SVector(
            method((s)->ξ(s)*mx(s, t[i])*ΔL/(s2-s1), s1, s2),
            method((s)->ξ(s)*my(s, t[i])*ΔL/(s2-s1), s1, s2),
            method((s)->ξ(s)*mz(s, t[i])*ΔL/(s2-s1), s1, s2)
        ) for i in eachindex(t)
    ]

    # integrate to get f1_follower for each time step
    f1_follower = [
        SVector(
            method((s)->(1-ξ(s))*fx_follower(s, t[i])*ΔL/(s2-s1), s1, s2),
            method((s)->(1-ξ(s))*fy_follower(s, t[i])*ΔL/(s2-s1), s1, s2),
            method((s)->(1-ξ(s))*fz_follower(s, t[i])*ΔL/(s2-s1), s1, s2)
        ) for i in eachindex(t)
    ]

    # integrate to get f2_follower for each time step
    f2_follower = [
        SVector(
            method((s)->ξ(s)*fx_follower(s, t[i])*ΔL/(s2-s1), s1, s2),
            method((s)->ξ(s)*fy_follower(s, t[i])*ΔL/(s2-s1), s1, s2),
            method((s)->ξ(s)*fz_follower(s, t[i])*ΔL/(s2-s1), s1, s2)
        ) for i in eachindex(t)
    ]

    # integrate to get m1_follower for each time step
    m1_follower = [
        SVector(
            method((s)->(1-ξ(s))*mx_follower(s, t[i])*ΔL/(s2-s1), s1, s2),
            method((s)->(1-ξ(s))*my_follower(s, t[i])*ΔL/(s2-s1), s1, s2),
            method((s)->(1-ξ(s))*mz_follower(s, t[i])*ΔL/(s2-s1), s1, s2)
        ) for i in eachindex(t)
    ]

    # integrate to get m2_follower for each time step
    m2_follower = [
        SVector(
            method((s)->ξ(s)*mx_follower(s, t[i])*ΔL/(s2-s1), s1, s2),
            method((s)->ξ(s)*my_follower(s, t[i])*ΔL/(s2-s1), s1, s2),
            method((s)->ξ(s)*mz_follower(s, t[i])*ΔL/(s2-s1), s1, s2)
        ) for i in eachindex(t)
    ]

    return DistributedLoads(f1, f2, m1, m2, f1_follower, f2_follower, m1_follower, m2_follower)
end
