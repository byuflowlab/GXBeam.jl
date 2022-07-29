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