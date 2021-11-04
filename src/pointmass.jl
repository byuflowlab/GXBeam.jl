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

PointMass(mass::AbstractMatrix{T}) where T = PointMass{T}(SMatrix{6,6,T,36}(mass))

"""
    point_mass_linear_momentum(point_mass, V, Ω)

Calculate the linear momentum of one or more masses which are rigidly connected to the 
center of an element given the element's linear and angular velocity.  Return this property 
in the deformed point mass reference frame.
"""
@inline function point_mass_linear_momentum(point_mass, V, Ω)
    massp = point_mass.mass
    massp11 = massp[SVector{3}(1:3), SVector{3}(1:3)]
    massp12 = massp[SVector{3}(1:3), SVector{3}(4:6)]
    return massp11*V + massp12*Ω
end

"""
    point_mass_angular_momentum(point_mass, V, Ω)

Calculate the angular momentum of one or more masses which are rigidly connected to the 
center of an element give the element's linear and angular velocity.  Return this property 
in the deformed point mass refeference frame.
"""
@inline function point_mass_angular_momentum(point_mass, V, Ω)
    massp = point_mass.mass
    massp21 = massp[SVector{3}(4:6), SVector{3}(1:3)]
    massp22 = massp[SVector{3}(4:6), SVector{3}(4:6)]
    return massp21*V + massp22*Ω
end

function static_point_mass_properties(point_mass)

    mass_p = point_mass.mass
    massp11 = mass_p[SVector{3}(1:3), SVector{3}(1:3)]
    massp12 = mass_p[SVector{3}(1:3), SVector{3}(4:6)]
    massp21 = mass_p[SVector{3}(4:6), SVector{3}(1:3)]
    massp22 = mass_p[SVector{3}(4:6), SVector{3}(4:6)]

    return massp11, massp12, massp21, massp22
end

function steady_state_point_mass_properties(point_mass, Cab, V, Ω)

    massp11, massp12, massp21, massp22 = static_point_mass_properties(point_mass)

    # velocities in point mass reference frame
    Vp = Cab*V
    Ωp = Cab*Ω

    # point mass linear and angular momentum
    Pp = massp11*Vp + massp12*Ωp
    Hp = massp21*Vp + massp22*Ωp

    return massp11, massp12, massp21, massp22, Vp, Ωp, Pp, Hp
end

function dynamic_point_mass_properties(point_mass, Ct, Ctdot, Cab, V, Ω, Vdot, Ωdot)

    massp11, massp12, massp21, massp22, Vp, Ωp, Pp, Hp = steady_state_point_mass_properties(point_mass, Cab, V, Ω)

    # accelerations in point mass reference frame
    Vpdot = Cab*Vdot
    Ωpdot = Cab*Ωdot

    # linear and angular momentum derivatives
    Ppdot = massp11*Vpdot + massp12*Ωpdot
    Hpdot = massp21*Vpdot + massp22*Ωpdot

    # calculate \dot{\bar{C^T*P}} and \dot{\bar{C^T*H}}
    CtPpdot = Ctdot*Pp + Ct*Ppdot
    CtHpdot = Ctdot*Hp + Ct*Hpdot

    return massp11, massp12, massp21, massp22, Vp, Ωp, Pp, Hp, Vpdot, Ωpdot, Ppdot, Hpdot, 
        CtPpdot, CtHpdot
end

@inline function static_point_mass_loads(massp11, massp12, C, Ct, gvec)

    # gravitational loads
    Fp = Ct*massp11*C*gvec
    Mp = -Ct*massp12*C*gvec

    return Fp, Mp
end

@inline function steady_state_point_mass_loads(massp11, massp12, massp21, massp22, C, Ct, 
    u, Vp, Pp, Hp, ω, a, α, gvec)

    # linear and angular acceleration in deformed point mass reference frame
    ap = C*(a - gvec + cross(α, u))
    αp = C*α

    # force and moment due to accelerating reference frame
    Fp = -Ct*(massp11*ap + massp12*αp) 
    Mp = -Ct*(massp21*ap + massp22*αp)

    # force and moment due to point mass linear and angular momentum
    Fp -= cross(ω, Ct*Pp)
    Mp -= cross(ω, Ct*Hp) + Ct*cross(Vp, Pp)

    # return result
    return Fp, Mp
end

@inline function dynamic_point_mass_loads(massp11, massp12, massp21, massp22, 
    C, Ct, Ctdot, u, Vp, Pp, Hp, Ppdot, Hpdot, ω, a, α, gvec)

    # steady state point mass loads
    Fp, Mp = steady_state_point_mass_loads(massp11, massp12, massp21, massp22, C, Ct, 
        u, Vp, Pp, Hp, ω, a, α, gvec)

    # force and moment due to local frame acceleration
    Fp -= Ctdot*Pp + Ct*Ppdot
    Mp -= Ctdot*Hp + Ct*Hpdot

    # return result
    return Fp, Mp
end

@inline function static_point_mass_jacobian(massp11, massp12, C, Ct, gvec, C_θ1, C_θ2, C_θ3, 
    Ct_θ1, Ct_θ2, Ct_θ3)

    # calculate force and moment due to gravitational forces on point masses
    Fp_θ = mul3(Ct_θ1, Ct_θ2, Ct_θ3, massp11*C*gvec) + Ct*massp11*mul3(C_θ1, C_θ2, C_θ3, gvec)
    Mp_θ = -mul3(Ct_θ1, Ct_θ2, Ct_θ3, massp12*C*gvec) - Ct*massp12*mul3(C_θ1, C_θ2, C_θ3, gvec)

    # return result
    return Fp_θ, Mp_θ
end

@inline function steady_state_point_mass_jacobian(massp11, massp12, massp21, massp22, C, Ct, 
    Cab, u, Vp, Pp, Hp, ω, a, α, gvec, C_θ1, C_θ2, C_θ3, Ct_θ1, Ct_θ2, Ct_θ3)

    # linear and angular acceleration in deformed point mass reference frame
    ae = C*(a - gvec + cross(α, u))
    αe = C*α
    
    # force and moment due to accelerating reference frame
    Fp_u = -Ct*massp11*C*tilde(α)
    Mp_u = -Ct*massp21*C*tilde(α)

    Fp_θ = -mul3(Ct_θ1, Ct_θ2, Ct_θ3, massp11*ae + massp12*αe) - 
        Ct*massp11*mul3(C_θ1, C_θ2, C_θ3, a - gvec + cross(α, u)) -
        Ct*massp12*mul3(C_θ1, C_θ2, C_θ3, α)
    Mp_θ = -mul3(Ct_θ1, Ct_θ2, Ct_θ3, massp21*ae + massp22*αe) - 
        Ct*massp21*mul3(C_θ1, C_θ2, C_θ3, a - gvec + cross(α, u)) -
        Ct*massp22*mul3(C_θ1, C_θ2, C_θ3, α)

    # force and moment due to point mass linear and angular momentum
    Fp_θ -= tilde(ω)*mul3(Ct_θ1, Ct_θ2, Ct_θ3, Pp)
    Fp_V = -tilde(ω)*Ct*massp11*Cab
    Fp_Ω = -tilde(ω)*Ct*massp12*Cab

    Mp_θ -= tilde(ω)*mul3(Ct_θ1, Ct_θ2, Ct_θ3, Hp) + mul3(Ct_θ1, Ct_θ2, Ct_θ3, cross(Vp, Pp))
    Mp_V = -tilde(ω)*Ct*massp21*Cab - Ct*tilde(Vp)*massp11*Cab + Ct*tilde(Pp)*Cab
    Mp_Ω = -tilde(ω)*Ct*massp22*Cab - Ct*tilde(Vp)*massp12*Cab

    return Fp_u, Mp_u, Fp_θ, Mp_θ, Fp_V, Mp_V, Fp_Ω, Mp_Ω
end

@inline function initial_condition_point_mass_jacobian(massp11, massp12, massp21, 
    massp22, Ct, Ctdot, Cab, Vp, Pp, ω)

    # add loads due to local acceleration
    Fp_Vdot = -Ct*massp11*Cab
    Fp_Ωdot = -Ct*massp12*Cab
    Fp_V = -tilde(ω)*Ct*massp11*Cab - Ctdot*massp11*Cab
    Fp_Ω = -tilde(ω)*Ct*massp12*Cab - Ctdot*massp12*Cab

    Mp_Vdot = -Ct*massp21*Cab
    Mp_Ωdot = -Ct*massp22*Cab
    Mp_V = -tilde(ω)*Ct*massp21*Cab - Ct*tilde(Vp)*massp11*Cab + Ct*tilde(Pp)*Cab - Ctdot*massp21*Cab
    Mp_Ω = -tilde(ω)*Ct*massp22*Cab - Ct*tilde(Vp)*massp12*Cab - Ctdot*massp22*Cab

    # return result
    return Fp_Vdot, Fp_Ωdot, Fp_V, Fp_Ω, Mp_Vdot, Mp_Ωdot, Mp_V, Mp_Ω
end

@inline function newmark_point_mass_jacobian(massp11, massp12, massp21, massp22, u, Vp, Pp, Hp, 
    Ppdot, Hpdot, C, Ct, Ctdot, Cab, ω, a, α, gvec, dt, C_θ1, C_θ2, C_θ3, Ct_θ1, Ct_θ2, 
    Ct_θ3, Ctdot_θ1, Ctdot_θ2, Ctdot_θ3, Ctdot_θdot1, Ctdot_θdot2, Ctdot_θdot3)

    Fp_u, Mp_u, Fp_θ, Mp_θ, Fp_V, Mp_V, Fp_Ω, Mp_Ω = steady_state_point_mass_jacobian(massp11, massp12, 
        massp21, massp22, C, Ct, Cab, u, Vp, Pp, Hp, ω, a, α, gvec, C_θ1, C_θ2, C_θ3, Ct_θ1, Ct_θ2, Ct_θ3)

    # force and moment due to local frame acceleration
    Fp_θ -= mul3(Ctdot_θ1, Ctdot_θ2, Ctdot_θ3, Pp) + 2/dt*mul3(Ctdot_θdot1, Ctdot_θdot2, Ctdot_θdot3, Pp) + mul3(Ct_θ1, Ct_θ2, Ct_θ3, Ppdot)
    Mp_θ -= mul3(Ctdot_θ1, Ctdot_θ2, Ctdot_θ3, Hp) + 2/dt*mul3(Ctdot_θdot1, Ctdot_θdot2, Ctdot_θdot3, Hp) + mul3(Ct_θ1, Ct_θ2, Ct_θ3, Hpdot)

    Fp_V -= Ctdot*massp11*Cab + 2/dt*Ct*massp11*Cab
    Mp_V -= Ctdot*massp21*Cab + 2/dt*Ct*massp21*Cab

    Fp_Ω -= Ctdot*massp12*Cab + 2/dt*Ct*massp12*Cab
    Mp_Ω -= Ctdot*massp22*Cab + 2/dt*Ct*massp22*Cab    

    return Fp_u, Mp_u, Fp_θ, Mp_θ, Fp_V, Mp_V, Fp_Ω, Mp_Ω
end

@inline function dynamic_point_mass_jacobian(massp11, massp12, massp21, massp22, C, Ct, 
    Ctdot, Cab, u, Vp, Pp, Hp, Ppdot, Hpdot, ω, a, α, gvec, C_θ1, C_θ2, C_θ3, Ct_θ1, Ct_θ2, 
    Ct_θ3, Ctdot_θ1, Ctdot_θ2, Ctdot_θ3)

    Fp_u, Mp_u, Fp_θ, Mp_θ, Fp_V, Mp_V, Fp_Ω, Mp_Ω = steady_state_point_mass_jacobian(massp11, massp12, 
        massp21, massp22, C, Ct, Cab, u, Vp, Pp, Hp, ω, a, α, gvec, C_θ1, C_θ2, C_θ3, Ct_θ1, 
        Ct_θ2, Ct_θ3)
    
    # force and moment due to local frame acceleration
    Fp_θ -= mul3(Ctdot_θ1, Ctdot_θ2, Ctdot_θ3, Pp) + mul3(Ct_θ1, Ct_θ2, Ct_θ3, Ppdot)
    Mp_θ -= mul3(Ctdot_θ1, Ctdot_θ2, Ctdot_θ3, Hp) + mul3(Ct_θ1, Ct_θ2, Ct_θ3, Hpdot)

    Fp_V -= Ctdot*massp11*Cab
    Mp_V -= Ctdot*massp21*Cab

    Fp_Ω -= Ctdot*massp12*Cab
    Mp_Ω -= Ctdot*massp22*Cab

    return Fp_u, Mp_u, Fp_θ, Mp_θ, Fp_V, Mp_V, Fp_Ω, Mp_Ω
end

@inline function point_mass_rate_jacobian_properties(point_mass, Cab, V, Ω)

    massp11, massp12, massp21, massp22 = static_point_mass_properties(point_mass)
    
    # velocities in point mass reference frame
    Vp = Cab*V
    Ωp = Cab*Ω
    
    # point mass linear and angular momentum
    Pp = massp11*Vp + massp12*Ωp
    Hp = massp21*Vp + massp22*Ωp
    
    return massp11, massp12, massp21, massp22, Pp, Hp
end

@inline function point_mass_rate_jacobian(massp11, massp12, massp21, massp22, Ct, Cab, Pp, 
    Hp, Ctdot_θdot1, Ctdot_θdot2, Ctdot_θdot3)

    # force and moment due to local frame acceleration
    Fp_θdot = -mul3(Ctdot_θdot1, Ctdot_θdot2, Ctdot_θdot3, Pp)
    Mp_θdot = -mul3(Ctdot_θdot1, Ctdot_θdot2, Ctdot_θdot3, Hp)

    Fp_Vdot = -Ct*massp11*Cab
    Mp_Vdot = -Ct*massp21*Cab

    Fp_Ωdot = -Ct*massp12*Cab
    Mp_Ωdot = -Ct*massp22*Cab

    return Fp_θdot, Fp_Vdot, Fp_Ωdot, Mp_θdot, Mp_Vdot, Mp_Ωdot
end