# unit vectors and identity matrix
const e1 = SVector(1, 0, 0)
const e2 = SVector(0, 1, 0)
const e3 = SVector(0, 0, 1)
const I3 = @SMatrix [1 0 0; 0 1 0; 0 0 1]
const Z3 = @SMatrix zeros(3,3)

"""
    tilde(x)

Construct the cross product operator matrix
"""
@inline tilde(x) = @SMatrix [0 -x[3] x[2]; x[3] 0 -x[1]; -x[2] x[1] 0]

"""
    wiener_milenkovic(c)

Construct a Wiener-Milenkovic transformation matrix, given the three Wiener-Milenkovic
parameters in `c`.

Note that the corresponding rotation matrix is the transpose of this
transformation matrix.
"""
@inline function wiener_milenkovic(c)

    c0 = 2 - c'*c/8

    C = 1/(4-c0)^2*(@SMatrix [
        c0^2 + c[1]^2 - c[2]^2 - c[3]^2  2*(c[1]*c[2] + c0*c[3])          2*(c[1]*c[3]-c0*c[2]);
        2*(c[1]*c[2] - c0*c[3])          c0^2 - c[1]^2 + c[2]^2 - c[3]^2  2*(c[2]*c[3] + c0*c[1]);
        2*(c[1]*c[3] + c0*c[2])          2*(c[2]*c[3] - c0*c[1])          c0^2 - c[1]^2 - c[2]^2 + c[3]^2
        ]
    )

    return C
end

"""
    transform_properties(K, T)

Applies the transformation `T` to the stiffness or mass matrix `K`
"""
transform_properties(X, T) = [T' Z3; Z3 T'] * X * [T Z3; Z3 T]

"""
    rotation_parameter_scaling(θ)

Extract a scaling parameter which may be multiplied by the angular parameters
to yield the Wiener-Milenkovic rotation parameters.  Use of this scaling
parameter allows deflections greater than 360 degrees.
"""
@inline function rotation_parameter_scaling(θ)

    θ = real(θ)

    θnorm = sqrt(θ'*θ)

    scaling = 1

    m = div(div(θnorm, 4) + 1, 2) # number of times to add/subtract

    if θnorm > 4
        scaling -= 8*m/θnorm
    end

    return scaling
end

@inline function rotation_parameter_scaling_θ(θ)
    
    return rotation_parameter_scaling_θ(rotation_parameter_scaling(θ), θ)
end

@inline function rotation_parameter_scaling_θ(scaling, θ)
    
    return ifelse(isone(scaling), zero(θ), (1 - scaling)*θ/(θ'*θ))
end

@inline function rotation_parameter_scaling_θ_θ(θ)
    
    return rotation_parameter_scaling_θ_θ(rotation_parameter_scaling(θ), θ)
end

@inline function rotation_parameter_scaling_θ_θ(scaling, θ)
    
    return ifelse(isone(scaling), (@SMatrix zeros(eltype(θ),3,3)), (1 - scaling)/(θ'*θ)*(I - 3*θ*θ'/(θ'*θ)))
end
        
"""
    get_C(θ)

Returns the transformation matrix `C` given the three angular parameters in `θ`.
"""
@inline function get_C(θ)

    scaling = rotation_parameter_scaling(θ)

    c = scaling*θ

    return wiener_milenkovic(c)
end

"""
    get_C_θ([C, ] θ)

Calculate the derivative of the Wiener-Milenkovic transformation matrix `C` with
respect to each of the rotation parameters in `θ`.
"""
@inline get_C_θ(c) = get_C_θ(get_C(c), c)

@inline function get_C_θ(C, θ)

    scaling = rotation_parameter_scaling(θ)
    scaling_θ = rotation_parameter_scaling_θ(scaling, θ)

    c = scaling*θ
    c_θ = scaling_θ*θ' + scaling*I3

    c0 = 2 - c'*c/8
    c0_c = -c/4

    tmp = 1/(4-c0)^2
    tmp_c = 2/(4-c0)^3*c0_c

    C_c1 = tmp_c[1]*C/tmp + tmp*(@SMatrix [
        2*c0*c0_c[1] + 2*c[1]    2*(c[2] + c0_c[1]*c[3]) 2*(c[3] - c0_c[1]*c[2]);
         2*(c[2] - c0_c[1]*c[3])  2*c0*c0_c[1] - 2*c[1]   2*(c0_c[1]*c[1] + c0);
         2*(c[3] + c0_c[1]*c[2]) -2*(c0_c[1]*c[1] + c0)   2*c0*c0_c[1] - 2*c[1]
        ]
    )

    C_c2 = tmp_c[2]*C/tmp + tmp*(@SMatrix [
        2*c0*c0_c[2] - 2*c[2]   2*(c[1] + c0_c[2]*c[3]) -2*(c0_c[2]*c[2] + c0);
        2*(c[1] - c0_c[2]*c[3]) 2*c0*c0_c[2] + 2*c[2]    2*(c[3] + c0_c[2]*c[1]);
         2*(c0_c[2]*c[2] + c0)   2*(c[3] - c0_c[2]*c[1])  2*c0*c0_c[2] - 2*c[2]
        ]
    )

    C_c3 = tmp_c[3]*C/tmp + tmp*(@SMatrix [
         2*c0*c0_c[3] - 2*c[3]   2*(c0_c[3]*c[3] + c0)    2*(c[1] - c0_c[3]*c[2]);
        -2*(c0_c[3]*c[3] + c0)   2*c0*c0_c[3] - 2*c[3]    2*(c[2] + c0_c[3]*c[1]);
         2*(c[1] + c0_c[3]*c[2]) 2*(c[2] - c0_c[3]*c[1])  2*c0*c0_c[3] + 2*c[3]
        ]
    )

    C_θ1 = C_c1*c_θ[1,1] + C_c2*c_θ[2,1] + C_c3*c_θ[3,1]
    C_θ2 = C_c1*c_θ[1,2] + C_c2*c_θ[2,2] + C_c3*c_θ[3,2]
    C_θ3 = C_c1*c_θ[1,3] + C_c2*c_θ[2,3] + C_c3*c_θ[3,3]

    return C_θ1, C_θ2, C_θ3
end

"""
    get_Q(θ)

Calculate the matrix Q as defined in the paper "Geometrically nonlinear analysis
of composite beams using Wiener-Milenković parameters" by Qi Wang and Wenbin Yu
given the rotational parameters in `θ`.
"""
@inline function get_Q(θ)

    scaling = rotation_parameter_scaling(θ)

    c = scaling*θ

    c0 = 2 - c'*c/8

    return 1/(4-c0)^2*((4 - 1/4*c'*c)*I - 2*tilde(c) + 1/2*c*c')
end

"""
    get_Q_θ(θ)
    get_Q_θ(Q, θ)

Calculate the derivative of the matrix `Q` with respect to each of the rotation
parameters in `θ`.
"""
@inline get_Q_θ(θ) = get_Q_θ(get_Q(θ), θ)

@inline function get_Q_θ(Q, θ)

    scaling = rotation_parameter_scaling(θ)
    scaling_θ = rotation_parameter_scaling_θ(scaling, θ)

    c = scaling*θ
    c_θ = scaling_θ*θ' + scaling*I3

    c0 = 2 - c'*c/8
    c0_c = -c/4

    tmp = 1/(4-c0)^2
    tmp_c = 2/(4-c0)^3*c0_c

    Q_c1 = tmp_c[1]*Q/tmp + tmp*(-c[1]/2*I - 2*tilde(e1) + (e1*c' + c*e1')/2)

    Q_c2 = tmp_c[2]*Q/tmp + tmp*(-c[2]/2*I - 2*tilde(e2) + (e2*c' + c*e2')/2)

    Q_c3 = tmp_c[3]*Q/tmp + tmp*(-c[3]/2*I - 2*tilde(e3) + (e3*c' + c*e3')/2)

    Q_θ1 = Q_c1*c_θ[1,1] + Q_c2*c_θ[2,1] + Q_c3*c_θ[3,1]
    Q_θ2 = Q_c1*c_θ[1,2] + Q_c2*c_θ[2,2] + Q_c3*c_θ[3,2]
    Q_θ3 = Q_c1*c_θ[1,3] + Q_c2*c_θ[2,3] + Q_c3*c_θ[3,3]

    return Q_θ1, Q_θ2, Q_θ3
end

"""
    get_ΔQ(θ, Δθ [, Q])

Calculate the matrix `ΔQ` for structural damping calculations
"""
@inline function get_ΔQ(θ, Δθ, Q=get_Q(θ))

    scaling = rotation_parameter_scaling(θ)
    scaling_θ = rotation_parameter_scaling_θ(scaling, θ)

    c = scaling*θ
    c_θ = scaling_θ*θ' + scaling*I3

    c0 = 2 - c'*c/8
    c0_θ = -c'/4*c_θ

    tmp1 = 1/(4-c0)^2
    tmp1_θ = 2/(4-c0)^3*c0_θ

    tmp2 = -Δθ*c' + 4*tilde(Δθ) + c'*Δθ*I3 + c*Δθ'

    ΔQ = Q*Δθ*tmp1_θ/tmp1 + 1/2*tmp1*tmp2*(scaling*I3 + scaling_θ*θ')

    return ΔQ
end

"""
    get_ΔQ_θ(θ, Δθ, [Q, Q_θ1, Q_θ2, Q_θ3])

Calculate the derivative of the matrix `ΔQ` with respect to each of the rotation
parameters in `θ`.
"""
get_ΔQ_θ(θ, Δθ, Q=get_Q(θ)) = get_ΔQ_θ(θ, Δθ, Q, get_Q_θ(Q, θ)...)

@inline function get_ΔQ_θ(θ, Δθ, Q, Q_θ1, Q_θ2, Q_θ3)

    # calculate scaling factor
    scaling = GXBeam.rotation_parameter_scaling(θ)
    scaling_θ = GXBeam.rotation_parameter_scaling_θ(scaling, θ)
    scaling_θ_θ = GXBeam.rotation_parameter_scaling_θ_θ(scaling, θ)

    # scale rotation parameters
    c = scaling*θ    
    c_θ = scaling_θ*θ' + scaling*I3

    # calculate c0 constant
    c0 = 2 - c'*c/8
    c0_θ = -c'/4*c_θ
    c0_θ_θ = -I/4*c_θ'*c_θ - 1/2*c*scaling_θ' - 1/4*c'*θ*scaling_θ_θ

    # calculate tmp1 constant
    tmp1 = 1/(4-c0)^2    
    tmp1_θ = 2/(4-c0)^3*c0_θ
    tmp1_θ_θ = 6/(4-c0)^4*c0_θ'*c0_θ + 2/(4-c0)^3*c0_θ_θ

    # calculate tmp2 matrix
    tmp2 = -Δθ*c' + 4*tilde(Δθ) + c'*Δθ*I3 + c*Δθ'
    tmp2_θ1 = -Δθ*c_θ[1,:]' + c_θ[1,:]'*Δθ*I3 + c_θ[1,:]*Δθ'
    tmp2_θ2 = -Δθ*c_θ[2,:]' + c_θ[2,:]'*Δθ*I3 + c_θ[2,:]*Δθ'
    tmp2_θ3 = -Δθ*c_θ[3,:]' + c_θ[3,:]'*Δθ*I3 + c_θ[3,:]*Δθ'

    # calculate ΔQ
    ΔQ_θ1 = (Q_θ1*Δθ*tmp1_θ)/tmp1 + (Q*Δθ*tmp1_θ_θ[1,:]')/tmp1 - (Q*Δθ*tmp1_θ)/tmp1^2*tmp1_θ[1] +
        1/2*tmp1_θ[1]*tmp2*(scaling*I3 + scaling_θ*θ') + 1/2*tmp1*tmp2_θ1*(scaling*I3 + scaling_θ*θ') +
        1/2*tmp1*tmp2*(scaling_θ[1]*I3 + scaling_θ_θ[1,:]*θ' + scaling_θ*e1')
    ΔQ_θ2 = (Q_θ2*Δθ*tmp1_θ)/tmp1 + (Q*Δθ*tmp1_θ_θ[2,:]')/tmp1 - (Q*Δθ*tmp1_θ)/tmp1^2*tmp1_θ[2] +
        1/2*tmp1_θ[2]*tmp2*(scaling*I3 + scaling_θ*θ') + 1/2*tmp1*tmp2_θ2*(scaling*I3 + scaling_θ*θ') +
        1/2*tmp1*tmp2*(scaling_θ[2]*I3 + scaling_θ_θ[2,:]*θ' + scaling_θ*e2')
    ΔQ_θ3 = (Q_θ3*Δθ*tmp1_θ)/tmp1 + (Q*Δθ*tmp1_θ_θ[3,:]')/tmp1 - (Q*Δθ*tmp1_θ)/tmp1^2*tmp1_θ[3] +
        1/2*tmp1_θ[3]*tmp2*(scaling*I3 + scaling_θ*θ') + 1/2*tmp1*tmp2_θ3*(scaling*I3 + scaling_θ*θ') +
        1/2*tmp1*tmp2*(scaling_θ[3]*I3 + scaling_θ_θ[3,:]*θ' + scaling_θ*e3')

    return ΔQ_θ1, ΔQ_θ2, ΔQ_θ3
end

"""
    get_Qinv(θ)

Calculate the matrix inverse `Qinv` as defined in the paper "Geometrically
nonlinear analysis of composite beams using Wiener-Milenković parameters" by
Qi Wang and Wenbin Yu given the rotational parameters in `θ`.
"""
@inline function get_Qinv(θ)

    scaling = rotation_parameter_scaling(θ)

    c = scaling*θ

    Qinv = (1 - 1/16*c'*c)*I3 + 1/2*tilde(c) + 1/8*c*c'

    return Qinv
end

"""
    get_Qinv_θ(θ)

Calculate the derivative of the matrix inverse `Qinv` with respect to each of the
rotation parameters in `θ`.
"""
@inline function get_Qinv_θ(θ)

    scaling = rotation_parameter_scaling(θ)
    scaling_θ = rotation_parameter_scaling_θ(scaling, θ)

    c = scaling*θ
    c_θ = scaling_θ*θ' + scaling*I3

    Qinv_c1 = -c[1]/8*I3 + 1/2*tilde(e1) + 1/8*(e1*c' + c*e1')
    Qinv_c2 = -c[2]/8*I3 + 1/2*tilde(e2) + 1/8*(e2*c' + c*e2')
    Qinv_c3 = -c[3]/8*I3 + 1/2*tilde(e3) + 1/8*(e3*c' + c*e3')

    Qinv_θ1 = Qinv_c1*c_θ[1,1] + Qinv_c2*c_θ[2,1] + Qinv_c3*c_θ[3,1]
    Qinv_θ2 = Qinv_c1*c_θ[1,2] + Qinv_c2*c_θ[2,2] + Qinv_c3*c_θ[3,2]
    Qinv_θ3 = Qinv_c1*c_θ[1,3] + Qinv_c2*c_θ[2,3] + Qinv_c3*c_θ[3,3]

    return Qinv_θ1, Qinv_θ2, Qinv_θ3
end

"""
    mul3(A_1, A_2, A_3, b)

Return the product of a 3x3x3 tensor represented by `A_1`, `A_2`, and `A_3` with
the vector `b`.
"""
@inline mul3(A_1, A_2, A_3, b) = hcat(A_1*b, A_2*b, A_3*b)

"""
    gauss_quadrature(f, a, b)

Default gauss-quadrature function used for integrating distributed loads.
"""
@inline function gauss_quadrature(f, a, b)
    h = b - a
    c = (a + b)/2
    x = h/2*GAUSS_NODES .+ c
    return h/2*GAUSS_WEIGHTS'*f.(x)
end

"""
    rotate(xyz, r, theta)

Rotate the vectors in `xyz` about point `r` using the Wiener-Milenkovic
parameters in `theta`.
"""
rotate(xyz, r, theta) = rotate!(copy(xyz), r, theta)

"""
    rotate!(xyz, r, theta)

Pre-allocated version of [`rotate`](@ref)
"""
function rotate!(xyz, r, theta)
    # reshape cross section points (if necessary)
    rxyz = reshape(xyz, 3, :)
    # convert inputs to static arrays
    r = SVector{3}(r)
    theta = SVector{3}(theta)
    # create rotation matrix
    Ct = wiener_milenkovic(theta)'
    # rotate each point
    for ipt = 1:size(rxyz, 2)
        p = SVector(rxyz[1,ipt], rxyz[2,ipt], rxyz[3,ipt])
        rxyz[:,ipt] .= Ct*(p - r) + r
    end
    # return the result
    return xyz
end

"""
    translate(xyz, u)

Translate the points in `xyz` by the displacements in `u`.
"""
translate(xyz, u) = translate!(copy(xyz), u)

"""
    translate!(xyz, u)

Pre-allocated version of [`translate`](@ref)
"""
function translate!(xyz, u)
    # reshape cross section points (if necessary)
    rxyz = reshape(xyz, 3, :)
    # convert inputs to static arrays
    u = SVector{3}(u)
    # translate each point
    for ipt = 1:size(rxyz, 2)
        p = SVector(rxyz[1,ipt], rxyz[2,ipt], rxyz[3,ipt])
        rxyz[:,ipt] .= p .+ u
    end
    # return the result
    return xyz
end

"""
    deform_cross_section(xyz, r, u, theta)

Rotate the points in `xyz` (of shape (3, :)) about point `r` using
the Wiener-Milenkovic parameters in `theta`, then translate the points by the
displacements in `u`.
"""
deform_cross_section(xyz, r, u, theta) = deform_cross_section!(copy(xyz), r, u, theta)

"""
    deform_cross_section!(xyz, r, u, theta)

Pre-allocated version of [`deform_cross_section`](@ref)
"""
deform_cross_section!(xyz, r, u, theta) = translate!(rotate!(xyz, r, theta), u)