"""
    body_frame_displacement(x)

Extract the linear and angular displacement of the body frame from the state vector.
"""
function body_frame_displacement(x)

    ub = SVector(x[1], x[2], x[3])
    θb = SVector(x[4], x[5], x[6])

    return ub, θb
end

"""
    body_frame_velocity(x)

Extract the linear and angular velocity of the body frame from the state vector.
"""
function body_frame_velocity(x)

    vb = SVector(x[7], x[8], x[9])
    ωb = SVector(x[10], x[11], x[12])

    return vb, ωb
end

"""
    body_frame_acceleration(x)

Extract the linear and angular acceleration of the body frame from the state vector.
"""
function body_frame_acceleration(x)

    ab = SVector(x[13], x[14], x[15])
    αb = SVector(x[16], x[17], x[18])

    return ab, αb
end

"""
    steady_state_body_residual!(resid, x, icol_accel, ub_p, θb_p, vb_p, ωb_p, ab_p, αb_p)

Calculate and insert the residual entries corresponding to the motion of the body frame 
for a steady state analysis into the system residual vector.
"""
function steady_state_body_residual!(resid, x, icol_accel, ub_p, θb_p, vb_p, ωb_p, 
    ab_p, αb_p)
    
    # update prescribed accelerations (with accelerations for satisying prescribed conditions)
    ab_p, αb_p = prescribed_body_frame_acceleration(x, icol_accel, ab_p, αb_p)

    # extract body states
    ub, θb = body_frame_displacement(x)
    vb, ωb = body_frame_velocity(x)
    ab, αb = body_frame_acceleration(x)

    # construct residuals
    ru = ub - ub_p
    rθ = θb - θb_p
    rv = vb - vb_p
    rω = ωb - ωb_p
    ra = ab - ab_p
    rα = αb - αb_p

    # insert residuals into the residual vector
    resid[1:3] = ru
    resid[4:6] = rθ
    resid[7:9] = rv
    resid[10:12] = rω
    resid[13:15] = ra
    resid[16:18] = rα
    
    return resid
end

"""
    initial_condition_body_residual!(resid, x, icol_accel, ub_p, θb_p, vb_p, ωb_p, ab_p, αb_p)

Calculate and insert the residual entries corresponding to the motion of the body frame 
for the initialization of a time domain analysis into the system residual vector.
"""
const initial_condition_body_residual! = steady_state_body_residual!

"""
    newmark_body_residual!(resid, x, icol_accel, ab_p, αb_p, 
        ubdot_init, θbdot_init, vbdot_init, ωbdot_init)

Calculate and insert the residual entries corresponding to the motion of the body frame 
for a newmark-scheme time marching analysis into the system residual vector.
"""
function newmark_body_residual!(resid, x, icol_accel, ab_p, αb_p, ubdot_init, θbdot_init, 
    vbdot_init, ωbdot_init, dt)
    
    # update prescribed accelerations (with accelerations for satisying prescribed conditions)
    ab_p, αb_p = prescribed_body_frame_acceleration(x, icol_accel, ab_p, αb_p)

    # extract body states
    ub, θb = body_frame_displacement(x)
    vb, ωb = body_frame_velocity(x)
    ab, αb = body_frame_acceleration(x)

    # extract body rates
    ubdot = 2/dt*ub - SVector{3}(ubdot_init)
    θbdot = 2/dt*θb - SVector{3}(θbdot_init)
    vbdot = 2/dt*vb - SVector{3}(vbdot_init)
    ωbdot = 2/dt*ωb - SVector{3}(ωbdot_init)

    # rotation parameter matrices
    C = get_C(θb)
    Qinv = get_Qinv(θb)

    # construct residuals
    ru = vb - ubdot
    rθ = Qinv*C*ωb - θbdot
    rv = ab - vbdot
    rω = αb - ωbdot
    ra = ab - ab_p
    rα = αb - αb_p

    # insert residuals into the residual vector
    resid[1:3] = ru
    resid[4:6] = rθ
    resid[7:9] = rv
    resid[10:12] = rω
    resid[13:15] = ra
    resid[16:18] = rα
    
    return resid
end

"""
    dynamic_body_residual!(resid, dx, x, icol_accel, ab_p, αb_p)

Calculate and insert the residual entries corresponding to the motion of the body frame 
for a dynamic analysis into the system residual vector.
"""
function dynamic_body_residual!(resid, dx, x, icol_accel, ab_p, αb_p)
    
    # update prescribed accelerations (with accelerations for satisying prescribed conditions)
    ab_p, αb_p = prescribed_body_frame_acceleration(x, icol_accel, ab_p, αb_p)

    # extract body states
    ub, θb = body_frame_displacement(x)
    vb, ωb = body_frame_velocity(x)
    ab, αb = body_frame_acceleration(x)

    # extract body rates
    ubdot, θbdot = body_frame_displacement(dx)
    vbdot, ωbdot = body_frame_velocity(dx)

    # rotation parameter matrices
    C = get_C(θb)
    Qinv = get_Qinv(θb)

    # construct residuals
    ru = vb - ubdot
    rθ = Qinv*C*ωb - θbdot
    rv = ab - vbdot
    rω = αb - ωbdot
    ra = ab - ab_p
    rα = αb - αb_p

    # insert residuals into the residual vector
    resid[1:3] = ru
    resid[4:6] = rθ
    resid[7:9] = rv
    resid[10:12] = rω
    resid[13:15] = ra
    resid[16:18] = rα
    
    return resid
end

"""
    steady_state_body_jacobian!(jacob, x, icol_accel, ub_p, θb_p, vb_p, ωb_p, ab_p, αb_p)

Calculate and insert the jacobian entries corresponding to the motion of the body frame 
for a steady state analysis into the system jacobian matrix.
"""
function steady_state_body_jacobian!(jacob, x, icol_accel, ub_p, θb_p, vb_p, ωb_p, ab_p, αb_p)

    for i = 1:18
        jacob[i,i] = 1 
    end

    icol_accel[1] > 0 && setindex!(jacob, -1, 13, icol_accel[1])
    icol_accel[2] > 0 && setindex!(jacob, -1, 14, icol_accel[2])
    icol_accel[3] > 0 && setindex!(jacob, -1, 15, icol_accel[3])
    icol_accel[4] > 0 && setindex!(jacob, -1, 16, icol_accel[4])
    icol_accel[5] > 0 && setindex!(jacob, -1, 17, icol_accel[5])
    icol_accel[6] > 0 && setindex!(jacob, -1, 18, icol_accel[6])
    
    return jacob
end

"""
    initial_condition_body_jacobian!(jacob, x, icol_accel, ub_p, θb_p, vb_p, ωb_p, ab_p, αb_p)

Calculate and insert the jacobian entries corresponding to the motion of the body frame 
for the initialization of a time domain analysis into the system jacobian vector.
"""
const initial_condition_body_jacobian! = steady_state_body_jacobian!

"""
    newmark_body_jacobian!(jacob, x, icol_accel, ab_p, αb_p, 
        ubdot_init, θbdot_init, vbdot_init, ωbdot_init)

Calculate and insert the jacobian entries corresponding to the motion of the body frame 
for a steady state analysis into the system jacobian matrix.
"""
function newmark_body_jacobian!(jacob, x, icol_accel, ab_p, αb_p, 
    ubdot_init, θbdot_init, vbdot_init, ωbdot_init, dt)

    # update prescribed accelerations (with accelerations for satisying prescribed conditions)
    ab_p, αb_p = prescribed_body_frame_acceleration(x, icol_accel, ab_p, αb_p)

    # extract body states
    ub, θb = body_frame_displacement(x)
    vb, ωb = body_frame_velocity(x)
    ab, αb = body_frame_acceleration(x)

    # extract body rates
    ubdot_ub = 2/dt*I3
    θbdot_θb = 2/dt*I3
    vbdot_vb = 2/dt*I3
    ωbdot_ωb = 2/dt*I3

    # rotation parameter matrices
    C = get_C(θb)
    Qinv = get_Qinv(θb)
    C_θ1, C_θ2, C_θ3 = get_C_θ(C, θb)
    Qinv_θ1, Qinv_θ2, Qinv_θ3 = get_Qinv_θ(θb)

    # construct residuals
    ru_ub = -ubdot_ub
    ru_vb = I3
    
    rθ_θb = mul3(Qinv_θ1, Qinv_θ2, Qinv_θ3, C*ωb) + Qinv*mul3(C_θ1, C_θ2, C_θ3, ωb) - θbdot_θb
    rθ_ωb = Qinv*C

    rv_vb = -vbdot_vb
    rv_ab = I3

    rω_ωb = -ωbdot_ωb
    rω_αb = I3

    ra_ab = I3
    rα_αb = I3

    # insert jacobian entries into the jacobian matrix

    jacob[1:3,1:3] = ru_ub
    jacob[1:3,7:9] = ru_vb

    jacob[4:6,4:6] = rθ_θb
    jacob[4:6,10:12] = rθ_ωb

    jacob[7:9,7:9] = rv_vb
    jacob[7:9,13:15] = rv_ab

    jacob[10:12,10:12] = rω_ωb
    jacob[10:12,16:18] = rω_αb

    jacob[13:15,13:15] = ra_ab

    jacob[16:18,16:18] = rα_αb

    # add influence of prescribed acceleration state variables
    icol_accel[1] > 0 && setindex!(jacob, -1, 13, icol_accel[1])
    icol_accel[2] > 0 && setindex!(jacob, -1, 14, icol_accel[2])
    icol_accel[3] > 0 && setindex!(jacob, -1, 15, icol_accel[3])

    icol_accel[4] > 0 && setindex!(jacob, -1, 16, icol_accel[4])
    icol_accel[5] > 0 && setindex!(jacob, -1, 17, icol_accel[5])
    icol_accel[6] > 0 && setindex!(jacob, -1, 18, icol_accel[6])

    return jacob
end

function dynamic_body_jacobian!(jacob, dx, x, icol_accel, ab_p, αb_p)

    # update prescribed accelerations (with accelerations for satisying prescribed conditions)
    ab_p, αb_p = prescribed_body_frame_acceleration(x, icol_accel, ab_p, αb_p)

    # extract body states
    ub, θb = body_frame_displacement(x)
    vb, ωb = body_frame_velocity(x)
    ab, αb = body_frame_acceleration(x)

    # rotation parameter matrices
    C = get_C(θb)
    Qinv = get_Qinv(θb)
    C_θ1, C_θ2, C_θ3 = get_C_θ(C, θb)
    Qinv_θ1, Qinv_θ2, Qinv_θ3 = get_Qinv_θ(θb)

    # construct residuals
    ru_vb = I3
    
    rθ_θb = mul3(Qinv_θ1, Qinv_θ2, Qinv_θ3, C*ωb) + Qinv*mul3(C_θ1, C_θ2, C_θ3, ωb)
    rθ_ωb = Qinv*C

    rv_ab = I3

    rω_αb = I3

    ra_ab = I3
    rα_αb = I3

    # insert jacobian entries into the jacobian matrix

    jacob[1:3,7:9] = ru_vb

    jacob[4:6,4:6] = rθ_θb
    jacob[4:6,10:12] = rθ_ωb

    jacob[7:9,13:15] = rv_ab

    jacob[10:12,16:18] = rω_αb

    jacob[13:15,13:15] = ra_ab

    jacob[16:18,16:18] = rα_αb

    # add influence of prescribed acceleration state variables
    icol_accel[1] > 0 && setindex!(jacob, -1, 13, icol_accel[1])
    icol_accel[2] > 0 && setindex!(jacob, -1, 14, icol_accel[2])
    icol_accel[3] > 0 && setindex!(jacob, -1, 15, icol_accel[3])

    icol_accel[4] > 0 && setindex!(jacob, -1, 16, icol_accel[4])
    icol_accel[5] > 0 && setindex!(jacob, -1, 17, icol_accel[5])
    icol_accel[6] > 0 && setindex!(jacob, -1, 18, icol_accel[6])

    return jacob
end

function mass_matrix_body_jacobian!(jacob)
    
    # define body residuals
    jacob[1:3,1:3] = -I3
    jacob[4:6,4:6] = -I3
    jacob[7:9,7:9] = -I3
    jacob[10:12,10:12] = -I3
    
    return jacob
end
