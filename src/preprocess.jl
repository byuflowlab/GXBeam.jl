"""
    set_state!([x,] system::StaticSystem, prescribed_conditions; kwargs...)

Set the state variables in `system` (or in the vector `x`) to the provided values.

# Keyword Arguments
- `u`: Vector containing the linear displacement of each point.
- `theta`: Vector containing the angular displacement of each point.
- `F`: Vector containing the externally applied forces acting on each point
- `M`: Vector containing the externally applied moments acting on each point
- `Fi`: Vector containing internal forces for each beam element
- `Mi`: Vector containing internal moments for each beam element
"""
set_state!(system::StaticSystem, prescribed_conditions; kwargs...)

"""
    set_state!([x,] system::DynamicSystem, prescribed_conditions; kwargs...)

Set the state variables in `system` (or in the vector `x`) to the provided values.

# Keyword Arguments
 - `u`: Vector containing the linear displacement of each point.
 - `theta`: Vector containing the angular displacement of each point.
 - `V`: Vector containing the linear velocity of each point.
 - `Omega` Vector containing the angular velocity of each point
 - `F`: Vector containing the externally applied forces acting on each point
 - `M`: Vector containing the externally applied moments acting on each point
 - `Fi`: Vector containing internal forces for each beam element (in the deformed 
        element frame)
 - `Mi`: Vector containing internal moments for each beam element (in the deformed 
        element frame)
 - `linear_displacement`: Linear displacement of the body frame.
 - `angular_displacement`: Angular displacement of the body frame (Wiener Milenković)
 - `linear_velocity`: Linear velocity of the body frame
 - `angular_velocity`: Angular velocity of the body frame
 - `linear_acceleration`: Linear acceleration of the body frame
 - `angular_acceleration`: Angular acceleration of the body frame
"""
set_state!(system::DynamicSystem, prescribed_conditions; kwargs...)

"""
    set_state!([x,] system::ExpandedSystem, prescribed_conditions; kwargs...)

Set the state variables in `system` (or in the vector `x`) to the provided values.

# Keyword Arguments
 - `u`: Vector containing the linear displacement of each point.
 - `theta`: Vector containing the angular displacement of each point.
 - `F`: Vector containing the externally applied forces acting on each point
 - `M`: Vector containing the externally applied moments acting on each point
 - `F1`: Vector containing resultant forces at the start of each beam element (in the 
    deformed element frame)
 - `M1`: Vector containing resultant moments at the start of each beam element (in the 
    deformed element frame)
 - `F2`: Vector containing resultant forces at the end of each beam element (in the 
    deformed element frame)
 - `M2`: Vector containing resultant moments at the end of each beam element (in the 
    deformed element frame)
 - `V_p`: Vector containing the linear velocity of each point in a deformed 
    point reference frame.
 - `Omega_p` Vector containing the angular velocity of each point in a deformed 
    point reference frame.
 - `V_e`: Vector containing the linear velocity of each beam element in the deformed
    beam element reference frame.
 - `Omega_e` Vector containing the angular velocity of each beam element in the deformed
    beam element reference frame.
 - `linear_displacement`: Linear displacement of the body frame.
 - `angular_displacement`: Angular displacement of the body frame (Wiener Milenković)
 - `linear_velocity`: Linear velocity of the body frame
 - `angular_velocity`: Angular velocity of the body frame
 - `linear_acceleration`: Linear acceleration of the body frame
 - `angular_acceleration`: Angular acceleration of the body frame
"""
set_state!(system::ExpandedSystem, prescribed_conditions; kwargs...)

function set_state!(system, prescribed_conditions; kwargs...)
    x = set_state!(system.x, system, prescribed_conditions; kwargs...)
    return system
end

function set_state!(x, system, prescribed_conditions; u = nothing, theta = nothing, 
    V = nothing, Omega = nothing, F = nothing, M = nothing, Fi = nothing, Mi = nothing,
    F1 = nothing, M1 = nothing, F2 = nothing, M2 = nothing, 
    V_p = nothing, Omega_p = nothing, V_e = nothing, Omega_e = nothing,
    linear_displacement=nothing, angular_displacement=nothing, 
    linear_velocity=nothing, angular_velocity=nothing,
    linear_acceleration=nothing, angular_acceleration=nothing) 

    if !isnothing(u)
        for ipoint = 1:length(u)
            set_linear_displacement!(x, system, prescribed_conditions, u[ipoint], ipoint)
        end
    end

    if !isnothing(theta)
        for ipoint = 1:length(theta)
            set_angular_displacement!(x, system, prescribed_conditions, theta[ipoint], ipoint)
        end
    end

    if !isnothing(F)
        for ipoint = 1:length(F)
            set_external_forces!(x, system, prescribed_conditions, F[ipoint], ipoint)
        end
    end

    if !isnothing(M)
        for ipoint = 1:length(M)
            set_external_moments!(x, system, prescribed_conditions, M[ipoint], ipoint)
        end
    end

    if !isnothing(V)
        for ipoint = 1:length(V)
            set_linear_velocity!(x, system, V[ipoint], ipoint)
        end
    end

    if !isnothing(Omega)
        for ipoint = 1:length(Omega)
            set_angular_velocity!(x, system, Omega[ipoint], ipoint)
        end
    end

    if !isnothing(Fi)
        for ielem = 1:length(Fi)
            set_internal_forces!(x, system, Fi[ielem], ielem)
        end
    end

    if !isnothing(Mi)
        for ielem = 1:length(Mi)
            set_internal_moments!(x, system, Mi[ielem], ielem)
        end
    end

    if !isnothing(F1)
        for ielem = 1:length(F1)
            set_start_forces!(x, system, F1[ielem], ielem)
        end
    end

    if !isnothing(M1)
        for ielem = 1:length(M1)
            set_start_moments!(x, system, M1[ielem], ielem)
        end
    end

    if !isnothing(F2)
        for ielem = 1:length(F2)
            set_end_forces!(x, system, F2[ielem], ielem)
        end
    end

    if !isnothing(M2)
        for ielem = 1:length(M2)
            set_end_moments!(x, system, M2[ielem], ielem)
        end
    end

    if !isnothing(V_p)
        for ipoint = 1:length(V_p)
            set_point_linear_velocity!(x, system, V_p[ipoint], ipoint)
        end
    end

    if !isnothing(Omega_p)
        for ipoint = 1:length(Omega_p)
            set_point_angular_velocity!(x, system, Omega_p[ipoint], ipoint)
        end
    end

    if !isnothing(V_e)
        for ielem = 1:length(V_e)
            set_element_linear_velocity!(x, system, V_e[ielem], ielem)
        end
    end

    if !isnothing(Omega_e)
        for ielem = 1:length(Omega_e)
            set_element_angular_velocity!(x, system, Omega_e[ielem], ielem)
        end
    end

    if !isnothing(linear_displacement)
        set_body_linear_displacement!(x, system, linear_displacement)
    end

    if !isnothing(angular_displacement)
        set_body_angular_displacement!(x, system, angular_displacement)
    end

    if !isnothing(linear_velocity)
        set_body_linear_velocity!(x, system, linear_velocity)
    end

    if !isnothing(angular_velocity)
        set_body_angular_velocity!(x, system, angular_velocity)
    end

    if !isnothing(linear_acceleration)
        set_body_linear_acceleration!(x, system, prescribed_conditions, linear_acceleration)
    end

    if !isnothing(angular_acceleration)
        set_body_angular_acceleration!(x, system, prescribed_conditions, angular_acceleration)
    end

    return x
end

"""
    set_linear_displacement!([x,] system, prescribed_conditions, u, ipoint)

Set the state variables in `system` (or in the vector `x`) corresponding to the
linear deflection of point `ipoint` to the provided values.
"""
function set_linear_displacement!(system, prescribed_conditions, u, ipoint)
    set_linear_displacement!(system.x, system, prescribed_conditions, u, ipoint)
    return system
end

function set_linear_displacement!(x, system, prescribed_conditions, u, ipoint)
    
    icol = system.indices.icol_point[ipoint]
    
    prescribed = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(system.t)

    if haskey(prescribed, ipoint)
        !prescribed[ipoint].pd[1] && setindex!(x, u[1], icol)
        !prescribed[ipoint].pd[2] && setindex!(x, u[2], icol+1)
        !prescribed[ipoint].pd[3] && setindex!(x, u[3], icol+2)
    else
        x[icol  ] = u[1]
        x[icol+1] = u[2]
        x[icol+2] = u[3]
    end

    return x
end

"""
    set_angular_displacement!([x,] system, prescribed_conditions, theta, ipoint)

Set the state variables in `system` (or in the vector `x`) corresponding to the
angular deflection of point `ipoint` to the provided values.
"""
function set_angular_displacement!(system, prescribed_conditions, theta, ipoint)
    set_angular_displacement!(system.x, system, prescribed_conditions, theta, ipoint)
    return system
end

function set_angular_displacement!(x, system, prescribed_conditions, theta, ipoint)

    icol = system.indices.icol_point[ipoint]
    
    prescribed = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(system.t)

    if haskey(prescribed, ipoint)
        !prescribed[ipoint].pd[4] && setindex!(x, theta[1], icol+3)
        !prescribed[ipoint].pd[5] && setindex!(x, theta[2], icol+4)
        !prescribed[ipoint].pd[6] && setindex!(x, theta[3], icol+5)
    else
        x[icol+3] = theta[1]
        x[icol+4] = theta[2]
        x[icol+5] = theta[3]
    end

    return x
end

"""
    set_linear_velocity!([x,] system, V, ipoint)

Set the state variables in `system` (or in the vector `x`) corresponding to the
linear velocity of point `ipoint` to the provided values.
"""
function set_linear_velocity!(system, V, ipoint)
    set_linear_velocity!(system.x, system, V, ipoint)
    return system
end

function set_linear_velocity!(x, system, V, ipoint)

    icol = system.indices.icol_point[ipoint]

    x[icol+6] = V[1]
    x[icol+7] = V[2]
    x[icol+8] = V[3]

    return x
end

"""
    set_angular_velocity!([x,] system, Omega, ipoint)

Set the state variables in `system` (or in the vector `x`) corresponding to the
angular velocity of point `ipoint` to the provided values.
"""
function set_angular_velocity!(system, Omega, ipoint)
    set_angular_velocity!(system.x, system, Omega, ipoint)
    return system
end

function set_angular_velocity!(x, system, Omega, ipoint)

    icol = system.indices.icol_point[ipoint]

    x[icol+9] = Omega[1]
    x[icol+10] = Omega[2]
    x[icol+11] = Omega[3]

    return x
end

"""
    set_external_forces!([x,] system, prescribed_conditions, F, ipoint)

Set the state variables in `system` (or in the vector `x`) corresponding to the
external forces applied at point `ipoint` to the provided values.
"""
function set_external_forces!(system, prescribed_conditions, F, ipoint)
    set_external_forces!(system.x, system, prescribed_conditions, F, ipoint)
    return system
end

function set_external_forces!(x, system, prescribed_conditions, F, ipoint)
    
    icol = system.indices.icol_point[ipoint]

    force_scaling = system.force_scaling
    
    prescribed = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(system.t)

    if haskey(prescribed, ipoint)
        !prescribed[ipoint].pl[1] && setindex!(x, F[1] / force_scaling, icol)
        !prescribed[ipoint].pl[2] && setindex!(x, F[2] / force_scaling, icol+1)
        !prescribed[ipoint].pl[3] && setindex!(x, F[3] / force_scaling, icol+2)
    else
        x[icol  ] = F[1] / force_scaling
        x[icol+1] = F[2] / force_scaling
        x[icol+2] = F[3] / force_scaling
    end

    return x
end

"""
    set_external_moments!([x,] system, prescribed_conditions, M, ipoint)

Set the state variables in `system` (or in the vector `x`) corresponding to the
external moments applied at point `ipoint` to the provided values.
"""
function set_external_moments!(system, prescribed_conditions, M, ipoint)
    set_external_moments!(system.x, system, prescribed_conditions, M, ipoint)
    return system
end

function set_external_moments!(x, system, prescribed_conditions, M, ipoint)

    icol = system.indices.icol_point[ipoint]
    
    force_scaling = system.force_scaling

    prescribed = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(system.t)

    if haskey(prescribed, ipoint)
        !prescribed[ipoint].pl[4] && setindex!(x, M[1] / force_scaling, icol+3)
        !prescribed[ipoint].pl[5] && setindex!(x, M[2] / force_scaling, icol+4)
        !prescribed[ipoint].pl[6] && setindex!(x, M[3] / force_scaling, icol+5)
    else
        x[icol+3] = M[1] / force_scaling
        x[icol+4] = M[2] / force_scaling
        x[icol+5] = M[3] / force_scaling
    end

    return x
end

"""
    set_internal_forces!([x,] system::Union{StaticSystem,DynamicSystem}, Fi, ielem)

Set the state variables in `system` (or in the vector `x`) corresponding to the
internal forces of element `ielem` to the provided values.
"""
function set_internal_forces!(system::Union{StaticSystem,DynamicSystem}, Fi, ielem)
    set_internal_forces!(system.x, system, Fi, ielem)
    return system
end

function set_internal_forces!(x, system::Union{StaticSystem,DynamicSystem}, Fi, ielem)

    icol = system.indices.icol_elem[ielem]

    force_scaling = system.force_scaling

    x[icol  ] = Fi[1] / force_scaling
    x[icol+1] = Fi[2] / force_scaling
    x[icol+2] = Fi[3] / force_scaling

    return x
end

"""
    set_internal_moments!([x,] system::Union{StaticSystem, DynamicSystem}, Mi, ielem)

Set the state variables in `system` (or in the vector `x`) corresponding to the
internal moments of element `ielem` to the provided values.
"""
function set_internal_moments!(system::Union{StaticSystem, DynamicSystem}, Mi, ielem)
    set_internal_moments!(system.x, system, Mi, ielem)
    return system
end

function set_internal_moments!(x, system::Union{StaticSystem, DynamicSystem}, Mi, ielem)

    icol = system.indices.icol_elem[ielem]

    force_scaling = system.force_scaling

    x[icol+3] = Mi[1] / force_scaling
    x[icol+4] = Mi[2] / force_scaling
    x[icol+5] = Mi[3] / force_scaling

    return x
end

"""
    set_start_forces!([x,] system, F1, ielem)

Set the state variables in `system` (or in the vector `x`) corresponding to the
resultant forces at the start of element `ielem` to the provided values.
"""
function set_start_forces!(system::ExpandedSystem, F1, ielem)
    set_start_forces!(system.x, system, F1, ielem)
    return system
end

function set_start_forces!(x, system::ExpandedSystem, F1, ielem)

    icol = system.indices.icol_elem[ielem]

    force_scaling = system.force_scaling

    x[icol] = F1[1] / force_scaling
    x[icol+1] = F1[2] / force_scaling
    x[icol+2] = F1[3] / force_scaling

    return x
end

"""
    set_start_moments!([x,] system, M1, ielem)

Set the state variables in `system` (or in the vector `x`) corresponding to the
resultant moments at the start of element `ielem` to the provided values.
"""
function set_start_moments!(system::ExpandedSystem, M1, ielem)
    set_start_moments!(system.x, system, M1, ielem)
    return system
end

function set_start_moments!(x, system::ExpandedSystem, M1, ielem)

    icol = system.indices.icol_elem[ielem]

    force_scaling = system.force_scaling

    x[icol+3] = M1[1] / force_scaling
    x[icol+4] = M1[2] / force_scaling
    x[icol+5] = M1[3] / force_scaling

    return x
end

"""
    set_end_forces!([x,] system, F2, ielem)

Set the state variables in `system` (or in the vector `x`) corresponding to the
resultant forces at the end of element `ielem` to the provided values.
"""
function set_end_forces!(system::ExpandedSystem, F2, ielem)
    set_end_forces!(system.x, system, F2, ielem)
    return system
end

function set_end_forces!(x, system::ExpandedSystem, F2, ielem)

    icol = system.indices.icol_elem[ielem]

    force_scaling = system.force_scaling

    x[icol+6] = F2[1] / force_scaling
    x[icol+7] = F2[2] / force_scaling
    x[icol+8] = F2[3] / force_scaling

    return x
end

"""
    set_end_moments!([x,] system, M2, ielem)

Set the state variables in `system` (or in the vector `x`) corresponding to the
resultant moments at the end of element `ielem` to the provided values.
"""
function set_end_moments!(system::ExpandedSystem, M2, ielem)
    set_end_moments!(system.x, system, M2, ielem)
    return system
end

function set_end_moments!(x, system::ExpandedSystem, M2, ielem)

    icol = system.indices.icol_elem[ielem]

    force_scaling = system.force_scaling

    x[icol+9] = M2[1] / force_scaling
    x[icol+10] = M2[2] / force_scaling
    x[icol+11] = M2[3] / force_scaling

    return x
end

"""
    set_point_linear_velocity!([x,] system::ExpandedSystem, V, ipoint)

Set the state variables in `system` (or in the vector `x`) corresponding to the
linear velocity of point `ipoint` to the provided values.
"""
function set_point_linear_velocity!(system::ExpandedSystem, V, ipoint)
    set_point_linear_velocity!(system.x, system::ExpandedSystem, V, ipoint)
    return system
end

function set_point_linear_velocity!(x, system::ExpandedSystem, V, ipoint)

    icol = system.indices.icol_point[ipoint]

    x[icol+6] = V[1]
    x[icol+7] = V[2]
    x[icol+8] = V[3]

    return x
end

"""
    set_point_angular_velocity!([x,] system::ExpandedSystem, Omega, ipoint)

Set the state variables in `system` (or in the vector `x`) corresponding to the
angular velocity of point `ipoint` to the provided values.
"""
function set_point_angular_velocity!(system::ExpandedSystem, Omega, ipoint)
    set_point_angular_velocity!(system.x, system, Omega, ipoint)
    return system
end

function set_point_angular_velocity!(x, system::ExpandedSystem, Omega, ipoint)

    icol = system.indices.icol_point[ipoint]

    x[icol+9] = Omega[1]
    x[icol+10] = Omega[2]
    x[icol+11] = Omega[3]

    return x
end

"""
    set_element_linear_velocity!([x,] system::ExpandedSystem, V, ielem)

Set the state variables in `system` (or in the vector `x`) corresponding to the
linear velocity of beam element `ielem` to the provided values.
"""
function set_element_linear_velocity!(system::ExpandedSystem, V, ielem)
    set_element_linear_velocity!(system.x, system::ExpandedSystem, V, ielem)
    return system
end

function set_element_linear_velocity!(x, system::ExpandedSystem, V, ielem)

    icol = system.indices.icol_elem[ielem]

    x[icol+12] = V[1]
    x[icol+13] = V[2]
    x[icol+14] = V[3]

    return x
end

"""
    set_element_angular_velocity!([x,] system::ExpandedSystem, Omega, ielem)

Set the state variables in `system` (or in the vector `x`) corresponding to the
angular velocity of element `ielem` to the provided values.
"""
function set_element_angular_velocity!(system::ExpandedSystem, Omega, ielem)
    set_element_angular_velocity!(system.x, system, Omega, ielem)
    return system
end

function set_element_angular_velocity!(x, system::ExpandedSystem, Omega, ielem)

    icol = system.indices.icol_elem[ielem]

    x[icol+15] = Omega[1]
    x[icol+16] = Omega[2]
    x[icol+17] = Omega[3]

    return x
end

"""
    set_body_linear_displacement!([x,] system, linear_displacement)

Set the state variables in `system` (or in the vector `x`) corresponding to the
linear displacement of the body frame to the provided values.
"""
function set_body_linear_displacement!(system, linear_displacement)
    set_body_linear_displacement!(system.x, system, linear_displacement)
    return system
end

function set_body_linear_displacement!(x, system, linear_displacement)
        
    x[1:3] .= linear_displacement

    return x
end

"""
    set_body_angular_displacement!([x,] system, angular_displacement)

Set the state variables in `system` (or in the vector `x`) corresponding to the
angular deflection of point `ipoint` to the provided values.
"""
function set_body_angular_displacement!(system, angular_displacement)
    set_body_angular_displacement!(system.x, system, angular_displacement)
    return system
end

function set_body_angular_displacement!(x, system, angular_displacement)

    x[4:6] .= angular_displacement

    return x
end

"""
    set_body_linear_velocity!([x,] system, linear_velocity)

Set the state variables in `system` (or in the vector `x`) corresponding to the
linear velocity of the body frame to the provided values.
"""
function set_body_linear_velocity!(system, linear_velocity)
    set_body_linear_velocity!(system.x, system, linear_velocity)
    return system
end

function set_body_linear_velocity!(x, system, linear_velocity)
        
    x[7:9] .= linear_velocity

    return x
end

"""
    set_body_angular_velocity!([x,] system, angular_velocity)

Set the state variables in `system` (or in the vector `x`) corresponding to the
angular deflection of point `ipoint` to the provided values.
"""
function set_body_angular_velocity!(system, angular_velocity)
    set_body_angular_velocity!(system.x, system, angular_velocity)
    return system
end

function set_body_angular_velocity!(x, system, angular_velocity)

    x[10:12] .= angular_velocity

    return x
end

"""
    set_body_linear_acceleration!([x,] system, prescribed_conditions, linear_acceleration)

Set the state variables in `system` (or in the vector `x`) corresponding to the
linear acceleration of the body frame to the provided values.
"""
function set_body_linear_acceleration!(system, prescribed_conditions, linear_acceleration)
    set_body_linear_acceleration!(system.x, system, prescribed_conditions, linear_acceleration)
    return system
end

function set_body_linear_acceleration!(x, system, prescribed_conditions, linear_acceleration)
        
    prescribed = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(system.t)

    icol = body_frame_acceleration_indices(system, prescribed)

    !iszero(icol[1]) && setindex!(x, linear_acceleration[1], icol[1])
    !iszero(icol[2]) && setindex!(x, linear_acceleration[2], icol[2])
    !iszero(icol[3]) && setindex!(x, linear_acceleration[3], icol[3])

    return x
end

"""
    set_body_angular_acceleration!([x,] system, prescribed_conditions, angular_acceleration)

Set the state variables in `system` (or in the vector `x`) corresponding to the
angular deflection of point `ipoint` to the provided values.
"""
function set_body_angular_acceleration!(system, prescribed_conditions, angular_acceleration)
    set_body_angular_acceleration!(system.x, system, prescribed_conditions, angular_acceleration)
    return system
end

function set_body_angular_acceleration!(x, system, prescribed_conditions, angular_acceleration)

    prescribed = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(system.t)

    icol = body_frame_acceleration_indices(system, prescribed)

    !iszero(icol[4]) && setindex!(x, angular_acceleration[1], icol[4])
    !iszero(icol[5]) && setindex!(x, angular_acceleration[2], icol[5])
    !iszero(icol[6]) && setindex!(x, angular_acceleration[3], icol[6])

    return x
end
