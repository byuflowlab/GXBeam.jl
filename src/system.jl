"""
    System{TF, TV<:AbstractVector{TF}, TM<:AbstractMatrix{TF}}

Contains the system state, residual vector, and jacobian matrices as well as
pointers to be able to access their contents.  Also contains additional storage
needed for time domain simulations.

# Fields:
 - `static`: Flag indicating whether system matrices are only valid for static analyses
 - `x`: State vector
 - `r`: Residual vector
 - `K`: System jacobian matrix with respect to the state variables
 - `M`: System jacobian matrix with respect to the time derivative of the state variables
 - `force_scaling`: Scaling for state variables corresponding to forces/moments
 - `irow_point`: Row index of first equilibrium equation for each point
 - `irow_elem`: Row index of first equation for just this beam element
 - `irow_elem1`: Row index of first equation for the left side of each beam
 - `irow_elem2`: Row index of first equation for the right side of each beam
 - `icol_point`: Row/Column index of first state variable for each point
 - `icol_elem`: Row/Column index of first state variable for each beam element
 - `udot`: Time derivative of state variable `u` for each beam element
 - `θdot`: Time derivative of state variable `θ` for each beam element
 - `Fdot`: Time derivative of state variable `F` for each beam element
 - `Mdot`: Time derivative of state variable `M` for each beam element
 - `Vdot`: Time derivative of state variable `V` for each beam element
 - `Ωdot`: Time derivative of state variable `Ω` for each beam element
 - `t`: Current system time
"""
mutable struct System{TF, TV<:AbstractVector{TF}, TM<:AbstractMatrix{TF}}
    static::Bool
    x::TV
    r::TV
    K::TM
    M::TM
    force_scaling::TF
    irow_point::Vector{Int}
    irow_elem::Vector{Int}
    irow_elem1::Vector{Int}
    irow_elem2::Vector{Int}
    icol_point::Vector{Int}
    icol_elem::Vector{Int}
    udot::Vector{SVector{3,TF}}
    θdot::Vector{SVector{3,TF}}
    Fdot::Vector{SVector{3,TF}}
    Mdot::Vector{SVector{3,TF}}
    Vdot::Vector{SVector{3,TF}}
    Ωdot::Vector{SVector{3,TF}}
    t::TF
end
Base.eltype(::System{TF, TV, TM}) where {TF, TV, TM} = TF

"""
    System([TF=eltype(assembly),] assembly, static; kwargs...)

Initialize an object of type `System` which stores the system state.

# Arguments:
 - `TF:`(optional) Floating point type, defaults to the floating point type of `assembly`
 - `assembly`: Assembly of rigidly connected nonlinear beam elements
 - `static`: Flag indicating whether the system corresponds to a static system.

# Keyword Arguments
 - `prescribed_points`: Point indices corresponding to the points whose equations
    and state variables should be included in the system of equations.  By default,
    all point indices are included in the system of equations.
 - `force_scaling`: Factor used to scale system forces/moments internally.  If
    not specified, a suitable default will be chosen based on the entries of the
    compliance matrix.

Note that points with prescribed conditions must be included in the system of
equations.
"""
function System(assembly, static; kwargs...)

    return System(eltype(assembly), assembly, static; kwargs...)
end

function System(TF, assembly, static;
    prescribed_points = 1:length(assembly.points),
    force_scaling = default_force_scaling(assembly)
    )

    # system dimensions
    npoint = length(assembly.points)
    nelem = length(assembly.elements)

    # initialize system pointers
    N, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem =
        system_indices(assembly.start, assembly.stop, static; prescribed_points)

    # initialize system matrices
    x = zeros(TF, N)
    r = zeros(TF, N)
    K = spzeros(TF, N, N)
    M = spzeros(TF, N, N)

    # initialize storage for time domain simulations
    udot = [@SVector zeros(TF, 3) for i = 1:nelem]
    θdot = [@SVector zeros(TF, 3) for i = 1:nelem]
    Fdot = [@SVector zeros(TF, 3) for i = 1:nelem]
    Mdot = [@SVector zeros(TF, 3) for i = 1:nelem]
    Vdot = [@SVector zeros(TF, 3) for i = 1:nelem]
    Ωdot = [@SVector zeros(TF, 3) for i = 1:nelem]

    # initialize current time
    t = 0.0

    # set system types
    TV = promote_type(typeof(x), typeof(r))
    TM = promote_type(typeof(K), typeof(M))

    return System{TF, TV, TM}(static, x, r, K, M, force_scaling,
        irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem, udot, θdot,
        Fdot, Mdot, Vdot, Ωdot, t)
end

function default_force_scaling(assembly)

    TF = eltype(assembly)

    # Count and sum all nonzero entries
    force_sum = 0.0
    N_entries = 0
    for elem in assembly.elements
        for val in elem.compliance
            if abs(val) > eps(TF)
                force_sum += abs(val)
                N_entries += 1
            end
        end
    end

    # set force scaling based on nonzero compliance matrix entries
    if N_entries == 0
        force_scaling = 1.0
    else
        force_scaling = nextpow(2.0, N_entries/force_sum/100)
    end

    return force_scaling
end

"""
    system_indices(start, stop, static; kwargs...)

Return indices for accessing the equations and state variables associated with
each point and beam element in a system given its connectivity.

# Arguments:
 - `start`: Vector containing the point indices where each beam element starts
 - `stop`: Vector containing the point indices where each beam element stops
 - `static`: Flag indicating whether the analysis is static (rather than dynamic).
    Defaults to `false`.

# Keyword Arguments:
 - `prescribed_points`: Point indices corresponding to the points whose equations
    and state variables should be included in the system of equations.  By default,
    all point indices are included in the system of equations.

# Return Arguments:
 - `N`: total number of equations and unknowns in the system
 - `irow_point`: Row index of the first equation corresponding to each point
 - `irow_elem`: Row index of the first equation corresponding to each beam element
 - `irow_elem1`: Row index of the first equation corresponding to the start of
        each beam element
 - `irow_elem2`: Row index of the first equation corresponding to the end of
        each beam element
 - `icol_point`: Column index of the first state variable corresponding to each
        point
 - `icol_elem`: Column index of the first state variable corresponding to each
        beam element

Negative indices indicate that the equations and/or state variables associated
with the point/beam element have been omitted from the system of equations.
"""
function system_indices(start, stop, static;
    prescribed_points = 1:max(maximum(start), maximum(stop)))

    npoint = max(maximum(start), maximum(stop))
    nelem = length(start)

    keep = [i in prescribed_points for i = 1:npoint]

    add_necessary_points!(keep, start, stop)

    # indicates whether a point has associated equations and state variables
    assigned = fill(false, npoint)

    # initialize pointers for equations
    irow_point = Vector{Int}(undef, npoint)
    irow_elem = Vector{Int}(undef, nelem)
    irow_elem1 = Vector{Int}(undef, nelem)
    irow_elem2 = Vector{Int}(undef, nelem)

    # initialize pointers for state variables
    icol_point = Vector{Int}(undef, npoint)
    icol_elem = Vector{Int}(undef, nelem)

    irow = 1
    icol = 1
    for ielem = 1:nelem
        # add state variables/equations for the start of the beam element
        ipt = start[ielem]

        # check if the point has associated equations and state variables
        if !assigned[ipt]
            # add equations and state variables for this point
            assigned[ipt] = true

            # add 6 equilibrium equations + 6 compatability equations
            irow_point[ipt] = irow
            irow_elem1[ielem] = irow
            irow += 12

            # add 6 state variables for this point (if necessary)
            if keep[ipt]
                icol_point[ipt] = icol
                icol += 6
            else
                icol_point[ipt] = -1
            end
        else
            # add additional equations for this point

            # add compatibility equations for this point (if necessary)
            if keep[ipt]
                irow_elem1[ielem] = irow
                irow += 6
            else
                irow_elem1[ielem] = -1
            end
        end

        # add state variables/equations for the beam element

        # add 12 state variables for this element
        icol_elem[ielem] = icol
        icol += 12

        # add an additional 6 state variables and equations for this element (if necessary)
        if static
            irow_elem[ielem] = -1
        else
            irow_elem[ielem] = irow
            irow += 6
            icol += 6
        end

        # add state variables/equations for the end of the beam element
        ipt = stop[ielem]

        # check if the point has associated equations and state variables
        if !assigned[ipt]
            # add equations and state variables for this point
            assigned[ipt] = true

            # add 6 equilibrium equations + 6 compatability equations
            irow_point[ipt] = irow
            irow_elem2[ielem] = irow
            irow += 12

            # add 6 state variables for this point (if necessary)
            if keep[ipt]
                icol_point[ipt] = icol
                icol += 6
            else
                icol_point[ipt] = -1
            end
        else
            # add additional equations for this point

            # add compatibility equations for this point (if necessary)
            if keep[ipt]
                irow_elem2[ielem] = irow
                irow += 6
            else
                irow_elem2[ielem] = -1
            end
        end
    end

    # number of state variables/equations
    N = irow - 1

    return N, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem
end

function add_necessary_points!(keep, start, stop)

    npoint = length(keep)

    for ipt = 1:npoint
        # calculate number of connections to this point
        ncon = count(x -> x == ipt, start) + count(x -> x == ipt, stop)

        if ncon != 2
            # this point must be preserved in the system of equations
            if keep[ipt] != true
                # only mutate if necessary
                keep[ipt] = true
            end
        end
    end

    return keep
end

"""
    system_state(system)

Return a vector containing the state variables of `system`.
"""
system_state(system) = system.x

"""
    set_state!([x,] system, prescribed_conditions; kwargs...)

Set the state variables in `system` (or in the vector `x`) to the provided values.
If values are not provided for a given keyword argument, then the state variables
corresponding to the keyword argument are not updated.

# Keyword Arguments
 - `u_e`: Vector containing initial displacements (in the body frame) for each beam element.
 - `theta_e`: Vector containing rotation variables (in the body frame) for each beam element.
 - `F_e`: Vector containing resultant forces (in the deformed local beam element frame) for each beam element.
 - `M_e`: Vector containing resultant moments (in the deformed local beam element frame) for each beam element.
 - `V_e`: Vector containing linear velocity (in the deformed local beam element frame) for each beam element.
 - `Ω_e`: Vector containing angular velocity (in the deformed local beam element frame) for each beam element.
 - `u_p`: Vector containing initial displacements (in the body frame) for each point.
 - `theta_p`: Vector containing rotation variables (in the body frame) for each point.
 - `F_p`: Vector containing externally applied forces (in the body frame) for each point.
 - `M_p`: Vector containing externally applied moments (in the body frame) for each point.
"""
set_state!

function set_state!(system, prescribed_conditions; kwargs...)
    x = set_state!(system.x, system, prescribed_conditions; kwargs...)
    return system
end

function set_state!(x, system, prescribed_conditions; kwargs...)
    return set_state!(x, system.icol_elem, system.icol_point; kwargs...)
end

function set_state!(x, icol_elem, icol_point, force_scaling,
    prescribed_conditions; u_e = nothing, theta_e = nothing,
    F_e = nothing, M_e = nothing, V_e = nothing, Ω_e = nothing,
    u_p = nothing, theta_p = nothing, F_p = nothing, M_p = nothing)

    nelem = length(icol_elem)
    npoint = length(icol_point)

    for ielem = 1:nelem
        icol = icol_elem[ielem]
        if !isnothing(u_e)
            set_element_deflection!(x, icol, u_e)
        end
        if !isnothing(theta_e)
            set_element_rotation!(x, icol, theta_e)
        end
        if !isnothing(F_e)
            set_element_forces!(x, icol, F_e, force_scaling)
        end
        if !isnothing(M_e)
            set_element_moments!(x, icol, M_e, force_scaling)
        end
        if !isnothing(V_e)
            set_element_linear_velocity!(x, icol, V_e)
        end
        if !isnothing(Ω_e)
            set_element_angular_velocity!(x, icol, Ω_e)
        end
    end

    for ipoint = 1:npoint
        icol = icol_point[ipoint]
        if !isnothing(u_p)
            set_point_deflections!(x, icol, u_p, prescribed_conditions)
        end
        if !isnothing(theta_p)
            set_point_rotations!(x, icol, theta_p, prescribed_conditions)
        end
        if !isnothing(F_p)
            set_point_forces!(x, icol, F_p, force_scaling, prescribed_conditions)
        end
        if !isnothing(M_p)
            set_point_moments!(x, icol, M_p, force_scaling, prescribed_conditions)
        end
    end

    return x
end

"""
    set_element_deflection!([x,] system, u_e, ielem)

Set the state variables in `system` (or in the vector `x`) corresponding to the
deflections of element `ielem` to the provided values.
"""
function set_element_deflection!(system::System, u_e, ielem)
    set_element_deflection!(system.x, system, u_e, ielem)
    return system
end

function set_element_deflection!(x, system::System, u_e, ielem)
    icol = system.icol_elem[ielem]
    set_element_deflection!(x, icol, u_e)
    return x
end

function set_element_deflection!(x, icol, u_e)
    x[icol:icol+2] .= u_e
    return x
end

"""
    set_element_rotation!([x,] system, θ_b, ielem)

Set the state variables in `system` (or in the vector `x`) corresponding to the
rotations of element `ielem` to the provided values.
"""
function set_element_rotation!(system::System, θ_e, ielem)
    set_element_rotation!(system.x, system, θ_e, ielem)
    return system
end

function set_element_rotation!(x, system::System, θ_e, ielem)
    icol = system.icol_elem[ielem]
    set_element_rotation!(x, icol, θ_e)
    return x
end

function set_element_rotation!(x, icol, θ_e)
    x[icol+3:icol+5] .= θ_e
    return x
end

"""
    set_element_forces!([x,] system, F_e, ielem)

Set the state variables in `system` (or in the vector `x`) corresponding to the
resultant forces of element `ielem` to the provided values.
"""
function set_element_forces!(system::System, F_e, ielem)
    set_element_forces!(system.x, system, F_e, ielem)
    return system
end

function set_element_forces!(x, system::System, F_e, ielem)
    icol = system.icol_elem[ielem]
    force_scaling = system.force_scaling
    set_element_forces!(x, icol, F_e, force_scaling)
    return x
end

function set_element_forces!(x, icol, F_e, force_scaling)
    x[icol+6:icol+8] .= F_e ./ force_scaling
    return x
end

"""
    set_element_moments!([x,] system, u_e, ielem)

Set the state variables in `system` (or in the vector `x`) corresponding to the
resultant moments of element `ielem` to the provided values.
"""
function set_element_moments!(system::System, M_e, ielem)
    set_element_moments!(system.x, system, M_e, ielem)
    return system
end

function set_element_moments!(x, system::System, M_e, ielem)
    icol = system.icol_elem[ielem]
    force_scaling = system.force_scaling
    set_element_moments!(x, icol, M_e, force_scaling)
    return x
end

function set_element_moments!(x, icol, M_e, force_scaling)
    x[icol+9:icol+11] .= M_e ./ force_scaling
    return x
end

"""
    set_element_linear_velocity!([x,] system, V_e, ielem)

Set the state variables in `system` (or in the vector `x`) corresponding to the
linear velocity of element `ielem` to the provided values.
"""
function set_element_linear_velocity!(system::System, V_e, ielem)
    set_element_linear_velocity!(x, system, V_e, ielem)
    return system
end

function set_element_linear_velocity!(x, system::System, V_e, ielem)
    @assert !system.static
    icol = system.icol_elem[ielem]
    set_element_linear_velocity!(x, icol, V_e)
    return x
end

function set_element_linear_velocity!(x, icol, V_e)
    x[icol+12:icol+14] .= V_e
    return x
end

"""
   set_element_angular_velocity!([x,] system, Ω_e, ielem)

Set the state variables in `system` (or in the vector `x`) corresponding to the
angular velocity of element `ielem` to the provided values.
"""
function set_element_angular_velocity!(system::System, Ω_e, ielem)
    set_element_angular_velocity!(x, system, Ω_e, ielem)
    return system
end

function set_element_angular_velocity!(x, system::System, Ω_e, ielem)
    @assert !system.static
    icol = system.icol_elem[ielem]
    set_element_angular_velocity!(x, icol, Ω_e)
    return x
end

function set_element_angular_velocity!(x, icol, Ω_e)
    x[icol+15:icol+17] .= Ω_e
    return x
end

"""
    set_point_deflections!([x,] system, u_e, ipoint, prescribed_conditions)

Set the state variables in `system` (or in the vector `x`) corresponding to the
deflections (or externally applied forces) at point `ipoint` to the provided values.
"""
set_point_deflections!

function set_point_deflections!(system::System, u_p, ipoint, prescribed_conditions)
    x = set_point_deflections!(system.x, system, u_p, ipoint, prescribed_conditions)
    return system
end

function set_point_deflections!(x, system::System, u_p, ipoint, prescribed_conditions)

    icol = system.icol_point[ipoint]

    if ipoint in keys(prescribed_conditions)
        prescribed_forces = prescribed_conditions[ipoint].force
    else
        prescribed_forces = @SVector ones(Bool, 6)
    end

    set_point_deflections!(x, icol, u_p, prescribed_forces)

    return x
end

function set_point_deflections!(x, icol, u_p, prescribed_forces)

    if icol <= 0
        # point variables are not state variables
        return x
    end

    # point variables are state variables
    for k = 1:3
        if prescribed_forces[k]
            # applied force is prescribed, deflection is a state variable
            x[icol+k-1] = u_p[k]
        end
    end

    return x
end

"""
    set_point_rotations!([x,] system, θ_e, ipoint, prescribed_conditions)

Set the state variables in `system` (or in the vector `x`) corresponding to the
rotations (or externally applied moments) at point `ipoint` to the provided values.
"""
set_point_rotations!

function set_point_rotations!(system::System, θ_p, ipoint, prescribed_conditions)
    set_point_rotations!(system.x, system, θ_p, ipoint, prescribed_conditions)
    return system
end

function set_point_rotations!(x, system::System, θ_p, ipoint, prescribed_conditions)

    icol = system.icol_point[ipoint]

    if ipoint in keys(prescribed_conditions)
        prescribed_forces = prescribed_conditions[ipoint].force
    else
        prescribed_forces = @SVector ones(Bool, 6)
    end

    set_point_rotations!(x, icol, θ_p, prescribed_forces)

    return x
end

function set_point_rotations!(x, icol, θ_p, prescribed_forces)

    if icol <= 0
        # point variables are not state variables
        return x
    end

    # point variables are state variables
    for k = 1:3
        if prescribed_forces[3+k]
            # applied moment is prescribed, rotation is a state variable
            x[icol+k+2] = θ_p[k]
        end
    end

    return x
end

"""
    set_point_forces!([x,] system, F_p, ipoint, prescribed_conditions)

Set the state variables in `system` (or in the vector `x`) corresponding to the
external forces applied on point `ipoint` to the provided values.
"""
set_point_forces!

function set_point_forces!(system::System, F_p, ipoint, prescribed_conditions)
    set_point_forces!(system.x, system, F_p, ipoint, prescribed_conditions)
    return system
end

function set_point_forces!(x, system::System, F_p, ipoint, prescribed_conditions)

    icol = icol_point[ipoint]

    if ipoint in keys(prescribed_conditions)
        prescribed_forces = prescribed_conditions[ipoint].force
    else
        prescribed_forces = @SVector ones(Bool, 6)
    end

    set_point_forces!(x, icol, F_p, prescribed_forces, force_scaling)

    return x
end

function set_point_forces!(x, icol, F_p, prescribed_forces, force_scaling)

    if icol <= 0
        # point variables are not state variables
        return x
    end

    # point variables are state variables
    for k = 1:3
        # prescribed conditions exist for this point
        if !prescribed_forces[k]
            # deflections are prescribed, applied force is a state variable
            x[icol+k-1] = F_p[k] / force_scaling
        end
    end

    return x
end


"""
    set_point_moments!([x,] system, M_p, ipoint, prescribed_conditions)

Set the state variables in `system` (or in the vector `x`) corresponding to the
external moments applied on point `ipoint` to the provided values.
"""
function set_point_moments!(system::System, M_p, ipoint, prescribed_conditions)
    set_point_moments!(system.x, system, M_p, ipoint, prescribed_conditions)
    return system
end

function set_point_moments!(x, system::System, M_p, ipoint, prescribed_conditions)

    icol = icol_point[ipoint]

    if ipoint in keys(prescribed_conditions)
        prescribed_forces = prescribed_conditions[ipoint].force
    else
        prescribed_forces = @SVector ones(Bool, 6)
    end

    set_point_moments!(x, icol, M_p, prescribed_forces, force_scaling)

    return x
end

function set_point_moments!(x, icol, M_p, prescribed_forces, force_scaling)

    if icol <= 0
        # point variables are not state variables
        return x
    end

    # point variables are state variables
    for k = 1:3
        # prescribed conditions exist for this point
        if !prescribed_forces[3+k]
            # rotations are prescribed, applied moment is a state variable
            x[icol+k+2] = M_p[k] / force_scaling
        end
    end

    return x
end

"""
    reset_state!(system)

Reset the state variables in `system` (stored in `system.x`).
"""
function reset_state!(system)
    system.x .= 0
    return system
end

"""
    get_sparsity(system, assembly)

Return a matrix indicating the sparsity structure of the jacobian matrix.
"""
function get_sparsity(system, assembly)

    N = length(system.x)
    irow_point = system.irow_point
    irow_elem = system.irow_elem
    irow_elem1 = system.irow_elem1
    irow_elem2 = system.irow_elem2
    icol_point = system.icol_point
    icol_elem = system.icol_elem

    sparsity = spzeros(Bool, N, N)

    # --- add constributions from beam element state variables --- #

    nelem = length(icol_elem)

    for ielem = 1:nelem
        # indices for this beam
        icol = icol_elem[ielem]
        irow_e = irow_elem[ielem]
        irow_e1 = irow_elem1[ielem]
        irow_p1 = irow_point[assembly.start[ielem]]
        irow_e2 = irow_elem2[ielem]
        irow_p2 = irow_point[assembly.stop[ielem]]

        # --- static, beam element jacobian contributions --- #

        # left endpoint equilibrium equations
        sparsity[irow_p1:irow_p1+2, icol:icol+8] .= true
        sparsity[irow_p1+3:irow_p1+5, icol+3:icol+11] .= true
        # left point compatability equations
        irow = ifelse(irow_e1 == irow_p1 || irow_e1 <= 0, irow_p1+6, irow_e1)
        sparsity[irow:irow+2, icol:icol+11] .= true
        sparsity[irow+3:irow+5, icol+3:icol+11] .= true
        # right endpoint equilibrium equations
        sparsity[irow_p2:irow_p2+2, icol:icol+8] .= true
        sparsity[irow_p2+3:irow_p2+5, icol+3:icol+11] .= true
        # right point compatability equations
        irow = ifelse(irow_e2 == irow_p2 || irow_e2 <= 0, irow_p2 + 6, irow_e2)
        sparsity[irow:irow+2, icol:icol+11] .= true
        sparsity[irow+3:irow+5, icol+3:icol+11] .= true

        # --- dynamic, beam element jacobian contributions --- #
        if !system.static
            # left endpoint equilibrium equations
            sparsity[irow_p1:irow_p1+2, icol+12:icol+14] .= true
            sparsity[irow_p1+3:irow_p1+5, icol+12:icol+17] .= true
            # right endpoint equilibrium equations
            sparsity[irow_p2:irow_p2+2, icol+12:icol+14] .= true
            sparsity[irow_p2+3:irow_p2+5, icol+12:icol+17] .= true
            # element residual equations
            sparsity[irow_e:irow_e+2, icol:icol+5] .= true
            sparsity[irow_e:irow_e+2, icol+12:icol+17] .= true
            sparsity[irow_e+3:irow_e+5, icol+3:icol+5] .= true
            sparsity[irow_e+3:irow_e+5, icol+12:icol+17] .= true
        end
    end

    # --- add constributions from point state variables --- #

    npoint = length(icol_point)

    for ipoint = 1:npoint

        # skip if the unknowns have been eliminated from the system of equations
        if icol_point[ipoint] <= 0
            continue
        end

        icol = icol_point[ipoint]
        irow_p = irow_point[ipoint]

        for ielem = 1:nelem
            # check left side of beam
            if ipoint == assembly.start[ielem]
                # add jacobian entries for the beam endpoint
                irow_e = irow_elem1[ielem]
                if irow_e == irow_p
                    for i = 1:6
                        for j = 4:6
                            sparsity[irow_e+i-1, icol+j-1] = true
                        end
                        sparsity[irow_e+i+5, icol+i-1] = true
                        sparsity[irow_e+i-1, icol+i-1] = true
                    end
                else
                    for i = 1:6
                        sparsity[irow_e+i-1, icol+i-1] = true
                    end
                end
            end
            # check right side of beam
            if ipoint == assembly.stop[ielem]
                # add jacobian entries for the beam endpoint
                irow_e = irow_elem2[ielem]
                if irow_e == irow_p
                    for i = 1:6
                        for j = 4:6
                            sparsity[irow_e+i-1, icol+j-1] = true
                        end
                        sparsity[irow_e+i+5, icol+i-1] = true
                        sparsity[irow_e+i-1, icol+i-1] = true
                    end
                else
                    for i = 1:6
                        sparsity[irow_e+i-1, icol+i-1] = true
                    end
                end
            end
        end
    end

    return sparsity
end

"""
    static_system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads, point_masses, gvec,
        force_scaling, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem)

Populate the system residual vector `resid` for a static analysis

# Arguments
 - `resid`: system residual vector
 - `x`: current states of the system
 - `assembly`: assembly of nonlinear beam elements
 - `prescribed_conditions`: dictionary of prescribed conditions
 - `distributed_loads`: dictionary of distributed loads
 - `point_masses`: dictionary of point masses 
 - `gvec`: gravity vector
 - `force_scaling`: scaling parameter for forces/moments
 - `irow_point`: row index of first equilibrium equation for each point
 - `irow_elem`: row index of first equation for just this beam element
 - `irow_elem1`: row index of first equation for the left side of each beam
 - `irow_elem2`: row index of first equation for the right side of each beam
 - `icol_point`: column index of first state variable for each point
 - `icol_elem`: column index of first state variable for each beam element
"""
function static_system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads, point_masses, gvec,
    force_scaling, irow_point, irow_elem1, irow_elem2, icol_point, icol_elem)

    npoint = length(assembly.points)
    nelem = length(assembly.elements)

    # add contributions to residual equations from the elements
    for ielem = 1:nelem

        # get pointers for element
        icol = icol_elem[ielem]
        irow_e1 = irow_elem1[ielem]
        irow_p1 = irow_point[assembly.start[ielem]]
        irow_e2 = irow_elem2[ielem]
        irow_p2 = irow_point[assembly.stop[ielem]]

        static_element_residual!(resid, x, ielem, assembly.elements[ielem],
            distributed_loads, point_masses, gvec, force_scaling, icol, irow_e1,
            irow_p1, irow_e2, irow_p2)
    end

    # add contributions to the residual equations from the prescribed point conditions
    for ipoint = 1:npoint

        # skip if the unknowns have been eliminated from the system of equations
        if icol_point[ipoint] <= 0
            continue
        end

        icol = icol_point[ipoint]
        irow_p = irow_point[ipoint]

        point_residual!(resid, x, ipoint, assembly, prescribed_conditions,
            force_scaling, icol, irow_p, irow_elem1, irow_elem2)
    end

    return resid
end

"""
    steady_state_system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads, point_masses, gvec,
        force_scaling, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem,
        x0, v0, ω0, a0, α0)

Populate the system residual vector `resid` for a steady state analysis.

# Arguments
 - `resid`: system residual vector
 - `x`: current states of the system
 - `assembly`: assembly of nonlinear beam elements
 - `prescribed_conditions`: dictionary of prescribed conditions
 - `distributed_loads`: dictionary of distributed loads
 - `point_masses`: dictionary of point masses 
 - `gvec`: gravity vector
 - `force_scaling`: scaling parameter for forces/moments
 - `irow_point`: row index of first equilibrium equation for each point
 - `irow_elem`: row index of first equation for just this beam element
 - `irow_elem1`: row index of first equation for the left side of each beam
 - `irow_elem2`: row index of first equation for the right side of each beam
 - `icol_point`: column index of first state variable for each point
 - `icol_elem`: column index of first state variable for each beam element
 - `x0`: body frame origin
 - `v0`: body frame linear velocity
 - `ω0`: body frame angular velocity
 - `a0`: body frame linear acceleration
 - `α0`: body frame angular acceleration
"""
function steady_state_system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads, point_masses, gvec,
    force_scaling, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem,
    x0, v0, ω0, a0, α0)

    npoint = length(assembly.points)
    nelem = length(assembly.elements)

    # add contributions to residual equations from the elements
    for ielem = 1:nelem

        # get pointers for element
        icol = icol_elem[ielem]
        irow_e = irow_elem[ielem]
        irow_e1 = irow_elem1[ielem]
        irow_p1 = irow_point[assembly.start[ielem]]
        irow_e2 = irow_elem2[ielem]
        irow_p2 = irow_point[assembly.stop[ielem]]

        steady_state_element_residual!(resid, x, ielem, assembly.elements[ielem],
             distributed_loads, point_masses, gvec, force_scaling, icol, irow_e, irow_e1,
            irow_p1, irow_e2, irow_p2, x0, v0, ω0, a0, α0)
    end

    # add contributions to the residual equations from the prescribed point conditions
    for ipoint = 1:npoint

        # skip if the unknowns have been eliminated from the system of equations
        if icol_point[ipoint] <= 0
            continue
        end

        icol = icol_point[ipoint]
        irow_p = irow_point[ipoint]

        point_residual!(resid, x, ipoint, assembly, prescribed_conditions,
            force_scaling, icol, irow_p, irow_elem1, irow_elem2)
    end

    return resid
end

"""
    initial_condition_system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads, point_masses, gvec,
        force_scaling, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem,
        x0, v0, ω0, a0, α0, u, θ, udot, θdot, Fdot, Mdot)

Populate the system residual vector `resid` for an initial conditions analysis.

 # Arguments
 - `resid`: system residual vector
 - `x`: current states of the system
 - `assembly`: assembly of nonlinear beam elements
 - `prescribed_conditions`: dictionary of prescribed conditions
 - `distributed_loads`: dictionary of distributed loads
 - `point_masses`: dictionary of point masses 
 - `gvec`: gravity vector
 - `force_scaling`: scaling parameter for forces/moments
 - `irow_point`: row index of first equilibrium equation for each point
 - `irow_elem`: row index of first equation for just this beam element
 - `irow_elem1`: row index of first equation for the left side of each beam
 - `irow_elem2`: row index of first equation for the right side of each beam
 - `icol_point`: column index of first state variable for each point
 - `icol_elem`: column index of first state variable for each beam element
 - `x0`: body frame origin
 - `v0`: body frame linear velocity
 - `ω0`: body frame angular velocity
 - `a0`: body frame linear acceleration
 - `α0`: body frame angular acceleration
 - `u`: initial linear deflections for each beam element
 - `θ`: initial angular deflections for each beam element
 - `udot`: initial linear deflection rates for each beam element
 - `θdot`: initial angular deflection rates for each beam element
 - `Fdot`: initial elastic force rate for each beam element
 - `Mdot`: initial elastic moment rate for each beam element
"""
function initial_condition_system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads, point_masses, gvec,
    force_scaling, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem,
    x0, v0, ω0, a0, α0, u, θ, udot, θdot, Fdot, Mdot)

    npoint = length(assembly.points)
    nelem = length(assembly.elements)

    # add contributions to residual equations from the elements
    for ielem = 1:nelem

        # get pointers for element
        icol = icol_elem[ielem]
        irow_e = irow_elem[ielem]
        irow_e1 = irow_elem1[ielem]
        irow_p1 = irow_point[assembly.start[ielem]]
        irow_e2 = irow_elem2[ielem]
        irow_p2 = irow_point[assembly.stop[ielem]]

        initial_condition_element_residual!(resid, x, ielem, assembly.elements[ielem],
            distributed_loads, point_masses, gvec, force_scaling, icol, irow_e, irow_e1,
            irow_p1, irow_e2, irow_p2, x0, v0, ω0, a0, α0,
            u[ielem], θ[ielem], udot[ielem], θdot[ielem], Fdot[ielem], Mdot[ielem])
    end

    # add contributions to the residual equations from the prescribed point conditions
    for ipoint = 1:npoint

        # skip if the unknowns have been eliminated from the system of equations
        if icol_point[ipoint] <= 0
            continue
        end

        icol = icol_point[ipoint]
        irow_p = irow_point[ipoint]

        point_residual!(resid, x, ipoint, assembly, prescribed_conditions,
            force_scaling, icol, irow_p, irow_elem1, irow_elem2)
    end

    return resid
end

"""
    newmark_system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads, point_masses, gvec,
        force_scaling, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem,
        x0, v0, ω0, a0, α0, udot_init, θdot_init, Fdot_init, Mdot_init, Vdot_init, Ωdot_init, dt)

Populate the system residual vector `resid` for a Newmark scheme time-marching analysis.

# Arguments
 - `resid`: system residual vector
 - `x`: current states of the system
 - `assembly`: assembly of nonlinear beam elements
 - `prescribed_conditions`: dictionary of prescribed conditions
 - `distributed_loads`: dictionary of distributed loads
 - `point_masses`: dictionary of point masses 
 - `gvec`: gravity vector
 - `force_scaling`: scaling parameter for forces/moments
 - `irow_point`: row index of first equilibrium equation for each point
 - `irow_elem`: row index of first equation for just this beam element
 - `irow_elem1`: row index of first equation for the left side of each beam
 - `irow_elem2`: row index of first equation for the right side of each beam
 - `icol_point`: column index of first state variable for each point
 - `icol_elem`: column index of first state variable for each beam element
 - `x0`: body frame origin
 - `v0`: body frame linear velocity
 - `ω0`: body frame angular velocity
 - `a0`: body frame linear acceleration
 - `α0`: body frame angular acceleration
 - `udot_init`: `2/dt*u + udot` for each beam element from the previous time step
 - `θdot_init`: `2/dt*θ + θdot` for each beam element from the previous time step
 - `Fdot_init`: `2/dt*F + Fdot` for each beam element from the previous time step
 - `Mdot_init`: `2/dt*M + Mdot` for each beam element from the previous time step
 - `Vdot_init`: `2/dt*V + Vdot` for each beam element from the previous time step
 - `Ωdot_init`: `2/dt*Ω + Ωdot` for each beam element from the previous time step
 - `dt`: time step size
"""
function newmark_system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads, point_masses, gvec,
    force_scaling, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem,
    x0, v0, ω0, a0, α0, udot_init, θdot_init, Fdot_init, Mdot_init, Vdot_init, Ωdot_init, dt)

    nelem = length(assembly.elements)
    npoint = length(assembly.points)

    # add contributions to residual equations from the beam elements
    for ielem = 1:nelem

        # get pointers for element
        icol = icol_elem[ielem]
        irow_e = irow_elem[ielem]
        irow_e1 = irow_elem1[ielem]
        irow_p1 = irow_point[assembly.start[ielem]]
        irow_e2 = irow_elem2[ielem]
        irow_p2 = irow_point[assembly.stop[ielem]]

        newmark_element_residual!(resid, x, ielem, assembly.elements[ielem],
            distributed_loads, point_masses, gvec, force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2,
            x0, v0, ω0, a0, α0,
            udot_init[ielem], θdot_init[ielem], Fdot_init[ielem], Mdot_init[ielem],
            Vdot_init[ielem], Ωdot_init[ielem], dt)
    end

    # add contributions to the residual equations from the prescribed point conditions
    for ipoint = 1:npoint

        # skip if the unknowns have been eliminated from the system of equations
        if icol_point[ipoint] <= 0
            continue
        end

        icol = icol_point[ipoint]
        irow_p = irow_point[ipoint]

        point_residual!(resid, x, ipoint, assembly, prescribed_conditions,
            force_scaling, icol, irow_p, irow_elem1, irow_elem2)
    end

    return resid
end

"""
    dynamic_system_residual!(resid, dx, x, assembly, prescribed_conditions, distributed_loads, point_masses, gvec,
        force_scaling, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem,
        x0, v0, ω0, a0, α0)

Populate the system residual vector `resid` for a general dynamic system analysis.

# Arguments
 - `resid`: system residual vector
 - `dx`: current state rates of the system
 - `x`: current states of the system
 - `assembly`: assembly of nonlinear beam elements
 - `prescribed_conditions`: dictionary of prescribed conditions
 - `distributed_loads`: dictionary of distributed loads
 - `point_masses`: dictionary of point masses 
 - `gvec`: gravity vector
 - `force_scaling`: scaling parameter for forces/moments
 - `irow_point`: row index of first equilibrium equation for each point
 - `irow_elem`: row index of first equation for just this beam element
 - `irow_elem1`: row index of first equation for the left side of each beam
 - `irow_elem2`: row index of first equation for the right side of each beam
 - `icol_point`: column index of first state variable for each point
 - `icol_elem`: column index of first state variable for each beam element
 - `x0`: body frame origin
 - `v0`: body frame linear velocity
 - `ω0`: body frame angular velocity
 - `a0`: body frame linear acceleration
 - `α0`: body frame angular acceleration
"""
function dynamic_system_residual!(resid, dx, x, assembly, prescribed_conditions, distributed_loads, point_masses, gvec,
    force_scaling, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem,
    x0, v0, ω0, a0, α0)

    nelem = length(assembly.elements)
    npoint = length(assembly.points)

    # add contributions to residual equations from the beam elements
    for ielem = 1:nelem

        # get pointers for element
        icol = icol_elem[ielem]
        irow_e = irow_elem[ielem]
        irow_e1 = irow_elem1[ielem]
        irow_p1 = irow_point[assembly.start[ielem]]
        irow_e2 = irow_elem2[ielem]
        irow_p2 = irow_point[assembly.stop[ielem]]

        # set state rates for element
        udot = SVector(dx[icol], dx[icol+1], dx[icol+2])
        θdot = SVector(dx[icol+3], dx[icol+4], dx[icol+5])
        Fdot = SVector(dx[icol+6], dx[icol+7], dx[icol+8]) .* force_scaling
        Mdot = SVector(dx[icol+9], dx[icol+10], dx[icol+11]) .* force_scaling
        Vdot = SVector(dx[icol+12], dx[icol+13], dx[icol+14])
        Ωdot = SVector(dx[icol+15], dx[icol+16], dx[icol+17])

        dynamic_element_residual!(resid, x, ielem, assembly.elements[ielem],
             distributed_loads, point_masses, gvec, force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2,
             x0, v0, ω0, a0, α0, udot, θdot, Fdot, Mdot, Vdot, Ωdot)

    end

    # add contributions to the residual equations from the prescribed point conditions
    for ipoint = 1:npoint

        # skip if the unknowns have been eliminated from the system of equations
        if icol_point[ipoint] <= 0
            continue
        end

        icol = icol_point[ipoint]
        irow_p = irow_point[ipoint]

        point_residual!(resid, x, ipoint, assembly, prescribed_conditions,
            force_scaling, icol, irow_p, irow_elem1, irow_elem2)
    end

    return resid
end

"""
    static_system_jacobian!(jacob, x, assembly,
        prescribed_conditions, distributed_loads, point_masses, gvec, force_scaling,
        irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem)

Populate the system jacobian matrix `jacob` for a static analysis.

# Arguments
 - `jacob`: system jacobian matrix
 - `x`: current states of the system
 - `assembly`: assembly of nonlinear beam elements
 - `prescribed_conditions`: dictionary of prescribed conditions
 - `distributed_loads`: dictionary of distributed loads
 - `point_masses`: dictionary of point masses 
 - `gvec`: gravity vector
 - `force_scaling`: scaling parameter for forces/moments
 - `irow_point`: row index of first equilibrium equation for each point
 - `irow_elem1`: row index of first equation for the left side of each beam
 - `irow_elem2`: row index of first equation for the right side of each beam
 - `icol_point`: column index of first state variable for each point
 - `icol_elem`: column index of first state variable for each beam element
"""
@inline function static_system_jacobian!(jacob, x, assembly,
    prescribed_conditions, distributed_loads, point_masses, gvec, force_scaling,
    irow_point, irow_elem1, irow_elem2, icol_point, icol_elem)

    jacob .= 0

    npoint = length(assembly.points)
    nelem = length(assembly.elements)

    # add contributions to residual equations from the elements
    for ielem = 1:nelem

        icol = icol_elem[ielem]
        irow_e1 = irow_elem1[ielem]
        irow_p1 = irow_point[assembly.start[ielem]]
        irow_e2 = irow_elem2[ielem]
        irow_p2 = irow_point[assembly.stop[ielem]]

        jacob = static_element_jacobian!(jacob, x, ielem, assembly.elements[ielem],
            distributed_loads, point_masses, gvec, force_scaling, icol, irow_e1,
            irow_p1, irow_e2, irow_p2)
    end

    # add contributions to the system jacobian matrix from the prescribed point conditions
    for ipoint = 1:npoint

        # skip if the unknowns have been eliminated from the system of equations
        if icol_point[ipoint] <= 0
            continue
        end

        icol = icol_point[ipoint]
        irow_p = irow_point[ipoint]

        jacob = point_jacobian!(jacob, x, ipoint, assembly, prescribed_conditions,
            force_scaling, icol, irow_p, irow_elem1, irow_elem2)

    end

    return jacob
end

"""
    steady_state_system_jacobian!(jacob, x, assembly,
        prescribed_conditions, distributed_loads, point_masses, gvec, force_scaling,
        irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem,
        x0, v0, ω0, a0, α0)

Populate the system jacobian matrix `jacob` for a steady state analysis.

# Arguments
 - `jacob`: system jacobian matrix
 - `x`: current states of the system
 - `assembly`: assembly of nonlinear beam elements
 - `prescribed_conditions`: dictionary of prescribed conditions
 - `distributed_loads`: dictionary of distributed loads
 - `point_masses`: dictionary of point masses 
 - `gvec`: gravity vector
 - `force_scaling`: scaling parameter for forces/moments
 - `irow_point`: row index of first equilibrium equation for each point
 - `irow_elem`: row index of first equation for just this beam element
 - `irow_elem1`: row index of first equation for the left side of each beam
 - `irow_elem2`: row index of first equation for the right side of each beam
 - `icol_point`: column index of first state variable for each point
 - `icol_elem`: column index of first state variable for each beam element
 - `x0`: body frame origin
 - `v0`: body frame linear velocity
 - `ω0`: body frame angular velocity
 - `a0`: body frame linear acceleration
 - `α0`: body frame angular acceleration
"""
@inline function steady_state_system_jacobian!(jacob, x, assembly,
    prescribed_conditions, distributed_loads, point_masses, gvec, force_scaling,
    irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem,
    x0, v0, ω0, a0, α0)

    jacob .= 0

    npoint = length(assembly.points)
    nelem = length(assembly.elements)

    # add contributions to residual equations from the elements
    for ielem = 1:nelem

        icol = icol_elem[ielem]
        irow_e = irow_elem[ielem]
        irow_e1 = irow_elem1[ielem]
        irow_p1 = irow_point[assembly.start[ielem]]
        irow_e2 = irow_elem2[ielem]
        irow_p2 = irow_point[assembly.stop[ielem]]

        steady_state_element_jacobian!(jacob, x, ielem, assembly.elements[ielem],
            distributed_loads, point_masses, gvec, force_scaling, icol, irow_e, irow_e1,
            irow_p1, irow_e2, irow_p2, x0, v0, ω0, a0, α0)
    end

    # add contributions to the system jacobian matrix from the prescribed point conditions
    for ipoint = 1:npoint

        # skip if the unknowns have been eliminated from the system of equations
        if icol_point[ipoint] <= 0
            continue
        end

        icol = icol_point[ipoint]
        irow_p = irow_point[ipoint]

        point_jacobian!(jacob, x, ipoint, assembly, prescribed_conditions,
             force_scaling, icol, irow_p, irow_elem1, irow_elem2)
    end

    return jacob
end

"""
    initial_condition_system_jacobian!(jacob, x, assembly,
        prescribed_conditions, distributed_loads, point_masses, gvec, force_scaling,
        irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem,
        x0, v0, ω0, a0, α0, u, θ, udot, θdot, Fdot, Mdot)

Populate the system jacobian matrix `jacob` for an initial conditions analysis.

# Arguments
 - `jacob`: system jacobian matrix
 - `x`: current states of the system
 - `assembly`: assembly of nonlinear beam elements
 - `prescribed_conditions`: dictionary of prescribed conditions
 - `distributed_loads`: dictionary of distributed loads
 - `point_masses`: dictionary of point masses 
 - `gvec`: gravity vector
 - `force_scaling`: scaling parameter for forces/moments
 - `irow_point`: row index of first equilibrium equation for each point
 - `irow_elem`: row index of first equation for just this beam element
 - `irow_elem1`: row index of first equation for the left side of each beam
 - `irow_elem2`: row index of first equation for the right side of each beam
 - `icol_point`: column index of first state variable for each point
 - `icol_elem`: column index of first state variable for each beam element
 - `x0`: body frame origin
 - `v0`: body frame linear velocity
 - `ω0`: body frame angular velocity
 - `a0`: body frame linear acceleration
 - `α0`: body frame angular acceleration
 - `u`: initial linear deflections for each beam element
 - `θ`: initial angular deflections for each beam element
 - `udot`: initial linear deflection rates for each beam element
 - `θdot`: initial angular deflection rates for each beam element
 - `Fdot`: initial resultant force rates for each beam element
 - `Mdot`: initial resultant moment rates for each beam element
"""
@inline function initial_condition_system_jacobian!(jacob, x, assembly,
    prescribed_conditions, distributed_loads, point_masses, gvec, force_scaling,
    irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem,
    x0, v0, ω0, a0, α0, u, θ, udot, θdot, Fdot, Mdot)

    jacob .= 0

    npoint = length(assembly.points)
    nelem = length(assembly.elements)

    # add contributions to residual equations from the elements
    for ielem = 1:nelem

        icol = icol_elem[ielem]
        irow_e = irow_elem[ielem]
        irow_e1 = irow_elem1[ielem]
        irow_p1 = irow_point[assembly.start[ielem]]
        irow_e2 = irow_elem2[ielem]
        irow_p2 = irow_point[assembly.stop[ielem]]

        initial_condition_element_jacobian!(jacob, x, ielem, assembly.elements[ielem],
            distributed_loads, point_masses, gvec, force_scaling, icol, irow_e, irow_e1,
            irow_p1, irow_e2, irow_p2, x0, v0, ω0, a0, α0,
            u[ielem], θ[ielem], udot[ielem], θdot[ielem], Fdot[ielem], Mdot[ielem])
    end

    # add contributions to the system jacobian matrix from the prescribed point conditions
    for ipoint = 1:npoint

        # skip if the unknowns have been eliminated from the system of equations
        if icol_point[ipoint] <= 0
            continue
        end

        icol = icol_point[ipoint]
        irow_p = irow_point[ipoint]

        point_jacobian!(jacob, x, ipoint, assembly, prescribed_conditions,
            force_scaling, icol, irow_p, irow_elem1, irow_elem2)
    end

    return jacob
end

"""
    newmark_system_jacobian!(jacob, x, assembly, prescribed_conditions, distributed_loads, point_masses, gvec,
        force_scaling, irow_point, irow_elem, irow_elem1, irow_elem2,
        icol_point, icol_elem, x0, v0, ω0, a0, α0, udot_init, θdot_init, Fdot_init, 
        Mdot_init, Vdot_init, Ωdot_init, dt)

Populate the system jacobian matrix `jacob` for a Newmark scheme time-marching analysis.

# Arguments
 - `jacob`: system jacobian matrix
 - `x`: current states of the system
 - `assembly`: assembly of nonlinear beam elements
 - `prescribed_conditions`: dictionary of prescribed conditions
 - `distributed_loads`: dictionary of distributed loads
 - `point_masses`: dictionary of point masses 
 - `gvec`: gravity vector
 - `force_scaling`: scaling parameter for forces/moments
 - `irow_point`: row index of first equilibrium equation for each point
 - `irow_elem`: row index of first equation for just this beam element
 - `irow_elem1`: row index of first equation for the left side of each beam
 - `irow_elem2`: row index of first equation for the right side of each beam
 - `icol_point`: column index of first state variable for each point
 - `icol_elem`: column index of first state variable for each beam element
 - `x0`: body frame origin
 - `v0`: body frame linear velocity
 - `ω0`: body frame angular velocity
 - `a0`: body frame linear acceleration
 - `α0`: body frame angular acceleration
 - `udot_init`: `2/dt*u + udot` for each beam element from the previous time step
 - `θdot_init`: `2/dt*θ + θdot` for each beam element from the previous time step
 - `Fdot_init`: `2/dt*F + Fdot` for each beam element from the previous time step
 - `Mdot_init`: `2/dt*M + Mdot` for each beam element from the previous time step
 - `Vdot_init`: `2/dt*V + Vdot` for each beam element from the previous time step
 - `Ωdot_init`: `2/dt*Ω + Ωdot` for each beam element from the previous time step
 - `dt`: time step size
"""
@inline function newmark_system_jacobian!(jacob, x, assembly,
    prescribed_conditions, distributed_loads, point_masses, gvec, force_scaling,
    irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem,
    x0, v0, ω0, a0, α0, udot_init, θdot_init, Fdot_init, Mdot_init, Vdot_init, Ωdot_init, dt)

    jacob .= 0

    npoint = length(assembly.points)
    nelem = length(assembly.elements)

    # add contributions to residual equations from the elements
    for ielem = 1:nelem

        icol = icol_elem[ielem]
        irow_e = irow_elem[ielem]
        irow_e1 = irow_elem1[ielem]
        irow_p1 = irow_point[assembly.start[ielem]]
        irow_e2 = irow_elem2[ielem]
        irow_p2 = irow_point[assembly.stop[ielem]]

        newmark_element_jacobian!(jacob, x, ielem, assembly.elements[ielem],
            distributed_loads, point_masses, gvec, force_scaling, icol, irow_e, irow_e1,
            irow_p1, irow_e2, irow_p2, x0, v0, ω0, a0, α0,
            udot_init[ielem], θdot_init[ielem], Fdot_init[ielem], Mdot_init[ielem],
            Vdot_init[ielem], Ωdot_init[ielem], dt)
    end

    # add contributions to the system jacobian matrix from the prescribed point conditions
    for ipoint = 1:npoint

        # skip if the unknowns have been eliminated from the system of equations
        if icol_point[ipoint] <= 0
            continue
        end

        icol = icol_point[ipoint]
        irow_p = irow_point[ipoint]

        point_jacobian!(jacob, x, ipoint, assembly, prescribed_conditions,
            force_scaling, icol, irow_p, irow_elem1, irow_elem2)
    end

    # # zero out near-zero values ( < eps() )
    # jacob = droptol!(jacob, eps(eltype(jacob)))

    return jacob
end

"""
    dynamic_system_jacobian!(jacob, dx, x, assembly,
        prescribed_conditions, distributed_loads, point_masses, gvec, force_scaling,
        irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem,
        x0, v0, ω0, a0, α0)

Populate the jacobian matrix `jacob` for a general dynamic analysis.

# Arguments
 - `jacob`: system jacobian matrix
 - `dx`: current state rates of the system
 - `x`: current states of the system
 - `assembly`: assembly of nonlinear beam elements
 - `prescribed_conditions`: dictionary of prescribed conditions
 - `distributed_loads`: dictionary of distributed loads
 - `point_masses`: dictionary of point masses 
 - `gvec`: gravity vector
 - `force_scaling`: scaling parameter for forces/moments
 - `irow_point`: row index of first equilibrium equation for each point
 - `irow_elem`: row index of first equation for just this beam element
 - `irow_elem1`: row index of first equation for the left side of each beam
 - `irow_elem2`: row index of first equation for the right side of each beam
 - `icol_point`: column index of first state variable for each point
 - `icol_elem`: column index of first state variable for each beam element
 - `x0`: body frame origin
 - `v0`: body frame linear velocity
 - `ω0`: body frame angular velocity
 - `a0`: body frame linear acceleration
 - `α0`: body frame angular acceleration
"""
@inline function dynamic_system_jacobian!(jacob, dx, x, assembly,
    prescribed_conditions, distributed_loads, point_masses, gvec, force_scaling,
    irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem,
    x0, v0, ω0, a0, α0)

    jacob .= 0

    npoint = length(assembly.points)
    nelem = length(assembly.elements)

    # add contributions to residual equations from the elements
    for ielem = 1:nelem

        icol = icol_elem[ielem]
        irow_e = irow_elem[ielem]
        irow_e1 = irow_elem1[ielem]
        irow_p1 = irow_point[assembly.start[ielem]]
        irow_e2 = irow_elem2[ielem]
        irow_p2 = irow_point[assembly.stop[ielem]]

        # set state rates for element
        udot = SVector(dx[icol], dx[icol+1], dx[icol+2])
        θdot = SVector(dx[icol+3], dx[icol+4], dx[icol+5])
        Fdot = SVector(dx[icol+6], dx[icol+7], dx[icol+8]) .* force_scaling
        Mdot = SVector(dx[icol+9], dx[icol+10], dx[icol+11]) .* force_scaling
        Vdot = SVector(dx[icol+12], dx[icol+13], dx[icol+14])
        Ωdot = SVector(dx[icol+15], dx[icol+16], dx[icol+17])

        dynamic_element_jacobian!(jacob, x, ielem, assembly.elements[ielem],
            distributed_loads, point_masses, gvec, force_scaling, icol, irow_e, irow_e1,
            irow_p1, irow_e2, irow_p2, x0, v0, ω0, a0, α0,
            udot, θdot, Fdot, Mdot, Vdot, Ωdot)
    end

    # add contributions to the system jacobian matrix from the prescribed point conditions
    for ipoint = 1:npoint

        # skip if the unknowns have been eliminated from the system of equations
        if icol_point[ipoint] <= 0
            continue
        end

        icol = icol_point[ipoint]
        irow_p = irow_point[ipoint]

        point_jacobian!(jacob, x, ipoint, assembly, prescribed_conditions,
            force_scaling, icol, irow_p, irow_elem1, irow_elem2)
    end

    # # zero out near-zero values ( < eps() )
    # jacob = droptol!(jacob, eps(eltype(jacob)))

    return jacob
end

"""
    system_mass_matrix!(jacob, x, assembly, point_masses, force_scaling,
        irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem)

Populate the system mass matrix for a general dynamic analysis

# Arguments
 - `jacob`: jacobian of the residuals with respect to the state rates
 - `x`: current states of the system
 - `assembly`: assembly of nonlinear beam elements
 - `point_masses`: dictionary of point masses 
 - `force_scaling`: scaling parameter for forces/moments
 - `irow_point`: row index of the first equilibrium equation for each point
 - `irow_elem`: row index of the first linear/angular velocity residual for each element
 - `irow_elem1`: row index of the first equation for the left side of each beam
 - `irow_elem2`: row index of the first equation for the right side of each beam
 - `icol_point`: column index of the first state variable for each point
 - `icol_elem`: column index of the first state variable for each beam element
"""
function system_mass_matrix!(jacob, x, assembly, point_masses, force_scaling, 
    irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem)

    jacob .= 0

    npoint = length(assembly.points)
    nelem = length(assembly.elements)

    # add contributions to "mass matrix" from the beam elements
    for ielem = 1:nelem

        icol = icol_elem[ielem]
        irow_e = irow_elem[ielem]
        irow_e1 = irow_elem1[ielem]
        irow_p1 = irow_point[assembly.start[ielem]]
        irow_e2 = irow_elem2[ielem]
        irow_p2 = irow_point[assembly.stop[ielem]]

        element_mass_matrix!(jacob, x, ielem, assembly.elements[ielem], point_masses, 
            force_scaling, icol, irow_e, irow_p1, irow_p2)
    end

    # no contributions to "mass matrix" from point state variables

    # # zero out near-zero values ( < eps() )
    # jacob = droptol!(jacob, eps(eltype(jacob)))

    return jacob
end

"""
    system_mass_matrix!(jacob, gamma, x, dx, assembly, point_masses, force_scaling,
        irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem)

Add the system mass matrix to `jacob`, scaled by the scaling parameter `gamma`.
"""
function system_mass_matrix!(jacob, gamma, x, assembly, point_masses, force_scaling, 
    irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem)

    npoint = length(assembly.points)
    nelem = length(assembly.elements)

    # add contributions to "mass matrix" from the beam elements
    for ielem = 1:nelem

        icol = icol_elem[ielem]
        irow_e = irow_elem[ielem]
        irow_e1 = irow_elem1[ielem]
        irow_p1 = irow_point[assembly.start[ielem]]
        irow_e2 = irow_elem2[ielem]
        irow_p2 = irow_point[assembly.stop[ielem]]

        element_mass_matrix!(jacob, gamma, x, ielem, assembly.elements[ielem], point_masses,
            force_scaling, icol, irow_e, irow_p1, irow_p2)
    end

    # no contributions to "mass matrix" from point state variables

    return jacob
end
