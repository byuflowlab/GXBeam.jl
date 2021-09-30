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
 - `mass_scaling`: Scaling for state variables corresponding to masses/inertias
 - `irow_point`: Row index of first equilibrium equation for each point
 - `irow_elem`: Row index of first equation for just this beam element
 - `irow_elem1`: Row index of first equation for the left side of each beam
 - `irow_elem2`: Row index of first equation for the right side of each beam
 - `icol_point`: Row/Column index of first state variable for each point
 - `icol_elem`: Row/Column index of first state variable for each beam element
 - `udot`: `2/dt*u + udot` for each beam element from the previous time step
 - `θdot`: `2/dt*θ + θdot` for each beam element from the previous time step
 - `Pdot`: `2/dt*C'*Cab*P + C'*Cab*Pdot` for each beam element from the previous time step
 - `Hdot`: `2/dt*C'*Cab*H + C'*Cab*Hdot` for each beam element from the previous time step
 - `t`: Current system time
"""
mutable struct System{TF, TV<:AbstractVector{TF}, TM<:AbstractMatrix{TF}}
    static::Bool
    x::TV
    r::TV
    K::TM
    M::TM
    force_scaling::TF
    mass_scaling::TF
    irow_point::Vector{Int}
    irow_elem::Vector{Int}
    irow_elem1::Vector{Int}
    irow_elem2::Vector{Int}
    icol_point::Vector{Int}
    icol_elem::Vector{Int}
    udot::Vector{SVector{3,TF}}
    θdot::Vector{SVector{3,TF}}
    Pdot::Vector{SVector{3,TF}}
    Hdot::Vector{SVector{3,TF}}
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
 - `mass_scaling`: Factor used to scale system mass/inertia internally.  If not
    specified, a suitable default will be chosen based on the entries of the
    mass matrix.

Note that points with prescribed conditions must be included in the system of
equations.
"""
function System(assembly, static; kwargs...)

    return System(eltype(assembly), assembly, static; kwargs...)
end

function System(TF, assembly, static;
    prescribed_points = 1:length(assembly.points),
    force_scaling = default_force_scaling(assembly),
    mass_scaling = default_mass_scaling(assembly)
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
    Pdot = [@SVector zeros(TF, 3) for i = 1:nelem]
    Hdot = [@SVector zeros(TF, 3) for i = 1:nelem]

    # initialize current time
    t = 0.0

    # set system types
    TV = promote_type(typeof(x), typeof(r))
    TM = promote_type(typeof(K), typeof(M))

    return System{TF, TV, TM}(static, x, r, K, M, force_scaling, mass_scaling,
        irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem, udot, θdot,
        Pdot, Hdot, t)
end

function default_force_scaling(assembly)

    TF = eltype(assembly)

    # create vector of all compliance entries
    compliance_entries = vcat([vcat(elem.C11..., elem.C12..., elem.C22...) for elem in assembly.elements]...)

    # get all nonzero entries
    compliance_nonzero_indices = findall(xi -> abs(xi) > eps(TF), compliance_entries)

    # set force scaling based on nonzero compliance matrix entries
    if isempty(compliance_nonzero_indices)
        force_scaling = 1.0
    else
        nonzero_compliance_entries = compliance_entries[compliance_nonzero_indices]
        force_scaling = nextpow(2.0, length(nonzero_compliance_entries)/sum(nonzero_compliance_entries)/100)
    end

    return force_scaling
end

function default_mass_scaling(assembly)

    TF = eltype(assembly)

    # create vector of all inverse mass matrix entries
    minv_entries = vcat([vcat(elem.minv11..., elem.minv12..., elem.minv22...) for elem in assembly.elements]...)

    # get all nonzero entries
    minv_nonzero_indices = findall(xi -> abs(xi) > eps(TF), minv_entries)

    # set mass scaling based on nonzero inverse mass matrix entries
    if isempty(minv_nonzero_indices)
        mass_scaling = 1.0
    else
        nonzero_minv_entries = minv_entries[minv_nonzero_indices]
        mass_scaling = nextpow(2.0, length(nonzero_minv_entries)/sum(nonzero_minv_entries))
    end

    return mass_scaling
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
 - `u_e`: Vector containing initial displacements (in the global frame) for each beam element.
 - `theta_e`: Vector containing rotation variables (in the global frame) for each beam element.
 - `F_e`: Vector containing resultant forces (in the beam coordinate frame) for each beam element.
 - `M_e`: Vector containing resultant moments (in the beam coordinate frame) for each beam element.
 - `P_e`: Vector containing linear momenta (in the beam coordinate frame) for each beam element.
 - `H_e`: Vector containing angular momenta (in the beam coordinate frame) for each beam element.
 - `u_p`: Vector containing initial displacements (in the global frame) for each point.
 - `theta_p`: Vector containing rotation variables (in the global frame) for each point.
 - `F_p`: Vector containing externally applied forces (in the global frame) for each point.
 - `M_p`: Vector containing externally applied moments (in the global frame) for each point.
"""
set_state!

function set_state!(system, prescribed_conditions; kwargs...)
    x = set_state!(system.x, system, prescribed_conditions; kwargs...)
    return system
end

function set_state!(x, system, prescribed_conditions; kwargs...)
    return set_state!(x, system.icol_elem, system.icol_point; kwargs...)
end

function set_state!(x, icol_elem, icol_point, force_scaling, mass_scaling,
    prescribed_conditions; u_e = nothing, theta_e = nothing,
    F_e = nothing, M_e = nothing, P_e = nothing, H_e = nothing,
    u_p = nothing, theta_p = nothing, F_p = nothing, M_p = nothing)

    nelem = length(icol_elem)
    npoint = length(icol_point)

    for ielem = 1:nelem
        icol = icol_elem[ielem]
        if !isnothing(u_e)
            set_element_deflections!(x, icol, u_e)
        end
        if !isnothing(theta_e)
            set_element_rotations!(x, icol, theta_e)
        end
        if !isnothing(F_e)
            set_element_forces!(x, icol, F_e, force_scaling)
        end
        if !isnothing(M_e)
            set_element_moments!(x, icol, M_e, force_scaling)
        end
        if !isnothing(P_e)
            set_element_linear_momenta!(x, icol, P_e, mass_scaling)
        end
        if !isnothing(H_e)
            set_element_angular_momenta!(x, icol, H_e, mass_scaling)
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
    set_element_deflections!([x,] system, u_e, ielem)

Set the state variables in `system` (or in the vector `x`) corresponding to the
deflections of element `ielem` to the provided values.
"""
function set_element_deflections!(system::System, u_e, ielem)
    set_element_deflections!(system.x, system, u_e, ielem)
    return system
end

function set_element_deflections!(x, system::System, u_e, ielem)
    icol = system.icol_elem[ielem]
    set_element_deflections!(x, icol, u_e)
    return x
end

function set_element_deflections!(x, icol, u_e)
    x[icol:icol+2] .= u_e
    return x
end

"""
    set_element_rotations!([x,] system, θ_b, ielem)

Set the state variables in `system` (or in the vector `x`) corresponding to the
rotations of element `ielem` to the provided values.
"""
function set_element_rotations!(system::System, θ_e, ielem)
    set_element_rotations!(system.x, system, θ_e, ielem)
    return system
end

function set_element_rotations!(x, system::System, θ_e, ielem)
    icol = system.icol_elem[ielem]
    set_element_rotations!(x, icol, θ_e)
    return x
end

function set_element_rotations!(x, icol, θ_e)
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
    set_element_linear_momenta!([x,] system, u_e, ielem)

Set the state variables in `system` (or in the vector `x`) corresponding to the
linear_momenta of element `ielem` to the provided values.
"""
function set_element_linear_momenta!(system::System, P_e, ielem)
    set_element_linear_momenta!(x, system, P_e, ielem)
    return system
end

function set_element_linear_momenta!(x, system::System, P_e, ielem)
    @assert !system.static
    icol = system.icol_elem[ielem]
    mass_scaling = system.mass_scaling
    set_element_linear_momenta!(x, icol, P_e, mass_scaling)
    return x
end

function set_element_linear_momenta!(x, icol, P_e, mass_scaling)
    x[icol+12:icol+14] .= P_e ./ mass_scaling
    return x
end

"""
   set_element_angular_momenta!([x,] system, u_e, ielem)

Set the state variables in `system` (or in the vector `x`) corresponding to the
angular_momenta of element `ielem` to the provided values.
"""
function set_element_angular_momenta!(system::System, H_e, ielem)
    set_element_angular_momenta!(x, system, H_e, ielem)
    return system
end

function set_element_angular_momenta!(x, system::System, H_e, ielem)
    @assert !system.static
    icol = system.icol_elem[ielem]
    mass_scaling = system.mass_scaling
    set_element_angular_momenta!(x, icol, H_e, mass_scaling)
    return x
end

function set_element_angular_momenta!(x, icol, H_e, mass_scaling)
    x[icol+15:icol+17] .= H_e ./ mass_scaling
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
    system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads, gvec,
        force_scaling, mass_scaling, irow_point, irow_elem, irow_elem1, irow_elem2,
        icol_point, icol_elem)
    system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads, gvec,
        force_scaling, mass_scaling, irow_point, irow_elem, irow_elem1, irow_elem2,
        icol_point, icol_elem, x0, v0, ω0)
    system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads, gvec,
        force_scaling, mass_scaling, irow_point, irow_elem, irow_elem1, irow_elem2,
        icol_point, icol_elem, x0, v0, ω0, u, θ, udot, θdot)
    system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads, gvec,
        force_scaling, mass_scaling, irow_point, irow_elem, irow_elem1, irow_elem2,
        icol_point, icol_elem, x0, v0, ω0, udot_init, θdot_init, CtCabPdot_init,
        CtCabHdot_init, dt)

Populate the residual vector `resid` with the results of the residual equations
for the system.

There are four implementations corresponding to the following analysis types:
 - Static
 - Dynamic - Steady State
 - Dynamic - Initial Step (for initializing time domain simulations)
 - Dynamic - Time Marching

 See "GEBT: A general-purpose nonlinear analysis tool for composite beams" by
 Wenbin Yu and "Geometrically nonlinear analysis of composite beams using
 Wiener-Milenković parameters" by Qi Wang and Wenbin Yu.

# Arguments
 - `resid`: System residual vector
 - `x`: Current state variables of the system
 - `assembly`: Assembly of rigidly connected nonlinear beam elements
 - `prescribed_conditions`: Dictionary of prescribed conditions at all time steps
 - `distributed_loads`: Dictionary of distributed loads at all time steps
 - `gvec`: Gravity vector
 - `force_scaling`: Scaling parameter for forces/moments
 - `mass_scaling`: Scaling parameter for masses/inertias
 - `irow_point`: Row index of first equilibrium equation for each point
 - `irow_elem`: Row index of first equation for just this beam element
 - `irow_elem1`: Row index of first equation for the left side of each beam
 - `irow_elem2`: Row index of first equation for the right side of each beam
 - `icol_point`: Column index of first state variable for each point
 - `icol_elem`: Column index of first state variable for each beam element

# Additional Arguments for Dynamic Analyses
 - `x0`: Global frame origin (for the current time step)
 - `v0`: Global frame linear velocity (for the current time step)
 - `ω0`: Global frame angular velocity (for the current time step)

# Additional Arguments for Initial Step Analyses
 - `u`: deflection variables for each beam element
 - `θ`: rotation variables for each beam element
 - `udot`: time derivative of u for each beam element
 - `θdot`: time derivative of θ for each beam element

# Additional Arguments for Time Marching Analyses
 - `udot_init`: `-2/dt*u - udot` for each beam element from the previous time step
 - `θdot_init`: `-2/dt*θ - θdot` for each beam element from the previous time step
 - `CtCabPdot_init`: `-2/dt*C'*Cab*P - C'*Cab*Pdot` for each beam element from the previous time step
 - `CtCabHdot_init`: `-2/dt*C'*Cab*H - C'*Cab*Hdot` for each beam element from the previous time step
 - `dt`: Time step size. If set to `nothing`, `udot_init`,
"""
system_residual!

# static
function system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads, gvec,
    force_scaling, mass_scaling, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem)

    return static_system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads, gvec,
        force_scaling, mass_scaling, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem)
end

# dynamic - steady state
function system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads, gvec,
    force_scaling, mass_scaling, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem,
    x0, v0, ω0)

    return steady_state_system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads, gvec,
        force_scaling, mass_scaling, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem,
        x0, v0, ω0)
end

# dynamic - initial step
function system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads, gvec,
    force_scaling, mass_scaling, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem,
    x0, v0, ω0, u, θ, udot, θdot)

    return initial_step_system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads, gvec,
    force_scaling, mass_scaling, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem,
    x0, v0, ω0, u, θ, udot, θdot)
end

# dynamic - newmark scheme time marching
function system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads, gvec,
    force_scaling, mass_scaling, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem,
    x0, v0, ω0, udot_init, θdot_init, CtCabPdot_init, CtCabHdot_init, dt)

    # dynamic - newmark scheme time marching
    return newmark_system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads, gvec,
        force_scaling, mass_scaling, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem,
        x0, v0, ω0, udot_init, θdot_init, CtCabPdot_init, CtCabHdot_init, dt)
end

# static
function static_system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads, gvec,
    force_scaling, mass_scaling, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem)

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

        resid = static_element_residual!(resid, x, ielem, assembly.elements[ielem],
             distributed_loads, gvec, force_scaling, mass_scaling, icol, irow_e, irow_e1,
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

        resid = point_residual!(resid, x, ipoint, assembly, prescribed_conditions,
            force_scaling, icol, irow_p, irow_elem1, irow_elem2)
    end

    return resid
end

# dynamic - steady state
function steady_state_system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads, gvec,
    force_scaling, mass_scaling, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem,
    x0, v0, ω0)

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

        resid = steady_state_element_residual!(resid, x, ielem, assembly.elements[ielem],
             distributed_loads, gvec, force_scaling, mass_scaling, icol, irow_e, irow_e1,
            irow_p1, irow_e2, irow_p2, x0, v0, ω0)
    end

    # add contributions to the residual equations from the prescribed point conditions
    for ipoint = 1:npoint

        # skip if the unknowns have been eliminated from the system of equations
        if icol_point[ipoint] <= 0
            continue
        end

        icol = icol_point[ipoint]
        irow_p = irow_point[ipoint]

        resid = point_residual!(resid, x, ipoint, assembly, prescribed_conditions,
            force_scaling, icol, irow_p, irow_elem1, irow_elem2)
    end

    return resid
end

# dynamic - initial step
function initial_step_system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads, gvec,
    force_scaling, mass_scaling, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem,
    x0, v0, ω0, u, θ, udot, θdot)

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

        resid = initial_step_element_residual!(resid, x, ielem, assembly.elements[ielem],
             distributed_loads, gvec, force_scaling, mass_scaling, icol, irow_e, irow_e1,
            irow_p1, irow_e2, irow_p2, x0, v0, ω0,
            u[ielem], θ[ielem], udot[ielem], θdot[ielem])
    end

    # add contributions to the residual equations from the prescribed point conditions
    for ipoint = 1:npoint

        # skip if the unknowns have been eliminated from the system of equations
        if icol_point[ipoint] <= 0
            continue
        end

        icol = icol_point[ipoint]
        irow_p = irow_point[ipoint]

        resid = point_residual!(resid, x, ipoint, assembly, prescribed_conditions,
            force_scaling, icol, irow_p, irow_elem1, irow_elem2)
    end

    return resid
end

# dynamic - newmark scheme time marching
function newmark_system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads, gvec,
    force_scaling, mass_scaling, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem,
    x0, v0, ω0, udot_init, θdot_init, CtCabPdot_init, CtCabHdot_init, dt)

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

        resid = newmark_element_residual!(resid, x, ielem, assembly.elements[ielem],
             distributed_loads, gvec, force_scaling, mass_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2,
             x0, v0, ω0,
             udot_init[ielem], θdot_init[ielem],
             CtCabPdot_init[ielem], CtCabHdot_init[ielem], dt)
    end

    # add contributions to the residual equations from the prescribed point conditions
    for ipoint = 1:npoint

        # skip if the unknowns have been eliminated from the system of equations
        if icol_point[ipoint] <= 0
            continue
        end

        icol = icol_point[ipoint]
        irow_p = irow_point[ipoint]

        resid = point_residual!(resid, x, ipoint, assembly, prescribed_conditions,
            force_scaling, icol, irow_p, irow_elem1, irow_elem2)
    end

    return resid
end

# dynamic - general
function dynamic_system_residual!(resid, x, dx, assembly, prescribed_conditions, distributed_loads, gvec,
    force_scaling, mass_scaling, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem,
    x0, v0, ω0)

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
        Pdot = SVector(dx[icol+12], dx[icol+13], dx[icol+14]) .* mass_scaling
        Hdot = SVector(dx[icol+15], dx[icol+16], dx[icol+17]) .* mass_scaling

        resid = dynamic_element_residual!(resid, x, ielem, assembly.elements[ielem],
             distributed_loads, gvec, force_scaling, mass_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2,
             x0, v0, ω0, udot, θdot, Pdot, Hdot)

    end

    # add contributions to the residual equations from the prescribed point conditions
    for ipoint = 1:npoint

        # skip if the unknowns have been eliminated from the system of equations
        if icol_point[ipoint] <= 0
            continue
        end

        icol = icol_point[ipoint]
        irow_p = irow_point[ipoint]

        resid = point_residual!(resid, x, ipoint, assembly, prescribed_conditions,
            force_scaling, icol, irow_p, irow_elem1, irow_elem2)
    end

    return resid
end

"""
    system_jacobian!(jacob, x, assembly, prescribed_conditions, distributed_loads, gvec,
        force_scaling, mass_scaling, irow_point, irow_elem, irow_elem1, irow_elem2,
        icol_point, icol_elem)
    system_jacobian!(jacob, x, assembly, prescribed_conditions, distributed_loads, gvec,
        force_scaling, mass_scaling, irow_point, irow_elem, irow_elem1, irow_elem2,
        icol_point, icol_elem, x0, v0, ω0)
    system_jacobian!(jacob, x, assembly, prescribed_conditions, distributed_loads, gvec,
        force_scaling, mass_scaling, irow_point, irow_elem, irow_elem1, irow_elem2,
        icol_point, icol_elem, x0, v0, ω0, u, θ, udot, θdot)
    system_jacobian!(jacob, x, assembly, prescribed_conditions, distributed_loads, gvec,
        force_scaling, mass_scaling, irow_point, irow_elem, irow_elem1, irow_elem2,
        icol_point, icol_elem, x0, v0, ω0, udot_init, θdot_init, CtCabPdot_init,
        CtCabHdot_init, dt)

Populate the jacobian matrix `jacob` with the jacobian of the residual vector
with respect to the state variables.

There are four implementations corresponding to the following analysis types:
 - Static
 - Dynamic - Steady State
 - Dynamic - Initial Step (for initializing time domain simulations)
 - Dynamic - Time Marching

 See "GEBT: A general-purpose nonlinear analysis tool for composite beams" by
 Wenbin Yu and "Geometrically nonlinear analysis of composite beams using
 Wiener-Milenković parameters" by Qi Wang and Wenbin Yu.

# Arguments
 - `jacob`: Jacobian matrix
 - `x`: Vector containing current state variables of the system
 - `assembly`: Assembly of rigidly connected nonlinear beam elements
 - `prescribed_conditions`: Dictionary of prescribed conditions at all time steps
 - `distributed_loads`: Dictionary of distributed loads at all time steps
 - `gvec`: Gravity vector
 - `force_scaling`: Scaling parameter for forces/moments
 - `mass_scaling`: Scaling parameter for masses/inertia
 - `irow_point`: Row index of first equilibrium equation for each point
 - `irow_elem`: Row index of first equation for just this beam element
 - `irow_elem1`: Row index of first equation for the left side of each beam
 - `irow_elem2`: Row index of first equation for the right side of each beam
 - `icol_point`: Column index of first state variable for each point
 - `icol_elem`: Column index of first state variable for each beam element

# Additional Arguments for Dynamic Analyses
 - `x0`: Global frame origin (for the current time step)
 - `v0`: Global frame linear velocity (for the current time step)
 - `ω0`: Global frame angular velocity (for the current time step)

# Additional Arguments for Initial Step Analyses
 - `u`: deflection variables for each beam element
 - `θ`: rotation variables for each beam element
 - `udot`: time derivative of u for each beam element
 - `θdot`: time derivative of θ for each beam element

# Additional Arguments for Time Marching Analyses
 - `udot`: `-2/dt*u - udot` for each beam element from the previous time step
 - `θdot_init`: `-2/dt*θ - θdot` for each beam element from the previous time step
 - `CtCabPdot`: `-2/dt*C'*Cab*P - C'*Cab*Pdot` for each beam element from the previous time step
 - `CtCabHdot`: `-2/dt*C'*Cab*H - C'*Cab*Hdot` for each beam element from the previous time step
 - `dt`: time step size
"""
system_jacobian!

# static
@inline function system_jacobian!(jacob, x, assembly, prescribed_conditions, distributed_loads, gvec,
    force_scaling, mass_scaling, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem)

    return static_system_jacobian!(jacob, x, assembly, prescribed_conditions, distributed_loads, gvec,
        force_scaling, mass_scaling, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem)
end

# dynamic - steady state
@inline function system_jacobian!(jacob, x, assembly, prescribed_conditions, distributed_loads, gvec,
    force_scaling, mass_scaling, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem,
    x0, v0, ω0)

    return steady_state_system_jacobian!(jacob, x, assembly, prescribed_conditions, distributed_loads, gvec,
        force_scaling, mass_scaling, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem,
        x0, v0, ω0)
end

# dynamic - initial step
@inline function system_jacobian!(jacob, x, assembly, prescribed_conditions, distributed_loads, gvec,
    force_scaling, mass_scaling, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem,
    x0, v0, ω0, u, θ, udot, θdot)

    return initial_step_system_jacobian!(jacob, x, assembly, prescribed_conditions, distributed_loads, gvec,
        force_scaling, mass_scaling, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem,
        x0, v0, ω0, u, θ, udot, θdot)
end

# dynamic - newmark scheme time marching
@inline function system_jacobian!(jacob, x, assembly, prescribed_conditions,
    distributed_loads, gvec, force_scaling, mass_scaling,
    irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem,
    x0, v0, ω0, udot_init, θdot_init, CtCabPdot_init, CtCabHdot_init, dt)

    return newmark_system_jacobian!(jacob, x, assembly, prescribed_conditions,
        distributed_loads, gvec, force_scaling, mass_scaling,
        irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem,
        x0, v0, ω0, udot_init, θdot_init, CtCabPdot_init, CtCabHdot_init, dt)
end

# static
@inline function static_system_jacobian!(jacob, x, assembly,
    prescribed_conditions, distributed_loads, gvec, force_scaling, mass_scaling,
    irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem)

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

        jacob = element_jacobian!(jacob, x, ielem, assembly.elements[ielem],
            distributed_loads, gvec, force_scaling, mass_scaling, icol, irow_e, irow_e1,
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

# dynamic - steady state
@inline function steady_state_system_jacobian!(jacob, x, assembly,
    prescribed_conditions, distributed_loads, gvec, force_scaling, mass_scaling,
    irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem,
    x0, v0, ω0)

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

        jacob = element_jacobian!(jacob, x, ielem, assembly.elements[ielem],
            distributed_loads, gvec, force_scaling, mass_scaling, icol, irow_e, irow_e1,
            irow_p1, irow_e2, irow_p2, x0, v0, ω0)
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

# dynamic - initial step
@inline function initial_step_system_jacobian!(jacob, x, assembly,
    prescribed_conditions, distributed_loads, gvec, force_scaling, mass_scaling,
    irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem,
    x0, v0, ω0, u, θ, udot, θdot)

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

        jacob = element_jacobian!(jacob, x, ielem, assembly.elements[ielem],
            distributed_loads, gvec, force_scaling, mass_scaling, icol, irow_e, irow_e1,
            irow_p1, irow_e2, irow_p2, x0, v0, ω0,
            u[ielem], θ[ielem], udot[ielem], θdot[ielem])
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

# dynamic - newmark scheme time marching
@inline function newmark_system_jacobian!(jacob, x, assembly,
    prescribed_conditions, distributed_loads, gvec, force_scaling, mass_scaling,
    irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem,
    x0, v0, ω0, udot_init, θdot_init, CtCabPdot_init, CtCabHdot_init, dt)

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

        jacob = element_jacobian!(jacob, x, ielem, assembly.elements[ielem],
            distributed_loads, gvec, force_scaling, mass_scaling, icol, irow_e, irow_e1,
            irow_p1, irow_e2, irow_p2, x0, v0, ω0,
            udot_init[ielem], θdot_init[ielem],
            CtCabPdot_init[ielem], CtCabHdot_init[ielem], dt)
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

    # # zero out near-zero values ( < eps() )
    # jacob = droptol!(jacob, eps(eltype(jacob)))

    return jacob
end

# dynamic - general
@inline function dynamic_system_jacobian!(jacob, x, dx, assembly,
    prescribed_conditions, distributed_loads, gvec, force_scaling, mass_scaling,
    irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem,
    x0, v0, ω0)

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
        Pdot = SVector(dx[icol+12], dx[icol+13], dx[icol+14]) .* mass_scaling
        Hdot = SVector(dx[icol+15], dx[icol+16], dx[icol+17]) .* mass_scaling

        jacob = dynamic_element_jacobian!(jacob, x, ielem, assembly.elements[ielem],
            distributed_loads, gvec, force_scaling, mass_scaling, icol, irow_e, irow_e1,
            irow_p1, irow_e2, irow_p2, x0, v0, ω0,
            udot, θdot, Pdot, Hdot)
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

    # # zero out near-zero values ( < eps() )
    # jacob = droptol!(jacob, eps(eltype(jacob)))

    return jacob
end

"""
    system_mass_matrix!(jacob, x, assembly, force_scaling, mass_scaling,
        irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem)

Populate the system "mass matrix", the jacobian of the residual vector
with respect to the time derivatives of the state variables.

See "GEBT: A general-purpose nonlinear analysis tool for composite beams" by
Wenbin Yu and "Geometrically nonlinear analysis of composite beams using
Wiener-Milenković parameters" by Qi Wang and Wenbin Yu.

# Arguments
 - `jacob`: Jacobian matrix
 - `x`: Vector containing current state variables of the system
 - `assembly`: Assembly of rigidly connected nonlinear beam elements
 - `force_scaling`: Scaling parameter for forces
 - `mass_scaling`: Scaling parameter for masses
 - `irow_point`: Row index of first equilibrium equation for each point
 - `irow_elem`: Row index of first equation for just this beam element
 - `irow_elem1`: Row index of first equation for the left side of each beam
 - `irow_elem2`: Row index of first equation for the right side of each beam
 - `icol_point`: Column index of first state variable for each point
 - `icol_elem`: Column index of first state variable for each beam element
"""
function system_mass_matrix!(jacob, x, assembly, force_scaling, mass_scaling, irow_point,
    irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem)

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

        element_mass_matrix!(jacob, x, assembly.elements[ielem], force_scaling,
            mass_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2)
    end

    # no contributions to "mass matrix" from point state variables

    # # zero out near-zero values ( < eps() )
    # jacob = droptol!(jacob, eps(eltype(jacob)))

    return jacob
end

"""
    system_mass_matrix!(jacob, gamma, x, dx, assembly, force_scaling, mass_scaling,
        irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem)

Add the system mass matrix to `jacob`, scaled by the scaling parameter `gamma`.
"""
function system_mass_matrix!(jacob, gamma, x, assembly, force_scaling, mass_scaling, irow_point, irow_elem,
    irow_elem1, irow_elem2, icol_point, icol_elem)

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

        element_mass_matrix!(jacob, gamma, x, assembly.elements[ielem],
            force_scaling, mass_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2)
    end

    # no contributions to "mass matrix" from point state variables

    return jacob
end
