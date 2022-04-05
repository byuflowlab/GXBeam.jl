"""
    SystemIndices

Structure for holding indices for accessing the state variables and equations associated 
with each point and beam element in a system.
"""
struct SystemIndices
    irow_point::Vector{Int}
    irow_elem::Vector{Int}
    icol_point::Vector{Int}
    icol_elem::Vector{Int}
end

"""
    SystemIndices(start, stop)

Define indices for accessing the state variables and equations associated with each point 
and beam element in an assembly using the connectivity of each beam element.
"""
function SystemIndices(start, stop, static)

    # number of points
    np = max(maximum(start), maximum(stop))
    
    # number of elements
    ne = length(start)

    # keep track of whether state variables have been assigned to each point
    assigned = fill(false, np)

    # initialize pointers
    irow_point = Vector{Int}(undef, np)
    irow_elem = Vector{Int}(undef, ne)
    icol_point = Vector{Int}(undef, np)
    icol_elem = Vector{Int}(undef, ne)

    # define pointers for state variables and equations
    irow = 1
    icol = 1
    for ielem = 1:ne

        # add state variables and equations for the start of the beam element
        ipt = start[ielem]
        if !assigned[ipt]

            assigned[ipt] = true

            # add point state variables
            icol_point[ipt] = icol
            icol += 6

            # add equilibrium equations
            irow_point[ipt] = irow
            irow += 6

            if !static
                # additional states and equations for dynamic simulations
                icol += 6
                irow += 6
            end

        end

        # add beam state variables
        icol_elem[ielem] = icol
        icol += 6

        # add compatability equations
        irow_elem[ielem] = irow
        irow += 6

        # add state variables and equations for the end of the beam element
        ipt = stop[ielem]

        if !assigned[ipt]

            assigned[ipt] = true

            # add point state variables
            icol_point[ipt] = icol
            icol += 6

            # add equilibrium equations
            irow_point[ipt] = irow
            irow += 6

            if !static
                # additional states and equations for dynamic simulations
                icol += 6
                irow += 6
            end

        end
    end

    return SystemIndices(irow_point, irow_elem, icol_point, icol_elem)
end

"""
    System{TF, TV<:AbstractVector{TF}, TM<:AbstractMatrix{TF}}

Contains the system state, residual vector, and jacobian matrices as well as
pointers to be able to access their contents.  Also contains additional storage
needed for time domain simulations.

# Fields:
 - `x`: State vector
 - `r`: Residual vector
 - `K`: System jacobian matrix with respect to the state variables
 - `M`: System jacobian matrix with respect to the time derivative of the state variables
 - `force_scaling`: Scaling for state variables corresponding to forces/moments
 - `static_indices`: Indices for indexing into the state and residual vectors.
 - `dynamic_indices`: Indices for indexing into the state and residual vectors.
 - `udot`: Time derivative of state variable `u` for each point
 - `θdot`: Time derivative of state variable `θ` for each point
 - `Vdot`: Time derivative of state variable `V` for each point
 - `Ωdot`: Time derivative of state variable `Ω` for each point
 - `t`: Current system time
"""
mutable struct System{TF, TV<:AbstractVector{TF}, TM<:AbstractMatrix{TF}}
    x::TV
    r::TV
    K::TM
    M::TM
    force_scaling::TF
    static_indices::SystemIndices
    dynamic_indices::SystemIndices
    udot::Vector{SVector{3,TF}}
    θdot::Vector{SVector{3,TF}}
    Vdot::Vector{SVector{3,TF}}
    Ωdot::Vector{SVector{3,TF}}
    t::TF
end
Base.eltype(::System{TF, TV, TM}) where {TF, TV, TM} = TF

"""
    System([TF=eltype(assembly),] assembly; kwargs...)

Initialize an object of type `System` which stores the system state.

# Arguments:
 - `TF:`(optional) Floating point type, defaults to the floating point type of `assembly`
 - `assembly`: Assembly of rigidly connected nonlinear beam elements

# Keyword Arguments
 - `prescribed_points`: Point indices at which prescribed conditions may be applied.  By 
    default, prescribed conditions may be applied at any point.  Using a limited number of 
    points reduces computational expenses.
 - `force_scaling`: Factor used to scale system forces/moments internally.  If
    not specified, a suitable default will be chosen based on the entries of the
    compliance matrix.
"""
function System(assembly; kwargs...)

    return System(eltype(assembly), assembly; kwargs...)
end

function System(TF, assembly; force_scaling = default_force_scaling(assembly))


    # system dimensions
    np = length(assembly.points)
    ne = length(assembly.elements)
    nx = 12*np + 6*ne

    # initialize system pointers
    static_indices = SystemIndices(assembly.start, assembly.stop, true)
    dynamic_indices = SystemIndices(assembly.start, assembly.stop, false)

    # initialize system matrices
    x = zeros(TF, nx)
    r = zeros(TF, nx)
    K = spzeros(TF, nx, nx)
    M = spzeros(TF, nx, nx)

    # initialize storage for time domain simulations
    udot = [@SVector zeros(TF, 3) for i = 1:np]
    θdot = [@SVector zeros(TF, 3) for i = 1:np]
    Vdot = [@SVector zeros(TF, 3) for i = 1:np]
    Ωdot = [@SVector zeros(TF, 3) for i = 1:np]

    # initialize current time
    t = 0.0

    # set system types
    TV = promote_type(typeof(x), typeof(r))
    TM = promote_type(typeof(K), typeof(M))

    return System{TF, TV, TM}(x, r, K, M, force_scaling, static_indices, dynamic_indices, 
        udot, θdot, Vdot, Ωdot, t)
end

function default_force_scaling(assembly)

    TF = eltype(assembly)

    nsum = 0
    csum = 1.0
    for elem in assembly.elements
        for val in elem.compliance
            csum += val
            if eps(TF) < abs(val)
                nsum += 1
            end
        end
    end

    force_scaling = nextpow(2.0, nsum/csum)

    return force_scaling
end

"""
    system_state(system)

Return a vector containing the state variables of `system`.
"""
system_state(system) = system.x

"""
    reset_state!(system)

Reset the state variables in `system` (stored in `system.x`).
"""
function reset_state!(system)
    system.x .= 0
    return system
end

"""
    get_static_state(system, x=system.x)

Return the state vector `x` for a static system
"""
function get_static_state(system, x=system.x)

    np = length(system.static_indices.icol_point)
    ne = length(system.static_indices.icol_elem)

    xs = zeros(eltype(x), 6*np+6*ne)

    for (is, id) in zip(system.static_indices.icol_point, system.dynamic_indices.icol_point)
        xs[is:is+5] .= x[id:id+5]
    end

    for (is, id) in zip(system.static_indices.icol_elem, system.dynamic_indices.icol_elem)
        xs[is:is+5] .= x[id:id+5]
    end

    return xs
end

"""
    set_static_state!(system, x)

Set the static state variables in `system` to the values in the state vector `x`
"""
function set_static_state!(system, x)

    system.x .= 0

    for (is, id) in zip(system.static_indices.icol_point, system.dynamic_indices.icol_point)
        system.x[id:id+5] .= x[is:is+5]
    end

    for (is, id) in zip(system.static_indices.icol_elem, system.dynamic_indices.icol_elem)
        system.x[id:id+5] .= x[is:is+5]
    end

    return system
end

@inline function static_system_residual!(resid, x, indices, force_scaling, 
    assembly, prescribed_conditions, distributed_loads, point_masses, gravity)

    for ipoint = 1:length(assembly.points)
        static_point_residual!(resid, x, indices, force_scaling, assembly, ipoint, 
            prescribed_conditions, point_masses, gravity)
    end

    for ielem = 1:length(assembly.elements)
        static_element_residual!(resid, x, indices, force_scaling, assembly, ielem, 
            prescribed_conditions, distributed_loads, gravity)
    end

    return resid
end

@inline function steady_state_system_residual!(resid, x, indices, force_scaling, 
    assembly, prescribed_conditions, distributed_loads, point_masses, gravity, 
    x0, v0, ω0, a0, α0)

    for ipoint = 1:length(assembly.points)
        steady_state_point_residual!(resid, x, indices, force_scaling, assembly, ipoint, 
            prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)
    end

    for ielem = 1:length(assembly.elements)
        steady_state_element_residual!(resid, x, indices, force_scaling, assembly, ielem, 
            prescribed_conditions, distributed_loads, gravity, x0, v0, ω0, a0, α0)
    end

    return resid
end

@inline function initial_condition_system_residual!(resid, x, indices, force_scaling, structural_damping, 
    assembly, prescribed_conditions, distributed_loads, point_masses, gravity, x0, v0, ω0, a0, α0,
    u0, θ0, udot0, θdot0)
       
    for ipoint = 1:length(assembly.points)
        initial_condition_point_residual!(resid, x, indices, force_scaling, assembly, ipoint, 
            prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0, u0, θ0, udot0, θdot0)
    end
    
    for ielem = 1:length(assembly.elements)
        initial_condition_element_residual!(resid, x, indices, force_scaling, assembly, ielem,
            prescribed_conditions, distributed_loads, gravity, x0, v0, ω0, a0, α0, u0, θ0, udot0, θdot0)
    end
    
    return resid
end

@inline function newmark_system_residual!(resid, x, indices, force_scaling, structural_damping, 
    assembly, prescribed_conditions, distributed_loads, point_masses, gravity, x0, v0, ω0, a0, α0,
    udot_init, θdot_init, Vdot_init, Ωdot_init, dt)
       
    for ipoint = 1:length(assembly.points)
        newmark_point_residual!(resid, x, indices, force_scaling, assembly, ipoint, 
            prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0, 
            udot_init, θdot_init, Vdot_init, Ωdot_init, dt)
    end
    
    for ielem = 1:length(assembly.elements)
        newmark_element_residual!(resid, x, indices, force_scaling, assembly, ielem,
            prescribed_conditions, distributed_loads, gravity, x0, v0, ω0, a0, α0, 
            udot_init, θdot_init, Vdot_init, Ωdot_init, dt)
    end
    
    return resid
end

@inline function dynamic_system_residual!(resid, dx, x, indices, force_scaling, structural_damping, 
    assembly, prescribed_conditions, distributed_loads, point_masses, gravity, x0, v0, ω0, a0, α0)
       
    for ipoint = 1:length(assembly.points)
        dynamic_point_residual!(resid, dx, x, indices, force_scaling, assembly, ipoint, 
            prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)
    end
    
    for ielem = 1:length(assembly.elements)
        dynamic_element_residual!(resid, dx, x, indices, force_scaling, assembly, ielem,
            prescribed_conditions, distributed_loads, gravity, x0, v0, ω0, a0, α0)
    end
    
    return resid
end

@inline function static_system_jacobian!(jacob, x, indices, force_scaling, 
    assembly, prescribed_conditions, distributed_loads, point_masses, gravity)

    jacob .= 0
    
    for ipoint = 1:length(assembly.points)
        static_point_jacobian!(jacob, x, indices, force_scaling, assembly, ipoint, 
            prescribed_conditions, point_masses, gravity)
    end
    
    for ielem = 1:length(assembly.elements)
        static_element_jacobian!(jacob, x, indices, force_scaling, assembly, ielem, 
            prescribed_conditions, distributed_loads, gravity)
    end
    
    return jacob
end

@inline function steady_state_system_jacobian!(jacob, x, indices, force_scaling, 
    assembly, prescribed_conditions, distributed_loads, point_masses, gravity, 
    x0, v0, ω0, a0, α0)

    jacob .= 0
    
    for ipoint = 1:length(assembly.points)
        steady_state_point_jacobian!(jacob, x, indices, force_scaling, assembly, ipoint, 
            prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)
    end
    
    for ielem = 1:length(assembly.elements)
        steady_state_element_jacobian!(jacob, x, indices, force_scaling, assembly, ielem, 
            prescribed_conditions, distributed_loads, gravity, x0, v0, ω0, a0, α0)
    end
    
    return jacob
end

@inline function initial_condition_system_jacobian!(jacob, x, indices, force_scaling, structural_damping, 
    assembly, prescribed_conditions, distributed_loads, point_masses, gravity, x0, v0, ω0, a0, α0,
    u0, θ0, udot0, θdot0)
    
    jacob .= 0
    
    for ipoint = 1:length(assembly.points)
        initial_condition_point_jacobian!(jacob, x, indices, force_scaling, assembly, ipoint, 
            prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0, u0, θ0, udot0, θdot0)
    end
    
    for ielem = 1:length(assembly.elements)
        initial_condition_element_jacobian!(jacob, x, indices, force_scaling, assembly, ielem,
            prescribed_conditions, distributed_loads, gravity, x0, v0, ω0, a0, α0, u0, θ0, udot0, θdot0)
    end
    
    return jacob
end

@inline function newmark_system_jacobian!(jacob, x, indices, force_scaling, structural_damping, 
    assembly, prescribed_conditions, distributed_loads, point_masses, gravity, x0, v0, ω0, a0, α0,
    udot_init, θdot_init, Vdot_init, Ωdot_init, dt)
    
    jacob .= 0
    
    for ipoint = 1:length(assembly.points)
        newmark_point_jacobian!(jacob, x, indices, force_scaling, assembly, ipoint, 
            prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0, 
            udot_init, θdot_init, Vdot_init, Ωdot_init, dt)
    end
    
    for ielem = 1:length(assembly.elements)
        newmark_element_jacobian!(jacob, x, indices, force_scaling, assembly, ielem,
            prescribed_conditions, distributed_loads, gravity, x0, v0, ω0, a0, α0, 
            udot_init, θdot_init, Vdot_init, Ωdot_init, dt)
    end
    
    return jacob
end

@inline function dynamic_system_jacobian!(jacob, dx, x, indices, force_scaling, structural_damping, 
    assembly, prescribed_conditions, distributed_loads, point_masses, gravity, x0, v0, ω0, a0, α0)
    
    jacob .= 0
    
    for ipoint = 1:length(assembly.points)
        dynamic_point_jacobian!(jacob, dx, x, indices, force_scaling, assembly, ipoint, 
            prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)
    end
    
    for ielem = 1:length(assembly.elements)
        dynamic_element_jacobian!(jacob, dx, x, indices, force_scaling, assembly, ielem,
            prescribed_conditions, distributed_loads, gravity, x0, v0, ω0, a0, α0)
    end
    
    return jacob
end

function system_mass_matrix!(jacob, x, indices, force_scaling, structural_damping, 
    assembly, prescribed_conditions, point_masses)

    jacob .= 0

    gamma = 1

    system_mass_matrix!(jacob, gamma, x, indices, force_scaling, structural_damping, 
        assembly, prescribed_conditions, point_masses)

    return jacob
end

function system_mass_matrix!(jacob, gamma, x, indices, force_scaling, structural_damping, 
    assembly, prescribed_conditions, point_masses)

    for ipoint = 1:length(assembly.points)
        point_mass_matrix!(jacob, gamma, x, indices, force_scaling, assembly, ipoint, 
            prescribed_conditions, point_masses)
    end
    
    for ielem = 1:length(assembly.elements)
        element_mass_matrix!(jacob, gamma, x, indices, force_scaling, assembly, ielem, 
            prescribed_conditions)
    end

    return jacob
end