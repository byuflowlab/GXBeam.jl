"""
    SystemIndices

Structure for holding indices for accessing the state variables and equations associated 
with each point and beam element in a system.
"""
struct SystemIndices
    nstates::Int
    irow_point::Vector{Int}
    irow_elem::Vector{Int}
    icol_point::Vector{Int}
    icol_elem::Vector{Int}
end

"""
    SystemIndices(start, stop, case)

Define indices for accessing the state variables and equations associated with each point 
and beam element in an assembly using the connectivity of each beam element.
"""
function SystemIndices(start, stop; static=false, expanded=false)

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

        if expanded
            # add equilibrium equations
            icol += 6
            irow += 6

            if !static
                # additional states and equations for dynamic simulations
                icol += 6
                irow += 6
            end
        end

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

    nstates = icol - 1

    return SystemIndices(nstates, irow_point, irow_elem, icol_point, icol_elem)
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
 - `static_indices`: Indices for indexing into the state and residual vectors of a static system.
 - `dynamic_indices`: Indices for indexing into the state and residual vectors of a dynamic system.
 - `expanded_indices`: Indices for indexing into the state and residual vectors of an expanded system.
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
    expanded_indices::SystemIndices
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
 - `force_scaling`: Factor used to scale system forces/moments internally.  If
    not specified, a suitable default will be chosen based on the entries of the
    beam element compliance matrices.
"""
function System(assembly; kwargs...)

    return System(eltype(assembly), assembly; kwargs...)
end

function System(TF, assembly; force_scaling = default_force_scaling(assembly))

    # initialize system pointers
    static_indices = SystemIndices(assembly.start, assembly.stop, static=true, expanded=false)
    dynamic_indices = SystemIndices(assembly.start, assembly.stop, static=false, expanded=false)
    expanded_indices = SystemIndices(assembly.start, assembly.stop, static=false, expanded=true)

    # initialize system matrices
    nx = dynamic_indices.nstates
    x = zeros(TF, nx)
    r = zeros(TF, nx)
    K = spzeros(TF, nx, nx)
    M = spzeros(TF, nx, nx)

    # initialize storage for time domain simulations
    np = length(assembly.points)
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
        expanded_indices, udot, θdot, Vdot, Ωdot, t)
end

function default_force_scaling(assembly)

    TF = eltype(assembly)

    nsum = 0
    csum = zero(TF)
    for elem in assembly.elements
        for val in elem.compliance
            csum += abs(val)
            if eps(TF) < abs(val)
                nsum += 1
            end
        end
    end

    force_scaling = nextpow(2.0, nsum/csum/100)

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
    set_static_state!(system, xs)

Set the static state variables in `system` to the values in the static state vector `xs`
"""
function set_static_state!(system, xs)

    system.x .= 0

    for (is, id) in zip(system.static_indices.icol_point, system.dynamic_indices.icol_point)
        system.x[id:id+5] .= xs[is:is+5]
    end

    for (is, id) in zip(system.static_indices.icol_elem, system.dynamic_indices.icol_elem)
        system.x[id:id+5] .= xs[is:is+5]
    end

    return system
end

"""
    get_expanded_state(system, assembly, prescribed_conditions, x=system.x)

Return the state vector `x` for a constant mass matrix system
"""
function get_expanded_state(system, assembly, prescribed_conditions, x = system.x)
    
    xe = zeros(eltype(x), system.expanded_indices.nstates)

    force_scaling = 1.0

    pcond = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(0)

    for ipoint = 1:length(assembly.points)

        icol = system.expanded_indices.icol_point[ipoint]

        # linear and angular displacements
        u, θ = point_displacement(x, ipoint, system.dynamic_indices.icol_point, pcond)
        
        # external forces and moments
        F, M = point_loads(x, ipoint, system.dynamic_indices.icol_point, force_scaling, pcond)

        # linear and angular velocities
        V, Ω = point_velocities(x, ipoint, system.dynamic_indices.icol_point)

        # rotation matrices
        C = get_C(θ)

        # rotate external forces for the point into the deformed frame
        F = C*F
        M = C*M

        # rotate linear and angular velocities into the deformed frame
        V = C*V
        Ω = C*Ω

        # set load and deflection state variables
        if haskey(pcond, ipoint)
            !pcond[ipoint].isforce[1] ? setindex!(xe, F[1], icol)   : setindex!(xe, u[1], icol)  
            !pcond[ipoint].isforce[2] ? setindex!(xe, F[2], icol+1) : setindex!(xe, u[2], icol+1)
            !pcond[ipoint].isforce[3] ? setindex!(xe, F[3], icol+2) : setindex!(xe, u[3], icol+2)
            !pcond[ipoint].isforce[4] ? setindex!(xe, M[1], icol+3) : setindex!(xe, u[1], icol+3)  
            !pcond[ipoint].isforce[5] ? setindex!(xe, M[2], icol+4) : setindex!(xe, u[2], icol+4)
            !pcond[ipoint].isforce[6] ? setindex!(xe, M[3], icol+5) : setindex!(xe, u[3], icol+5)
        else
            xe[icol  ] = F[1]
            xe[icol+1] = F[2]
            xe[icol+2] = F[3]
            xe[icol+3] = M[1]
            xe[icol+4] = M[2]
            xe[icol+5] = M[3]
        end
        
        # set velocity state variables
        xe[icol+6:icol+8] .= V
        xe[icol+9:icol+11] .= Ω
    end

    for ielem = 1:length(assembly.elements)

        icol = system.expanded_indices.icol_elem[ielem]

        # unpack element parameters
        @unpack L, Cab, compliance, mass, mu = assembly.elements[ielem]

        # linear and angular displacement
        u1, θ1 = point_displacement(x, assembly.start[ielem], system.dynamic_indices.icol_point, pcond)
        u2, θ2 = point_displacement(x, assembly.stop[ielem], system.dynamic_indices.icol_point, pcond)
        u = (u1 + u2)/2
        θ = (θ1 + θ2)/2

        # rotation parameter matrices
        C = get_C(θ)
        CtCab = C'*Cab
        
        # linear and angular velocity
        V1, Ω1 = point_velocities(x, assembly.start[ielem], system.dynamic_indices.icol_point)
        V2, Ω2 = point_velocities(x, assembly.stop[ielem], system.dynamic_indices.icol_point)
        V = (V1 + V2)/2
        Ω = (Ω1 + Ω2)/2

        # forces and moments
        F, M = element_loads(x, ielem, system.dynamic_indices.icol_elem, force_scaling)

        # rotate linear and angular velocity into the deformed frame
        V = CtCab'*V
        Ω = CtCab'*Ω

        # set endpoint loads (TODO: initialize this correctly)
        xe[icol:icol+2] .= F
        xe[icol+3:icol+5] .= M
        xe[icol+6:icol+8] .= F
        xe[icol+9:icol+11] .= M

        # set new element velocities
        xe[icol+12:icol+14] .= V
        xe[icol+15:icol+17] .= Ω
    end

    return xe
end

"""
    set_expanded_state!(system, assembly, prescribed_conditions, xe)

Set the state variables in `system` to the values in the expanded state vector `xe`
"""
function set_expanded_state!(system, assembly, prescribed_conditions, xe)

    for ipoint = 1:length(assembly.points)
        # extract point state variables
        u, θ = point_displacement(xe, ipoint, system.expanded_indices.icol_point, prescribed_conditions)
        CV, CΩ = point_velocities(xe, ipoint, system.expanded_indices.icol_point)
        C = get_C(θ)

        # indices for this point
        id = system.dynamic_indices.icol_point[ipoint]
        ie = system.expanded_indices.icol_point[ipoint]

        # set point state variables
        system.x[id:id+5] .= xe[ie:ie+5]
        system.x[id+6:id+8] .= C'*CV
        system.x[id+9:id+11] .= C'*CΩ
    end

    for ielem = 1:length(assembly.elements)
        # indices for this element
        id = system.dynamic_indices.icol_elem[ielem]
        ie = system.expanded_indices.icol_elem[ielem]

        # set element loads
        system.x[id:id+2] .= (xe[ie:ie+2] .+ xe[ie+6:ie+8]) ./ 2
        system.x[id+3:id+5] .= (xe[ie+3:ie+5] .+ xe[ie+9:ie+11]) ./ 2
    end

    return system
end

"""
    set_state_variables!([x,] system, prescribed_conditions; kwargs...)

Set the state variables in `system` (or in the vector `x`) to the provided values.

# Keyword Arguments
- `u`: Vector containing the linear displacement of each point.
- `theta`: Vector containing the angular displacement of each point.
- `V`: Vector containing the linear velocity of each point.
- `Omega` Vector containing the angular velocity of each point
- `F`: Vector containing the externally applied forces acting on each point
- `M`: Vector containing the externally applied moments acting on each point
- `Fi`: Vector containing internal forces for each beam element
- `Mi`: Vector containing internal moments for each beam element
"""
set_state_variables!

function set_state_variables!(system, prescribed_conditions; kwargs...)
    x = set_state_variables!(system.x, system, prescribed_conditions; kwargs...)
    return system
end

function set_state_variables!(x, system, prescribed_conditions; u = nothing, theta = nothing, 
    V = nothing, Omega = nothing, F = nothing, M = nothing, Fi = nothing, Mi = nothing) 

    if !isnothing(u)
        for ipoint = 1:length(u)
            set_linear_deflection!(x, system, prescribed_conditions, u[ipoint], ipoint)
        end
    end

    if !isnothing(theta)
        for ipoint = 1:length(theta)
            set_angular_deflection!(x, system, prescribed_conditions, theta[ipoint], ipoint)
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

    if !isnothing(F)
        for ipoint = 1:length(F)
            set_external_forces!(x, system, F[ipoint], ipoint)
        end
    end

    if !isnothing(M)
        for ipoint = 1:length(M)
            set_external_moments!(x, system, M[ipoint], ipoint)
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

    return x
end

"""
    set_linear_deflection!([x,] system, prescribed_conditions, u, ipoint)

Set the state variables in `system` (or in the vector `x`) corresponding to the
linear deflection of point `ipoint` to the provided values.
"""
function set_linear_deflection!(system::System, prescribed_conditions, u, ipoint)
    set_linear_deflection!(system.x, system, prescribed_conditions, u, ipoint)
    return system
end

function set_linear_deflection!(x, system::System, prescribed_conditions, u, ipoint)
    
    icol = system.dynamic_indices.icol_point[ipoint]
    
    prescribed = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(t)

    if haskey(prescribed, ipoint)
        prescribed[ipoint].isforce[1] && setindex!(x, u[1], icol)
        prescribed[ipoint].isforce[2] && setindex!(x, u[2], icol+1)
        prescribed[ipoint].isforce[3] && setindex!(x, u[3], icol+2)
    else
        x[icol  ] = u[1]
        x[icol+1] = u[2]
        x[icol+2] = u[3]
    end

    return x
end

"""
    set_angular_deflection!([x,] system, prescribed_conditions, theta, ipoint)

Set the state variables in `system` (or in the vector `x`) corresponding to the
angular deflection of point `ipoint` to the provided values.
"""
function set_angular_deflection!(system::System, prescribed_conditions, theta, ipoint)
    set_angular_deflection!(system.x, system, prescribed_conditions, theta, ipoint)
    return system
end

function set_angular_deflection!(x, system::System, prescribed_conditions, theta, ipoint)

    icol = system.dynamic_indices.icol_point[ipoint]
    
    prescribed = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(t)

    if haskey(prescribed, ipoint)
        prescribed[ipoint].isforce[4] && setindex!(x, theta[1], icol+3)
        prescribed[ipoint].isforce[5] && setindex!(x, theta[2], icol+4)
        prescribed[ipoint].isforce[6] && setindex!(x, theta[3], icol+5)
    else
        x[icol+3] = theta[1]
        x[icol+4] = theta[2]
        x[icol+5] = theta[3]
    end

    return x
end

"""
    set_linear_velocity([x,] system, V, ipoint)

Set the state variables in `system` (or in the vector `x`) corresponding to the
linear velocity of point `ipoint` to the provided values.
"""
function set_linear_velocity(system::System, V, ipoint)
    set_linear_velocity(system.x, system, V, ipoint)
    return system
end

function set_linear_velocity(x, system::System, V, ipoint)

    icol = system.dynamic_indices.icol_point[ipoint]

    x[icol  ] = V[1]
    x[icol+1] = V[2]
    x[icol+2] = V[3]

    return x
end

"""
    set_angular_velocity!([x,] system, Omega, ipoint)

Set the state variables in `system` (or in the vector `x`) corresponding to the
angular velocity of point `ipoint` to the provided values.
"""
function set_angular_velocity!(system::System, Omega, ipoint)
    set_angular_velocity!(system.x, system, Omega, ipoint)
    return system
end

function set_angular_velocity!(x, system::System, Omega, ipoint)

    icol = system.dynamic_indices.icol_elem[ipoint]

    x[icol+3] = Omega[1]
    x[icol+4] = Omega[2]
    x[icol+5] = Omega[3]

    return x
end


"""
    set_external_forces!([x,] system, prescribed_conditions, F, ipoint)

Set the state variables in `system` (or in the vector `x`) corresponding to the
external forces applied at point `ipoint` to the provided values.
"""
function set_external_forces!(system::System, prescribed_conditions, F, ipoint)
    set_external_forces!(system.x, system, prescribed_conditions, F, ipoint)
    return system
end

function set_external_forces!(x, system::System, prescribed_conditions, F, ipoint)
    
    icol = system.dynamic_indices.icol_point[ipoint]
    
    prescribed = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(t)

    if haskey(prescribed, ipoint)
        !prescribed[ipoint].isforce[1] && setindex!(x, F[1], icol)
        !prescribed[ipoint].isforce[2] && setindex!(x, F[2], icol+1)
        !prescribed[ipoint].isforce[3] && setindex!(x, F[3], icol+2)
    else
        x[icol  ] = F[1]
        x[icol+1] = F[2]
        x[icol+2] = F[3]
    end

    return x
end

"""
    set_external_moments([x,] system, prescribed_conditions, M, ipoint)

Set the state variables in `system` (or in the vector `x`) corresponding to the
external moments applied at point `ipoint` to the provided values.
"""
function set_external_moments(system::System, prescribed_conditions, M, ipoint)
    set_external_moments(system.x, system, prescribed_conditions, M, ipoint)
    return system
end

function set_external_moments(x, system::System, prescribed_conditions, M, ipoint)

    icol = system.dynamic_indices.icol_point[ipoint]
    
    prescribed = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(t)

    if haskey(prescribed, ipoint)
        !prescribed[ipoint].isforce[4] && setindex!(x, M[1], icol+3)
        !prescribed[ipoint].isforce[5] && setindex!(x, M[2], icol+4)
        !prescribed[ipoint].isforce[6] && setindex!(x, M[3], icol+5)
    else
        x[icol+3] = M[1]
        x[icol+4] = M[2]
        x[icol+5] = M[3]
    end

    return x
end

"""
    set_internal_forces!([x,] system, Fi, ielem)

Set the state variables in `system` (or in the vector `x`) corresponding to the
internal forces of element `ielem` to the provided values.
"""
function set_internal_forces!(system::System, Fi, ielem)
    set_internal_forces!(system.x, system, Fi, ielem)
    return system
end

function set_internal_forces!(x, system::System, Fi, ielem)

    icol = system.dynamic_indices.icol_elem[ielem]

    force_scaling = system.force_scaling

    x[icol  ] = Fi[1] / force_scaling
    x[icol+1] = Fi[2] / force_scaling
    x[icol+2] = Fi[3] / force_scaling

    return x
end

"""
    set_internal_moments!([x,] system, Mi, ielem)

Set the state variables in `system` (or in the vector `x`) corresponding to the
internal moments of element `ielem` to the provided values.
"""
function set_internal_moments!(system::System, Mi, ielem)
    set_internal_moments!(system.x, system, Mi, ielem)
    return system
end

function set_internal_moments!(x, system::System, Mi, ielem)

    icol = system.dynamic_indices.icol_elem[ielem]

    force_scaling = system.force_scaling

    x[icol+3] = Mi[1] / force_scaling
    x[icol+4] = Mi[2] / force_scaling
    x[icol+5] = Mi[3] / force_scaling

    return x
end

"""
    static_system_residual!(resid, x, indices, force_scaling, 
        assembly, prescribed_conditions, distributed_loads, point_masses, gravity)

Populate the system residual vector `resid` for a static analysis
"""
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

"""
    steady_state_system_residual!(resid, x, indices, force_scaling, structural_damping, 
        assembly, prescribed_conditions, distributed_loads, point_masses, gravity, 
        x0, v0, ω0, a0, α0)

Populate the system residual vector `resid` for a steady state analysis
"""
@inline function steady_state_system_residual!(resid, x, indices, force_scaling, structural_damping,
    assembly, prescribed_conditions, distributed_loads, point_masses, gravity, 
    x0, v0, ω0, a0, α0)

    for ipoint = 1:length(assembly.points)
        steady_state_point_residual!(resid, x, indices, force_scaling, assembly, ipoint, 
            prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)
    end

    for ielem = 1:length(assembly.elements)
        steady_state_element_residual!(resid, x, indices, force_scaling, structural_damping, 
            assembly, ielem, prescribed_conditions, distributed_loads, gravity, x0, v0, ω0, a0, α0)
    end

    return resid
end

"""
    initial_condition_system_residual!(resid, x, indices, force_scaling, structural_damping, 
        assembly, prescribed_conditions, distributed_loads, point_masses, gravity, 
        x0, v0, ω0, a0, α0, u0, θ0, udot0, θdot0)

Populate the system residual vector `resid` for the initialization of a time domain 
simulation.
"""
@inline function initial_condition_system_residual!(resid, x, indices, force_scaling, structural_damping, 
    assembly, prescribed_conditions, distributed_loads, point_masses, gravity, x0, v0, ω0, a0, α0,
    u0, θ0, udot0, θdot0)

    for ipoint = 1:length(assembly.points)
        initial_condition_point_residual!(resid, x, indices, force_scaling, assembly, ipoint, 
            prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0, u0, θ0, udot0, θdot0)
    end
    
    for ielem = 1:length(assembly.elements)
        initial_condition_element_residual!(resid, x, indices, force_scaling, structural_damping, 
            assembly, ielem, prescribed_conditions, distributed_loads, gravity, 
            x0, v0, ω0, a0, α0, u0, θ0, udot0, θdot0)
    end
    
    # NOTE: If a point and the elements connected to it are massless, then Vdot and Ωdot for
    # the point do not appear in the system of equations and thus cannot be found during the 
    # time domain initialization.  Therefore when a point and the elements connected to it 
    # are massless we replace the force equilibrium equations for the point with the linear 
    # and angular velocity constraints ``\dot{V}=0`` and ``\dot{\Omega}=0``.
    
    for ipoint = 1:length(assembly.points)
        # check if Vdot and Ωdot for this point are used 
        hasmass = haskey(point_masses, ipoint)
        for ielem = 1:length(assembly.elements)
            if ipoint == assembly.start[ielem] || ipoint == assembly.stop[ielem]
                hasmass = hasmass || (!iszero(assembly.elements[ielem].L) && !iszero(assembly.elements[ielem].mass))
            end
        end

        # replace equilibrium constraints with linear and angular velocity constraints
        if !hasmass
            irow = indices.irow_point[ipoint]

            Vdot, Ωdot = point_displacement_rates(x, ipoint, indices.icol_point, prescribed_conditions)

            if haskey(prescribed_conditions, ipoint)
                prescribed_conditions[ipoint].isforce[1] && setindex!(resid, Vdot[1], irow) 
                prescribed_conditions[ipoint].isforce[2] && setindex!(resid, Vdot[2], irow+1) 
                prescribed_conditions[ipoint].isforce[3] && setindex!(resid, Vdot[3], irow+2)  
                prescribed_conditions[ipoint].isforce[4] && setindex!(resid, Ωdot[1], irow+3)  
                prescribed_conditions[ipoint].isforce[5] && setindex!(resid, Ωdot[2], irow+4)  
                prescribed_conditions[ipoint].isforce[6] && setindex!(resid, Ωdot[3], irow+5)   
            else
                resid[irow] = Vdot[1]
                resid[irow+1] = Vdot[2]
                resid[irow+2] = Vdot[3]
                resid[irow+3] = Ωdot[1]
                resid[irow+4] = Ωdot[2]
                resid[irow+5] = Ωdot[3] 
            end
        end

    end

    return resid
end

"""
    newmark_system_residual!(resid, x, indices, force_scaling, structural_damping, 
        assembly, prescribed_conditions, distributed_loads, point_masses, gravity, 
        x0, v0, ω0, a0, α0, udot_init, θdot_init, Vdot_init, Ωdot_init, dt)

Populate the system residual vector `resid` for a Newmark scheme time marching analysis.
"""
@inline function newmark_system_residual!(resid, x, indices, force_scaling, structural_damping, 
    assembly, prescribed_conditions, distributed_loads, point_masses, gravity, x0, v0, ω0, a0, α0,
    udot_init, θdot_init, Vdot_init, Ωdot_init, dt)
       
    for ipoint = 1:length(assembly.points)
        newmark_point_residual!(resid, x, indices, force_scaling, assembly, ipoint, 
            prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0, 
            udot_init, θdot_init, Vdot_init, Ωdot_init, dt)
    end
    
    for ielem = 1:length(assembly.elements)
        newmark_element_residual!(resid, x, indices, force_scaling, structural_damping, 
            assembly, ielem, prescribed_conditions, distributed_loads, gravity, 
            x0, v0, ω0, a0, α0, Vdot_init, Ωdot_init, dt)
    end
    
    return resid
end

"""
    dynamic_system_residual!(resid, dx, x, indices, force_scaling, structural_damping, 
        assembly, prescribed_conditions, distributed_loads, point_masses, gravity, 
        x0, v0, ω0, a0, α0)

Populate the system residual vector `resid` for a general dynamic analysis.
"""
@inline function dynamic_system_residual!(resid, dx, x, indices, force_scaling, structural_damping, 
    assembly, prescribed_conditions, distributed_loads, point_masses, gravity, x0, v0, ω0, a0, α0)

    for ipoint = 1:length(assembly.points)
        dynamic_point_residual!(resid, dx, x, indices, force_scaling, assembly, ipoint, 
            prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)
    end
    
    for ielem = 1:length(assembly.elements)
        dynamic_element_residual!(resid, dx, x, indices, force_scaling, structural_damping, 
            assembly, ielem, prescribed_conditions, distributed_loads, gravity, 
            x0, v0, ω0, a0, α0)
    end
    
    return resid
end

"""
    expanded_system_residual!(resid, x, indices, force_scaling, structural_damping, 
        assembly, prescribed_conditions, distributed_loads, point_masses, gravity, 
        x0, v0, ω0, a0, α0)

Populate the system residual vector `resid` for a constant mass matrix system.
"""
@inline function expanded_system_residual!(resid, x, indices, force_scaling, structural_damping, 
    assembly, prescribed_conditions, distributed_loads, point_masses, gravity, x0, v0, ω0, a0, α0)

    for ipoint = 1:length(assembly.points)
        expanded_point_residual!(resid, x, indices, force_scaling, assembly, ipoint, 
            prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)
    end
    
    for ielem = 1:length(assembly.elements)
        expanded_element_residual!(resid, x, indices, force_scaling, structural_damping, 
            assembly, ielem, prescribed_conditions, distributed_loads, gravity, 
            x0, v0, ω0, a0, α0)
    end
    
    return resid
end

"""
    static_system_jacobian!(jacob, x, indices, force_scaling, 
        assembly, prescribed_conditions, distributed_loads, point_masses, gravity)

Populate the system jacobian matrix `jacob` for a static analysis
"""
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

"""
    steady_state_system_jacobian!(jacob, x, indices, force_scaling, structural_damping,
        assembly, prescribed_conditions, distributed_loads, point_masses, gravity, 
        x0, v0, ω0, a0, α0)

Populate the system jacobian matrix `jacob` for a steady-state analysis
"""
@inline function steady_state_system_jacobian!(jacob, x, indices, force_scaling, structural_damping,
    assembly, prescribed_conditions, distributed_loads, point_masses, gravity, 
    x0, v0, ω0, a0, α0)

    jacob .= 0
    
    for ipoint = 1:length(assembly.points)
        steady_state_point_jacobian!(jacob, x, indices, force_scaling, assembly, ipoint, 
            prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)
    end
    
    for ielem = 1:length(assembly.elements)
        steady_state_element_jacobian!(jacob, x, indices, force_scaling, structural_damping, assembly, ielem, 
            prescribed_conditions, distributed_loads, gravity, x0, v0, ω0, a0, α0)
    end
    
    return jacob
end

"""
    initial_condition_system_jacobian!(jacob, x, indices, force_scaling, structural_damping, 
        assembly, prescribed_conditions, distributed_loads, point_masses, gravity, x0, v0, ω0, a0, α0,
        u0, θ0, udot0, θdot0)

Populate the system jacobian matrix `jacob` for the initialization of a time domain 
simulation.
"""
@inline function initial_condition_system_jacobian!(jacob, x, indices, force_scaling, structural_damping, 
    assembly, prescribed_conditions, distributed_loads, point_masses, gravity, x0, v0, ω0, a0, α0,
    u0, θ0, udot0, θdot0)
    
    jacob .= 0
    
    for ipoint = 1:length(assembly.points)
        initial_condition_point_jacobian!(jacob, x, indices, force_scaling, assembly, ipoint, 
            prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0, u0, θ0, udot0, θdot0)
    end
    
    for ielem = 1:length(assembly.elements)
        initial_condition_element_jacobian!(jacob, x, indices, force_scaling, structural_damping, 
            assembly, ielem, prescribed_conditions, distributed_loads, gravity, 
            x0, v0, ω0, a0, α0, u0, θ0, udot0, θdot0)
    end

    # NOTE: If a point and the elements connected to it are massless, then Vdot and Ωdot for
    # the point do not appear in the system of equations and thus cannot be found during the 
    # time domain initialization.  Therefore when a point and the elements connected to it 
    # are massless we replace the force equilibrium equations for the point with the linear 
    # and angular velocity constraints ``\dot{V}=0`` and ``\dot{\Omega}=0``.
    
    for ipoint = 1:length(assembly.points)
        # check if Vdot and Ωdot for this point are used 
        hasmass = haskey(point_masses, ipoint)
        for ielem = 1:length(assembly.elements)
            if ipoint == assembly.start[ielem] || ipoint == assembly.stop[ielem]
                hasmass = hasmass || (!iszero(assembly.elements[ielem].L) && !iszero(assembly.elements[ielem].mass))
            end
        end

        # replace equilibrium constraints with linear and angular velocity constraints
        if !hasmass
            irow = indices.irow_point[ipoint]
            icol = indices.icol_point[ipoint]

            if haskey(prescribed_conditions, ipoint)
                if prescribed_conditions[ipoint].isforce[1]
                    jacob[irow, :] .= 0
                    jacob[irow, icol] = 1
                end 
                if prescribed_conditions[ipoint].isforce[2] 
                    jacob[irow+1, :] .= 0
                    jacob[irow+1, icol+1] = 1
                end 
                if prescribed_conditions[ipoint].isforce[3] 
                    jacob[irow+2, :] .= 0
                    jacob[irow+2, icol+2] = 1
                end 
                if prescribed_conditions[ipoint].isforce[4] 
                    jacob[irow+3, :] .= 0
                    jacob[irow+3, icol+3] = 1
                end 
                if prescribed_conditions[ipoint].isforce[5] 
                    jacob[irow+4, :] .= 0
                    jacob[irow+4, icol+4] = 1
                end 
                if prescribed_conditions[ipoint].isforce[6] 
                    jacob[irow+5, :] .= 0
                    jacob[irow+5, icol+5] = 1
                end 
            else
                jacob[irow, :] .= 0
                jacob[irow, icol] = 1

                jacob[irow+1, :] .= 0
                jacob[irow+1, icol+1] = 1

                jacob[irow+2, :] .= 0
                jacob[irow+2, icol+2] = 1

                jacob[irow+3, :] .= 0
                jacob[irow+3, icol+3] = 1

                jacob[irow+4, :] .= 0
                jacob[irow+4, icol+4] = 1

                jacob[irow+5, :] .= 0
                jacob[irow+5, icol+5] = 1
            end
        end

    end
    
    return jacob
end

"""
    newmark_system_jacobian!(jacob, x, indices, force_scaling, structural_damping, 
        assembly, prescribed_conditions, distributed_loads, point_masses, gravity, 
        x0, v0, ω0, a0, α0, udot_init, θdot_init, Vdot_init, Ωdot_init, dt)

Populate the system jacobian matrix `jacob` for a Newmark scheme time marching analysis.
"""
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
        newmark_element_jacobian!(jacob, x, indices, force_scaling, structural_damping, 
            assembly, ielem, prescribed_conditions, distributed_loads, gravity, 
            x0, v0, ω0, a0, α0, Vdot_init, Ωdot_init, dt)
    end
    
    return jacob
end

"""
    dynamic_system_jacobian!(jacob, dx, x, indices, force_scaling, structural_damping, 
        assembly, prescribed_conditions, distributed_loads, point_masses, gravity, 
        x0, v0, ω0, a0, α0)

Populate the system jacobian matrix `jacob` for a general dynamic analysis.
"""
@inline function dynamic_system_jacobian!(jacob, dx, x, indices, force_scaling, structural_damping, 
    assembly, prescribed_conditions, distributed_loads, point_masses, gravity, x0, v0, ω0, a0, α0)
    
    jacob .= 0
    
    for ipoint = 1:length(assembly.points)
        dynamic_point_jacobian!(jacob, dx, x, indices, force_scaling, assembly, ipoint, 
            prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)
    end
    
    for ielem = 1:length(assembly.elements)
        dynamic_element_jacobian!(jacob, dx, x, indices, force_scaling, structural_damping, 
            assembly, ielem, prescribed_conditions, distributed_loads, gravity, 
            x0, v0, ω0, a0, α0)
    end
    
    return jacob
end

"""
    expanded_system_jacobian!(jacob, x, indices, force_scaling, structural_damping, 
        assembly, prescribed_conditions, distributed_loads, point_masses, gravity, 
        x0, v0, ω0, a0, α0)

Populate the system jacobian matrix `jacob` for a general expanded analysis.
"""
@inline function expanded_system_jacobian!(jacob, x, indices, force_scaling, structural_damping, 
    assembly, prescribed_conditions, distributed_loads, point_masses, gravity, x0, v0, ω0, a0, α0)
    
    jacob .= 0
    
    for ipoint = 1:length(assembly.points)
        expanded_point_jacobian!(jacob, x, indices, force_scaling, assembly, ipoint, 
            prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)
    end
    
    for ielem = 1:length(assembly.elements)
        expanded_element_jacobian!(jacob, x, indices, force_scaling, structural_damping, 
            assembly, ielem, prescribed_conditions, distributed_loads, gravity, 
            x0, v0, ω0, a0, α0)
    end
    
    return jacob
end

"""
    system_mass_matrix!(jacob, x, indices, force_scaling,  assembly, prescribed_conditions, 
        point_masses)

Calculate the jacobian of the residual expressions with respect to the state rates.
"""
function system_mass_matrix!(jacob, x, indices, force_scaling, assembly, 
    prescribed_conditions, point_masses)

    jacob .= 0

    gamma = 1

    system_mass_matrix!(jacob, gamma, x, indices, force_scaling,  assembly, 
        prescribed_conditions, point_masses)

    return jacob
end

"""
    system_mass_matrix!(jacob, gamma, x, indices, force_scaling, assembly, 
        prescribed_conditions, point_masses)

Calculate the jacobian of the residual expressions with respect to the state rates and 
add the result multiplied by `gamma` to `jacob`.
"""
function system_mass_matrix!(jacob, gamma, x, indices, force_scaling, assembly, 
    prescribed_conditions, point_masses)

    for ipoint = 1:length(assembly.points)
        mass_matrix_point_jacobian!(jacob, gamma, x, indices, force_scaling, assembly, ipoint, 
            prescribed_conditions, point_masses)
    end
    
    for ielem = 1:length(assembly.elements)
        mass_matrix_element_jacobian!(jacob, gamma, x, indices, force_scaling, assembly, ielem, 
            prescribed_conditions)
    end

    return jacob
end

"""
    expanded_system_mass_matrix(system, assembly, prescribed_conditions, point_masses)

Calculate the jacobian of the residual expressions with respect to the state rates.
"""
function expanded_system_mass_matrix(system, assembly, prescribed_conditions, point_masses)

    @unpack expanded_indices, force_scaling = system

    TF = eltype(system)
    nx = expanded_indices.nstates
    jacob = spzeros(TF, nx, nx)
    gamma = -1
    pcond = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(0)
    pmass = typeof(point_masses) <: AbstractDict ? point_masses : point_masses(0)

    expanded_system_mass_matrix!(jacob, gamma, indices, force_scaling, assembly, pcond, pmass) 

    return jacob
end

"""
    expanded_system_mass_matrix!(jacob, indices, force_scaling,  assembly, prescribed_conditions, 
        point_masses)

Calculate the jacobian of the residual expressions with respect to the state rates.
"""
function expanded_system_mass_matrix!(jacob, indices, force_scaling, assembly, 
    prescribed_conditions, point_masses)

    jacob .= 0

    gamma = 1

    expanded_system_mass_matrix!(jacob, gamma, indices, force_scaling, assembly, 
        prescribed_conditions, point_masses)

    return jacob
end

"""
    expanded_system_mass_matrix!(jacob, gamma, indices, force_scaling, assembly, 
        prescribed_conditions, point_masses)

Calculate the jacobian of the residual expressions with respect to the state rates and 
add the result multiplied by `gamma` to `jacob`.
"""
function expanded_system_mass_matrix!(jacob, gamma, indices, force_scaling, assembly, 
    prescribed_conditions, point_masses)

    for ipoint = 1:length(assembly.points)
        expanded_mass_matrix_point_jacobian!(jacob, gamma, indices, force_scaling, assembly, 
            ipoint, prescribed_conditions, point_masses)
    end
    
    for ielem = 1:length(assembly.elements)
        expanded_mass_matrix_element_jacobian!(jacob, gamma, indices, force_scaling, assembly, 
            ielem, prescribed_conditions)
    end

    return jacob
end