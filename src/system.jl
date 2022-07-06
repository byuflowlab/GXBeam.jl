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

    # add other states and equations
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

    # total number of state variables
    nstates = icol - 1

    return SystemIndices(nstates, irow_point, irow_elem, icol_point, icol_elem)
end

"""
    body_frame_acceleration_indices(system, prescribed_conditions)

Return the state vector indices corresponding to the components of the body-frame linear
and angular acceleration.  Return `0` for a component's index when accelerations are 
prescribed.
"""
function body_frame_acceleration_indices(system, prescribed_conditions)

    # Use only one point's state variables for linear and angular acceleration states
    @assert count(p -> p.second.pl[1] & p.second.pd[1], prescribed_conditions) <= 1 &&
        count(p -> p.second.pl[2] & p.second.pd[2], prescribed_conditions) <= 1 &&
        count(p -> p.second.pl[3] & p.second.pd[3], prescribed_conditions) <= 1 &&
        count(p -> p.second.pl[4] & p.second.pd[4], prescribed_conditions) <= 1 &&
        count(p -> p.second.pl[5] & p.second.pd[5], prescribed_conditions) <= 1 &&
        count(p -> p.second.pl[6] & p.second.pd[6], prescribed_conditions) <= 1 "Forces "*
        "and displacements corresponding to the same degree of freedom cannot be "
        "prescribed at the same time at more than one point"

    # point index for linear and angular acceleration state variables
    ipoint_accel = @SVector [findfirst(p -> p.pl[i] & p.pd[i], prescribed_conditions) for i = 1:6]
    ipoint_accel = replace(ipoint_accel, nothing => 0)
    
    # column index for linear and angular acceleration state variables
    icol = SVector(
        iszero(ipoint_accel[1]) ? 0 : system.indices.icol_point[ipoint_accel[1]],
        iszero(ipoint_accel[2]) ? 0 : system.indices.icol_point[ipoint_accel[2]]+1,
        iszero(ipoint_accel[3]) ? 0 : system.indices.icol_point[ipoint_accel[3]]+2,
        iszero(ipoint_accel[4]) ? 0 : system.indices.icol_point[ipoint_accel[4]]+3,
        iszero(ipoint_accel[5]) ? 0 : system.indices.icol_point[ipoint_accel[5]]+4,
        iszero(ipoint_accel[6]) ? 0 : system.indices.icol_point[ipoint_accel[6]]+5,
    )

    return icol
end

"""
    body_frame_acceleration(x, icol, a0, α0)

Extract the linear and angular acceleration of the body frame from the state variable vector
or prescribed linear and angular acceleration.
"""
function body_frame_acceleration(x, icol, a0, α0)
    
    # adjust body-frame linear and angular acceleration to account for state variables
    a0 = SVector(
        iszero(icol[1]) ? a0[1] : x[icol[1]],
        iszero(icol[2]) ? a0[2] : x[icol[2]],
        iszero(icol[3]) ? a0[3] : x[icol[3]],
    )

    α0 = SVector(
        iszero(icol[4]) ? α0[1] : x[icol[4]],
        iszero(icol[5]) ? α0[2] : x[icol[5]],
        iszero(icol[6]) ? α0[3] : x[icol[6]],
    )
 
    return a0, α0
end

"""
    AbstractSystem

Supertype for types which contain the system state, residual vector, and jacobian matrix. 
"""
abstract type AbstractSystem end

"""
    StaticSystem{TF, TV<:AbstractVector{TF}, TM<:AbstractMatrix{TF}} <: AbstractSystem

Contains the system state, residual vector, and jacobian matrix for a static system. 
"""
mutable struct StaticSystem{TF, TV<:AbstractVector{TF}, TM<:AbstractMatrix{TF}} <: AbstractSystem
    x::TV
    r::TV
    K::TM
    indices::SystemIndices
    force_scaling::TF
    t::TF
end
Base.eltype(::StaticSystem{TF, TV, TM}) where {TF, TV, TM} = TF

"""
    StaticSystem([TF=eltype(assembly),] assembly; kwargs...)

Initialize an object of type [`StaticSystem`](@ref).

# Arguments:
 - `TF:`(optional) Floating point type, defaults to the floating point type of `assembly`
 - `assembly`: Assembly of rigidly connected nonlinear beam elements

# Keyword Arguments
 - `force_scaling`: Factor used to scale system forces/moments internally.  If
    not specified, a suitable default will be chosen based on the entries of the
    beam element compliance matrices.
"""
function StaticSystem(assembly; kwargs...)

    return StaticSystem(eltype(assembly), assembly; kwargs...)
end

function StaticSystem(TF, assembly; force_scaling = default_force_scaling(assembly))

    # initialize system pointers
    indices = SystemIndices(assembly.start, assembly.stop, static=true, expanded=false)

    # initialize system states
    x = zeros(TF, indices.nstates)
    r = zeros(TF, indices.nstates)
    K = spzeros(TF, indices.nstates, indices.nstates)

    # initialize current time
    t = 0.0

    x, r = promote(x, r)

    return StaticSystem{TF, Vector{TF}, SparseMatrixCSC{TF, Int64}}(x, r, K, indices, force_scaling, t)
end

"""
    DynamicSystem{TF, TV<:AbstractVector{TF}, TM<:AbstractMatrix{TF}} <: AbstractSystem

Contains the system state, residual vector, and jacobian matrix for a dynamic system. 
"""
mutable struct DynamicSystem{TF, TV<:AbstractVector{TF}, TM<:AbstractMatrix{TF}} <: AbstractSystem
    x::TV
    r::TV
    K::TM
    M::TM
    indices::SystemIndices
    force_scaling::TF
    udot::Vector{SVector{3,TF}}
    θdot::Vector{SVector{3,TF}}
    Vdot::Vector{SVector{3,TF}}
    Ωdot::Vector{SVector{3,TF}}
    t::TF
end
Base.eltype(::DynamicSystem{TF, TV, TM}) where {TF, TV, TM} = TF

"""
    DynamicSystem([TF=eltype(assembly),] assembly; kwargs...)

Initialize an object of type [`DynamicSystem`](@ref).

# Arguments:
 - `TF:`(optional) Floating point type, defaults to the floating point type of `assembly`
 - `assembly`: Assembly of rigidly connected nonlinear beam elements

# Keyword Arguments
 - `force_scaling`: Factor used to scale system forces/moments internally.  If
    not specified, a suitable default will be chosen based on the entries of the
    beam element compliance matrices.
"""
function DynamicSystem(assembly; kwargs...)

    return DynamicSystem(eltype(assembly), assembly; kwargs...)
end

function DynamicSystem(TF, assembly; force_scaling = default_force_scaling(assembly))

    # initialize system pointers
    indices = SystemIndices(assembly.start, assembly.stop, static=false, expanded=false)

    # initialize system states
    x = zeros(TF, indices.nstates)
    r = zeros(TF, indices.nstates)
    K = spzeros(TF, indices.nstates, indices.nstates)
    M = spzeros(TF, indices.nstates, indices.nstates)

    # initialize storage for a Newmark-Scheme time marching analysis
    udot = [@SVector zeros(TF, 3) for point in assembly.points]
    θdot = [@SVector zeros(TF, 3) for point in assembly.points]
    Vdot = [@SVector zeros(TF, 3) for point in assembly.points]
    Ωdot = [@SVector zeros(TF, 3) for point in assembly.points]

    # initialize current time
    t = 0.0

    return DynamicSystem{TF, Vector{TF}, SparseMatrixCSC{TF, Int64}}(x, r, K, M, indices, 
        force_scaling, udot, θdot, Vdot, Ωdot, t)
end

"""
    ExpandedSystem{TF, TV<:AbstractVector{TF}, TM<:AbstractMatrix{TF}} <: AbstractSystem

Contains the system state, residual vector, and jacobian matrix for a constant mass matrix 
system. 
"""
mutable struct ExpandedSystem{TF, TV<:AbstractVector{TF}, TM<:AbstractMatrix{TF}} <: AbstractSystem
    x::TV
    r::TV
    K::TM
    M::TM
    indices::SystemIndices
    force_scaling::TF
    t::TF
end
Base.eltype(::ExpandedSystem{TF, TV, TM}) where {TF, TV, TM} = TF

"""
    ExpandedSystem([TF=eltype(assembly),] assembly; kwargs...)

Initialize an object of type [`ExpandedSystem`](@ref).

# Arguments:
 - `TF:`(optional) Floating point type, defaults to the floating point type of `assembly`
 - `assembly`: Assembly of rigidly connected nonlinear beam elements

# Keyword Arguments
 - `force_scaling`: Factor used to scale system forces/moments internally.  If
    not specified, a suitable default will be chosen based on the entries of the
    beam element compliance matrices.
"""
function ExpandedSystem(assembly; kwargs...)

    return ExpandedSystem(eltype(assembly), assembly; kwargs...)
end

function ExpandedSystem(TF, assembly; force_scaling = default_force_scaling(assembly))

    # initialize system pointers
    indices = SystemIndices(assembly.start, assembly.stop, static=false, expanded=true)

    # initialize system states
    x = zeros(TF, indices.nstates)
    r = zeros(TF, indices.nstates)
    K = spzeros(TF, indices.nstates, indices.nstates)
    M = spzeros(TF, indices.nstates, indices.nstates)

    # initialize current time
    t = 0.0

    return ExpandedSystem{TF, Vector{TF}, SparseMatrixCSC{TF, Int64}}(x, r, K, M, indices, 
        force_scaling, t)
end

# default system is a DynamicSystem
const System = DynamicSystem

"""
    default_force_scaling(assembly)

Defines a suitable default force scaling factor based on the nonzero elements of the 
compliance matrices in `assembly`.
"""
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

    force_scaling = iszero(nsum) ? 1.0 : nextpow(2.0, nsum/csum/100)

    return force_scaling
end

"""
    reset_state!(system)

Sets the state variables in `system` to zero.
"""
function reset_state!(system)
    system.x .= 0
    return system
end

"""
    copy_state!(system1, system2, assembly; kwargs...)

Copy the state variables from `system2` into `system1`

# General Keyword Arguments
 - `prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}()`:
        A dictionary with keys corresponding to the points at
        which prescribed conditions are applied and values of type
        [`PrescribedConditions`](@ref) which describe the prescribed conditions
        at those points.  If time varying, this input may be provided as a
        function of time.
 - `distributed_loads = Dict{Int,DistributedLoads{Float64}}()`: A dictionary
        with keys corresponding to the elements to which distributed loads are
        applied and values of type [`DistributedLoads`](@ref) which describe
        the distributed loads on those elements.  If time varying, this input may
        be provided as a function of time.
 - `point_masses = Dict{Int,PointMass{Float64}}()`: A dictionary with keys 
        corresponding to the points to which point masses are attached and values 
        of type [`PointMass`](@ref) which contain the properties of the attached 
        point masses.  If time varying, this input may be provided as a function of time.
 - `origin = zeros(3)`: Body frame origin vector. If time varying, this input
        may be provided as a function of time.
 - `linear_velocity = zeros(3)`: Body frame linear velocity vector. If time
        varying, this vector may be provided as a function of time.
 - `angular_velocity = zeros(3)`: Body frame angular velocity vector. If time
        varying, this vector may be provided as a function of time.
 - `linear_acceleration = zeros(3)`: Body frame linear acceleration vector. If time
        varying, this vector may be provided as a function of time.
 - `angular_acceleration = zeros(3)`: Body frame angular acceleration vector. If time
        varying, this vector may be provided as a function of time.
 - `gravity = [0,0,0]`: Gravity vector.  If time varying, this input may be provided as a 
       function of time. 
 - `t`: Current time or time vector. Defaults to the current time stored in `system2`  
 
 # Control Flag Keyword Arguments
 - `reset_state = true`: Flag indicating whether the system state variables should be 
        set to zero prior to copying over the new state variables.
"""
function copy_state!(system1, system2, assembly; 
    # general keyword arguments
    prescribed_conditions=Dict{Int,PrescribedConditions{Float64}}(),
    distributed_loads=Dict{Int,DistributedLoads{Float64}}(),
    point_masses=Dict{Int,PointMass{Float64}}(),
    origin=SVector(0,0,0),
    linear_velocity=SVector(nothing, nothing, nothing),
    angular_velocity=SVector(nothing, nothing, nothing),
    linear_acceleration=SVector(nothing, nothing, nothing),
    angular_acceleration=SVector(nothing, nothing, nothing),
    gravity=SVector(0,0,0),
    time=system2.t,
    # control flag keyword arguments
    structural_damping=true,
    reset_state=true,
    )

    # reset state, if specified
    if reset_state
        reset_state!(system1)
    end

    # get current time
    t = first(time)

    # update stored time
    system1.t = t

    # current parameters
    pcond = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(t)
    dload = typeof(distributed_loads) <: AbstractDict ? distributed_loads : distributed_loads(t)
    gvec = typeof(gravity) <: AbstractVector ? SVector{3}(gravity) : SVector{3}(gravity(t))
    x0 = typeof(origin) <: AbstractVector ? SVector{3}(origin) : SVector{3}(origin(t))
    v0 = typeof(linear_velocity) <: AbstractVector ? SVector{3}(linear_velocity) : SVector{3}(linear_velocity(t))
    ω0 = typeof(angular_velocity) <: AbstractVector ? SVector{3}(angular_velocity) : SVector{3}(angular_velocity(t))
    a0 = typeof(linear_acceleration) <: AbstractVector ? SVector{3}(linear_acceleration) : SVector{3}(linear_acceleration(t))
    α0 = typeof(angular_acceleration) <: AbstractVector ? SVector{3}(angular_acceleration) : SVector{3}(angular_acceleration(t))

    # extract state vectors for each system
    x1 = system1.x
    x2 = system2.x

    # extract point state variable pointers for each system
    icol_point1 = system1.indices.icol_point
    icol_point2 = system2.indices.icol_point

    # extract element state variable pointers for each system
    icol_elem1 =  system1.indices.icol_elem
    icol_elem2 =  system2.indices.icol_elem

    # extract force scaling parameter for each system
    force_scaling1 = system1.force_scaling
    force_scaling2 = system2.force_scaling

    # copy over rigid body accelerations
    icol_accel = body_frame_acceleration_indices(system2, pcond)
    !iszero(icol_accel[1]) && setindex!(x1, x2[icol_accel[1]], icol_accel[1])
    !iszero(icol_accel[2]) && setindex!(x1, x2[icol_accel[2]], icol_accel[2])
    !iszero(icol_accel[3]) && setindex!(x1, x2[icol_accel[3]], icol_accel[3])
    !iszero(icol_accel[4]) && setindex!(x1, x2[icol_accel[4]], icol_accel[4])
    !iszero(icol_accel[5]) && setindex!(x1, x2[icol_accel[5]], icol_accel[5])
    !iszero(icol_accel[6]) && setindex!(x1, x2[icol_accel[6]], icol_accel[6])

    # copy over state variables for each point
    for ipoint = 1:length(assembly.points)
        
        # linear and angular displacement
        u, θ = point_displacement(x2, ipoint, icol_point2, pcond)

        # external forces and moments
        F, M = point_loads(x2, ipoint, icol_point2, force_scaling2, pcond)

        # transformation matrix to the deformed frame
        C = get_C(θ)

        # linear and angular velocity
        if typeof(system2) <: StaticSystem
            V = Ω = @SVector zeros(3)
        elseif typeof(system2) <: DynamicSystem
            V, Ω = point_velocities(x2, ipoint, icol_point2)
        elseif typeof(system2) <: ExpandedSystem
            CV, CΩ = point_velocities(x2, ipoint, icol_point2)
            V = C'*CV
            Ω = C'*CΩ
        end

        # copy over new state variables
        icol = icol_point1[ipoint]

        # displacement and external load state variables
        if haskey(pcond, ipoint)
            # linear and angular displacement
            !pcond[ipoint].pd[1] && setindex!(x1, u[1], icol)
            !pcond[ipoint].pd[2] && setindex!(x1, u[2], icol+1)
            !pcond[ipoint].pd[3] && setindex!(x1, u[3], icol+2)
            !pcond[ipoint].pd[4] && setindex!(x1, θ[1], icol+3)
            !pcond[ipoint].pd[5] && setindex!(x1, θ[2], icol+4)
            !pcond[ipoint].pd[6] && setindex!(x1, θ[3], icol+5)
            # external forces and moments
            !pcond[ipoint].pl[1] && setindex!(x1, F[1] / force_scaling1, icol)
            !pcond[ipoint].pl[2] && setindex!(x1, F[2] / force_scaling1, icol+1)
            !pcond[ipoint].pl[3] && setindex!(x1, F[3] / force_scaling1, icol+2)
            !pcond[ipoint].pl[4] && setindex!(x1, M[1] / force_scaling1, icol+3)
            !pcond[ipoint].pl[5] && setindex!(x1, M[2] / force_scaling1, icol+4)
            !pcond[ipoint].pl[6] && setindex!(x1, M[3] / force_scaling1, icol+5)
        else
            # linear and angular displacement
            x1[icol] = u[1]
            x1[icol+1] = u[2]
            x1[icol+2] = u[3]
            x1[icol+3] = θ[1]
            x1[icol+4] = θ[2]
            x1[icol+5] = θ[3]
        end
        
        # linear and angular velocity state variables
        if typeof(system1) <: DynamicSystem
            x1[icol+6:icol+8] .= V
            x1[icol+9:icol+11] .= Ω
        elseif typeof(system1) <: ExpandedSystem
            x1[icol+6:icol+8] .= C*V
            x1[icol+9:icol+11] .= C*Ω
        end
    end

    # copy over state variables for each element
    for ielem = 1:length(assembly.elements)
        
        # resultant forces and moments
        if typeof(system1) <: ExpandedSystem

            if typeof(system2) <: StaticSystem

                # compute static element properties
                properties = static_element_properties(x2, system2.indices, force_scaling2, 
                    assembly, ielem, pcond, gvec)

                # compute static element resultants
                resultants = static_element_resultants(properties, dload, ielem)

                # unpack element resultants
                @unpack F1, M1, F2, M2 = resultants

                # rotate element resultants into the deformed element frame
                CtCab = properties.CtCab
                
                F1 = CtCab'*F1
                F2 = CtCab'*F2
                M1 = CtCab'*M1
                M2 = CtCab'*M2

            elseif typeof(system2) <: DynamicSystem

                # compute steady state element properties
                properties = steady_state_element_properties(x2, system2.indices, icol_accel, 
                    force_scaling2, structural_damping, assembly, ielem, pcond, gvec, 
                    x0, v0, ω0, a0, α0)

                @unpack C, Cab, CtCab, mass11, mass12, mass21, mass22, V1, V2, Ω1, Ω2, V, Ω, ω = properties

                # linear and angular velocity rates
                V1dot = system2.Vdot[assembly.start[ielem]]
                Ω1dot = system2.Ωdot[assembly.start[ielem]]
            
                V2dot = system2.Vdot[assembly.stop[ielem]]
                Ω2dot = system2.Ωdot[assembly.stop[ielem]]
            
                Vdot = (V1dot + V2dot)/2
                Ωdot = (Ω1dot + Ω2dot)/2

                # linear and angular momentum rates
                CtCabdot = C'*tilde(Ω - ω)*Cab
                
                Pdot = CtCab*mass11*CtCab'*Vdot + CtCab*mass12*CtCab'*Ωdot +
                    CtCab*mass11*CtCabdot'*V + CtCab*mass12*CtCabdot'*Ω +
                    CtCabdot*mass11*CtCab'*V + CtCabdot*mass12*CtCab'*Ω
                
                Hdot = CtCab*mass21*CtCab'*Vdot + CtCab*mass22*CtCab'*Ωdot +
                    CtCab*mass21*CtCabdot'*V + CtCab*mass22*CtCabdot'*Ω +
                    CtCabdot*mass21*CtCab'*V + CtCabdot*mass22*CtCab'*Ω

                # combine steady state and dynamic element properties
                properties = (; properties..., CtCabdot, Vdot, Ωdot, Pdot, Hdot) 

                # compute dynamic element resultants
                resultants = dynamic_element_resultants(properties, dload, ielem)

                # unpack element resultants
                @unpack F1, M1, F2, M2 = resultants

                # rotate element resultants into the deformed element frame
                CtCab = properties.CtCab
                F1 = CtCab'*F1
                F2 = CtCab'*F2
                M1 = CtCab'*M1
                M2 = CtCab'*M2

            elseif typeof(system2) <: ExpandedSystem
                F1, M1, F2, M2 = expanded_element_loads(x2, ielem, icol_elem2, force_scaling2)
            end

        else

            # internal forces and moments
            if typeof(system2) <: StaticSystem || typeof(system2) <: DynamicSystem
                F, M = element_loads(x2, ielem, icol_elem2, force_scaling2)
            elseif typeof(system2) <: ExpandedSystem
                F1, M1, F2, M2 = expanded_element_loads(x2, ielem, icol_elem2, force_scaling2)
                F = (F1 + F2)/2
                M = (M1 + M2)/2
            end

        end

        # copy over new state variables
        icol = icol_elem1[ielem]

        # set internal load state variables
        if typeof(system1) <: StaticSystem || typeof(system1) <: DynamicSystem

            # copy over internal forces and moments
            x1[icol] = F[1] / force_scaling1
            x1[icol+1] = F[2] / force_scaling1
            x1[icol+2] = F[3] / force_scaling1
            x1[icol+3] = M[1] / force_scaling1
            x1[icol+4] = M[2] / force_scaling1
            x1[icol+5] = M[3] / force_scaling1

        elseif typeof(system1) <: ExpandedSystem

            x1[icol] = F1[1] / force_scaling1
            x1[icol+1] = F1[2] / force_scaling1
            x1[icol+2] = F1[3] / force_scaling1
            x1[icol+3] = M1[1] / force_scaling1
            x1[icol+4] = M1[2] / force_scaling1
            x1[icol+5] = M1[3] / force_scaling1
            x1[icol+6] = F2[1] / force_scaling1
            x1[icol+7] = F2[2] / force_scaling1
            x1[icol+8] = F2[3] / force_scaling1
            x1[icol+9] = M2[1] / force_scaling1
            x1[icol+10] = M2[2] / force_scaling1
            x1[icol+11] = M2[3] / force_scaling1

        end

        # set element velocities
        if typeof(system1) <: ExpandedSystem

            # linear and angular velocity
            if typeof(system2) <: StaticSystem
                V = Ω = @SVector zeros(3)
            elseif typeof(system2) <: DynamicSystem
                @unpack CtCab, V, Ω, = properties
                V = CtCab'*V
                Ω = CtCab'*Ω
            elseif typeof(system2) <: ExpandedSystem
                V, Ω = expanded_element_velocities(x2, ielem, icol_point2)
            end

            x1[icol+12] = V[1]
            x1[icol+13] = V[2]
            x1[icol+14] = V[3]
            x1[icol+15] = Ω[1]
            x1[icol+16] = Ω[2]
            x1[icol+17] = Ω[3]

        end

    end

    return system1
end

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
- `linear_acceleration`: Linear acceleration of the body frame
- `angular_acceleration`: Angular acceleration of the body frame
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

    if !isnothing(linear_acceleration)
        set_linear_acceleration!(x, system, prescribed_conditions, linear_acceleration)
    end

    if !isnothing(angular_acceleration)
        set_angular_acceleration!(x, system, prescribed_conditions, angular_acceleration)
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
    
    prescribed = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(t)

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
    
    prescribed = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(t)

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
    
    prescribed = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(t)

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

    prescribed = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(t)

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
    set_linear_acceleration!([x,] system, prescribed_conditions, linear_acceleration)

Set the state variables in `system` (or in the vector `x`) corresponding to the
linear acceleration of the body frame to the provided values.
"""
function set_linear_acceleration!(system, prescribed_conditions, linear_acceleration)
    set_linear_acceleration!(system.x, system, prescribed_conditions, linear_acceleration)
    return system
end

function set_linear_acceleration!(x, system, prescribed_conditions, linear_acceleration)
        
    prescribed = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(t)

    icol = body_frame_acceleration_indices(system, prescribed)

    !iszero(icol[1]) && setindex!(x, linear_acceleration[1], icol[1])
    !iszero(icol[2]) && setindex!(x, linear_acceleration[2], icol[2])
    !iszero(icol[3]) && setindex!(x, linear_acceleration[3], icol[3])

    return x
end

"""
    set_angular_acceleration!([x,] system, prescribed_conditions, angular_acceleration)

Set the state variables in `system` (or in the vector `x`) corresponding to the
angular deflection of point `ipoint` to the provided values.
"""
function set_angular_acceleration!(system, prescribed_conditions, angular_acceleration)
    set_angular_acceleration!(system.x, system, prescribed_conditions, angular_acceleration)
    return system
end

function set_angular_acceleration!(x, system, prescribed_conditions, angular_acceleration)

    prescribed = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(t)

    icol = body_frame_acceleration_indices(system, prescribed)

    !iszero(icol[4]) && setindex!(x, angular_acceleration[1], icol[4])
    !iszero(icol[5]) && setindex!(x, angular_acceleration[2], icol[5])
    !iszero(icol[6]) && setindex!(x, angular_acceleration[3], icol[6])

    return x
end

"""
    static_system_residual!(resid, x, indices, force_scaling, 
        assembly, prescribed_conditions, distributed_loads, point_masses, gravity)

Populate the system residual vector `resid` for a static analysis
"""
function static_system_residual!(resid, x, indices, force_scaling, 
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
    steady_state_system_residual!(resid, x, indices, icol_accel, force_scaling, 
        structural_damping, assembly, prescribed_conditions, distributed_loads, 
        point_masses, gravity, x0, v0, ω0, a0, α0)

Populate the system residual vector `resid` for a steady state analysis
"""
function steady_state_system_residual!(resid, x, indices, icol_accel, force_scaling, 
    structural_damping, assembly, prescribed_conditions, distributed_loads, point_masses, 
    gravity, x0, v0, ω0, a0, α0)

    # contributions to the residual vector from points
    for ipoint = 1:length(assembly.points)
        steady_state_point_residual!(resid, x, indices, icol_accel, force_scaling, 
            assembly, ipoint, prescribed_conditions, point_masses, gravity, 
            x0, v0, ω0, a0, α0)
    end

    # contributions to the residual vector from elements
    for ielem = 1:length(assembly.elements)
        steady_state_element_residual!(resid, x, indices, icol_accel, force_scaling, 
            structural_damping, assembly, ielem, prescribed_conditions, 
            distributed_loads, gravity, x0, v0, ω0, a0, α0)
    end

    return resid
end

"""
    initial_condition_system_residual!(resid, x, indices, rate_vars, icol_accel, 
        force_scaling, structural_damping, assembly, prescribed_conditions, 
        distributed_loads, point_masses, gravity, x0, v0, ω0, a0, α0, 
        u0, theta0, V0, Omega0, Vdot0, Omegadot0)

Populate the system residual vector `resid` for the initialization of a time domain 
simulation.
"""
function initial_condition_system_residual!(resid, x, indices, rate_vars, icol_accel, 
    force_scaling, structural_damping, assembly, prescribed_conditions, distributed_loads, 
    point_masses, gravity, x0, v0, ω0, a0, α0, u0, θ0, V0, Ω0, Vdot0, Ωdot0)

    # contributions to the residual vector from points
    for ipoint = 1:length(assembly.points)
        initial_condition_point_residual!(resid, x, indices, rate_vars, icol_accel, 
            force_scaling, assembly, ipoint, prescribed_conditions, point_masses, gravity, 
            x0, v0, ω0, a0, α0, u0, θ0, V0, Ω0, Vdot0, Ωdot0)
    end
    
    # contributions to the residual vector from elements
    for ielem = 1:length(assembly.elements)
        initial_condition_element_residual!(resid, x, indices, rate_vars, icol_accel, 
            force_scaling, structural_damping, assembly, ielem, prescribed_conditions, 
            distributed_loads, gravity, x0, v0, ω0, a0, α0, u0, θ0, V0, Ω0, Vdot0, Ωdot0)
    end

    return resid
end

"""
    newmark_system_residual!(resid, x, indices, icol_accel, force_scaling, 
        structural_damping, assembly, prescribed_conditions, distributed_loads, 
        point_masses, gravity, x0, v0, ω0, a0, α0, udot_init, θdot_init, 
        Vdot_init, Ωdot_init, dt)

Populate the system residual vector `resid` for a Newmark scheme time marching analysis.
"""
function newmark_system_residual!(resid, x, indices, icol_accel, force_scaling, 
    structural_damping, assembly, prescribed_conditions, distributed_loads, point_masses, 
    gravity, x0, v0, ω0, a0, α0, udot_init, θdot_init, Vdot_init, Ωdot_init, dt)
    
    # contributions to the residual vector from points
    for ipoint = 1:length(assembly.points)
        newmark_point_residual!(resid, x, indices, icol_accel, force_scaling, assembly, 
            ipoint, prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0, 
            udot_init, θdot_init, Vdot_init, Ωdot_init, dt)
    end
    
    # contributions to the residual vector from elements
    for ielem = 1:length(assembly.elements)
        newmark_element_residual!(resid, x, indices, icol_accel, force_scaling, 
            structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, 
            gravity, x0, v0, ω0, a0, α0, Vdot_init, Ωdot_init, dt)
    end
    
    return resid
end

"""
    dynamic_system_residual!(resid, dx, x, indices, icol_accel, force_scaling, 
        structural_damping, assembly, prescribed_conditions, distributed_loads, 
        point_masses, gravity, x0, v0, ω0, a0, α0)

Populate the system residual vector `resid` for a general dynamic analysis.
"""
function dynamic_system_residual!(resid, dx, x, indices, icol_accel, force_scaling, 
    structural_damping, assembly, prescribed_conditions, distributed_loads, point_masses, 
    gravity, x0, v0, ω0, a0, α0)

    # contributions to the residual vector from points
    for ipoint = 1:length(assembly.points)
        dynamic_point_residual!(resid, dx, x, indices, icol_accel, force_scaling, assembly, 
            ipoint, prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)
    end
    
    # contributions to the residual vector from elements
    for ielem = 1:length(assembly.elements)
        dynamic_element_residual!(resid, dx, x, indices, icol_accel, force_scaling, 
            structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, 
            gravity, x0, v0, ω0, a0, α0)
    end
    
    return resid
end

"""
    expanded_system_residual!(resid, x, indices, icol_accel, force_scaling, 
        structural_damping, assembly, prescribed_conditions, distributed_loads, 
        point_masses, gravity, x0, v0, ω0, a0, α0)

Populate the system residual vector `resid` for a constant mass matrix system.
"""
function expanded_system_residual!(resid, x, indices, icol_accel, force_scaling, 
    structural_damping, assembly, prescribed_conditions, distributed_loads, point_masses, 
    gravity, x0, v0, ω0, a0, α0)

    # point residuals
    for ipoint = 1:length(assembly.points)
        expanded_point_residual!(resid, x, indices, icol_accel, force_scaling, assembly, 
            ipoint, prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)
    end
    
    # element residuals
    for ielem = 1:length(assembly.elements)
        expanded_element_residual!(resid, x, indices, icol_accel, force_scaling, 
            structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, 
            gravity, x0, v0, ω0, a0, α0)
    end
    
    return resid
end

"""
    static_system_jacobian!(jacob, x, indices, force_scaling, 
        assembly, prescribed_conditions, distributed_loads, point_masses, gravity)

Populate the system jacobian matrix `jacob` for a static analysis
"""
function static_system_jacobian!(jacob, x, indices, force_scaling, 
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
    steady_state_system_jacobian!(jacob, x, indices, icol_accel, force_scaling, 
        structural_damping, assembly, prescribed_conditions, distributed_loads, 
        point_masses, gravity, x0, v0, ω0, a0, α0)

Populate the system jacobian matrix `jacob` for a steady-state analysis
"""
function steady_state_system_jacobian!(jacob, x, indices, icol_accel, force_scaling, 
    structural_damping, assembly, prescribed_conditions, distributed_loads, 
    point_masses, gravity, x0, v0, ω0, a0, α0)

    jacob .= 0
    
    for ipoint = 1:length(assembly.points)
        steady_state_point_jacobian!(jacob, x, indices, icol_accel, force_scaling, assembly, 
            ipoint, prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)
    end
    
    for ielem = 1:length(assembly.elements)
        steady_state_element_jacobian!(jacob, x, indices, icol_accel, force_scaling, 
            structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, 
            gravity, x0, v0, ω0, a0, α0)
    end
    
    return jacob
end

"""
    initial_condition_system_jacobian!(jacob, x, indices, rate_vars, icol_accel,
        force_scaling, structural_damping, assembly, prescribed_conditions, 
        distributed_loads, point_masses, gravity, x0, v0, ω0, a0, α0, u0, θ0, V0, Ω0, 
        Vdot0, Ωdot0)

Populate the system jacobian matrix `jacob` for the initialization of a time domain 
simulation.
"""
function initial_condition_system_jacobian!(jacob, x, indices, rate_vars, icol_accel, 
    force_scaling, structural_damping, assembly, prescribed_conditions, distributed_loads, 
    point_masses, gravity, x0, v0, ω0, a0, α0, u0, θ0, V0, Ω0, Vdot0, Ωdot0)
    
    jacob .= 0
    
    for ipoint = 1:length(assembly.points)
        initial_condition_point_jacobian!(jacob, x, indices, rate_vars, icol_accel, force_scaling, 
            assembly, ipoint, prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0, 
            u0, θ0, V0, Ω0, Vdot0, Ωdot0)
    end
    
    for ielem = 1:length(assembly.elements)
        initial_condition_element_jacobian!(jacob, x, indices, rate_vars, icol_accel, 
            force_scaling, structural_damping, assembly, ielem, prescribed_conditions, 
            distributed_loads, gravity, x0, v0, ω0, a0, α0, u0, θ0, V0, Ω0, Vdot0, Ωdot0)
    end

    return jacob
end

"""
    newmark_system_jacobian!(jacob, x, indices, icol_accel, force_scaling, 
        structural_damping, assembly, prescribed_conditions, distributed_loads, 
        point_masses, gravity, x0, v0, ω0, a0, α0, udot_init, θdot_init, Vdot_init, 
        Ωdot_init, dt)

Populate the system jacobian matrix `jacob` for a Newmark scheme time marching analysis.
"""
function newmark_system_jacobian!(jacob, x, indices, icol_accel, force_scaling, structural_damping, 
    assembly, prescribed_conditions, distributed_loads, point_masses, gravity, x0, v0, ω0, a0, α0,
    udot_init, θdot_init, Vdot_init, Ωdot_init, dt)
    
    jacob .= 0
    
    for ipoint = 1:length(assembly.points)
        newmark_point_jacobian!(jacob, x, indices, icol_accel, force_scaling, assembly, ipoint, 
            prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0, 
            udot_init, θdot_init, Vdot_init, Ωdot_init, dt)
    end
    
    for ielem = 1:length(assembly.elements)
        newmark_element_jacobian!(jacob, x, indices, icol_accel, force_scaling, structural_damping, 
            assembly, ielem, prescribed_conditions, distributed_loads, gravity, 
            x0, v0, ω0, a0, α0, Vdot_init, Ωdot_init, dt)
    end
    
    return jacob
end

"""
    dynamic_system_jacobian!(jacob, dx, x, indices, icol_accel, force_scaling, 
        structural_damping, assembly, prescribed_conditions, distributed_loads, 
        point_masses, gravity, x0, v0, ω0, a0, α0)

Populate the system jacobian matrix `jacob` for a general dynamic analysis.
"""
function dynamic_system_jacobian!(jacob, dx, x, indices, icol_accel, force_scaling, 
    structural_damping, assembly, prescribed_conditions, distributed_loads, point_masses, 
    gravity, x0, v0, ω0, a0, α0)
    
    jacob .= 0
    
    for ipoint = 1:length(assembly.points)
        dynamic_point_jacobian!(jacob, dx, x, indices, icol_accel, force_scaling, assembly, 
            ipoint, prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)
    end
    
    for ielem = 1:length(assembly.elements)
        dynamic_element_jacobian!(jacob, dx, x, indices, icol_accel, force_scaling, 
            structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, 
            gravity, x0, v0, ω0, a0, α0)
    end
    
    return jacob
end

"""
    expanded_system_jacobian!(jacob, x, indices, force_scaling, structural_damping, 
        assembly, prescribed_conditions, distributed_loads, point_masses, gravity, 
        x0, v0, ω0, a0, α0)

Populate the system jacobian matrix `jacob` for a general expanded analysis.
"""
function expanded_system_jacobian!(jacob, x, indices, icol_accel, force_scaling, 
    structural_damping, assembly, prescribed_conditions, distributed_loads, point_masses, 
    gravity, x0, v0, ω0, a0, α0)
    
    jacob .= 0
    
    for ipoint = 1:length(assembly.points)
        expanded_point_jacobian!(jacob, x, indices, icol_accel, force_scaling, assembly, 
            ipoint, prescribed_conditions, point_masses, gravity, x0, v0, ω0, a0, α0)
    end
    
    for ielem = 1:length(assembly.elements)
        expanded_element_jacobian!(jacob, x, indices, icol_accel, force_scaling, 
            structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, 
            gravity, x0, v0, ω0, a0, α0)
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
    expanded_system_mass_matrix(system, assembly;
        prescribed_conditions=Dict{Int, PrescribedConditions}(), 
        point_masses=Dict{Int, PointMass}())

Calculate the jacobian of the residual expressions with respect to the state rates for a 
constant mass matrix system.
"""
function expanded_system_mass_matrix(system, assembly;
    prescribed_conditions=Dict{Int, PrescribedConditions}(), 
    point_masses=Dict{Int, PointMass}())

    @unpack expanded_indices, force_scaling = system

    TF = eltype(system)
    nx = expanded_indices.nstates
    jacob = spzeros(TF, nx, nx)
    gamma = -1
    pcond = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(0)
    pmass = typeof(point_masses) <: AbstractDict ? point_masses : point_masses(0)

    expanded_system_mass_matrix!(jacob, gamma, expanded_indices, force_scaling, assembly, pcond, pmass) 

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

