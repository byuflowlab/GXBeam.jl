"""
    AbstractSystem

Supertype for types which contain the system state, residual vector, and jacobian matrix. 
"""
abstract type AbstractSystem end

"""
    SystemIndices

Structure for holding indices for accessing the state variables and equations associated 
with each point and beam element in a system.
"""
struct SystemIndices
    nstates::Int
    irow_point::Vector{Int}
    irow_elem::Vector{Int}
    icol_body::Vector{Int}
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
    icol_body = Vector{Int}(undef, 6)
    icol_point = Vector{Int}(undef, np)
    icol_elem = Vector{Int}(undef, ne)

    # default to no body acceleration states
    icol_body .= 0

    # define pointers for state variables and equations
    irow = 1
    icol = 1

    # add other states and equations
    for ielem = 1:ne

        # add state variables and equations for the start of the beam element
        ipt = start[ielem]
        if !assigned[ipt]

            assigned[ipt] = true

            # add point state variables: u/F, θ/M
            icol_point[ipt] = icol
            icol += 6

            # add equilibrium equations: ∑F=0, ∑M=0
            irow_point[ipt] = irow
            irow += 6

            if !static
                # add velocity states and residuals: V, Ω
                icol += 6
                irow += 6
            end

        end

        # add element state variables
        irow_elem[ielem] = irow
        icol_elem[ielem] = icol
        if expanded
            # states: F1, F2, M1, M2, V, Ω  
            # residuals: compatability (x1), velocity (x1), equilibrium (x1)  
            irow += 18
            icol += 18
        else
            # states: F, M
            # residuals: compatability (x1)
            irow += 6
            icol += 6
        end

        # add state variables and equations for the end of the beam element
        ipt = stop[ielem]

        if !assigned[ipt]

            assigned[ipt] = true

            # add point state variables: u/F, θ/M
            icol_point[ipt] = icol
            icol += 6

            # add equilibrium equations: ∑F=0, ∑M=0
            irow_point[ipt] = irow
            irow += 6

            if !static
                # add velocity states and residuals: V, Ω
                icol += 6
                irow += 6
            end

        end
    end

    # total number of state variables
    nstates = icol - 1

    return SystemIndices(nstates, irow_point, irow_elem, icol_body, icol_point, icol_elem)
end

"""
    update_body_acceleration_indices!(system, prescribed_conditions)

Updates the state variable indices corresponding to the body frame accelerations to 
correspond to the provided prescribed conditions.
"""
function update_body_acceleration_indices!(system, prescribed_conditions)

    for i = 1:6
        ipoint = findfirst(p -> p.pl[i] & p.pd[i], prescribed_conditions)
        if isnothing(ipoint)
            system.indices.icol_body[i] = 0
        else
            system.indices.icol_body[i] = system.indices.icol_point[ipoint]+i-1
        end
    end

    return system
end

"""
    body_accelerations(system, x=system.x; linear_acceleration=zeros(3), angular_acceleration=zeros(3))

Extract the linear and angular acceleration of the body frame from the state vector or 
prescribed accelerations.  Valid for the results of steady state and initial condition 
analyses.
"""
function body_accelerations(system::AbstractSystem, x=system.x; 
    linear_acceleration=(@SVector zeros(3)), 
    angular_acceleration=(@SVector zeros(3))
    )

    return body_accelerations(x, system.indices.icol_indices, linear_acceleration, angular_acceleration)
end

function body_accelerations(x, icol, ab=(@SVector zeros(3)), αb=(@SVector zeros(3)))

    ab = SVector(
        iszero(icol[1]) ? ab[1] : x[icol[1]],
        iszero(icol[2]) ? ab[2] : x[icol[2]],
        iszero(icol[3]) ? ab[3] : x[icol[3]],
    )

    αb = SVector(
        iszero(icol[4]) ? αb[1] : x[icol[4]],
        iszero(icol[5]) ? αb[2] : x[icol[5]],
        iszero(icol[6]) ? αb[3] : x[icol[6]],
    )

    return ab, αb
end

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
    indices = SystemIndices(assembly.start, assembly.stop; static=false, expanded=false)

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

    # initialize current body frame states and time
    t = zero(TF)

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
    indices = SystemIndices(assembly.start, assembly.stop; static=false, expanded=true)

    # initialize system states
    x = zeros(TF, indices.nstates)
    r = zeros(TF, indices.nstates)
    K = spzeros(TF, indices.nstates, indices.nstates)
    M = spzeros(TF, indices.nstates, indices.nstates)

    # initialize current time
    t = zero(TF)

    return ExpandedSystem{TF, Vector{TF}, SparseMatrixCSC{TF, Int64}}(x, r, K, M, indices, 
        force_scaling, t)
end

# default system is a DynamicSystem
const System = DynamicSystem

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
 - `linear_velocity = zeros(3)`: Prescribed body frame linear velocity vector. 
    This vector may also be provided as a function of time.
 - `angular_velocity = zeros(3)`: Prescribed body frame angular velocity vector. 
    This vector may also be provided as a function of time.
 - `linear_acceleration = zeros(3)`: Prescribed body frame linear acceleration vector. 
    This vector may also be provided as a function of time.
 - `angular_acceleration = zeros(3)`: Prescribed body frame angular acceleration vector. 
    This vector may also be provided as a function of time.
 - `gravity = [0,0,0]`: Gravity vector.  This vector may also be provided as a function 
    of time.
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
    linear_velocity=(@SVector zeros(3)),
    angular_velocity=(@SVector zeros(3)),
    linear_acceleration=(@SVector zeros(3)),
    angular_acceleration=(@SVector zeros(3)),
    gravity=(@SVector zeros(3)),
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
    vb_p = typeof(linear_velocity) <: AbstractVector ? SVector{3}(linear_velocity) : SVector{3}(linear_velocity(t))
    ωb_p = typeof(angular_velocity) <: AbstractVector ? SVector{3}(angular_velocity) : SVector{3}(angular_velocity(t))
    ab_p = typeof(linear_acceleration) <: AbstractVector ? SVector{3}(linear_acceleration) : SVector{3}(linear_acceleration(t))
    αb_p = typeof(angular_acceleration) <: AbstractVector ? SVector{3}(angular_acceleration) : SVector{3}(angular_acceleration(t))

    # update body frame acceleration indices for both systems
    update_body_acceleration_indices!(system1, pcond)
    update_body_acceleration_indices!(system2, pcond)

    # extract state vectors for each system
    x1 = system1.x
    x2 = system2.x

    # extract body frame state variable pointers for each system
    icol_body1 = system1.indices.icol_body
    icol_body2 = system2.indices.icol_body

    # extract point state variable pointers for each system
    icol_point1 = system1.indices.icol_point
    icol_point2 = system2.indices.icol_point

    # extract element state variable pointers for each system
    icol_elem1 =  system1.indices.icol_elem
    icol_elem2 =  system2.indices.icol_elem

    # extract force scaling parameter for each system
    force_scaling1 = system1.force_scaling
    force_scaling2 = system2.force_scaling

    # copy over body frame accelerations
    !iszero(icol_body1[1]) && setindex!(x1, x2[icol_body2[1]], icol_body1[1])
    !iszero(icol_body1[2]) && setindex!(x1, x2[icol_body2[2]], icol_body1[2])
    !iszero(icol_body1[3]) && setindex!(x1, x2[icol_body2[3]], icol_body1[3])
    !iszero(icol_body1[4]) && setindex!(x1, x2[icol_body2[4]], icol_body1[4])
    !iszero(icol_body1[5]) && setindex!(x1, x2[icol_body2[5]], icol_body1[5])
    !iszero(icol_body1[6]) && setindex!(x1, x2[icol_body2[6]], icol_body1[6])

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
                properties = steady_element_properties(x2, system2.indices, force_scaling2, 
                    structural_damping, assembly, ielem, pcond, gvec, vb_p, ωb_p, ab_p, αb_p)

                @unpack ωb, C, Cab, CtCab, mass11, mass12, mass21, mass22, V1, V2, Ω1, Ω2, V, Ω = properties

                # linear and angular velocity rates
                V1dot = system2.Vdot[assembly.start[ielem]]
                Ω1dot = system2.Ωdot[assembly.start[ielem]]
            
                V2dot = system2.Vdot[assembly.stop[ielem]]
                Ω2dot = system2.Ωdot[assembly.stop[ielem]]
            
                Vdot = (V1dot + V2dot)/2
                Ωdot = (Ω1dot + Ω2dot)/2

                # linear and angular momentum rates
                CtCabdot = tilde(Ω - ωb)*CtCab
                
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

        if typeof(system1) <: StaticSystem || typeof(system1) <: DynamicSystem

            # copy over internal forces and moments
            x1[icol] = F[1] / force_scaling1
            x1[icol+1] = F[2] / force_scaling1
            x1[icol+2] = F[3] / force_scaling1
            x1[icol+3] = M[1] / force_scaling1
            x1[icol+4] = M[2] / force_scaling1
            x1[icol+5] = M[3] / force_scaling1

        elseif typeof(system1) <: ExpandedSystem

            # copy over internal forces and moments
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
    static_system_residual!(resid, x, indices, two_dimensional, force_scaling, 
        assembly, prescribed_conditions, distributed_loads, point_masses, gravity)

Populate the system residual vector `resid` for a static analysis
"""
function static_system_residual!(resid, x, indices, two_dimensional, force_scaling, 
    assembly, prescribed_conditions, distributed_loads, point_masses, gravity)

    for ipoint = 1:length(assembly.points)
        static_point_residual!(resid, x, indices, force_scaling, assembly, ipoint, 
            prescribed_conditions, point_masses, gravity)
    end

    for ielem = 1:length(assembly.elements)
        static_element_residual!(resid, x, indices, force_scaling, assembly, ielem, 
            prescribed_conditions, distributed_loads, gravity)
    end

    if two_dimensional
        two_dimensional_residual!(resid, x)
    end

    return resid
end

"""
    steady_system_residual!(resid, x, indices, two_dimensional, force_scaling, 
        structural_damping, assembly, prescribed_conditions, distributed_loads, 
        point_masses, gravity, vb, ωb, ab, αb)

Populate the system residual vector `resid` for a steady state analysis
"""
function steady_system_residual!(resid, x, indices, two_dimensional, force_scaling, 
    structural_damping, assembly, prescribed_conditions, distributed_loads, point_masses, 
    gravity, linear_velocity, angular_velocity, linear_acceleration, angular_acceleration)

    # overwrite prescribed body frame accelerations (if necessary)
    linear_acceleration, angular_acceleration = body_accelerations(x, indices.icol_body, 
        linear_acceleration, angular_acceleration)

    # contributions to the residual vector from points
    for ipoint = 1:length(assembly.points)
        steady_point_residual!(resid, x, indices, force_scaling, 
            assembly, ipoint, prescribed_conditions, point_masses, gravity, 
            linear_velocity, angular_velocity, linear_acceleration, angular_acceleration)
    end

    # contributions to the residual vector from elements
    for ielem = 1:length(assembly.elements)
        steady_element_residual!(resid, x, indices, force_scaling, structural_damping, 
            assembly, ielem, prescribed_conditions, distributed_loads, gravity, 
            linear_velocity, angular_velocity, linear_acceleration, angular_acceleration)
    end

    # restrict the analysis to the x-y plane (if requested)
    if two_dimensional
        two_dimensional_residual!(resid, x)
    end

    return resid
end

"""
    initial_system_residual!(resid, x, indices, rate_vars1, rate_vars2, 
        two_dimensional, force_scaling, structural_damping, assembly, prescribed_conditions, 
        distributed_loads, point_masses, gravity, linear_velocity, angular_velocity, 
        linear_acceleration, angular_acceleration, u0, θ0, V0, Ω0, Vdot0, Ωdot0)

Populate the system residual vector `resid` for the initialization of a time domain 
simulation.
"""
function initial_system_residual!(resid, x, indices, rate_vars1, rate_vars2,
    two_dimensional, force_scaling, structural_damping, assembly, prescribed_conditions, 
    distributed_loads, point_masses, gravity, linear_velocity, angular_velocity, 
    linear_acceleration, angular_acceleration, u0, θ0, V0, Ω0, Vdot0, Ωdot0)

    # overwrite prescribed body frame accelerations (if necessary)
    linear_acceleration, angular_acceleration = body_accelerations(x, indices.icol_body, 
        linear_acceleration, angular_acceleration)

    # contributions to the residual vector from points
    for ipoint = 1:length(assembly.points)
        initial_point_residual!(resid, x, indices, rate_vars2, force_scaling, 
            assembly, ipoint, prescribed_conditions, point_masses, gravity, 
            linear_velocity, angular_velocity, linear_acceleration, angular_acceleration, 
            u0, θ0, V0, Ω0, Vdot0, Ωdot0)
    end
    
    # contributions to the residual vector from elements
    for ielem = 1:length(assembly.elements)
        initial_element_residual!(resid, x, indices, rate_vars2, force_scaling, 
            structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, 
            gravity, linear_velocity, angular_velocity, linear_acceleration, angular_acceleration, 
            u0, θ0, V0, Ω0, Vdot0, Ωdot0)
    end

    # replace equilibrium equations, if necessary
    for ipoint = 1:length(assembly.points)
        irow = indices.irow_point[ipoint]
        icol = indices.icol_point[ipoint]
        Vdot, Ωdot = initial_point_velocity_rates(x, ipoint, indices.icol_point, 
            prescribed_conditions, Vdot0, Ωdot0, rate_vars2)
        if haskey(prescribed_conditions, ipoint)
            pd = prescribed_conditions[ipoint].pd
            !pd[1] && !rate_vars1[icol+6] && rate_vars2[icol+6] && setindex!(resid, Vdot[1] - Vdot0[ipoint][1], irow)
            !pd[2] && !rate_vars1[icol+7] && rate_vars2[icol+7] && setindex!(resid, Vdot[2] - Vdot0[ipoint][2], irow+1)
            !pd[3] && !rate_vars1[icol+8] && rate_vars2[icol+8] && setindex!(resid, Vdot[3] - Vdot0[ipoint][3], irow+2)
            !pd[4] && !rate_vars1[icol+9] && rate_vars2[icol+9] && setindex!(resid, Ωdot[1] - Ωdot0[ipoint][1], irow+3)
            !pd[5] && !rate_vars1[icol+10] && rate_vars2[icol+10] && setindex!(resid, Ωdot[2] - Ωdot0[ipoint][2], irow+4)
            !pd[6] && !rate_vars1[icol+11] && rate_vars2[icol+11] && setindex!(resid, Ωdot[3] - Ωdot0[ipoint][3], irow+5)
        else
            !rate_vars1[icol+6] && rate_vars2[icol+6] && setindex!(resid, Vdot[1] - Vdot0[ipoint][1], irow)
            !rate_vars1[icol+7] && rate_vars2[icol+7] && setindex!(resid, Vdot[2] - Vdot0[ipoint][2], irow+1)
            !rate_vars1[icol+8] && rate_vars2[icol+8] && setindex!(resid, Vdot[3] - Vdot0[ipoint][3], irow+2)
            !rate_vars1[icol+9] && rate_vars2[icol+9] && setindex!(resid, Ωdot[1] - Ωdot0[ipoint][1], irow+3)
            !rate_vars1[icol+10] && rate_vars2[icol+10] && setindex!(resid, Ωdot[2] - Ωdot0[ipoint][2], irow+4)
            !rate_vars1[icol+11] && rate_vars2[icol+11] && setindex!(resid, Ωdot[3] - Ωdot0[ipoint][3], irow+5)
        end
    end

    # restrict the analysis to the x-y plane (if requested)
    if two_dimensional
        two_dimensional_residual!(resid, x)
    end

    return resid
end

"""
    newmark_system_residual!(resid, x, indices, two_dimensional, force_scaling, structural_damping, 
        assembly, prescribed_conditions, distributed_loads, point_masses, gravity,
        linear_velocity, angular_velocity, udot_init, θdot_init, Vdot_init, Ωdot_init, dt)

Populate the system residual vector `resid` for a Newmark scheme time marching analysis.
"""
function newmark_system_residual!(resid, x, indices, two_dimensional, force_scaling, structural_damping, 
    assembly, prescribed_conditions, distributed_loads, point_masses, gravity,
    linear_velocity, angular_velocity, udot_init, θdot_init, Vdot_init, Ωdot_init, dt)
    
    # contributions to the residual vector from points
    for ipoint = 1:length(assembly.points)
        newmark_point_residual!(resid, x, indices, force_scaling, assembly, ipoint, 
            prescribed_conditions, point_masses, gravity, linear_velocity, angular_velocity, 
            udot_init, θdot_init, Vdot_init, Ωdot_init, dt)
    end
    
    # contributions to the residual vector from elements
    for ielem = 1:length(assembly.elements)
        newmark_element_residual!(resid, x, indices, force_scaling, structural_damping, 
            assembly, ielem, prescribed_conditions, distributed_loads, gravity, 
            linear_velocity, angular_velocity, Vdot_init, Ωdot_init, dt)
    end
    
    # restrict the analysis to the x-y plane (if requested)
    if two_dimensional
        two_dimensional_residual!(resid, x)
    end

    return resid
end

"""
    dynamic_system_residual!(resid, dx, x, indices, two_dimensional, force_scaling, 
        structural_damping, assembly, prescribed_conditions, distributed_loads, 
        point_masses, gravity, linear_velocity, angular_velocity)

Populate the system residual vector `resid` for a general dynamic analysis.
"""
function dynamic_system_residual!(resid, dx, x, indices, two_dimensional, force_scaling, 
    structural_damping, assembly, prescribed_conditions, distributed_loads, point_masses, 
    gravity, linear_velocity, angular_velocity)

    # contributions to the residual vector from points
    for ipoint = 1:length(assembly.points)
        dynamic_point_residual!(resid, dx, x, indices, force_scaling, assembly, 
            ipoint, prescribed_conditions, point_masses, gravity, linear_velocity, angular_velocity)
    end
    
    # contributions to the residual vector from elements
    for ielem = 1:length(assembly.elements)
        dynamic_element_residual!(resid, dx, x, indices, force_scaling, structural_damping, 
            assembly, ielem, prescribed_conditions, distributed_loads, gravity, 
            linear_velocity, angular_velocity)
    end

    # restrict the analysis to the x-y plane (if requested)
    if two_dimensional
        two_dimensional_residual!(resid, x)
    end
    
    return resid
end

"""
    expanded_steady_system_residual!(resid, x, indices, two_dimensional, force_scaling, structural_damping, 
        assembly, prescribed_conditions, distributed_loads, point_masses, gravity, 
        linear_velocity, angular_velocity, linear_acceleration, angular_acceleration)

Populate the system residual vector `resid` for a constant mass matrix system.
"""
function expanded_steady_system_residual!(resid, x, indices, two_dimensional, force_scaling, structural_damping, 
    assembly, prescribed_conditions, distributed_loads, point_masses, gravity, 
    linear_velocity, angular_velocity, linear_acceleration, angular_acceleration)

    # overwrite prescribed body frame accelerations (if necessary)
    linear_acceleration, angular_acceleration = body_accelerations(x, indices.icol_body, 
        linear_acceleration, angular_acceleration)

    # point residuals
    for ipoint = 1:length(assembly.points)
        expanded_steady_point_residual!(resid, x, indices, force_scaling, assembly, ipoint, 
            prescribed_conditions, point_masses, gravity, linear_velocity, angular_velocity, 
            linear_acceleration, angular_acceleration)
    end
    
    # element residuals
    for ielem = 1:length(assembly.elements)
        expanded_steady_element_residual!(resid, x, indices, force_scaling, structural_damping, 
            assembly, ielem, prescribed_conditions, distributed_loads, gravity, 
            linear_velocity, angular_velocity, linear_acceleration, angular_acceleration)
    end

    # restrict the analysis to the x-y plane (if requested)
    if two_dimensional
        two_dimensional_residual!(resid, x)
    end
    
    return resid
end

"""
    expanded_dynamic_system_residual!(resid, dx, x, indices, two_dimensional, force_scaling, 
        structural_damping, assembly, prescribed_conditions, distributed_loads, 
        point_masses, gravity, linear_velocity, angular_velocity)

Populate the system residual vector `resid` for a constant mass matrix system.
"""
function expanded_dynamic_system_residual!(resid, dx, x, indices, two_dimensional, force_scaling, 
    structural_damping, assembly, prescribed_conditions, distributed_loads, point_masses, 
    gravity, linear_velocity, angular_velocity)

    # point residuals
    for ipoint = 1:length(assembly.points)
        expanded_dynamic_point_residual!(resid, dx, x, indices, force_scaling, assembly, ipoint, 
            prescribed_conditions, point_masses, gravity, linear_velocity, angular_velocity)
    end
    
    # element residuals
    for ielem = 1:length(assembly.elements)
        expanded_dynamic_element_residual!(resid, dx, x, indices, force_scaling, structural_damping, 
            assembly, ielem, prescribed_conditions, distributed_loads, gravity, 
            linear_velocity, angular_velocity)
    end
    
    # restrict the analysis to the x-y plane (if requested)
    if two_dimensional
        two_dimensional_residual!(resid, x)
    end

    return resid
end

"""
    static_system_jacobian!(jacob, x, indices, two_dimensional, force_scaling, 
        assembly, prescribed_conditions, distributed_loads, point_masses, gravity)

Populate the system jacobian matrix `jacob` for a static analysis
"""
function static_system_jacobian!(jacob, x, indices, two_dimensional, force_scaling, 
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

    if two_dimensional
        two_dimensional_jacobian!(jacob, x)
    end
    
    return jacob
end

"""
    steady_system_jacobian!(jacob, x, indices, two_dimensional, force_scaling, 
        structural_damping, assembly, prescribed_conditions, distributed_loads, 
        point_masses, gravity, linear_velocity, angular_velocity, linear_acceleration, 
        angular_acceleration)

Populate the system jacobian matrix `jacob` for a steady-state analysis
"""
function steady_system_jacobian!(jacob, x, indices, two_dimensional, force_scaling, structural_damping, 
    assembly, prescribed_conditions, distributed_loads, point_masses, gravity, 
    linear_velocity, angular_velocity, linear_acceleration, angular_acceleration)

    jacob .= 0
    
    # overwrite prescribed body frame accelerations (if necessary)
    linear_acceleration, angular_acceleration = body_accelerations(x, indices.icol_body, 
        linear_acceleration, angular_acceleration)

    for ipoint = 1:length(assembly.points)
        steady_point_jacobian!(jacob, x, indices, force_scaling, assembly, 
            ipoint, prescribed_conditions, point_masses, gravity, linear_velocity, 
            angular_velocity, linear_acceleration, angular_acceleration)
    end
    
    for ielem = 1:length(assembly.elements)
        steady_element_jacobian!(jacob, x, indices, force_scaling, structural_damping, 
            assembly, ielem, prescribed_conditions, distributed_loads, gravity, 
            linear_velocity, angular_velocity, linear_acceleration, angular_acceleration)
    end

    if two_dimensional
        two_dimensional_jacobian!(jacob, x)
    end
    
    return jacob
end

"""
    initial_system_jacobian!(jacob, x, indices, rate_vars1, rate_vars2, two_dimensional, force_scaling, 
        structural_damping, assembly, prescribed_conditions, distributed_loads, 
        point_masses, gravity, linear_velocity, angular_velocity, linear_acceleration, 
        angular_acceleration, u0, θ0, V0, Ω0, Vdot0, Ωdot0)

Populate the system jacobian matrix `jacob` for the initialization of a time domain 
simulation.
"""
function initial_system_jacobian!(jacob, x, indices, rate_vars1, rate_vars2, two_dimensional, force_scaling, 
    structural_damping, assembly, prescribed_conditions, distributed_loads, point_masses, 
    gravity, linear_velocity, angular_velocity, linear_acceleration, angular_acceleration,
    u0, θ0, V0, Ω0, Vdot0, Ωdot0)
    
    jacob .= 0

    # overwrite prescribed body frame accelerations (if necessary)
    linear_acceleration, angular_acceleration = body_accelerations(x, indices.icol_body, 
        linear_acceleration, angular_acceleration)
    
    for ipoint = 1:length(assembly.points)
        initial_point_jacobian!(jacob, x, indices, rate_vars2, force_scaling, 
            assembly, ipoint, prescribed_conditions, point_masses, gravity,  
            linear_velocity, angular_velocity, linear_acceleration, angular_acceleration, 
            u0, θ0, V0, Ω0, Vdot0, Ωdot0)
    end
    
    for ielem = 1:length(assembly.elements)
        initial_element_jacobian!(jacob, x, indices, rate_vars2, force_scaling, 
            structural_damping, assembly, ielem, prescribed_conditions, 
            distributed_loads, gravity, linear_velocity, angular_velocity, 
            linear_acceleration, angular_acceleration, u0, θ0, V0, Ω0, Vdot0, Ωdot0)
    end

    # replace equilibrium equations, if necessary
    for ipoint = 1:length(assembly.points)
        irow = indices.irow_point[ipoint]
        icol = indices.icol_point[ipoint]
        if haskey(prescribed_conditions, ipoint)
            pd = prescribed_conditions[ipoint].pd
            # displacements not prescribed, Vdot and Ωdot are arbitrary, F and M from compatability
            if !pd[1] && !rate_vars1[icol+6] && rate_vars2[icol+6]
                jacob[irow,:] .= 0
                jacob[irow,icol] = 1
            end
            if !pd[2] && !rate_vars1[icol+7] && rate_vars2[icol+7]
                jacob[irow+1,:] .= 0
                jacob[irow+1,icol+1] = 1
            end
            if !pd[3] && !rate_vars1[icol+8] && rate_vars2[icol+8]
                jacob[irow+2,:] .= 0
                jacob[irow+2,icol+2] = 1
            end
            if !pd[4] && !rate_vars1[icol+9] && rate_vars2[icol+9]
                jacob[irow+3,:] .= 0
                jacob[irow+3,icol+3] = 1
            end
            if !pd[5] && !rate_vars1[icol+10] && rate_vars2[icol+10]
                jacob[irow+4,:] .= 0
                jacob[irow+4,icol+4] = 1
            end
            if !pd[6] && !rate_vars1[icol+11] && rate_vars2[icol+11]
                jacob[irow+5,:] .= 0
                jacob[irow+5,icol+5] = 1
            end
        else
            if !rate_vars1[icol+6] && rate_vars2[icol+6]
                jacob[irow,:] .= 0
                jacob[irow,icol] = 1
            end
            if !rate_vars1[icol+7] && rate_vars2[icol+7]
                jacob[irow+1,:] .= 0
                jacob[irow+1,icol+1] = 1
            end
            if !rate_vars1[icol+8] && rate_vars2[icol+8]
                jacob[irow+2,:] .= 0
                jacob[irow+2,icol+2] = 1
            end
            if !rate_vars1[icol+9] && rate_vars2[icol+9]
                jacob[irow+3,:] .= 0
                jacob[irow+3,icol+3] = 1
            end
            if !rate_vars1[icol+10] && rate_vars2[icol+10]
                jacob[irow+4,:] .= 0
                jacob[irow+4,icol+4] = 1
            end
            if !rate_vars1[icol+11] && rate_vars2[icol+11]
                jacob[irow+5,:] .= 0
                jacob[irow+5,icol+5] = 1
            end
        end
    end

    if two_dimensional
        two_dimensional_jacobian!(jacob, x)
    end

    return jacob
end

"""
    newmark_system_jacobian!(jacob, x, indices, two_dimensional, force_scaling, structural_damping, 
        assembly, prescribed_conditions, distributed_loads, point_masses, gravity, 
        linear_velocity, angular_velocity, udot_init, θdot_init, Vdot_init, Ωdot_init, dt)

Populate the system jacobian matrix `jacob` for a Newmark scheme time marching analysis.
"""
function newmark_system_jacobian!(jacob, x, indices, two_dimensional, force_scaling, structural_damping, 
    assembly, prescribed_conditions, distributed_loads, point_masses, gravity, 
    linear_velocity, angular_velocity, udot_init, θdot_init, Vdot_init, Ωdot_init, dt)
    
    jacob .= 0

    for ipoint = 1:length(assembly.points)
        newmark_point_jacobian!(jacob, x, indices, force_scaling, assembly, ipoint, 
            prescribed_conditions, point_masses, gravity, linear_velocity, angular_velocity, 
            udot_init, θdot_init, Vdot_init, Ωdot_init, dt)
    end
    
    for ielem = 1:length(assembly.elements)
        newmark_element_jacobian!(jacob, x, indices, force_scaling, structural_damping, 
            assembly, ielem, prescribed_conditions, distributed_loads, gravity, 
            linear_velocity, angular_velocity, Vdot_init, Ωdot_init, dt)
    end

    if two_dimensional
        two_dimensional_jacobian!(jacob, x)
    end
    
    return jacob
end

"""
    dynamic_system_jacobian!(jacob, dx, x, indices, two_dimensional, force_scaling, 
        structural_damping, assembly, prescribed_conditions, distributed_loads, 
        point_masses, gravity, linear_velocity, angular_velocity)

Populate the system jacobian matrix `jacob` for a general dynamic analysis.
"""
function dynamic_system_jacobian!(jacob, dx, x, indices, two_dimensional, force_scaling, 
    structural_damping, assembly, prescribed_conditions, distributed_loads, point_masses, 
    gravity, linear_velocity, angular_velocity)
    
    jacob .= 0

    for ipoint = 1:length(assembly.points)
        dynamic_point_jacobian!(jacob, dx, x, indices, force_scaling, assembly, 
            ipoint, prescribed_conditions, point_masses, gravity,
            linear_velocity, angular_velocity)
    end
    
    for ielem = 1:length(assembly.elements)
        dynamic_element_jacobian!(jacob, dx, x, indices, force_scaling, structural_damping, 
            assembly, ielem, prescribed_conditions, distributed_loads, gravity, 
            linear_velocity, angular_velocity)
    end
    
    if two_dimensional
        two_dimensional_jacobian!(jacob, x)
    end

    return jacob
end

"""
    expanded_steady_system_jacobian!(jacob, x, indices, two_dimensional, force_scaling, structural_damping, 
        assembly, prescribed_conditions, distributed_loads, point_masses, gravity, 
        ub_p, θb_p, vb_p, ωb_p, ab_p, αb_p)

Populate the system jacobian matrix `jacob` for a general dynamic analysis with a 
constant mass matrix system.
"""
function expanded_steady_system_jacobian!(jacob, x, indices, two_dimensional, force_scaling, 
    structural_damping, assembly, prescribed_conditions, distributed_loads, point_masses, 
    gravity, linear_velocity, angular_velocity, linear_acceleration, angular_acceleration)
    
    jacob .= 0

    # overwrite prescribed body frame accelerations (if necessary)
    linear_acceleration, angular_acceleration = body_accelerations(x, indices.icol_body, 
        linear_acceleration, angular_acceleration)

    for ipoint = 1:length(assembly.points)
        expanded_steady_point_jacobian!(jacob, x, indices, force_scaling, 
            assembly, ipoint, prescribed_conditions, point_masses, gravity, 
            linear_velocity, angular_velocity, linear_acceleration, angular_acceleration)
    end
    
    for ielem = 1:length(assembly.elements)
        expanded_steady_element_jacobian!(jacob, x, indices, force_scaling, 
            structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, 
            gravity, linear_velocity, angular_velocity, linear_acceleration, angular_acceleration)
    end
   
    if two_dimensional
        two_dimensional_jacobian!(jacob, x)
    end

    return jacob
end

"""
    expanded_dynamic_system_jacobian!(jacob, dx, x, indices, two_dimensional, force_scaling, structural_damping, 
        assembly, prescribed_conditions, distributed_loads, point_masses, gravity, 
        linear_velocity, angular_velocity)

Populate the system jacobian matrix `jacob` for a general dynamic analysis with a 
constant mass matrix system.
"""
function expanded_dynamic_system_jacobian!(jacob, dx, x, indices, two_dimensional, force_scaling, 
    structural_damping, assembly, prescribed_conditions, distributed_loads, point_masses, 
    gravity, linear_velocity, angular_velocity)
    
    jacob .= 0

    for ipoint = 1:length(assembly.points)
        expanded_dynamic_point_jacobian!(jacob, dx, x, indices, force_scaling, assembly, 
            ipoint, prescribed_conditions, point_masses, gravity, linear_velocity, angular_velocity)
    end
    
    for ielem = 1:length(assembly.elements)
        expanded_dynamic_element_jacobian!(jacob, dx, x, indices, force_scaling, 
            structural_damping, assembly, ielem, prescribed_conditions, distributed_loads, 
            gravity, linear_velocity, angular_velocity)
    end
    
    if two_dimensional
        two_dimensional_jacobian!(jacob, x)
    end

    return jacob
end

"""
    system_mass_matrix!(jacob, x, indices, two_dimensional, force_scaling,  assembly, 
        prescribed_conditions, point_masses)

Calculate the jacobian of the residual expressions with respect to the state rates.
"""
function system_mass_matrix!(jacob, x, indices, two_dimensional, force_scaling, assembly, 
    prescribed_conditions, point_masses)

    jacob .= 0

    gamma = 1

    system_mass_matrix!(jacob, gamma, x, indices, two_dimensional, force_scaling,  assembly, 
        prescribed_conditions, point_masses)

    return jacob
end

"""
    system_mass_matrix!(jacob, gamma, x, indices, two_dimensional, force_scaling, assembly, 
        prescribed_conditions, point_masses)

Calculate the jacobian of the residual expressions with respect to the state rates and 
add the result multiplied by `gamma` to `jacob`.
"""
function system_mass_matrix!(jacob, gamma, x, indices, two_dimensional, force_scaling, assembly, 
    prescribed_conditions, point_masses)

    for ipoint = 1:length(assembly.points)
        mass_matrix_point_jacobian!(jacob, gamma, x, indices, two_dimensional, force_scaling, 
            assembly, ipoint, prescribed_conditions, point_masses)
    end
    
    for ielem = 1:length(assembly.elements)
        mass_matrix_element_jacobian!(jacob, gamma, x, indices, two_dimensional, force_scaling, 
            assembly, ielem, prescribed_conditions)
    end

    return jacob
end

"""
    expanded_system_mass_matrix(system, assembly;
        two_dimensional = false,
        prescribed_conditions=Dict{Int, PrescribedConditions}(), 
        point_masses=Dict{Int, PointMass}())

Calculate the jacobian of the residual expressions with respect to the state rates for a 
constant mass matrix system.
"""
function expanded_system_mass_matrix(system, assembly;
    two_dimensional=false,
    prescribed_conditions=Dict{Int, PrescribedConditions}(), 
    point_masses=Dict{Int, PointMass}())

    @unpack indices, force_scaling = system

    TF = eltype(system)
    nx = expanded_indices.nstates
    jacob = spzeros(TF, nx, nx)
    gamma = -1
    pcond = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(0)
    pmass = typeof(point_masses) <: AbstractDict ? point_masses : point_masses(0)

    expanded_system_mass_matrix!(jacob, gamma, indices, two_dimensional, force_scaling, assembly, pcond, pmass) 

    return jacob
end

"""
    expanded_system_mass_matrix!(jacob, indices, two_dimensional, force_scaling,  assembly, prescribed_conditions, 
        point_masses)

Calculate the jacobian of the residual expressions with respect to the state rates.
"""
function expanded_system_mass_matrix!(jacob, indices, two_dimensional, force_scaling, assembly, 
    prescribed_conditions, point_masses)

    jacob .= 0

    gamma = 1

    expanded_system_mass_matrix!(jacob, gamma, indices, two_dimensional, force_scaling, assembly, 
        prescribed_conditions, point_masses)

    return jacob
end

"""
    expanded_system_mass_matrix!(jacob, gamma, indices, two_dimensional, force_scaling, assembly, 
        prescribed_conditions, point_masses)

Calculate the jacobian of the residual expressions with respect to the state rates and 
add the result multiplied by `gamma` to `jacob`.
"""
function expanded_system_mass_matrix!(jacob, gamma, indices, two_dimensional, force_scaling, assembly, 
    prescribed_conditions, point_masses)

    for ipoint = 1:length(assembly.points)
        expanded_mass_matrix_point_jacobian!(jacob, gamma, indices, two_dimensional, 
            force_scaling, assembly, ipoint, prescribed_conditions, point_masses)
    end
    
    for ielem = 1:length(assembly.elements)
        expanded_mass_matrix_element_jacobian!(jacob, gamma, indices, two_dimensional, 
            force_scaling, assembly, ielem, prescribed_conditions)
    end

    return jacob
end

function two_dimensional_residual!(resid, x)

    for (irow, icol) in zip(1:6:length(x), 1:6:length(x))
        resid[irow+2] = x[icol+2] # constrain linear component in z-direction to be zero
        resid[irow+3] = x[icol+3] # constrain angular component in x-direction to be zero
        resid[irow+4] = x[icol+4] # constrain angular component in y-direction to be zero
    end

    return resid
end

function two_dimensional_jacobian!(jacob, x)

    for (irow, icol) in zip(1:6:length(x), 1:6:length(x))
        # constrain linear component in z-direction to be zero
        jacob[irow+2,:] .= 0 
        jacob[irow+2,icol+2] = 1 
        # constrain angular component in x-direction to be zero
        jacob[irow+3,:] .= 0 
        jacob[irow+3,icol+3] = 1 
        # constrain angular component in y-direction to be zero
        jacob[irow+4,:] .= 0 
        jacob[irow+4,icol+4] = 1 
    end

    return jacob
end