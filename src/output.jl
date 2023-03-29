"""
    PointState

Holds the state variables for a point

# Fields:
 - `u`: Linear displacement
 - `udot`: Linear displacement rate
 - `theta`: Angular displacement (Wiener-Milenkovic parameters)
 - `thetadot`: Angular displacement rate
 - `V`: Linear velocity
 - `Vdot`: Linear velocity rate
 - `Omega`: Angular velocity
 - `Omegadot`: Angular velocity rate
 - `F`: External forces
 - `M`: External moments
"""
struct PointState{TF}
    u::SVector{3, TF}
    udot::SVector{3, TF}
    theta::SVector{3, TF}
    thetadot::SVector{3, TF}
    V::SVector{3, TF}
    Vdot::SVector{3, TF}
    Omega::SVector{3, TF}
    Omegadot::SVector{3, TF}
    F::SVector{3, TF}
    M::SVector{3, TF}
end
Base.eltype(::PointState{TF}) where {TF} = TF
Base.eltype(::Type{PointState{TF}}) where {TF} = TF

PointState{TF}(p::PointState) where {TF} = PointState{TF}(p.u, p.udot, p.theta, p.thetadot,
    p.V, p.Vdot, p.Omega, p.Omegadot, p.F, p.M)
Base.convert(::Type{PointState{TF}}, p::PointState) where {TF} = PointState{TF}(p)

"""
    ElementState

Holds the state variables for an element

# Fields:
 - `u`: Linear displacement
 - `udot`: Linear displacement rate
 - `theta`: Angular displacement (Wiener-Milenkovic parameters)
 - `thetadot`: Angular displacement rate
 - `V`: Linear velocity
 - `Vdot`: Linear velocity rate
 - `Omega`: Angular velocity
 - `Omegadot`: Angular velocity rate
 - `Fi`: Internal forces
 - `Mi`: Internal moments
"""
struct ElementState{TF}
    u::SVector{3, TF}
    udot::SVector{3, TF}
    theta::SVector{3, TF}
    thetadot::SVector{3, TF}
    V::SVector{3, TF}
    Vdot::SVector{3, TF}
    Omega::SVector{3, TF}
    Omegadot::SVector{3, TF}
    Fi::SVector{3, TF}
    Mi::SVector{3, TF}
end
Base.eltype(::ElementState{TF}) where {TF} = TF
Base.eltype(::Type{ElementState{TF}}) where {TF} = TF

ElementState{TF}(p::ElementState) where {TF} = ElementState{TF}(p.u, p.udot, p.theta, p.thetadot,
    p.V, p.Vdot, p.Omega, p.Omegadot, p.Fi, p.Mi)
Base.convert(::Type{ElementState{TF}}, p::ElementState) where {TF} = ElementState{TF}(p)

"""
    AssemblyState{TF, TP<:AbstractVector{PointState{TF}}, TE<:AbstractVector{ElementState{TF}}}

Struct for storing state variables for the points and elements in an assembly.

# Fields:
 - `points::TP`: Array of [`PointState`](@ref) for each point in the assembly
 - `elements::TE`: Array of [`ElementState`](@ref) for each element in the assembly
"""
struct AssemblyState{TF, TP<:AbstractVector{PointState{TF}}, TE<:AbstractVector{ElementState{TF}}}
    points::TP
    elements::TE
end
Base.eltype(::AssemblyState{TF, TP, TE}) where {TF, TP, TE} = TF
Base.eltype(::Type{AssemblyState{TF, TP, TE}}) where {TF, TP, TE} = TF

function AssemblyState{TF, TP, TE}(a::AssemblyState) where {TF, TP, TE}
    return AssemblyState{TF, TP, TE}(a.points, a.elements)
end
Base.convert(::Type{AssemblyState{TF, TP, TE}}, a::AssemblyState) where {TF, TP, TE} = AssemblyState{TF, TP, TE}(a)

"""
    AssemblyState(system, assembly; prescribed_conditions = Dict())
    AssemblyState(x, system, assembly; prescribed_conditions = Dict())
    AssemblyState(dx, x, system, assembly; prescribed_conditions = Dict())

Post-process the system state given the state vector `x` and rate vector `dx`.  Return an
object of type `AssemblyState` that defines the state of the assembly for the time step.

If `prescribed_conditions` is not provided, all point state variables are assumed
to be displacements/rotations, rather than their actual identities as used in the
analysis.
"""
function AssemblyState(args...; kwargs...)

    points = extract_point_states(args...; kwargs...)

    elements = extract_element_states(args...; kwargs...)

    return AssemblyState(points, elements)
end

"""
    extract_point_state(system, assembly, ipoint; prescribed_conditions = Dict())
    extract_point_state(x, system, assembly, ipoint; prescribed_conditions = Dict())
    extract_point_state(dx, x, system, assembly, ipoint; prescribed_conditions = Dict())

Return the state variables corresponding to point `ipoint` (see [`PointState`](@ref))
given the state vector `x` and rate vector `dx`.

If `prescribed_conditions` is not provided, all point state variables are assumed
to be displacements/rotations, rather than their actual identities as used in the
analysis.
"""
extract_point_state

function extract_point_state(system::StaticSystem, assembly, ipoint; kwargs...)

    return extract_point_state(system.x, system, assembly, ipoint; kwargs...)
end

function extract_point_state(system::Union{DynamicSystem,ExpandedSystem}, assembly, ipoint; kwargs...)

    return extract_point_state(system.dx, system.x, system, assembly, ipoint; kwargs...)
end

function extract_point_state(x, system::Union{DynamicSystem,ExpandedSystem}, assembly, ipoint; kwargs...)

    dx = FillArrays.Fill(NaN, length(x))

    return extract_point_state(dx, x, system::DynamicSystem, assembly, ipoint; kwargs...)
end

function extract_point_state(x, system::StaticSystem, assembly, ipoint;
    prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}())

    @unpack force_scaling, indices, t = system

    # check state variable length
    @assert length(x) == length(system.x) "state vector length does not match with the provided system"

    # get current prescribed conditions
    pc = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(t)

    # extract point state variables
    u, θ = point_displacement(x, ipoint, indices.icol_point, pc)
    F, M = point_loads(x, ipoint, indices.icol_point, force_scaling, pc)
    V = Ω = @SVector zeros(3)

    # extract point rate variables
    udot = θdot = Vdot = Ωdot = @SVector zeros(3)

    # convert rotation parameter to Wiener-Milenkovic parameters
    scaling = rotation_parameter_scaling(θ)
    θ *= scaling

    # promote all variables to the same type
    u, udot, θ, θdot, V, Vdot, Ω, Ωdot, F, M = promote(u, udot, θ, θdot, V, Vdot, Ω, Ωdot, F, M)

    return PointState(u, udot, θ, θdot, V, Vdot, Ω, Ωdot, F, M)
end

function extract_point_state(dx, x, system::DynamicSystem, assembly, ipoint;
    prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}())

    @unpack force_scaling, indices, t = system

    # check state variable length
    @assert length(x) == length(system.x) "state vector length does not match with the provided system"

    # get current prescribed conditions
    pc = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(t)

    # extract point state variables
    u, θ = point_displacement(x, ipoint, indices.icol_point, pc)
    F, M = point_loads(x, ipoint, indices.icol_point, force_scaling, pc)
    V, Ω = point_velocities(x, ipoint, indices.icol_point)

    # extract point rate variables
    udot, θdot = point_displacement_rates(dx, ipoint, indices.icol_point, pc)
    Vdot, Ωdot = point_velocities(dx, ipoint, indices.icol_point)

    # convert rotation parameter to Wiener-Milenkovic parameters
    scaling = rotation_parameter_scaling(θ)
    θ *= scaling

    # promote all variables to the same type
    u, udot, θ, θdot, V, Vdot, Ω, Ωdot, F, M = promote(u, udot, θ, θdot, V, Vdot, Ω, Ωdot, F, M)

    return PointState(u, udot, θ, θdot, V, Vdot, Ω, Ωdot, F, M)
end

function extract_point_state(dx, x, system::ExpandedSystem, assembly, ipoint;
    prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}())

    # system variables
    @unpack force_scaling, indices, t = system

    # check state variable length
    @assert length(x) == length(system.x) "state vector length does not match with the provided system"

    # get current prescribed conditions
    pc = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(t)

    # extract point state variables
    u, θ = point_displacement(x, ipoint, indices.icol_point, pc)
    V, Ω = point_velocities(x, ipoint, indices.icol_point)
    F, M = point_loads(x, ipoint, indices.icol_point, force_scaling, pc)

    # extract point rate variables
    udot, θdot = point_displacement_rates(dx, ipoint, indices.icol_point, pc)
    Vdot, Ωdot = point_velocities(dx, ipoint, indices.icol_point)

    # rotate linear and angular velocities into the body frame
    C = get_C(θ)
    V = C'*V
    Ω = C'*Ω

    # rotate linear and angular velocity rates into the body frame
    Cdot = -C*tilde(C'*get_Q(θ)*θdot)
    Vdot = Cdot'*V + C'*Vdot
    Ωdot = Cdot'*Ω + C'*Ωdot

    # convert rotation parameters to Wiener-Milenkovic parameters
    k = rotation_parameter_scaling(θ)
    kdot = rotation_parameter_scaling_θ(k, θ)'*θdot

    θ = k*θ
    θdot = kdot*θ + k*θdot

    # promote all variables to the same type
    u, udot, θ, θdot, V, Vdot, Ω, Ωdot, F, M = promote(u, udot, θ, θdot, V, Vdot, Ω, Ωdot, F, M)

    return PointState(u, udot, θ, θdot, V, Vdot, Ω, Ωdot, F, M)
end

"""
    extract_point_states(system, assembly; prescribed_conditions = Dict())
    extract_point_states(x, system, assembly; prescribed_conditions = Dict())
    extract_point_states(dx, x, system, assembly; prescribed_conditions = Dict())

Return the state variables corresponding to each point (see [`PointState`](@ref))
given the state vector `x` and rate vector `dx`.

If `prescribed_conditions` is not provided, all point state variables are assumed
to be displacements/rotations, rather than their actual identities as used in the
analysis.
"""
extract_point_states

function extract_point_states(system::StaticSystem, assembly; kwargs...)

    return extract_point_states(system.x, system, assembly; kwargs...)
end

function extract_point_states(system::Union{DynamicSystem,ExpandedSystem}, assembly; kwargs...)

    return extract_point_states(system.dx, system.x, system, assembly; kwargs...)
end

function extract_point_states(x, system::Union{DynamicSystem,ExpandedSystem}, assembly; kwargs...)

    dx = FillArrays.Fill(NaN, length(x))

    return extract_point_states(dx, x, system, assembly; kwargs...)
end

function extract_point_states(x, system::StaticSystem, assembly; kwargs...)

    return [extract_point_state(x, system, assembly, ipoint; kwargs...) for ipoint in eachindex(assembly.points)]
end

function extract_point_states(dx, x, system::Union{DynamicSystem,ExpandedSystem}, assembly; kwargs...)

    return [extract_point_state(dx, x, system, assembly, ipoint; kwargs...) for ipoint in eachindex(assembly.points)]
end

"""
    extract_element_state(system, assembly, ielem; prescribed_conditions = Dict())
    extract_element_state(x, system, assembly, ielem; prescribed_conditions = Dict())
    extract_element_state(dx, x, system, assembly, ielem; prescribed_conditions = Dict())

Return the state variables corresponding to element `ielem` (see [`ElementState`](@ref))
given the state vector `x` and rate vector `dx`.
"""
extract_element_state

function extract_element_state(system::StaticSystem, assembly, ielem; kwargs...)

    return extract_element_state(system.x, system, assembly, ielem; kwargs...)
end

function extract_element_state(system::Union{DynamicSystem,ExpandedSystem}, assembly, ielem; kwargs...)

    return extract_element_state(system.dx, system.x, system, assembly, ielem; kwargs...)
end

function extract_element_state(x, system::Union{DynamicSystem,ExpandedSystem}, assembly, ielem; kwargs...)

    dx = FillArrays.Fill(NaN, length(x))

    return extract_element_state(dx, x, system::DynamicSystem, assembly, ielem; kwargs...)
end

function extract_element_state(x, system::StaticSystem, assembly, ielem;
    prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}())

    # system variables
    @unpack force_scaling, indices, t = system

    # check state variable length
    @assert length(x) == length(system.x) "state vector length does not match with the provided system"

    # current prescribed conditions
    pc = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(t)

    # linear and angular displacement
    u1, θ1 = point_displacement(x, assembly.start[ielem], indices.icol_point, pc)
    u2, θ2 = point_displacement(x, assembly.stop[ielem], indices.icol_point, pc)
    u = (u1 + u2)/2
    θ = (θ1 + θ2)/2

    # internal forces and moments
    F, M = element_loads(x, ielem, indices.icol_elem, force_scaling)

    # linear and angular velocity
    V = Ω = @SVector zeros(3)

    # extract point rate variables
    udot = θdot = Vdot = Ωdot = @SVector zeros(3)

    # convert rotation parameter to Wiener-Milenkovic parameters
    scaling = rotation_parameter_scaling(θ)
    θ *= scaling

    # promote all variables to the same type
    u, udot, θ, θdot, V, Vdot, Ω, Ωdot, F, M = promote(u, udot, θ, θdot, V, Vdot, Ω, Ωdot, F, M)

    return ElementState(u, udot, θ, θdot, V, Vdot, Ω, Ωdot, F, M)
end

function extract_element_state(dx, x, system::DynamicSystem, assembly, ielem;
    prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}())

    # system variables
    @unpack force_scaling, indices, t = system

    # check state variable length
    @assert length(x) == length(system.x) "state vector length does not match with the provided system"

    # current prescribed conditions
    pc = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(t)

    # linear and angular displacement
    u1, θ1 = point_displacement(x, assembly.start[ielem], indices.icol_point, pc)
    u2, θ2 = point_displacement(x, assembly.stop[ielem], indices.icol_point, pc)
    u = (u1 + u2)/2
    θ = (θ1 + θ2)/2

    # internal forces and moments
    F, M = element_loads(x, ielem, indices.icol_elem, force_scaling)

    # linear and angular velocity
    V1, Ω1 = point_velocities(x, assembly.start[ielem], indices.icol_point)
    V2, Ω2 = point_velocities(x, assembly.stop[ielem], indices.icol_point)
    V = (V1 + V2)/2
    Ω = (Ω1 + Ω2)/2

    # linear and angular displacement rates
    u1dot, θ1dot = point_displacement_rates(dx, assembly.start[ielem], indices.icol_point, pc)
    u2dot, θ2dot = point_displacement_rates(dx, assembly.stop[ielem], indices.icol_point, pc)

    udot = (u1dot + u2dot)/2
    θdot = (θ1dot + θ2dot)/2

    # linear and angular velocity rates
    V1dot, Ω1dot = point_velocities(dx, assembly.start[ielem], indices.icol_point)
    V2dot, Ω2dot = point_velocities(dx, assembly.stop[ielem], indices.icol_point)

    Vdot = (V1dot + V2dot)/2
    Ωdot = (Ω1dot + Ω2dot)/2

    # convert rotation parameter to Wiener-Milenkovic parameters
    scaling = rotation_parameter_scaling(θ)
    θ *= scaling

    # promote all variables to the same type
    u, udot, θ, θdot, V, Vdot, Ω, Ωdot, F, M = promote(u, udot, θ, θdot, V, Vdot, Ω, Ωdot, F, M)

    return ElementState(u, udot, θ, θdot, V, Vdot, Ω, Ωdot, F, M)
end

function extract_element_state(dx, x, system::ExpandedSystem, assembly, ielem;
    prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}())

    # system variables
    @unpack force_scaling, indices, t = system

    # check state variable length
    @assert length(x) == length(system.x) "state vector length does not match with the provided system"

    # current prescribed conditions
    pc = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(t)

    # linear and angular displacement
    u1, θ1 = point_displacement(x, assembly.start[ielem], indices.icol_point, pc)
    u2, θ2 = point_displacement(x, assembly.stop[ielem], indices.icol_point, pc)
    u = (u1 + u2)/2
    θ = (θ1 + θ2)/2

    # element forces and moments
    F1, M1, F2, M2 = expanded_element_loads(x, ielem, indices.icol_elem, force_scaling)
    F = (F1 + F2)/2
    M = (M1 + M2)/2

    # linear and angular velocity
    V, Ω = expanded_element_velocities(x, ielem, indices.icol_elem)

    # linear and angular displacement rates
    u1dot, θ1dot = point_displacement_rates(dx, assembly.start[ielem], indices.icol_point, pc)
    u2dot, θ2dot = point_displacement_rates(dx, assembly.stop[ielem], indices.icol_point, pc)

    udot = (u1dot + u2dot)/2
    θdot = (θ1dot + θ2dot)/2

    # linear and angular velocity rates
    Vdot, Ωdot = expanded_element_velocities(dx, ielem, indices.icol_elem)

    # rotate linear and angular velocities into the body frame
    C = get_C(θ)
    Cab = assembly.elements[ielem].Cab
    CtCab = C'*Cab
    V = CtCab*V
    Ω = CtCab*Ω

    # rotate linear and angular velocity rates into the body frame
    CtCabdot = tilde(C'*get_Q(θ)*θdot)*CtCab
    Vdot = CtCabdot*V + CtCab*Vdot
    Ωdot = CtCabdot*Ω + CtCab*Ωdot

    # convert rotation parameters to Wiener-Milenkovic parameters
    k = rotation_parameter_scaling(θ)
    kdot = rotation_parameter_scaling_θ(k, θ)'*θdot

    θ = k*θ
    θdot = kdot*θ + k*θdot

    # promote all variables to the same type
    u, udot, θ, θdot, V, Vdot, Ω, Ωdot, F, M = promote(u, udot, θ, θdot, V, Vdot, Ω, Ωdot, F, M)

    return ElementState(u, udot, θ, θdot, V, Vdot, Ω, Ωdot, F, M)
end

"""
    extract_element_states(system, assembly, x = system.x, dx = system.dx;
        prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}())

Return the state variables corresponding to each element (see [`ElementState`](@ref))
given the solution vector `x`.
"""
extract_element_states

function extract_element_states(system::StaticSystem, assembly; kwargs...)

    return extract_element_states(system.x, system, assembly; kwargs...)
end

function extract_element_states(system::Union{DynamicSystem,ExpandedSystem}, assembly; kwargs...)

    return extract_element_states(system.dx, system.x, system, assembly; kwargs...)
end

function extract_element_states(x, system::Union{DynamicSystem,ExpandedSystem}, assembly; kwargs...)

    dx = FillArrays.Fill(NaN, length(x))

    return extract_element_states(dx, x, system, assembly; kwargs...)
end

function extract_element_states(x, system::StaticSystem, assembly; kwargs...)

    return [extract_element_state(x, system, assembly, ielem; kwargs...) for ielem in eachindex(assembly.elements)]
end

function extract_element_states(dx, x, system::Union{DynamicSystem,ExpandedSystem}, assembly; kwargs...)

    return [extract_element_state(dx, x, system, assembly, ielem; kwargs...) for ielem in eachindex(assembly.elements)]
end

"""
    set_state!([x=system.x,] system, assembly, state; kwargs...)

Set the state variables in `x` to the values in `state`
"""
function set_state!(system::AbstractSystem, assembly::Assembly, state::AssemblyState; kwargs...)
    return set_state!(system.x, system, assembly, state; kwargs...)
end

# static system
function set_state!(x, system::StaticSystem, assembly::Assembly, state::AssemblyState;
    prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}(), kwargs...)

    for ipoint in eachindex(state.points)
        set_linear_displacement!(x, system, prescribed_conditions, state.points[ipoint].u, ipoint)
        set_angular_displacement!(x, system, prescribed_conditions, state.points[ipoint].theta, ipoint)
        set_external_forces!(x, system, prescribed_conditions, state.points[ipoint].F, ipoint)
        set_external_moments!(x, system, prescribed_conditions, state.points[ipoint].M, ipoint)
    end

    for ielem in eachindex(state.elements)
        set_internal_forces!(x, system, state.elements[ielem].Fi, ielem)
        set_internal_moments!(x, system, state.elements[ielem].Mi, ielem)
    end

    return x
end

# dynamic system
function set_state!(x, system::DynamicSystem, assembly::Assembly, state::AssemblyState;
    prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}(), kwargs...)

    for ipoint in eachindex(state.points)
        set_linear_displacement!(x, system, prescribed_conditions, state.points[ipoint].u, ipoint)
        set_angular_displacement!(x, system, prescribed_conditions, state.points[ipoint].theta, ipoint)
        set_linear_velocity!(x, system, state.points[ipoint].V, ipoint)
        set_angular_velocity!(x, system, state.points[ipoint].Omega, ipoint)
        set_external_forces!(x, system, prescribed_conditions, state.points[ipoint].F, ipoint)
        set_external_moments!(x, system, prescribed_conditions, state.points[ipoint].M, ipoint)
    end

    for ielem in eachindex(state.elements)
        set_internal_forces!(x, system, state.elements[ielem].Fi, ielem)
        set_internal_moments!(x, system, state.elements[ielem].Mi, ielem)
    end

    return x
end

# expanded system (need to perform some computations for this one)
function set_state!(x, system::ExpandedSystem, assembly::Assembly, state::AssemblyState;
    prescribed_conditions=Dict{Int,PrescribedConditions{Float64}}(),
    distributed_loads=Dict{Int,DistributedLoads{Float64}}(),
    point_masses=Dict{Int,PointMass{Float64}}(),
    linear_velocity=(@SVector zeros(3)),
    angular_velocity=(@SVector zeros(3)),
    gravity=(@SVector zeros(3)),
    structural_damping=true,
    )

    # update point state variables
    for ipoint in eachindex(state.points)
        # extract states
        u = state.points[ipoint].u
        θ = state.points[ipoint].theta
        V = state.points[ipoint].V
        Ω = state.points[ipoint].Omega
        F = state.points[ipoint].F
        M = state.points[ipoint].M
        # rotate linear and angular velocities into the local frame
        C = get_C(θ)
        V = C*V
        Ω = C*Ω
        # set states
        set_linear_displacement!(x, system, prescribed_conditions, u, ipoint)
        set_angular_displacement!(x, system, prescribed_conditions, θ, ipoint)
        set_linear_velocity!(x, system, V, ipoint)
        set_angular_velocity!(x, system, Ω, ipoint)
        set_external_forces!(x, system, prescribed_conditions, F, ipoint)
        set_external_moments!(x, system, prescribed_conditions, M, ipoint)
    end

    # update element state variables
    for ielem in eachindex(state.elements)
        # extract states
        u = state.elements[ielem].u
        θ = state.elements[ielem].theta
        Fi = state.elements[ielem].Fi
        Mi = state.elements[ielem].Mi
        V = state.elements[ielem].V
        Ω = state.elements[ielem].Omega
        # rotate linear and angular velocities into the local frame
        CtCab = get_C(θ)'*assembly.elements[ielem].Cab
        V = CtCab'*V
        Ω = CtCab'*Ω
        # set states
        set_start_forces!(x, system, Fi, ielem)
        set_start_moments!(x, system, Mi, ielem)
        set_end_forces!(x, system, Fi, ielem)
        set_end_moments!(x, system, Mi, ielem)
        set_element_linear_velocity!(x, system, V, ielem)
        set_element_angular_velocity!(x, system, Ω, ielem)
    end

    # correct element state variables
    for ielem in eachindex(state.elements)
        # compute steady state element properties
        properties = expanded_steady_element_properties(x, system.indices, system.force_scaling,
            structural_damping, assembly, ielem, prescribed_conditions, gravity,
            linear_velocity, angular_velocity)
        # compute dynamic properties
        @unpack mass11, mass12, mass21, mass22 = properties
        Vdot = state.elements[ielem].Vdot
        Ωdot = state.elements[ielem].Omegadot
        Pdot = mass11*Vdot + mass12*Ωdot
        Hdot = mass21*Vdot + mass22*Ωdot
        # augment steady properties with dynamic properties
        properties = (; properties..., Vdot, Ωdot, Pdot, Hdot)
        # compute equilibrium residuals (internal force gradient)
        rF, rM = expanded_element_equilibrium_residuals(properties, distributed_loads, ielem)
        # apply gradient to internal forces
        F1 = state.elements[ielem].Fi + rF/2
        M1 = state.elements[ielem].Mi + rM/2
        F2 = state.elements[ielem].Fi - rF/2
        M2 = state.elements[ielem].Mi - rM/2
        # correct internal forces
        set_start_forces!(x, system, F1, ielem)
        set_start_moments!(x, system, M1, ielem)
        set_end_forces!(x, system, F2, ielem)
        set_end_moments!(x, system, M2, ielem)
    end

    return x
end

"""
    set_rate!([x=system.dx,] system, assembly, state::AssemblyState; kwargs...)

Set the state variable rates in `dx` to the values in `state`
"""
function set_rate!(system::AbstractSystem, assembly::Assembly, state::AssemblyState; kwargs...)
    return set_rate!(system.x, system, assembly, state; kwargs...)
end

# dynamic system
function set_rate!(dx, system::DynamicSystem, assembly::Assembly, state::AssemblyState;
    prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}(), kwargs...)

    dx .= 0

    for ipoint in eachindex(state.points)
        set_linear_displacement_rate!(dx, system, prescribed_conditions, state.points[ipoint].udot, ipoint)
        set_angular_displacement_rate!(dx, system, prescribed_conditions, state.points[ipoint].thetadot, ipoint)
        set_linear_velocity_rate!(dx, system, state.points[ipoint].Vdot, ipoint)
        set_angular_velocity_rate!(dx, system, state.points[ipoint].Omegadot, ipoint)
    end

    return dx
end

# expanded system
function set_rate!(dx, system::ExpandedSystem, assembly::Assembly, state::AssemblyState;
    prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}(), kwargs...)

    dx .= 0

    for ipoint in eachindex(state.points)
        # extract states
        u = state.points[ipoint].u
        θ = state.points[ipoint].theta
        V = state.points[ipoint].V
        Ω = state.points[ipoint].Omega
        # extract rates
        udot = state.points[ipoint].udot
        θdot = state.points[ipoint].thetadot
        Vdot = state.points[ipoint].Vdot
        Ωdot = state.points[ipoint].Omegadot
        # rotate linear and angular velocity rates into the local frame
        C = get_C(θ)
        Q = get_Q(θ)
        Cdot = -C*tilde(C'*Q*θdot)
        Vdot = Cdot*V + C*Vdot
        Ωdot = Cdot*Ω + C*Ωdot
        # set rates
        set_linear_displacement_rate!(dx, system, prescribed_conditions, udot, ipoint)
        set_angular_displacement_rate!(dx, system, prescribed_conditions, θdot, ipoint)
        set_linear_velocity_rate!(dx, system, Vdot, ipoint)
        set_angular_velocity_rate!(dx, system, Ωdot, ipoint)
    end

    for ielem in eachindex(state.elements)
        # extract states
        θ = state.elements[ielem].theta
        V = state.elements[ielem].V
        Ω = state.elements[ielem].Omega
        # extract rates
        θdot = state.points[ielem].thetadot
        Vdot = state.elements[ielem].Vdot
        Ωdot = state.elements[ielem].Omegadot
        # rotate linear and angular velocity rates into the local frame
        C = get_C(θ)
        Cab = assembly.elements[ielem].Cab
        CtCab = C'*Cab
        CtCabdot = tilde(C'*get_Q(θ)*θdot)*CtCab
        Vdot = CtCabdot'*V + CtCab'*Vdot
        Ωdot = CtCabdot'*Ω + CtCab'*Ωdot
        # set rates
        set_element_linear_velocity_rate!(dx, system, Vdot, ielem)
        set_element_angular_velocity_rate!(dx, system, Ωdot, ielem)
    end

    return dx
end

"""
    set_initial_state!([x=system.x,] system, rate_vars, assembly, state; kwargs...)

Set the state variables in `x` for an initial condition analysis to the values in `state`
"""
function set_initial_state!(system::AbstractSystem, rate_vars, assembly::Assembly, state::AssemblyState; kwargs...)
    return set_initial_state!(system.x, system, rate_vars, assembly, state; kwargs...)
end

function set_initial_state!(x, system::DynamicSystem, rate_vars, assembly::Assembly, state::AssemblyState;
    prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}(), kwargs...)

    for ipoint in eachindex(state.points)
        set_external_forces!(x, system, prescribed_conditions, state.points[ipoint].F, ipoint)
        set_external_moments!(x, system, prescribed_conditions, state.points[ipoint].M, ipoint)
        set_initial_linear_displacement!(x, system, rate_vars, prescribed_conditions, state.points[ipoint].u, ipoint)
        set_initial_angular_displacement!(x, system, rate_vars, prescribed_conditions, state.points[ipoint].theta, ipoint)
        set_initial_linear_displacement_rate!(x, system, state.points[ipoint].udot, ipoint)
        set_initial_angular_displacement_rate!(x, system, state.points[ipoint].θdot, ipoint)
        set_initial_linear_velocity_rates!(x, system, rate_vars, prescribed_conditions, state.points[ipoint].Vdot)
        set_initial_angular_velocity_rates!(x, system, rate_vars, prescribed_conditions, state.points[ipoint].Omegadot)
    end

    for ielem in eachindex(state.elements)
        set_internal_forces!(x, system, state.elements[ielem].Fi, ielem)
        set_internal_moments!(x, system, state.elements[ielem].Mi, ielem)
    end

    return x
end