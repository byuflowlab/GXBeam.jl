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

"""
    AssemblyState(system, assembly, x = system.x, dx=system.dx; prescribed_conditions = Dict())

Post-process the system state given the solution vector `x`.  Return an object
of type `AssemblyState` that defines the state of the assembly for the time step.

If `prescribed_conditions` is not provided, all point state variables are assumed
to be displacements/rotations, rather than their actual identities as used in the
analysis.
"""
function AssemblyState(system::AbstractSystem, assembly, args...; prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}())

    points = extract_point_states(system, assembly, args...; prescribed_conditions)

    elements = extract_element_states(system, assembly, args...; prescribed_conditions)

    return AssemblyState(points, elements)
end

"""
    extract_point_state(system, assembly, ipoint, x = system.x;
        prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}())

Return the state variables corresponding to point `ipoint` (see [`PointState`](@ref))
given the solution vector `x`.

If `prescribed_conditions` is not provided, all point state variables are assumed
to be displacements/rotations, rather than their actual identities as used in the
analysis.
"""
extract_point_state

function extract_point_state(system::StaticSystem, assembly, ipoint, x = system.x;
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

function extract_point_state(system::DynamicSystem, assembly, ipoint, x = system.x, dx=FillArrays.Fill(NaN, length(x));
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

function extract_point_state(system::ExpandedSystem, assembly, ipoint, x = system.x, dx=FillArrays.Fill(NaN, length(x));
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
    extract_point_states(system, assembly, x = system.x, dx=system.dx;
        prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}())

Return the state variables corresponding to each point (see [`PointState`](@ref))
given the solution vector `x`.

If `prescribed_conditions` is not provided, all point state variables are assumed
to be displacements/rotations, rather than their actual identities as used in the
analysis.
"""
function extract_point_states(system, assembly, args...; kwargs...)

    return [extract_point_state(system, assembly, ipoint, args...; kwargs...) for ipoint in eachindex(assembly.points)]
end

"""
    extract_point_states!(points, system, assembly, x = system.x;
        prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}())

Pre-allocated version of [`extract_point_states`](@ref)
"""
function extract_point_states!(points, system, assembly, args...; kwargs...)

    for ipoint in eachindex(points)
        points[ipoint] = extract_point_state(system, assembly, ipoint, args...; kwargs...)
    end

    return points
end

"""
    extract_element_state(system, assembly, ielem, x = system.x;
        prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}())

Return the state variables corresponding to element `ielem` (see [`ElementState`](@ref))
given the solution vector `x`.
"""
extract_element_state

function extract_element_state(system::StaticSystem, assembly, ielem, x = system.x;
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

function extract_element_state(system::DynamicSystem, assembly, ielem, x = system.x, dx = FillArrays.Fill(NaN, length(x));
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

function extract_element_state(system::ExpandedSystem, assembly, ielem, x = system.x, dx=FillArrays.Fill(NaN, length(x));
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

    # rotate linear and angular veocity rates into the body frame
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
function extract_element_states(system, assembly, args...; kwargs...)

    return [extract_element_state(system, assembly, ielem, args...; kwargs...) for ielem in eachindex(assembly.elements)]
end

"""
    extract_element_states!(elements, system, assembly, x = system.x, dx=system.dx;
        prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}())

Pre-allocated version of [`extract_element_states`](@ref)
"""
function extract_element_states!(elements, system, assembly, args...; kwargs...)

    for ielem in eachindex(elements)
        elements[ielem] = extract_element_state(system, assembly, ielem, args...; kwargs...)
    end

    return elements
end

"""
    set_state!([x=system.x,] system, assembly, state::AssemblyState; kwargs...)

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
    structural_damping=true,
    prescribed_conditions=Dict{Int,PrescribedConditions{Float64}}(),
    distributed_loads=Dict{Int,DistributedLoads{Float64}}(),
    point_masses=Dict{Int,PointMass{Float64}}(),
    linear_velocity=(@SVector zeros(3)),
    angular_velocity=(@SVector zeros(3)),
    gravity=(@SVector zeros(3)),
    )

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

    for ielem in eachindex(state.elements)
        # compute dynamic element properties
        properties = dynamic_element_properties(assembly, state, ielem, structural_damping,
            prescribed_conditions, gravity, linear_velocity, angular_velocity)
        # compute dynamic element resultants
        F1, F2, M1, M2 = dynamic_element_resultants(properties, distributed_loads, ielem)
        # unpack element properties
        CtCab, V, Ω = properties.CtCab
        # rotate element resultants into the deformed element frame
        F1 = CtCab'*F1
        F2 = CtCab'*F2
        M1 = CtCab'*M1
        M2 = CtCab'*M2
        # rotate element velocities into the deformed element frame
        V = CtCab'*V
        Ω = CtCab'*Ω
        # set states
        set_start_forces!(x, system, F1, ielem)
        set_start_moments!(x, system, M1, ielem)
        set_end_forces!(x, system, F2, ielem)
        set_end_moments!(x, system, M2, ielem)
        set_element_linear_velocity!(x, system, V, ielem)
        set_element_angular_velocity!(x, system, Ω, ielem)
    end

    return x
end

function dynamic_element_properties(assembly, state, ielem, structural_damping,
    prescribed_conditions, gravity, linear_velocity, angular_velocity)

    # unpack element parameters
    @unpack L, Cab, compliance, mass = assembly.elements[ielem]

    # scale compliance and mass matrices by the element length
    compliance *= L
    mass *= L

    # compliance submatrices
    S11 = compliance[SVector{3}(1:3), SVector{3}(1:3)]
    S12 = compliance[SVector{3}(1:3), SVector{3}(4:6)]
    S21 = compliance[SVector{3}(4:6), SVector{3}(1:3)]
    S22 = compliance[SVector{3}(4:6), SVector{3}(4:6)]

    # mass submatrices
    mass11 = mass[SVector{3}(1:3), SVector{3}(1:3)]
    mass12 = mass[SVector{3}(1:3), SVector{3}(4:6)]
    mass21 = mass[SVector{3}(4:6), SVector{3}(1:3)]
    mass22 = mass[SVector{3}(4:6), SVector{3}(4:6)]

    # linear and angular displacement
    u = state.elements[ielem].u
    θ = state.elements[ielem].θ

    # rotation parameter matrices
    C = get_C(θ)
    CtCab = C'*Cab
    Q = get_Q(θ)

    # forces and moments
    F = state.elements[ielem].F
    M = state.elements[ielem].M

    # strain and curvature
    γ = S11*F + S12*M
    κ = S21*F + S22*M

    # body frame velocity (use prescribed values)
    vb, ωb = SVector{3}(linear_velocity), SVector{3}(angular_velocity)

    # gravitational loads
    gvec = SVector{3}(gravity)

    # linear and angular velocity
    V = state.elements[ielem].V
    Ω = state.elements[ielem].Omega

    # linear and angular momentum
    P = CtCab*mass11*CtCab'*V + CtCab*mass12*CtCab'*Ω
    H = CtCab*mass21*CtCab'*V + CtCab*mass22*CtCab'*Ω

    # linear and angular displacement rates
    udot = state.elements[ielem].udot
    θdot = state.elements[ielem].θdot

    # linear and angular velocity rates
    Vdot = state.elements[ielem].Vdot
    Ωdot = state.elements[ielem].Omegadot

    # linear and angular momentum rates
    CtCabdot = tilde(Ω - ωb)*CtCab

    Pdot = CtCab*mass11*CtCab'*Vdot + CtCab*mass12*CtCab'*Ωdot +
        CtCab*mass11*CtCabdot'*V + CtCab*mass12*CtCabdot'*Ω +
        CtCabdot*mass11*CtCab'*V + CtCabdot*mass12*CtCab'*Ω

    Hdot = CtCab*mass21*CtCab'*Vdot + CtCab*mass22*CtCab'*Ωdot +
        CtCab*mass21*CtCabdot'*V + CtCab*mass22*CtCabdot'*Ω +
        CtCabdot*mass21*CtCab'*V + CtCabdot*mass22*CtCab'*Ω

    if structural_damping

        # damping coefficients
        μ = assembly.elements[ielem].mu

        # damping submatrices
        μ11 = @SMatrix [μ[1] 0 0; 0 μ[2] 0; 0 0 μ[3]]
        μ22 = @SMatrix [μ[4] 0 0; 0 μ[5] 0; 0 0 μ[6]]

        # linear displacement
        u1 = state.points[assembly.start[ielem]].u
        u2 = state.points[assembly.stop[ielem]].u

        # angular displacement
        θ1 = state.points[assembly.start[ielem]].theta
        θ2 = state.points[assembly.stop[ielem]].theta

        # linear displacement rates
        udot1 = state.points[assembly.start[ielem]].udot
        udot2 = state.points[assembly.stop[ielem]].udot

        # angular displacement rates
        θdot1 = state.points[assembly.start[ielem]].thetadot
        θdot2 = state.points[assembly.stop[ielem]].thetadot

        # change in linear and angular displacement
        Δu = u2 - u1
        Δθ = θ2 - θ1

        # change in linear and angular displacement rates
        Δudot = udot2 - udot1
        Δθdot = θdot2 - θdot1

        # ΔQ matrix (see structural damping theory)
        ΔQ = get_ΔQ(θ, Δθ, Q)

        # strain rates
        γdot = -CtCab'*tilde(Ω - ωb)*Δu + CtCab'*Δudot - L*CtCab'*tilde(Ω - ωb)*Cab*e1
        κdot = Cab'*Q*Δθdot + Cab'*ΔQ*θdot

        # adjust strains to account for strain rates
        γ -= μ11*γdot
        κ -= μ22*κdot
    end

    # return element properties
    return (; L, C, Cab, CtCab, mass11, mass12, mass21, mass22, F, M, γ, κ, gvec,
        ωb, V, Ω, P, H, Pdot, Hdot)
end

function set_rate!(dx, system::DynamicSystem, state::AssemblyState, prescribed_conditions)

    dx .= 0

    for ipoint in eachindex(state.points)
        set_linear_displacement_rate!(dx, system, prescribed_conditions, state.points[ipoint].u, ipoint)
    end

    for ipoint in eachindex(state.points)
        set_angular_displacement_rate!(dx, system, prescribed_conditions, state.points[ipoint].theta, ipoint)
    end

    for ipoint in eachindex(state.points)
        set_linear_velocity_rate!(dx, system, state.points[ipoint].V, ipoint)
    end

    for ipoint in eachindex(state.points)
        set_angular_velocity_rate!(dx, system, state.points[ipoint].Omega, ipoint)
    end

    return x
end