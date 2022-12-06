"""
    static_analysis(assembly; kwargs...)

Perform a static analysis for the system of nonlinear beams contained in `assembly`.
Return the resulting system, the post-processed solution, and a convergence flag
indicating whether the iteration procedure converged.

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
 - `gravity = [0,0,0]`: Gravity vector.  If time varying, this input may be provided as a
        function of time.
 - `time = 0.0`: Current time or vector of times corresponding to each step. May be used
        in conjunction with time varying prescribed conditions, distributed loads, and
        body frame motion to gradually increase displacements and loads.

# Control Flag Keyword Arguments
 - `reset_state = true`: Flag indicating whether the system state variables should be
        set to zero prior to performing this analysis.
 - `initial_state = nothing`: Object of type [`AssemblyState`](@ref) which contains the
        initial state variables.  If not provided (or set to `nothing`), then the state
        variables stored in `system` (which default to zeros) will be used as the initial
        state variables.
 - `linear = false`: Flag indicating whether a linear analysis should be performed.
 - `two_dimensional = false`: Flag indicating whether to constrain results to the x-y plane
 - `show_trace = false`: Flag indicating whether to display the solution progress.

 # Linear Analysis Keyword Arguments
 - `update_linearization = false`: Flag indicating whether to update the linearization state
        variables for a linear analysis with the instantaneous state variables. If `false`,
        then the initial set of state variables will be used for the linearization.

 # Nonlinear Analysis Keyword Arguments (see [NLsolve.jl](https://github.com/JuliaNLSolvers/NLsolve.jl))
 - `method = :newton`: Solution method for nonlinear systems of equations
 - `linesearch = LineSearches.BackTracking(maxstep=1e6)`: Line search for solving nonlinear
        systems of equations
 - `ftol = 1e-9`: Tolerance for solving the nonlinear system of equations
 - `iterations = 1000`: Iteration limit when solving nonlinear systems of equations

# Sensitivity Analysis Keyword Arguments
 - `pfunc = (p, t) -> (;)`: Function which returns a named tuple with fields corresponding
        to updated versions of the arguments `assembly`, `prescribed_conditions`,
        `distributed_loads`, `point_masses`, and `gravity`. Only fields contained in the
        resulting named tuple will be overwritten.
 - `p`: Sensitivity parameters, as defined in conjunction with the keyword argument `pfunc`.
        While not necessary, using `pfunc` and `p` to define the arguments to this function
        allows automatic differentiation sensitivities to be computed more efficiently
"""
function static_analysis(assembly; kwargs...)

    system = StaticSystem(assembly)

    return static_analysis!(system, assembly; kwargs..., reset_state=false)
end

"""
    static_analysis!(system, assembly; kwargs...)

Pre-allocated version of [`static_analysis`](@ref).
"""
function static_analysis!(system::StaticSystem, assembly;
    # general keyword arguments
    prescribed_conditions=Dict{Int,PrescribedConditions{Float64}}(),
    distributed_loads=Dict{Int,DistributedLoads{Float64}}(),
    point_masses=Dict{Int,PointMass{Float64}}(),
    gravity=(@SVector zeros(3)),
    time=0.0,
    # control flag keyword arguments
    reset_state=true,
    initial_state=nothing,
    linear=false,
    two_dimensional=false,
    show_trace=false,
    # linear analysis keyword arguments
    update_linearization=false,
    # nonlinear analysis keyword arguments
    method=:newton,
    linesearch=LineSearches.BackTracking(maxstep=1e6),
    ftol=1e-9,
    iterations=1000,
    # sensitivity analysis keyword arguments
    pfunc = (p, t) -> (;),
    p = nothing,
    )

    # reset state, if specified
    if reset_state
        reset_state!(system)
    end

    # unpack force scaling parameter and system indices
    @unpack force_scaling, indices = system

    # initialize convergence flag
    converged = Ref(false)

    # package up keyword arguments corresponding to this analysis
    constants = (;
        # assembly, indices, control flags, parameter function, and current time
        assembly, indices, two_dimensional, force_scaling, pfunc, t=first(time),
        # pointers to the pre-allocated storage and the convergence flag
        x=system.x, resid=system.r, jacob=system.K, converged=converged,
        # default parameters
        prescribed_conditions, distributed_loads, point_masses, gravity,
        # nonlinear analysis keyword arguments
        show_trace, method, linesearch, ftol, iterations
    )

    # construct initial state vector
    x0 = static_state_vector(system, initial_state, p, constants)

    # current state vector is the initial state vector
    x = x0

    # also update the state variables in system.x
    dual_safe_copy!(system.x, x)

    # begin time stepping
    for t in time

        # print the current time
        if show_trace
            println("Solving for t=$t")
        end

        # update the system time
        system.t = t

        # update the stored time
        constants = (; constants..., t)

        # solve for the new set of state variables
        if linear
            if update_linearization
                # use current state vector for the linearization
                x = static_lsolve!(x, p, constants)
            else
                # use the initial state vector for the linearization
                x = static_lsolve!(x0, p, constants)
            end
        else
            x = ImplicitAD.implicit(static_nlsolve!, static_residual!, p, constants; drdy=static_drdy)
        end

        # NOTE: `x`, `r`, `K`, and `converged` are updated in `lsolve!`/`nlsolve!`

    end

    # post process state variables
    state = static_output!(system, x, p, constants)

    # NOTE: `system.x` is updated in `static_output!`

    # return the system, state, and the (de-referenced) convergence flag
    return system, state, converged[]
end

# combine constant and variable parameters for a static analysis
function static_parameters(p, constants)

    # unpack default parameters, parameter function, and current time
    @unpack assembly, prescribed_conditions, distributed_loads, point_masses, gravity,
        pfunc, t = constants

    # overwrite default assembly and parameters (if applicable)
    parameters = pfunc(p, t)
    assembly = get(parameters, :assembly, assembly)
    prescribed_conditions = get(parameters, :prescribed_conditions, prescribed_conditions)
    distributed_loads = get(parameters, :distributed_loads, distributed_loads)
    point_masses = get(parameters, :point_masses, point_masses)
    gravity = get(parameters, :gravity, gravity)

    # get parameters corresponding to this time step
    pcond = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(t)
    dload = typeof(distributed_loads) <: AbstractDict ? distributed_loads : distributed_loads(t)
    pmass = typeof(point_masses) <: AbstractDict ? point_masses : point_masses(t)
    gvec = typeof(gravity) <: AbstractVector ? SVector{3}(gravity) : SVector{3}(gravity(t))

    # return results
    return assembly, pcond, dload, pmass, gvec
end

# returns state vector for a static system corresponding to `state`
function static_state_vector(system, state, p, constants)

    # get initial state vector
    if isnothing(state)
        # initial state vector is equal to the system state vector
        x0 = copy(system.x)
    else
        # initialize new state vector storage
        if isnothing(p)
            x0 = similar(system.x, promote_type(eltype(system), eltype(state)))
        else
            x0 = similar(system.x, promote_type(eltype(system), eltype(state), eltype(p)))
        end
        # combine constants and parameters
        assembly, pcond, dload, pmass, gvec = static_parameters(p, constants)
        # set initial state variables in `x`
        set_state!(x0, system, state, pcond)
    end

    return x0
end

# residual function for a static analysis (in format expected by ImplicitAD)
function static_residual!(resid, x, p, constants)

    # unpack indices and control flags
    @unpack indices, two_dimensional, force_scaling = constants

    # combine constants and parameters
    assembly, pcond, dload, pmass, gvec = static_parameters(p, constants)

    # compute and return the residual
    return static_system_residual!(resid, x, indices, two_dimensional, force_scaling,
        assembly, pcond, dload, pmass, gvec)
end

# jacobian function for a static analysis
function static_jacobian!(jacob, x, p, constants)

    # unpack indices, control flags, and jacobian storage
    @unpack indices, two_dimensional, force_scaling = constants

    # combine constants and parameters
    assembly, pcond, dload, pmass, gvec = static_parameters(p, constants)

    # compute and return the jacobian
    return static_system_jacobian!(jacob, x, indices, two_dimensional, force_scaling,
        assembly, pcond, dload, pmass, gvec)
end

# jacobian function (in format expected by ImplicitAD)
static_drdy(residual, x, p, constants) = static_jacobian!(constants.jacob, x, p, constants)

# linear solve for a static analysis (for use with ImplicitAD)
static_lsolve!(x0, p, constants) = lsolve!(x0, p, constants, static_residual!, static_jacobian!)

# nonlinear solve for a static analysis (for use with ImplicitAD)
static_nlsolve!(p, constants) = nlsolve!(p, constants, static_residual!, static_jacobian!)

# returns post-processed state and rate variable vectors
function static_output!(system, x, p, constants)

    # get new assembly and parameters
    assembly, pcond, dload, pmass, gvec = static_parameters(p, constants)

    # update the state variables in `system`
    dual_safe_copy!(system.x, x)

    # return post-processed state
    return AssemblyState(x, system, assembly; prescribed_conditions=pcond)
end

"""
    steady_state_analysis(assembly; kwargs...)

Perform a steady-state analysis for the system of nonlinear beams contained in
`assembly`.  Return the resulting system and a flag indicating whether the
iteration procedure converged.

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
 - `linear_velocity = zeros(3)`: Prescribed linear velocity of the body frame.
        If time varying, this input may be provided as a function of time.
 - `angular_velocity = zeros(3)`: Prescribed angular velocity of the body frame.
        If time varying, this input may be provided as a function of time.
 - `linear_acceleration = zeros(3)`: Prescribed linear acceleration of the body frame.
        If time varying, this input may be provided as a function of time.
 - `angular_acceleration = zeros(3)`: Prescribed angular acceleration of the body frame.
        If time varying, this input may be provided as a function of time.
 - `gravity = [0,0,0]`: Gravity vector in the body frame.  If time varying, this input
        may be provided as a function of time.
 - `time = 0.0`: Current time or vector of times corresponding to each step. May be used
        in conjunction with time varying prescribed conditions, distributed loads, and
        body frame motion to gradually increase displacements and loads.

# Control Flag Keyword Arguments
 - `reset_state = true`: Flag indicating whether the system state variables should be
        set to zero prior to performing this analysis.
 - `initial_state = nothing`: Object of type [`AssemblyState`](@ref) which contains the
        initial state variables.  If not provided (or set to `nothing`), then the state
        variables stored in `system` (which default to zeros) will be used as the initial
        state variables.
 - `structural_damping = false`: Indicates whether to enable structural damping
 - `linear = false`: Flag indicating whether a linear analysis should be performed.
 - `two_dimensional = false`: Flag indicating whether to constrain results to the x-y plane
 - `show_trace = false`: Flag indicating whether to display the solution progress.

# Linear Analysis Keyword Arguments
 - `update_linearization = false`: Flag indicating whether to update the linearization state
        variables for a linear analysis with the instantaneous state variables. If `false`,
        then the initial set of state variables will be used for the linearization.

# Nonlinear Analysis Keyword Arguments
 - `method = :newton`: Method (as defined in NLsolve) to solve nonlinear system of equations
 - `linesearch = LineSearches.BackTracking(maxstep=1e6)`: Line search used to solve the
        nonlinear system of equations
 - `ftol = 1e-9`: tolerance for solving the nonlinear system of equations
 - `iterations = 1000`: maximum iterations for solving the nonlinear system of equations=

# Sensitivity Analysis Keyword Arguments
 - `pfunc = (p, t) -> (;)`: Function which returns a named tuple with fields corresponding
        to updated versions of the arguments `assembly`, `prescribed_conditions`,
        `distributed_loads`, `point_masses`, `linear_velocity`, `angular_velocity`,
        `linear_acceleration`, `angular_acceleration`, and `gravity`. Only fields contained
        in the resulting named tuple will be overwritten.
 - `p`: Sensitivity parameters, as defined in conjunction with the keyword argument `pfunc`.
        While not necessary, using `pfunc` and `p` to define the arguments to this function
        allows automatic differentiation sensitivities to be computed more efficiently
"""
function steady_state_analysis(assembly; constant_mass_matrix=false, kwargs...)

    if constant_mass_matrix
        system = ExpandedSystem(assembly)
    else
        system = DynamicSystem(assembly)
    end

    return steady_state_analysis!(system, assembly; kwargs..., reset_state=true)
end

"""
    steady_state_analysis!(system, assembly; kwargs...)

Pre-allocated version of [`steady_state_analysis`](@ref).
"""
function steady_state_analysis!(system::Union{DynamicSystem, ExpandedSystem}, assembly;
    # general keyword arguments
    prescribed_conditions=Dict{Int,PrescribedConditions{Float64}}(),
    distributed_loads=Dict{Int,DistributedLoads{Float64}}(),
    point_masses=Dict{Int,PointMass{Float64}}(),
    linear_velocity=(@SVector zeros(3)),
    angular_velocity=(@SVector zeros(3)),
    linear_acceleration=(@SVector zeros(3)),
    angular_acceleration=(@SVector zeros(3)),
    gravity=(@SVector zeros(3)),
    time=0.0,
    # control flag keyword arguments
    reset_state=true,
    initial_state=nothing,
    structural_damping=false,
    linear=false,
    constant_mass_matrix=typeof(system)<:ExpandedSystem,
    two_dimensional=false,
    show_trace=false,
    # linear analysis keyword arguments
    update_linearization=false,
    # nonlinear analysis keyword arguments
    method=:newton,
    linesearch=LineSearches.BackTracking(maxstep=1e6),
    ftol=1e-9,
    iterations=1000,
    # sensitivity analysis keyword arguments
    pfunc = (p, t) -> (;),
    p = nothing,
    )

    # check if provided system is consistent with the provided keyword arguments
    if constant_mass_matrix
        @assert typeof(system) <: ExpandedSystem
    else
        @assert typeof(system) <: DynamicSystem
    end

    # reset state, if specified
    if reset_state
        reset_state!(system)
    end

    # unpack force scaling parameter and system pointers
    @unpack force_scaling, indices = system

    # initialize convergence flag
    converged = Ref(false)

    # package up keyword arguments corresponding to this analysis
    constants = (;
        # assembly, indices, control flags, parameter function, and current time
        assembly, indices, structural_damping, two_dimensional, force_scaling, pfunc, t=first(time),
        # pointers to the pre-allocated storage and the convergence flag
        x=system.x, resid=system.r, jacob=system.K, converged=converged,
        # default parameters
        prescribed_conditions, distributed_loads, point_masses, gravity,
        linear_velocity, angular_velocity, linear_acceleration, angular_acceleration,
        # nonlinear analysis keyword arguments
        show_trace, method, linesearch, ftol, iterations
    )

    # construct initial state vector
    x0 = steady_state_vector(system, initial_state, p, constants)

    # current state vector is the initial state vector
    x = x0

    # also update the state variables in system.x
    dual_safe_copy!(system.x, x)

    # begin time stepping
    for t in time

        # print the current time
        if show_trace
            println("Solving for t=$t")
        end

        # update the system time
        system.t = t

        # update the stored time
        constants = (; constants..., t)

        # solve for the new set of state variables
        if linear
            if constant_mass_matrix
                if update_linearization
                    x = expanded_steady_lsolve!(x, p, constants)
                else
                    x = expanded_steady_lsolve!(x0, p, constants)
                end
            else
                if update_linearization
                    x = steady_lsolve!(x, p, constants)
                else
                    x = steady_lsolve!(x0, p, constants)
                end
            end
        else
            if constant_mass_matrix
                x = ImplicitAD.implicit(expanded_steady_nlsolve!, expanded_steady_residual!, p, constants; drdy=expanded_steady_drdy)
            else
                x = ImplicitAD.implicit(steady_nlsolve!, steady_residual!, p, constants; drdy=steady_drdy)
            end
        end

        # NOTE: `x`, `r`, `K`, and `converged` are updated in lsolve!/nlsolve!

    end

    # post process state
    if constant_mass_matrix
        state = expanded_steady_output!(system, x, p, constants)
    else
        state = steady_output!(system, x, p, constants)
    end

    # NOTE: `system.x` and `system.dx` is updated in `steady_output!`

    # return the system, state, and the (de-referenced) convergence flag
    return system, state, converged[]
end

# combines constant and variable parameters for a steady state analysis
function steady_parameters(p, constants)

    # unpack default parameters, parameter function, and current time
    @unpack assembly, prescribed_conditions, distributed_loads, point_masses, gravity,
        linear_velocity, angular_velocity, linear_acceleration, angular_acceleration,
        pfunc, t = constants

    # overwrite default assembly and parameters (if applicable)
    parameters = pfunc(p, t)
    assembly = get(parameters, :assembly, assembly)
    prescribed_conditions = get(parameters, :prescribed_conditions, prescribed_conditions)
    distributed_loads = get(parameters, :distributed_loads, distributed_loads)
    point_masses = get(parameters, :point_masses, point_masses)
    gravity = get(parameters, :gravity, gravity)
    linear_velocity = get(parameters, :linear_velocity, linear_velocity)
    angular_velocity = get(parameters, :angular_velocity, angular_velocity)
    linear_acceleration = get(parameters, :linear_acceleration, linear_acceleration)
    angular_acceleration = get(parameters, :angular_acceleration, angular_acceleration)

    # get parameters corresponding to this time step
    pcond = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(t)
    dload = typeof(distributed_loads) <: AbstractDict ? distributed_loads : distributed_loads(t)
    pmass = typeof(point_masses) <: AbstractDict ? point_masses : point_masses(t)
    gvec = typeof(gravity) <: AbstractVector ? SVector{3}(gravity) : SVector{3}(gravity(t))
    vb_p = typeof(linear_velocity) <: AbstractVector ? SVector{3}(linear_velocity) : SVector{3}(linear_velocity(t))
    ωb_p = typeof(angular_velocity) <: AbstractVector ? SVector{3}(angular_velocity) : SVector{3}(angular_velocity(t))
    ab_p = typeof(linear_acceleration) <: AbstractVector ? SVector{3}(linear_acceleration) : SVector{3}(linear_acceleration(t))
    αb_p = typeof(angular_acceleration) <: AbstractVector ? SVector{3}(angular_acceleration) : SVector{3}(angular_acceleration(t))

    # update body acceleration frame indices
    update_body_acceleration_indices!(constants.indices, pcond)

    return assembly, pcond, dload, pmass, gvec, vb_p, ωb_p, ab_p, αb_p
end

# returns state vector for a steady system corresponding to `state`
function steady_state_vector(system, state, p, constants)

    # get initial state vector
    if isnothing(state)
        # initial state vector is equal to the system state vector
        x0 = copy(system.x)
    else
        # initialize new state vector storage
        if isnothing(p)
            x0 = similar(system.x, promote_type(eltype(system), eltype(state)))
        else
            x0 = similar(system.x, promote_type(eltype(system), eltype(state), eltype(p)))
        end
        # combine constants and parameters
        assembly, pcond, dload, pmass, gvec, vb_p, ωb_p, ab_p, αb_p = steady_parameters(p, constants)
        # set initial state variables in `x`
        set_state!(x0, system, assembly, state; prescribed_conditions=pcond)
    end

    return x0
end

# residual function for a steady analysis (in format expected by ImplicitAD)
function steady_residual!(resid, x, p, constants)

    # unpack indices and control flags
    @unpack indices, structural_damping, two_dimensional, force_scaling, pfunc = constants

    # combine constants and parameters
    assembly, pcond, dload, pmass, gvec, vb_p, ωb_p, ab_p, αb_p = steady_parameters(p, constants)

    # compute and return the residual
    return steady_system_residual!(resid, x, indices, two_dimensional, force_scaling,
        structural_damping, assembly, pcond, dload, pmass, gvec, vb_p, ωb_p, ab_p, αb_p)
end

# jacobian function for a steady analysis
function steady_jacobian!(jacob, x, p, constants)

    # unpack indices, control flags, and jacobian storage
    @unpack indices, structural_damping, two_dimensional, force_scaling = constants

    # combine constants and parameters
    assembly, pcond, dload, pmass, gvec, vb_p, ωb_p, ab_p, αb_p = steady_parameters(p, constants)

    # compute and return the jacobian
    return steady_system_jacobian!(jacob, x, indices, two_dimensional, force_scaling,
        structural_damping, assembly, pcond, dload, pmass, gvec, vb_p, ωb_p, ab_p, αb_p)
end

# jacobian function (in format expected by ImplicitAD)
steady_drdy(residual, x, p, constants) = steady_jacobian!(constants.jacob, x, p, constants)

# defines the linear solver (for use with ImplicitAD)
steady_lsolve!(x0, p, constants) = lsolve!(x0, p, constants, steady_residual!, steady_jacobian!)

# defines the nonlinear solver (for use with ImplicitAD)
steady_nlsolve!(p, constants) = nlsolve!(p, constants, steady_residual!, steady_jacobian!)

# returns post-processed state and rate variable vectors
function steady_output!(system, x, p, constants)

    # unpack indices
    @unpack indices = constants

    # combine constants and parameters
    assembly, pcond, dload, pmass, gvec, vb_p, ωb_p, ab_p, αb_p = steady_parameters(p, constants)

    # initialize state rate vector (if necessary)
    dx = typeof(system.dx) <: typeof(x) ? system.dx .= 0 : zero(x)

    # extract body frame accelerations
    ab, αb = body_accelerations(x, indices.icol_body, ab_p, αb_p)

    # define differentiable variables in the state rate vector
    for ipoint in eachindex(assembly.points)
        # compute state rates
        Δx = assembly.points[ipoint]
        u, _ = point_displacement(x, ipoint, indices.icol_point, pcond)
        udot = @SVector zeros(3)
        θdot = @SVector zeros(3)
        Vdot = ab + cross(αb, Δx) + cross(αb, u)
        Ωdot = αb
        # insert result into the state rate vector
        icol = indices.icol_point[ipoint]
        udot_udot, θdot_θdot = point_displacement_jacobians(ipoint, pcond)
        dx[icol:icol+2] = udot_udot*udot
        dx[icol+3:icol+5] = θdot_θdot*θdot
        dx[icol+6:icol+8] = Vdot
        dx[icol+9:icol+11] = Ωdot
    end

    # update the state and rate variables in `system`
    dual_safe_copy!(system.x, x)
    dual_safe_copy!(system.dx, dx)

    # return result
    return AssemblyState(dx, x, system, assembly; prescribed_conditions=pcond)
end

# residual function for a steady analysis (in format expected by ImplicitAD)
function expanded_steady_residual!(resid, x, p, constants)

    # unpack indices and control flags
    @unpack indices, structural_damping, two_dimensional, force_scaling = constants

    # combine constants and parameters
    assembly, pcond, dload, pmass, gvec, vb_p, ωb_p, ab_p, αb_p = steady_parameters(p, constants)

    # compute and return the residual
    return expanded_steady_system_residual!(resid, x, indices, two_dimensional, force_scaling,
        structural_damping, assembly, pcond, dload, pmass, gvec, vb_p, ωb_p, ab_p, αb_p)
end

# jacobian function for a steady analysis
function expanded_steady_jacobian!(jacob, x, p, constants)

    # unpack indices, control flags, and jacobian storage
    @unpack indices, structural_damping, two_dimensional, force_scaling = constants

    # combine constants and parameters
    assembly, pcond, dload, pmass, gvec, vb_p, ωb_p, ab_p, αb_p = steady_parameters(p, constants)

    # compute and return the jacobian
    return expanded_steady_system_jacobian!(jacob, x, indices, two_dimensional, force_scaling,
        structural_damping, assembly, pcond, dload, pmass, gvec, vb_p, ωb_p, ab_p, αb_p)
end

# jacobian function (in format expected by ImplicitAD)
expanded_steady_drdy(residual, x, p, constants) = expanded_steady_jacobian!(constants.jacob, x, p, constants)

# defines the linear solver (for use with ImplicitAD)
expanded_steady_lsolve!(x0, p, constants) = lsolve!(x0, p, constants, expanded_steady_residual!, expanded_steady_jacobian!)

# defines the nonlinear solver (for use with ImplicitAD)
expanded_steady_nlsolve!(p, constants) = nlsolve!(p, constants, expanded_steady_residual!, expanded_steady_jacobian!)

# returns post-processed state and rate variable vectors
function expanded_steady_output!(system, x, p, constants)

    # unpack indices
    @unpack indices = constants

    # combine constants and parameters
    assembly, pcond, dload, pmass, gvec, vb_p, ωb_p, ab_p, αb_p = steady_parameters(p, constants)

    # initialize state rate vector (if necessary)
    dx = typeof(system.dx) <: typeof(x) ? system.dx .= 0 : zero(x)

    # extract body frame accelerations
    ab, αb = body_accelerations(x, indices.icol_body, ab_p, αb_p)

    # define differentiable variables for each point in the state rate vector
    for ipoint in eachindex(assembly.points)
        # compute state rates for point
        Δx = assembly.points[ipoint]
        u, _ = point_displacement(x, ipoint, indices.icol_point, pcond)
        udot = @SVector zeros(3)
        θdot = @SVector zeros(3)
        Vdot = ab + cross(αb, Δx) + cross(αb, u)
        Ωdot = αb
        # insert into state rate vector
        icol = indices.icol_point[ipoint]
        udot_udot, θdot_θdot = point_displacement_jacobians(ipoint, pcond)
        dx[icol:icol+2] = udot_udot*udot
        dx[icol+3:icol+5] = θdot_θdot*θdot
        dx[icol+6:icol+8] = Vdot
        dx[icol+9:icol+11] = Ωdot
    end

    # define differentiable variables for each element in the state rate vector
    for ielem in eachindex(assembly.elements)
        # compute state rates for element
        Δx = assembly.elements[ielem].x
        u1, _ = point_displacement(x, assembly.start[ielem], indices.icol_point, pcond)
        u2, _ = point_displacement(x, assembly.stop[ielem], indices.icol_point, pcond)
        u = (u1 + u2)/2
        Vdot = ab + cross(αb, Δx) + cross(αb, u)
        Ωdot = αb
        # insert into state rate vector
        icol = indices.icol_elem[ielem]
        dx[icol+12:icol+14] = Vdot
        dx[icol+15:icol+17] = Ωdot
    end

    # update the state and rate variables in `system`
    dual_safe_copy!(system.x, x)
    dual_safe_copy!(system.dx, dx)

    # return result
    return AssemblyState(dx, x, system, assembly; prescribed_conditions=pcond)
end

"""
    linearize!(system, assembly; kwargs...)

Return the state variables, jacobian matrix, and mass matrix of a linearized system using
the current system state vector.

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
 - `linear_velocity = zeros(3)`: Prescribed linear velocity of the body frame.
        If time varying, this input may be provided as a function of time.
 - `angular_velocity = zeros(3)`: Prescribed angular velocity of the body frame.
        If time varying, this input may be provided as a function of time.
 - `linear_acceleration = zeros(3)`: Prescribed linear acceleration of the body frame.
        If time varying, this input may be provided as a function of time.
 - `angular_acceleration = zeros(3)`: Prescribed angular acceleration of the body frame.
        If time varying, this input may be provided as a function of time.
 - `gravity = [0,0,0]`: Gravity vector in the inertial frame.  If time varying, this input
        may be provided as a function of time.
 - `time = 0.0`: Current time or vector of times corresponding to each step. May be used
        in conjunction with time varying prescribed conditions, distributed loads, and
        body frame motion to gradually increase displacements and loads.

# Control Flag Keyword Arguments
 - `reset_state = true`: Flag indicating whether the system state variables should be
         set to zero prior to performing this analysis.
 - `initial_state = nothing`: Object of type [`AssemblyState`](@ref) which contains the
         initial state variables.  If not provided (or set to `nothing`), then the state
         variables stored in `system` (which default to zeros) will be used as the initial
         state variables.
 - `structural_damping = false`: Indicates whether to enable structural damping
 - `linear = false`: Flag indicating whether a linear analysis should be performed.
 - `two_dimensional = false`: Flag indicating whether to constrain results to the x-y plane
 - `show_trace = false`: Flag indicating whether to display the solution progress.

 # Sensitivity Analysis Keyword Arguments
 - `pfunc = (p, t) -> (;)`: Function which returns a named tuple with fields corresponding
        to updated versions of the arguments `assembly`, `prescribed_conditions`,
        `distributed_loads`, `point_masses`, `linear_velocity`, `angular_velocity`,
        `linear_acceleration`, `angular_acceleration`, and `gravity`. Only fields contained
        in the resulting named tuple will be overwritten.
 - `p`: Sensitivity parameters, as defined in conjunction with the keyword argument `pfunc`.
        While not necessary, using `pfunc` and `p` to define the arguments to this function
        allows automatic differentiation sensitivities to be computed more efficiently
"""
function linearize!(system, assembly;
    # general keyword arguments
    prescribed_conditions=Dict{Int,PrescribedConditions{Float64}}(),
    distributed_loads=Dict{Int,DistributedLoads{Float64}}(),
    point_masses=Dict{Int,PointMass{Float64}}(),
    linear_velocity=(@SVector zeros(3)),
    angular_velocity=(@SVector zeros(3)),
    linear_acceleration=(@SVector zeros(3)),
    angular_acceleration=(@SVector zeros(3)),
    gravity=(@SVector zeros(3)),
    time=0.0,
    # control flag keyword arguments
    reset_state=false,
    initial_state=nothing,
    structural_damping=false,
    constant_mass_matrix=typeof(system) <: ExpandedSystem,
    two_dimensional=false,
    show_trace=false,
    # sensitivity analysis keyword arguments
    pfunc = (p, t) -> (;),
    p = nothing,
    )

    # check if provided system is consistent with provided keyword arguments
    if constant_mass_matrix
        @assert typeof(system) <: ExpandedSystem
    else
        @assert typeof(system) <: DynamicSystem
    end

    # reset state, if specified
    if reset_state
        reset_state!(system)
    end

    # unpack internal storage, force scaling parameter, and system pointers
    @unpack x, K, M, force_scaling, indices = system

    # package up keyword arguments corresponding to this analysis
    constants = (;
        # assembly, indices, control flags, parameter function, and current time
        assembly, indices, structural_damping, two_dimensional, force_scaling, pfunc, t=first(time),
        # default parameters
        prescribed_conditions, distributed_loads, point_masses, gravity,
        linear_velocity, angular_velocity, linear_acceleration, angular_acceleration,
    )

    # construct initial state vector
    x0 = steady_state_vector(system, initial_state, p, constants)

    # current state vector is the initial state vector
    x = x0

    # also update the state variables in system.x
    dual_safe_copy!(system.x, x)

    # extract the current time
    t = first(time)

    # print the current time
    if show_trace
        println("Linearizing at t=$t")
    end

    # update the system time
    system.t = t

    # update the stored time
    constants = (; constants..., t)

    # initialize storage for jacobian and mass matrices
    K = similar(system.K, eltype(x))
    M = similar(system.M, eltype(x))

    # solve for the system stiffness and mass matrices
    if constant_mass_matrix
        expanded_steady_jacobian!(K, x, p, constants)
        expanded_mass_matrix!(M, p, constants)
    else
        steady_jacobian!(K, x, p, constants)
        mass_matrix!(M, x, p, constants)
    end

    # return state and jacobian matrices
    return x, K, M
end

"""
    solve_eigensystem(x, K, M, nev)

Return the eigenvalues and eigenvectors of a linearized system.
"""
function solve_eigensystem(x, K, M, nev)

    # construct linear map
    T = eltype(x)
    nx = length(x)
    Kfact = lu(K)
    f! = (b, x) -> b .= ImplicitAD.implicit_linear(K,  M * x; Af=Kfact)
    fc! = (b, x) -> mul!(b, M', ImplicitAD.implicit_linear(K', x; Af=Kfact'))
    A = LinearMap{T}(f!, fc!, nx, nx; ismutating=true)

    # compute eigenvalues and eigenvectors
    λ, V = partialeigen(partialschur(A; nev=min(nx, nev), which=LM())[1])

    # sort eigenvalues by magnitude
    perm = sortperm(λ, by=(λ) -> (abs(λ), imag(λ)), rev=true)
    λ .= λ[perm]
    V .= V[:,perm]

    # eigenvalues are actually -1/λ, no modification necessary for eigenvectors
    λ .= -1 ./ λ

    return λ, V
end

"""
    left_eigenvectors(system, λ, V)
    left_eigenvectors(K, M, λ, V)

Compute the left eigenvector matrix `U` for the `system` using inverse power
iteration given the eigenvalues `λ` and the corresponding right eigenvector
matrix `V`.

The complex conjugate of each left eigenvector is stored in each row of the
matrix `U`

Left and right eigenvectors satisfy the following M-orthogonality condition:
 - u'*M*v = 1 if u and v correspond to the same eigenvalue
 - u'*M*v = 0 if u and v correspond to different eigenvalues
This means that U*M*V = I

This function assumes that `system` has not been modified since the eigenvalues
and right eigenvectors were computed.
"""
left_eigenvectors(system, λ, V) = left_eigenvectors(system.K, system.M, λ, V)

function left_eigenvectors(K, M, λ, V)

    # problem type and dimensions
    TC = eltype(V)
    nx = size(V,1)
    nev = size(V,2)

    # allocate storage
    U = rand(TC, nev, nx)
    u = Vector{TC}(undef, nx)
    tmp = Vector{TC}(undef, nx)

    # compute eigenvectors for each eigenvalue
    for iλ = 1:nev

        # factorize (K + λ*M)'
        KmλMfact = factorize(K' + λ[iλ]'*M')

        # initialize left eigenvector
        for i = 1:nx
            u[i] = U[iλ,i]
        end

        # perform a few iterations to converge the left eigenvector
        for ipass = 1:3
            # get updated u
            mul!(tmp, M, u)
            ldiv!(u, KmλMfact, tmp)
            # normalize u
            unorm = zero(TC)
            for i in axes(M, 1), j in axes(M, 2)
                unorm += conj(u[i])*M[i,j]*V[j,iλ]
            end
            rdiv!(u, conj(unorm))
        end

        # store conjugate of final eigenvector
        for i = 1:nx
            U[iλ,i] = conj(u[i])
        end
    end

    return U
end

function left_eigenvectors(K, M::SparseMatrixCSC, λ, V)

    # problem type and dimensions
    TC = eltype(V)
    nx = size(V,1)
    nev = size(V,2)

    # allocate storage
    U = rand(TC, nev, nx)
    u = Vector{TC}(undef, nx)
    tmp = Vector{TC}(undef, nx)

    # get entries in M
    iM, jM, valM = findnz(M)

    # compute eigenvectors for each eigenvalue
    for iλ = 1:nev

        # factorize (K + λ*M)'
        KmλMfact = factorize(K' + λ[iλ]'*M')

        # initialize left eigenvector
        for i = 1:nx
            u[i] = U[iλ,i]
        end

        # perform a few iterations to converge the left eigenvector
        for ipass = 1:3
            # get updated u
            mul!(tmp, M, u)
            ldiv!(u, KmλMfact, tmp)
            # normalize u
            unorm = zero(TC)
            for k in eachindex(valM)
                unorm += conj(u[iM[k]])*valM[k]*V[jM[k],iλ]
            end
            rdiv!(u, conj(unorm))
        end

        # store conjugate of final eigenvector
        for i = 1:nx
            U[iλ,i] = conj(u[i])
        end
    end

    return U
end

"""
    correlate_eigenmodes(C)

Return the permutation and the associated corruption index vector which associates
eigenmodes from the current iteration with those of the previous iteration given
the correlation matrix `C`.

The correlation matrix can take one of the following forms (in order of preference):
 - `C = U_p*M*V`
 - `C = U*M_p*V_p`
 - `C = V_p'*V`
 - `C = V'*V_p`
where `U` is a matrix of conjugated left eigenvectors, `M` is the system mass
matrix, `V` is a matrix of right eigenvectors, and `()_p` indicates a variable
from the previous iteration.

Note that the following two forms of the correlation matrix seem to be significantly
inferior to their counterparts listed above: `C = U*M*V_p` and `C = U_p*M_p*V`.
This is likely due to the way in which the left eigenvector matrix is calculated.

The corruption index is the largest magnitude in a given row of `C`
that was not chosen divided by the magnitude of the chosen eigenmode.  It is most
meaningful when using one of the forms of the correlation matrix that uses left
eigenvectors since correct eigenmodes will have magnitudes close to 1 and
incorrect eigenmodes will have magnitudes close to 0.

If the new mode number is already assigned, the next highest unassigned mode
number is used.  In this case a corruption index higher than 1 will be returned,
otherwise the values of the corruption index will always be bounded by 0 and 1.

See "New Mode Tracking Methods in Aeroelastic Analysis" by Eldred, Vankayya, and
Anderson.
"""
function correlate_eigenmodes(C)

    # get row permutation that puts maximum values on the diagonals
    nev = size(C, 1)
    tmp = Vector{real(eltype(C))}(undef, nev)
    perm = zeros(Int, nev)
    corruption = Vector{real(eltype(C))}(undef, nev)
    for i = 1:nev

        # rank each value in this row
        for j = 1:nev
            tmp[j] = abs(C[i,j])
        end
        ranked_modes = sortperm(tmp, rev=true)

        # choose the best fit that is not yet assigned
        i1 = ranked_modes[findfirst((x) -> !(x in perm), ranked_modes)]
        i2 = i1 == ranked_modes[1] ? ranked_modes[2] : ranked_modes[1]

        # assign best eigenmode fit, create corruption index
        perm[i] = i1
        corruption[i] = tmp[i2]/tmp[i1]
    end

    return perm, corruption
end

"""
    eigenvalue_analysis(assembly; kwargs...)

Compute the eigenvalues and eigenvectors of the system of nonlinear beams
contained in `assembly`.  Return the modified system, eigenvalues, eigenvectors,
and a convergence flag indicating whether the corresponding steady-state analysis
converged.

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
 - `linear_velocity = zeros(3)`: Prescribed linear velocity of the body frame.
        If time varying, this input may be provided as a function of time.
 - `angular_velocity = zeros(3)`: Prescribed angular velocity of the body frame.
        If time varying, this input may be provided as a function of time.
 - `linear_acceleration = zeros(3)`: Prescribed linear acceleration of the body frame.
        If time varying, this input may be provided as a function of time.
 - `angular_acceleration = zeros(3)`: Prescribed angular acceleration of the body frame.
        If time varying, this input may be provided as a function of time.
 - `gravity = [0,0,0]`: Gravity vector in the inertial frame.  If time varying, this input
        may be provided as a function of time.
 - `time = 0.0`: Current time or vector of times corresponding to each step. May be used
        in conjunction with time varying prescribed conditions, distributed loads, and
        body frame motion to gradually increase displacements and loads.

# Control Flag Keyword Arguments
 - `reset_state = true`: Flag indicating whether the system state variables should be
        set to zero prior to performing this analysis.
 - `initial_state = nothing`: Object of type [`AssemblyState`](@ref) which contains the
        initial state variables.  If not provided (or set to `nothing`), then the state
        variables stored in `system` (which default to zeros) will be used as the initial
        state variables.
 - `structural_damping = false`: Indicates whether to enable structural damping
 - `linear = false`: Flag indicating whether a linear analysis should be performed.
 - `two_dimensional = false`: Flag indicating whether to constrain results to the x-y plane
 - `show_trace = false`: Flag indicating whether to display the solution progress.

# Linear Analysis Keyword Arguments
 - `update_linearization = false`: Flag indicating whether to update the linearization state
        variables for a linear analysis with the instantaneous state variables. If `false`,
        then the initial set of state variables will be used for the linearization.

# Nonlinear Analysis Keyword Arguments
 - `method = :newton`: Method (as defined in NLsolve) to solve nonlinear system of equations
 - `linesearch = LineSearches.LineSearches.BackTracking(maxstep=1e6)`: Line search used to
        solve the nonlinear system of equations
 - `ftol = 1e-9`: tolerance for solving the nonlinear system of equations
 - `iterations = 1000`: maximum iterations for solving the nonlinear system of equations

# Sensitivity Analysis Keyword Arguments
- `pfunc = (p, t) -> (;)`: Function which returns a named tuple with fields corresponding
        to updated versions of the arguments `assembly`, `prescribed_conditions`,
        `distributed_loads`, `point_masses`, `linear_velocity`, `angular_velocity`,
        `linear_acceleration`, `angular_acceleration`, and `gravity`. Only fields contained
        in the resulting named tuple will be overwritten.
- `p`: Sensitivity parameters, as defined in conjunction with the keyword argument `pfunc`.
        While not necessary, using `pfunc` and `p` to define the arguments to this function
        allows automatic differentiation sensitivities to be computed more efficiently

# Eigenvalue Analysis Keyword Arguments
 - `nev = 6`: Number of eigenvalues to compute
 - `steady = reset_state && !linear`: Flag indicating whether the steady state
        solution should be found prior to performing the eigenvalue analysis.
 - `left = false`: Flag indicating whether to return left and right eigenvectors rather
        than just right eigenvectors.
 - `Uprev = nothing`: Previous left eigenvector matrix.  May be provided in order to
        reorder eigenvalues based on results from a previous iteration.
"""
function eigenvalue_analysis(assembly; constant_mass_matrix=false, kwargs...)

    if constant_mass_matrix
        system = ExpandedSystem(assembly)
    else
        system = DynamicSystem(assembly)
    end

    return eigenvalue_analysis!(system, assembly; kwargs..., reset_state=true)
end

"""
    eigenvalue_analysis!(system, assembly; kwargs...)

Pre-allocated version of `eigenvalue_analysis`.
"""
function eigenvalue_analysis!(system, assembly;
    # general keyword arguments
    prescribed_conditions=Dict{Int,PrescribedConditions{Float64}}(),
    distributed_loads=Dict{Int,DistributedLoads{Float64}}(),
    point_masses=Dict{Int,PointMass{Float64}}(),
    linear_velocity=(@SVector zeros(3)),
    angular_velocity=(@SVector zeros(3)),
    linear_acceleration=(@SVector zeros(3)),
    angular_acceleration=(@SVector zeros(3)),
    gravity=(@SVector zeros(3)),
    time=0.0,
    # control flag keyword arguments
    reset_state=true,
    initial_state=nothing,
    structural_damping=false,
    linear=false,
    constant_mass_matrix=typeof(system)<:ExpandedSystem,
    two_dimensional=false,
    show_trace=false,
    # linear analysis keyword arguments
    update_linearization=false,
    # nonlinear analysis keyword arguments
    method=:newton,
    linesearch=LineSearches.BackTracking(maxstep=1e6),
    ftol=1e-9,
    iterations=1000,
    # sensitivity analysis keyword arguments
    pfunc = (p, t) -> (;),
    p = nothing,
    # eigenvalue analysis keyword arguments
    nev = 6,
    steady_state=true,
    left=false,
    Uprev=nothing,
    )

    # check if provided system is consistent with provided keyword arguments
    if constant_mass_matrix
        @assert typeof(system) <: ExpandedSystem
    else
        @assert typeof(system) <: DynamicSystem
    end

    # reset state, if specified
    if reset_state
        reset_state!(system)
    end

    # set linearization state variables
    if steady_state
        # print status
        if show_trace
            println("Finding a Steady-State Solution")
        end
        # find the linearization state variables using a steady state analysis
        system, initial_state, converged = steady_state_analysis!(system, assembly;
            # general keyword arguments
            prescribed_conditions=prescribed_conditions,
            distributed_loads=distributed_loads,
            point_masses=point_masses,
            linear_velocity=linear_velocity,
            angular_velocity=angular_velocity,
            linear_acceleration=linear_acceleration,
            angular_acceleration=angular_acceleration,
            gravity=gravity,
            time=time,
            # control flag keyword arguments
            reset_state=reset_state,
            initial_state=initial_state,
            structural_damping=structural_damping,
            linear=linear,
            constant_mass_matrix=constant_mass_matrix,
            two_dimensional=two_dimensional,
            show_trace=show_trace,
            # linear analysis keyword arguments
            update_linearization=update_linearization,
            # nonlinear analysis keyword arguments
            method=method,
            linesearch=linesearch,
            ftol=ftol,
            iterations=iterations,
            # sensitivity analysis keyword arguments
            pfunc = pfunc,
            p = p,
            )
    else
        # use current state variables
        converged = true
    end

    # linearize the resulting system of equations
    x, K, M = linearize!(system, assembly;
        # general keyword arguments
        prescribed_conditions=prescribed_conditions,
        distributed_loads=distributed_loads,
        point_masses=point_masses,
        linear_velocity=linear_velocity,
        angular_velocity=angular_velocity,
        linear_acceleration=linear_acceleration,
        angular_acceleration=angular_acceleration,
        gravity=gravity,
        time=time,
        # control flag keyword arguments
        reset_state=false,
        initial_state=initial_state,
        structural_damping=structural_damping,
        constant_mass_matrix=typeof(system) <: ExpandedSystem,
        two_dimensional=two_dimensional,
        show_trace=show_trace,
        # sensitivity analysis keyword arguments
        pfunc=pfunc,
        p=p,
        )

    # solve the eigensystem
    λ, V = solve_eigensystem(x, K, M, nev)

    # correlate eigenmodes
    if !isnothing(Uprev)
        # construct correlation matrix
        C = Uprev*M*V
        # find correct permutation of eigenmodes
        perm, corruption = correlate_eigenmodes(C)
        # rearrange eigenmodes
        λ = λ[perm]
        V = V[:,perm]
    end

    if left
        # find the left eigenvector corresponding to each right eigenvector
        U = left_eigenvectors(K, M, λ, V)
        # return left and right eigenvectors
        return system, λ, U, V, converged
    end

    # return only right eigenvalues
    return system, λ, V, converged
end

"""
    initial_condition_analysis(assembly, t0; kwargs...)

Perform an analysis to obtain a consistent set of initial conditions.  Return the
resulting system and a flag indicating whether the iteration procedure converged.

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
 - `linear_velocity = zeros(3)`: Initial linear velocity of the body frame.
 - `angular_velocity = zeros(3)`: Initial angular velocity of the body frame.
 - `linear_acceleration = zeros(3)`: Initial linear acceleration of the body frame.
 - `angular_acceleration = zeros(3)`: Initial angular acceleration of the body frame.
 - `gravity = [0,0,0]`: Gravity vector in the inertial frame.

# Control Flag Keyword Arguments
 - `reset_state = true`: Flag indicating whether the system state variables should be
        set to zero prior to performing this analysis.
 - `initial_state = nothing`: Object of type [`AssemblyState`](@ref) which contains the
        initial state variables.  If not provided (or set to `nothing`), then the state
        variables stored in `system` (which default to zeros) will be used as the initial
        state variables.
 - `structural_damping = true`: Indicates whether to enable structural damping
 - `linear = false`: Flag indicating whether a linear analysis should be performed.
 - `two_dimensional = false`: Flag indicating whether to constrain results to the x-y plane
 - `steady_state=false`: Flag indicating whether to initialize by performing a steady state
        analysis.
 - `show_trace = false`: Flag indicating whether to display the solution progress.

 # Initial Condition Analysis Keyword Arguments
 - `u0 = fill(zeros(3), length(assembly.points))`: Initial linear displacement of
        each point **relative to the body frame**
 - `theta0 = fill(zeros(3), length(assembly.points))`: Initial angular displacement of
        each point **relative to the body frame** (using Wiener-Milenkovic Parameters)
 - `V0 = fill(zeros(3), length(assembly.points))`: Initial linear velocity of
        each point **relative to the body frame**
 - `Omega0 = fill(zeros(3), length(assembly.points))`: Initial angular velocity of
        each point **relative to the body frame**
 - `Vdot0 = fill(zeros(3), length(assembly.points))`: Initial linear acceleration of
        each point **relative to the body frame**
 - `Omegadot0 = fill(zeros(3), length(assembly.points))`: Initial angular acceleration of
        each point **relative to the body frame**

 # Linear Analysis Keyword Arguments
 - `update_linearization = false`: Flag indicating whether to update the linearization state
        variables for a linear analysis with the instantaneous state variables. If `false`,
        then the initial set of state variables will be used for the linearization.

# Nonlinear Analysis Keyword Arguments
 - `method = :newton`: Method (as defined in NLsolve) to solve nonlinear system of equations
 - `linesearch = LineSearches.BackTracking(maxstep=1e6)`: Line search used to solve the
        nonlinear system of equations
 - `ftol = 1e-9`: tolerance for solving the nonlinear system of equations
 - `iterations = 1000`: maximum iterations for solving the nonlinear system of equations

# Sensitivity Analysis Keyword Arguments
- `pfunc = (p, t) -> (;)`: Function which returns a named tuple with fields corresponding
        to updated versions of the arguments `assembly`, `prescribed_conditions`,
        `distributed_loads`, `point_masses`, `linear_velocity`, `angular_velocity`,
        `linear_acceleration`, `angular_acceleration`, `gravity`, `u0`, `theta0`, `V0`,
        `Omega0`, `Vdot0`, and `Omegadot0`. Only fields contained in the resulting named
        tuple will be overwritten.
- `p`: Sensitivity parameters, as defined in conjunction with the keyword argument `pfunc`.
        While not necessary, using `pfunc` and `p` to define the arguments to this function
        allows automatic differentiation sensitivities to be computed more efficiently
"""
function initial_condition_analysis(assembly, t0; constant_mass_matrix=false, kwargs...)

    if constant_mass_matrix
        system = ExpandedSystem(assembly)
    else
        system = DynamicSystem(assembly)
    end

    return initial_condition_analysis!(system, assembly, t0; kwargs..., reset_state=true)
end

"""
    initial_condition_analysis!(system, assembly, t0; kwargs...)

Pre-allocated version of `initial_condition_analysis`.
"""
function initial_condition_analysis!(system, assembly, t0;
    # general keyword arguments
    prescribed_conditions=Dict{Int,PrescribedConditions{Float64}}(),
    distributed_loads=Dict{Int,DistributedLoads{Float64}}(),
    point_masses=Dict{Int,PointMass{Float64}}(),
    linear_velocity=(@SVector zeros(3)),
    angular_velocity=(@SVector zeros(3)),
    linear_acceleration=(@SVector zeros(3)),
    angular_acceleration=(@SVector zeros(3)),
    gravity=(@SVector zeros(3)),
    # control flag keyword arguments
    reset_state=true,
    initial_state=nothing,
    structural_damping=true,
    linear=false,
    steady_state=false,
    constant_mass_matrix=typeof(system) <: ExpandedSystem,
    two_dimensional=false,
    show_trace=false,
    # initial condition analysis keyword arguments
    u0=fill((@SVector zeros(3)), length(assembly.points)),
    theta0=fill((@SVector zeros(3)), length(assembly.points)),
    V0=fill((@SVector zeros(3)), length(assembly.points)),
    Omega0=fill((@SVector zeros(3)), length(assembly.points)),
    Vdot0=fill((@SVector zeros(3)), length(assembly.points)),
    Omegadot0=fill((@SVector zeros(3)), length(assembly.points)),
    # linear analysis keyword arguments
    update_linearization=false,
    # nonlinear analysis keyword arguments
    method=:newton,
    linesearch=LineSearches.BackTracking(maxstep=1e6),
    ftol=1e-9,
    iterations=1000,
    # sensitivity analysis keyword arguments
    pfunc = (p, t) -> (;),
    p = nothing,
    )

    # perform steady state analysis instead (if requested)
    if steady_state
        return steady_state_analysis!(system, assembly;
            # general keyword arguments
            prescribed_conditions=prescribed_conditions,
            distributed_loads=distributed_loads,
            point_masses=point_masses,
            linear_velocity=linear_velocity,
            angular_velocity=angular_velocity,
            linear_acceleration=linear_acceleration,
            angular_acceleration=angular_acceleration,
            gravity=gravity,
            time=t0,
            # control flag keyword arguments
            reset_state=reset_state,
            initial_state=initial_state,
            structural_damping=structural_damping,
            linear=linear,
            constant_mass_matrix=constant_mass_matrix,
            two_dimensional=two_dimensional,
            show_trace=show_trace,
            # linear analysis keyword arguments
            update_linearization=update_linearization,
            # nonlinear analysis keyword arguments
            method=method,
            linesearch=linesearch,
            ftol=ftol,
            iterations=iterations,
            # sensitivity analysis keyword arguments
            pfunc = pfunc,
            p = p,
            )
    end

    # check if provided system is consistent with the provided keyword arguments
    if constant_mass_matrix
        @assert typeof(system) <: ExpandedSystem
    else
        @assert typeof(system) <: DynamicSystem
    end

    # reset state, if specified
    if reset_state
        reset_state!(system)
    end

    # save provided system
    original_system = system

    # construct dynamic system for this analysis
    if constant_mass_matrix
        system = DynamicSystem(assembly; force_scaling = system.force_scaling)
    end

    # unpack force scaling parameter and system pointers
    @unpack force_scaling, indices = system

    # initialize convergence flag
    converged = Ref(false)

    # package up keyword arguments corresponding to this analysis
    constants = (;
        # assembly, indices, control flags, parameter function, and current time
        assembly, indices, structural_damping, two_dimensional, force_scaling, pfunc, t=first(t0),
        # pointers to the pre-allocated storage and the convergence flag
        x=system.x, resid=system.r, jacob=system.K, converged=converged,
        # default parameters
        prescribed_conditions, distributed_loads, point_masses, gravity,
        linear_velocity, angular_velocity, linear_acceleration, angular_acceleration,
        u0, theta0, V0, Omega0, Vdot0, Omegadot0,
        # nonlinear analysis keyword arguments
        show_trace, method, linesearch, ftol, iterations
    )

    # solve for initial time
    t = first(t0)

    # print the current time
    if show_trace
        println("Solving for t=$t")
    end

    # update the system time
    system.t = t

    # update the stored time
    constants = (; constants..., t)

    # get new assembly and parameters
    assembly, pcond, dload, pmass, gvec, vb_p, ωb_p, ab_p, αb_p, u0, theta0, V0, Omega0,
        Vdot0, Omegadot0 = initial_parameters(p, constants)

    # extract initial state from the original system (if not provided)
    if isnothing(initial_state) && !reset_state
        initial_state = AssemblyState(original_system, assembly; prescribed_conditions=pcond)
    end

    # get floating point type
    TF = isnothing(p) ? eltype(system) : promote_type(eltype(system), eltype(p))

    # --- Determine whether Vdot and Ωdot may be found using the equilibrium equations --- #

    # NOTE: Our system of equations cannot be solved for Vdot and Ωdot if the corresponding
    # rows and columns of the mass matrix are zero.

    # Fortunately, though we cannot solve for these variables, their values are not used.

    # construct `rate_vars` vector
    M1 = similar(system.M, TF)
    M1 = system_mass_matrix!(M1, FillArrays.Zeros(system.x), indices, two_dimensional,
        force_scaling, assembly, pcond, pmass)
    rate_vars1 = .!(iszero.(sum(M1, dims=1)))

    # save results
    constants = (; constants..., rate_vars1)

    # --- Determine whether Fi and Mi may be found using the compatability equations --- #

    # NOTE: The compatability equations cannot be solved for Fi and Mi if the corresponding
    # rows and columns of the compliance matrix are zero.

    # Values for these variables must then be derived from the equilibrium equations.
    # Since the equilibrium equations cannot then be used to find Vdot and Ωdot,
    # compatible values for Vdot and Ωdot must be provided.

    # replace mass matrix with compliance matrix
    elements = copy(assembly.elements)
    for (i, e) in enumerate(assembly.elements)
        assembly.elements[i] = Element(e.L, e.x, e.compliance, e.compliance, e.Cab, e.mu)
    end

    # temporarily ignore the mass of all point masses
    pmass2 = Dict{Int,PointMass{Float64}}()

    # construct `rate_vars2` vector to test if the second case applies
    M2 = similar(system.M, TF)
    M2 = system_mass_matrix!(M2, FillArrays.Zeros(system.x), indices, two_dimensional,
        force_scaling, assembly, pcond, pmass2)
    rate_vars2 = .!(iszero.(sum(M2, dims=1)))

    # restore original element mass matrices
    for i in eachindex(assembly.elements)
        assembly.elements[i] = elements[i]
    end

    # save results
    constants = (; constants..., rate_vars2)

    # ---------------------------------------------------------------------------------- #

    # construct initial state vector
    x0 = initial_state_vector(system, initial_state, p, constants)

    # current state vector is the initial state vector
    x = x0

    # solve for the new set of state variables
    if linear
        if update_linearization
            x = initial_lsolve!(x, p, constants)
        else
            x = initial_lsolve!(x0, p, constants)
        end
    else
        x = ImplicitAD.implicit(initial_nlsolve!, initial_residual!, p, constants; drdy=initial_drdy)
    end

    # NOTE: `x`, `r`, `K`, and `converged` are updated in lsolve!/nlsolve!

    # postprocess the dual number state
    state = initial_output(system, x, p, constants)

    # copy the state and rate variables to the original system
    if eltype(original_system) <: eltype(state)
        # set the values of the original system directly
        set_state!(original_system, assembly, state; prescribed_conditions=pcond)
        set_rate!(original_system, assembly, state; prescribed_conditions=pcond)
    else
        # initialize storage
        dx = similar(original_system.dx, output_type)
        x = similar(original_system.x, output_type)
        # construct state and rate vectors
        set_rate!(dx, original_system, assembly, state; prescribed_conditions=pcond)
        set_state!(x, original_system, assembly, state; prescribed_conditions=pcond)
        # copy only the primal portion of the variables
        dual_safe_copy!(system.dx, dx)
        dual_safe_copy!(system.x, x)
    end

    # return the system, state, and the (de-referenced) convergence flag
    return original_system, state, converged[]
end

function initial_parameters(p, constants)

    # unpack default parameters, parameter function, and current time
    @unpack assembly, prescribed_conditions, distributed_loads, point_masses, gravity,
        linear_velocity, angular_velocity, linear_acceleration, angular_acceleration,
        pfunc, t = constants

    # also unpack initial conditions
    @unpack u0, theta0, V0, Omega0, Vdot0, Omegadot0 = constants

    # overwrite default assembly and parameters (if applicable)
    parameters = pfunc(p, t)
    assembly = get(parameters, :assembly, assembly)
    prescribed_conditions = get(parameters, :prescribed_conditions, prescribed_conditions)
    distributed_loads = get(parameters, :distributed_loads, distributed_loads)
    point_masses = get(parameters, :point_masses, point_masses)
    gravity = get(parameters, :gravity, gravity)
    linear_velocity = get(parameters, :linear_velocity, linear_velocity)
    angular_velocity = get(parameters, :angular_velocity, angular_velocity)
    linear_acceleration = get(parameters, :linear_acceleration, linear_acceleration)
    angular_acceleration = get(parameters, :angular_acceleration, angular_acceleration)

    # overwrite default initial conditions (if applicable)
    u0 = get(parameters, :u0, u0)
    theta0 = get(parameters, :theta0, theta0)
    V0 = get(parameters, :V0, V0)
    Omega0 = get(parameters, :Omega0, Omega0)
    Vdot0 = get(parameters, :Vdot0, Vdot0)
    Omegadot0 = get(parameters, :Omegadot0, Omegadot0)

    # get parameters corresponding to this time step
    pcond = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(t)
    dload = typeof(distributed_loads) <: AbstractDict ? distributed_loads : distributed_loads(t)
    pmass = typeof(point_masses) <: AbstractDict ? point_masses : point_masses(t)
    gvec = typeof(gravity) <: AbstractVector ? SVector{3}(gravity) : SVector{3}(gravity(t))
    vb_p = typeof(linear_velocity) <: AbstractVector ? SVector{3}(linear_velocity) : SVector{3}(linear_velocity(t))
    ωb_p = typeof(angular_velocity) <: AbstractVector ? SVector{3}(angular_velocity) : SVector{3}(angular_velocity(t))
    ab_p = typeof(linear_acceleration) <: AbstractVector ? SVector{3}(linear_acceleration) : SVector{3}(linear_acceleration(t))
    αb_p = typeof(angular_acceleration) <: AbstractVector ? SVector{3}(angular_acceleration) : SVector{3}(angular_acceleration(t))

    # update body acceleration frame indices
    update_body_acceleration_indices!(constants.indices, pcond)

    return assembly, pcond, dload, pmass, gvec, vb_p, ωb_p, ab_p, αb_p, u0, theta0,
        V0, Omega0, Vdot0, Omegadot0
end

# returns state vector for an initial condition system corresponding to `state`
function initial_state_vector(system, state, p, constants)

    # extract indices of differential variables
    rate_vars = constants.rate_vars2

    # get initial state vector
    if isnothing(state)
        # initial state vector is equal to the system state vector
        x0 = copy(system.x)
    else
        # initialize new state vector storage
        if isnothing(p)
            x0 = similar(system.x, promote_type(eltype(system), eltype(state)))
        else
            x0 = similar(system.x, promote_type(eltype(system), eltype(state), eltype(p)))
        end
        # combine constants and parameters
        assembly, pcond, dload, pmass, gvec, vb_p, ωb_p, ab_p, αb_p,
            u0, theta0, V0, Omega0, Vdot0, Omegadot0 = initial_parameters(p, constants)
        # set initial state variables in `x`
        set_initial_state!(x0, system, rate_vars, state, pcond)
    end

    return x0
end

# defines the initial residual (for use with ImplicitAD)
function initial_residual!(resid, x, p, constants)

    # unpack indices, control flags, and parameter function
    @unpack indices, structural_damping, two_dimensional, force_scaling, pfunc = constants

    @unpack rate_vars1, rate_vars2 = constants

    # combine constants and parameters
    assembly, pcond, dload, pmass, gvec, vb_p, ωb_p, ab_p, αb_p, u0, theta0, V0, Omega0,
        Vdot0, Omegadot0 = initial_parameters(p, constants)

    # update acceleration state variable indices
    update_body_acceleration_indices!(indices, pcond)

    # compute and return the residual
    return initial_system_residual!(resid, x, indices, rate_vars1, rate_vars2, two_dimensional,
        force_scaling, structural_damping, assembly, pcond, dload, pmass,
        gvec, vb_p, ωb_p, ab_p, αb_p, u0, theta0, V0, Omega0, Vdot0, Omegadot0)
end

# defines the initial jacobian (for use with ImplicitAD)
function initial_jacobian!(jacob, x, p, constants)

    # unpack indices, control flags, and parameter function
    @unpack indices, structural_damping, two_dimensional, force_scaling, pfunc = constants

    @unpack rate_vars1, rate_vars2 = constants

    # combine constants and parameters
    assembly, pcond, dload, pmass, gvec, vb_p, ωb_p, ab_p, αb_p, u0, theta0, V0, Omega0,
        Vdot0, Omegadot0 = initial_parameters(p, constants)

    # update acceleration state variable indices
    update_body_acceleration_indices!(indices, pcond)

    # compute and return the jacobian
    return initial_system_jacobian!(jacob, x, indices, rate_vars1, rate_vars2, two_dimensional,
        force_scaling, structural_damping, assembly, pcond, dload, pmass,
        gvec, vb_p, ωb_p, ab_p, αb_p, u0, theta0, V0, Omega0, Vdot0, Omegadot0)
end

# jacobian function (in format expected by ImplicitAD)
initial_drdy(residual, x, p, constants) = initial_jacobian!(constants.jacob, x, p, constants)

# defines the linear solver (for use with ImplicitAD)
initial_lsolve!(x0, p, constants) = lsolve!(x0, p, constants, initial_residual!, initial_jacobian!)

# defines the nonlinear solver (for use with ImplicitAD)
initial_nlsolve!(p, constants) = nlsolve!(p, constants, initial_residual!, initial_jacobian!)

# returns post-processed state and rate variable vectors
function initial_output(system, x, p, constants)

    # unpack indices
    @unpack indices, rate_vars1, rate_vars2, force_scaling = constants

    # combine constants and parameters
    assembly, pcond, dload, pmass, gvec, vb_p, ωb_p, ab_p, αb_p, u0, theta0, V0, Omega0,
        Vdot0, Omegadot0 = initial_parameters(p, constants)

    # initialize state rate vector (if necessary)
    dx = typeof(system.dx) <: typeof(x) ? system.dx : zero(x)

    # body velocities and accelerations
    vb, ωb = vb_p, ωb_p
    ab, αb = body_accelerations(x, indices.icol_body, ab_p, αb_p)

    # populate state and state rate vector
    for ipoint in eachindex(assembly.points)
        # states corresponding to this point
        u, θ = initial_point_displacement(x, ipoint, indices.icol_point, pcond, u0, theta0, rate_vars2)
        V, Ω = V0[ipoint], Omega0[ipoint]
        F, M = point_loads(x, ipoint, indices.icol_point, force_scaling, pcond)
        # rates corresponding to this point
        udot, θdot = initial_point_displacement_rates(x, ipoint, indices.icol_point)
        Vdot, Ωdot = initial_point_velocity_rates(x, ipoint, indices.icol_point, pcond,
            Vdot0, Omegadot0, rate_vars2)
        # adjust velocities and accelerations to account for rigid body motion
        Δx = assembly.points[ipoint]
        V += vb + cross(ωb, Δx + u)
        Ω += ωb
        Vdot += ab + cross(αb, Δx + u) + cross(ωb, udot)
        Ωdot += αb
        # insert result into the state vector
        set_linear_displacement!(x, system, pcond, u, ipoint)
        set_angular_displacement!(x, system, pcond, θ, ipoint)
        set_linear_velocity!(x, system, V, ipoint)
        set_angular_velocity!(x, system, Ω, ipoint)
        set_external_forces!(x, system, pcond, F, ipoint)
        set_external_moments!(x, system, pcond, M, ipoint)
        # insert result into the rate vector
        icol = indices.icol_point[ipoint]
        udot_udot, θdot_θdot = point_displacement_jacobians(ipoint, pcond)
        dx[icol:icol+2] = udot_udot*udot
        dx[icol+3:icol+5] = θdot_θdot*θdot
        dx[icol+6:icol+8] = Vdot
        dx[icol+9:icol+11] = Ωdot
    end

    # update the state and rate variables in `system`
    dual_safe_copy!(system.x, x)
    dual_safe_copy!(system.dx, dx)

    # return result
    state = AssemblyState(dx, x, system, assembly; prescribed_conditions=pcond)

    return state
end

"""
    time_domain_analysis(assembly, tvec; kwargs...)

Perform a time-domain analysis for the system of nonlinear beams contained in
`assembly` using the time vector `tvec`.  Return the final system, a post-processed
solution history, and a convergence flag indicating whether the iteration procedure
converged for every time step.

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
 - `linear_velocity = zeros(3)`: Prescribed linear velocity of the body frame.
 - `angular_velocity = zeros(3)`: Prescribed angular velocity of the body frame.
 - `linear_acceleration = zeros(3)`: Initial linear acceleration of the body frame.
 - `angular_acceleration = zeros(3)`: Initial angular acceleration of the body frame.
 - `gravity = [0,0,0]`: Gravity vector in the body frame.  If time varying, this input
        may be provided as a function of time.

 # Control Flag Keyword Arguments
 - `reset_state = true`: Flag indicating whether the system state variables should be
        set to zero prior to performing this analysis.
 - `initial_state = nothing`: Object of type `AssemblyState`, which defines the initial
        states and state rates corresponding to the analysis.  By default, this input is
        calculated using either `steady_state_analysis` or `initial_condition_analysis`.
 - `steady_state = false`: Flag indicating whether to compute the state variables
        corresponding to the keyword argument `initial_state` using `steady_state_analysis`
        (rather than `initial_condition_analysis`).
 - `structural_damping = true`: Flag indicating whether to enable structural damping
 - `linear = false`: Flag indicating whether a linear analysis should be performed.
 - `two_dimensional = false`: Flag indicating whether to constrain results to the x-y plane
 - `show_trace = false`: Flag indicating whether to display the solution progress.
 - `save = eachindex(tvec)`: Steps at which to save the time history

 # Initial Condition Analysis Arguments
 - `u0 = fill(zeros(3), length(assembly.points))`: Initial linear displacement of
        each point in the body frame
 - `theta0 = fill(zeros(3), length(assembly.points))`: Initial angular displacement of
        each point in the body frame (using Wiener-Milenkovic Parameters)
 - `V0 = fill(zeros(3), length(assembly.points))`: Initial linear velocity of
        each point in the body frame **excluding contributions from body frame motion**
 - `Omega0 = fill(zeros(3), length(assembly.points))`: Initial angular velocity of
        each point in the body frame **excluding contributions from body frame motion**
 - `Vdot0 = fill(zeros(3), length(assembly.points))`: Initial linear acceleration of
        each point in the body frame **excluding contributions from body frame motion**
 - `Omegadot0 = fill(zeros(3), length(assembly.points))`: Initial angular acceleration of
        each point in the body frame **excluding contributions from body frame motion**

 # Linear Analysis Keyword Arguments
 - `update_linearization = false`: Flag indicating whether to update the linearization state
        variables for a linear analysis with the instantaneous state variables. If `false`,
        then the initial set of state variables will be used for the linearization.

 # Nonlinear Analysis Keyword Arguments
 - `method = :newton`: Method (as defined in NLsolve) to solve nonlinear system of equations
 - `linesearch = LineSearches.BackTracking(maxstep=1e6)`: Line search used to solve
        nonlinear systems of equations
 - `ftol = 1e-9`: tolerance for solving the nonlinear system of equations
 - `iterations = 1000`: maximum iterations for solving the nonlinear system of equations

# Sensitivity Analysis Keyword Arguments
- `pfunc = (p, t) -> (;)`: Function which returns a named tuple with fields corresponding
        to updated versions of the arguments `assembly`, `prescribed_conditions`,
        `distributed_loads`, `point_masses`, `linear_velocity`, `angular_velocity`,
        `linear_acceleration`, `angular_acceleration`, `gravity`, `u0`, `theta0`, `V0`,
        `Omega0`, `Vdot0`, and `Omegadot0`. Only fields contained in the resulting named
        tuple will be overwritten.
- `p`: Sensitivity parameters, as defined in conjunction with the keyword argument `pfunc`.
        While not necessary, using `pfunc` and `p` to define the arguments to this function
        allows automatic differentiation sensitivities to be computed more efficiently
"""
function time_domain_analysis(assembly, tvec; kwargs...)

    system = DynamicSystem(assembly)

    return time_domain_analysis!(system, assembly, tvec; kwargs...)
end

"""
    time_domain_analysis!(system, assembly, tvec; kwargs...)

Pre-allocated version of [`time_domain_analysis`](@ref).
"""
function time_domain_analysis!(system::DynamicSystem, assembly, tvec;
    # general keyword arguments
    prescribed_conditions=Dict{Int,PrescribedConditions{Float64}}(),
    distributed_loads=Dict{Int,DistributedLoads{Float64}}(),
    point_masses=Dict{Int,PointMass{Float64}}(),
    linear_velocity=(@SVector zeros(3)),
    angular_velocity=(@SVector zeros(3)),
    linear_acceleration=(@SVector zeros(3)),
    angular_acceleration=(@SVector zeros(3)),
    gravity=(@SVector zeros(3)),
    # control flag keyword arguments
    reset_state=true,
    initialize=true,
    initial_state=nothing,
    steady_state=false,
    structural_damping=true,
    linear=false,
    two_dimensional=false,
    show_trace=false,
    save=eachindex(tvec),
    # initial condition analysis keyword arguments
    u0=fill((@SVector zeros(3)), length(assembly.points)),
    theta0=fill((@SVector zeros(3)), length(assembly.points)),
    V0=fill((@SVector zeros(3)), length(assembly.points)),
    Omega0=fill((@SVector zeros(3)), length(assembly.points)),
    Vdot0=fill((@SVector zeros(3)), length(assembly.points)),
    Omegadot0=fill((@SVector zeros(3)), length(assembly.points)),
    # linear analysis keyword arguments
    update_linearization=false,
    # nonlinear analysis keyword arguments
    method=:newton,
    linesearch=LineSearches.BackTracking(maxstep=1e6),
    ftol=1e-9,
    iterations=1000,
    # sensitivity analysis keyword arguments
    pfunc = (p, t) -> (;),
    p = nothing,
    )

    # --- Find the Initial States and Rates --- #

    if initialize
        # perform an analysis to find the initial state
        system, initial_state, converged = initial_condition_analysis!(system, assembly, first(tvec);
            # general keyword arguments
            prescribed_conditions=prescribed_conditions,
            distributed_loads=distributed_loads,
            point_masses=point_masses,
            linear_velocity=linear_velocity,
            angular_velocity=angular_velocity,
            linear_acceleration=linear_acceleration,
            angular_acceleration=angular_acceleration,
            gravity=gravity,
            # control flag keyword arguments
            reset_state=reset_state,
            initial_state=initial_state,
            steady_state=steady_state,
            structural_damping=structural_damping,
            linear=linear,
            two_dimensional=two_dimensional,
            show_trace=show_trace,
            # initial condition analysis keyword arguments
            u0=u0,
            theta0=theta0,
            V0=V0,
            Omega0=Omega0,
            Vdot0=Vdot0,
            Omegadot0=Omegadot0,
            # linear analysis keyword arguments
            update_linearization=update_linearization,
            # nonlinear analysis keyword arguments
            method=method,
            linesearch=linesearch,
            ftol=ftol,
            iterations=iterations,
            # sensitivity analysis keyword arguments
            pfunc=pfunc,
            p=p,
            )
    else
        # converged by default
        converged = true
    end

    # --- Initialize Analysis --- #

    # unpack force scaling parameter and system indices
    @unpack force_scaling, indices = system

    # current time
    t = first(tvec)
    dt = zero(t)

    # update the system time
    system.t = t

    # initialize temporary storage for state rate initialization terms
    udot_init = [(@SVector zeros(eltype(p), 3)) for i = 1:length(assembly.points)]
    θdot_init = [(@SVector zeros(eltype(p), 3)) for i = 1:length(assembly.points)]
    Vdot_init = [(@SVector zeros(eltype(p), 3)) for i = 1:length(assembly.points)]
    Ωdot_init = [(@SVector zeros(eltype(p), 3)) for i = 1:length(assembly.points)]

    # augment the parameter vector with space for the newmark-scheme initialization terms
    p_i = zeros(eltype(initial_state), 12*length(assembly.points))
    p_p = p
    p = isnothing(p) ? p_i : vcat(p_i, p_p)

    # initialize convergence flag
    converged = Ref(converged)

    # package up keyword arguments corresponding to this analysis
    constants = (;
        # assembly, indices, control flags, parameter function, current time, and time step
        assembly, indices, two_dimensional, structural_damping, force_scaling, pfunc, t, dt,
        # pointers to the pre-allocated storage and the convergence flag
        x=system.x, resid=system.r, jacob=system.K, converged=converged,
        # default parameters
        prescribed_conditions, distributed_loads, point_masses, gravity,
        linear_velocity, angular_velocity,
        # pointers to the pre-allocated storage for state rate initialization terms
        udot_init, θdot_init, Vdot_init, Ωdot_init,
        # nonlinear analysis keyword arguments
        show_trace, method, linesearch, ftol, iterations,
    )

    # --- Process the Initial States --- #

    # construct initial state and rate vector
    x0, dx0 = newmark_state_vector(system, initial_state, p, constants)

    # current state and rate vector is the initial state and rate vector
    x = x0
    dx = dx0

    # also update the state and rate variables in system.x
    dual_safe_copy!(system.x, x)
    dual_safe_copy!(system.dx, dx)

    # --- Perform Time-Domain Simulation --- #

    # initialize storage for each time step
    history = Vector{AssemblyState{eltype(p)}}(undef, length(save))
    isave = 1

    # add initial state to the solution history
    if isave in save
        history[isave] = initial_state
        isave += 1
    end

    # define the newmark-scheme initialization terms for the first iteration
    if length(tvec) > 1
        # step size corresponding to the next time step
        next_dt = tvec[2] - tvec[1]
        # save initialization terms to the augmented parameter vector
        for ipoint = 1:length(assembly.points)
            irate = 12 * (ipoint - 1)
            p[irate+1:irate+3] = 2/next_dt*initial_state.points[ipoint].u + initial_state.points[ipoint].udot
            p[irate+4:irate+6] = 2/next_dt*initial_state.points[ipoint].theta + initial_state.points[ipoint].thetadot
            p[irate+7:irate+9] = 2/next_dt*initial_state.points[ipoint].V + initial_state.points[ipoint].Vdot
            p[irate+10:irate+12] = 2/next_dt*initial_state.points[ipoint].Omega + initial_state.points[ipoint].Omegadot
        end
    end

    # begin time stepping
    for it in eachindex(tvec)[2:end]

        # current time
        t = tvec[it]

        # current time step size
        dt = tvec[it] - tvec[it-1]

        # print the current time
        if show_trace
            println("Solving for t=$t")
        end

        # update the system time
        system.t = t

        # update the stored time and time step
        constants = (; constants..., t, dt)

        # solve for the new set of state variables
        if linear
            if update_linearization
                x = newmark_lsolve!(x, p, constants)
            else
                x = newmark_lsolve!(x0, p, constants)
            end
        else
            x = ImplicitAD.implicit(newmark_nlsolve!, newmark_residual!, p, constants; drdy=newmark_drdy)
        end

        # add state to history
        if it in save
            history[isave] = newmark_output(system, x, p, constants)
            isave += 1
        end

        # stop early if unconverged
        if !converged[]
            # print error message
            if show_trace
                println("Solution failed to converge")
            end
            # trim the time history
            history = history[1:it]
            # exit time simulation
            break
        end

        # define the newmark-scheme initialization terms for the next iteration
        if length(tvec) > it
            # get updated set of parameters corresponding to the analysis
            assembly, pcond, dload, pmass, gvec, vb_p, ωb_p, udot_init, θdot_init, Vdot_init, Ωdot_init, dt = newmark_parameters(p, constants)
            # step size corresponding to the next time step
            next_dt = (tvec[it+1] - tvec[it])
            # save initialization terms to the augmented parameter vector
            for ipoint = 1:length(assembly.points)
                u, θ = point_displacement(x, ipoint, indices.icol_point, pcond)
                V, Ω = point_velocities(x, indices.icol_point[ipoint])
                # pointer to initialization terms for this point
                irate = 12 * (ipoint - 1)
                # extract current state rates
                udot = 2/dt*u - SVector{3}(p[irate+1], p[irate+2], p[irate+3])
                θdot = 2/dt*θ - SVector{3}(p[irate+4], p[irate+5], p[irate+6])
                Vdot = 2/dt*V - SVector{3}(p[irate+7], p[irate+8], p[irate+9])
                Ωdot = 2/dt*Ω - SVector{3}(p[irate+10], p[irate+11], p[irate+12])
                # save initialization terms
                p[irate+1:irate+3] = 2/next_dt*u + udot
                p[irate+4:irate+6] = 2/next_dt*θ + θdot
                p[irate+7:irate+9] = 2/next_dt*V + Vdot
                p[irate+10:irate+12] = 2/next_dt*Ω + Ωdot
            end
        end

    end

    return system, history, converged[]
end

# combines constant and variable parameters for a newmark-scheme time marching analysis
function newmark_parameters(p, constants)

    # unpack default parameters, parameter function, current time, and step size
    @unpack assembly, prescribed_conditions, distributed_loads, point_masses, gravity,
        linear_velocity, angular_velocity, pfunc, t, dt = constants

    # also unpack indices and temporary storage for rate variables
    @unpack indices, udot_init, θdot_init, Vdot_init, Ωdot_init = constants

    # number of state variables and points
    np = length(assembly.points)

    # separate initialization terms from the parameter vector
    p_i = view(p, 1:12*np) # initialization terms
    p_p = view(p, 12*np+1:length(p)) # provided parameter vector

    # extract state variable initialization terms
    for ipoint = 1:length(assembly.points)
        irate = 12 * (ipoint - 1)
        udot_init[ipoint] = SVector{3}(p_i[irate+1], p_i[irate+2], p_i[irate+3])
        θdot_init[ipoint] = SVector{3}(p_i[irate+4], p_i[irate+5], p_i[irate+6])
        Vdot_init[ipoint] = SVector{3}(p_i[irate+7], p_i[irate+8], p_i[irate+9])
        Ωdot_init[ipoint] = SVector{3}(p_i[irate+10], p_i[irate+11], p_i[irate+12])
    end

    # overwrite default assembly and parameters (if applicable)
    parameters = pfunc(p_p, t)
    assembly = get(parameters, :assembly, assembly)
    prescribed_conditions = get(parameters, :prescribed_conditions, prescribed_conditions)
    distributed_loads = get(parameters, :distributed_loads, distributed_loads)
    point_masses = get(parameters, :point_masses, point_masses)
    gravity = get(parameters, :gravity, gravity)
    linear_velocity = get(parameters, :linear_velocity, linear_velocity)
    angular_velocity = get(parameters, :angular_velocity, angular_velocity)

    # get parameters corresponding to this time step
    pcond = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(t)
    dload = typeof(distributed_loads) <: AbstractDict ? distributed_loads : distributed_loads(t)
    pmass = typeof(point_masses) <: AbstractDict ? point_masses : point_masses(t)
    gvec = typeof(gravity) <: AbstractVector ? SVector{3}(gravity) : SVector{3}(gravity(t))
    vb_p = typeof(linear_velocity) <: AbstractVector ? SVector{3}(linear_velocity) : SVector{3}(linear_velocity(t))
    ωb_p = typeof(angular_velocity) <: AbstractVector ? SVector{3}(angular_velocity) : SVector{3}(angular_velocity(t))

    # return results
    return assembly, pcond, dload, pmass, gvec, vb_p, ωb_p, udot_init, θdot_init, Vdot_init, Ωdot_init, dt
end

# returns state vector for a newmark system corresponding to `state`
function newmark_state_vector(system, state, p, constants)

    # get initial state vector
    if isnothing(state)
        # initial state vector is equal to the system state vector
        x0 = copy(system.x)
        dx0 = copy(system.dx)
    else
        # initialize new state vector storage
        if isnothing(p)
            x0 = similar(system.x, promote_type(eltype(system), eltype(state)))
            dx0 = similar(system.x, promote_type(eltype(system), eltype(state)))
        else
            x0 = similar(system.x, promote_type(eltype(system), eltype(state), eltype(p)))
            dx0 = similar(system.x, promote_type(eltype(system), eltype(state), eltype(p)))
        end
        # combine constants and parameters
        assembly, pcond, dload, pmass, gvec, vb_p, ωb_p,
            udot_init, θdot_init, Vdot_init, Ωdot_init, dt = newmark_parameters(p, constants)
        # set initial state variables in `x`
        set_state!(x0, system, assembly, state; prescribed_conditions=pcond)
        # set initial rate variables in `dx`
        set_rate!(dx0, system, assembly, state; prescribed_conditions=pcond)
    end

    return x0, dx0
end

# defines the newmark residual (for use with ImplicitAD)
function newmark_residual!(resid, x, p, constants)

    # unpack indices and control flags
    @unpack indices, structural_damping, two_dimensional, force_scaling = constants

    # combine constants and parameters
    assembly, pcond, dload, pmass, gvec, vb_p, ωb_p, udot_init, θdot_init, Vdot_init, Ωdot_init, dt = newmark_parameters(p, constants)

    # update acceleration state variable indices
    update_body_acceleration_indices!(indices, pcond)

    # compute and return the residual
    return newmark_system_residual!(resid, x, indices, two_dimensional, force_scaling,
        structural_damping, assembly, pcond, dload, pmass, gvec, vb_p, ωb_p,
        udot_init, θdot_init, Vdot_init, Ωdot_init, dt)
end

# defines the newmark jacobian (for use with ImplicitAD)
function newmark_jacobian!(jacob, x, p, constants)

    # unpack indices and control flags
    @unpack indices, structural_damping, two_dimensional, force_scaling = constants

    # combine constants and parameters
    assembly, pcond, dload, pmass, gvec, vb_p, ωb_p, udot_init, θdot_init, Vdot_init, Ωdot_init, dt = newmark_parameters(p, constants)

    # update acceleration state variable indices
    update_body_acceleration_indices!(indices, pcond)

    # compute and return the jacobian
    return newmark_system_jacobian!(jacob, x, indices, two_dimensional, force_scaling,
        structural_damping, assembly, pcond, dload, pmass, gvec, vb_p, ωb_p,
        udot_init, θdot_init, Vdot_init, Ωdot_init, dt)
end

# jacobian function (in format expected by ImplicitAD)
newmark_drdy(residual, x, p, constants) = newmark_jacobian!(constants.jacob, x, p, constants)

# defines the linear solver (for use with ImplicitAD)
newmark_lsolve!(x0, p, constants) = lsolve!(x0, p, constants, newmark_residual!, newmark_jacobian!)

# defines the nonlinear solver (for use with ImplicitAD)
newmark_nlsolve!(p, constants) = nlsolve!(p, constants, newmark_residual!, newmark_jacobian!)

# returns post-processed state and rate variable vectors
function newmark_output(system, x, p, constants)

    # unpack indices
    @unpack indices = constants

    # combine constants and parameters
    assembly, pcond, dload, pmass, gvec, vb_p, ωb_p, udot_init, θdot_init, Vdot_init, Ωdot_init, dt = newmark_parameters(p, constants)

    # update acceleration state variable indices
    update_body_acceleration_indices!(indices, pcond)

    # initialize state rate vector (if necessary)
    dx = typeof(system.dx) <: typeof(x) ? system.dx : zero(x)

    # populate state rate vector with differentiable variables
    for ipoint in eachindex(assembly.points)
        # compute state rates
        u, θ = point_displacement(x, ipoint, indices.icol_point, pcond)
        V, Ω = point_velocities(x, ipoint, indices.icol_point)
        udot = 2/dt*u - udot_init[ipoint]
        θdot = 2/dt*θ - θdot_init[ipoint]
        Vdot = 2/dt*V - Vdot_init[ipoint]
        Ωdot = 2/dt*Ω - Ωdot_init[ipoint]
        # insert into state rate vector
        icol = indices.icol_point[ipoint]
        udot_udot, θdot_θdot = point_displacement_jacobians(ipoint, pcond)
        dx[icol:icol+2] = udot_udot*udot
        dx[icol+3:icol+5] = θdot_θdot*θdot
        dx[icol+6:icol+8] = Vdot
        dx[icol+9:icol+11] = Ωdot
    end

    # update the state and rate variables in `system`
    dual_safe_copy!(system.x, x)
    dual_safe_copy!(system.dx, dx)

    # return result
    return AssemblyState(dx, x, system, assembly; prescribed_conditions=pcond)
end

# linear analysis function
function lsolve!(x0, p, constants, residual!, jacobian!)

    # unpack pre-allocated storage and the convergence flag
    @unpack resid, jacob, converged = constants

    # update the residual
    residual!(resid, x0, p, constants)

    # update the jacobian
    jacobian!(jacob, x0, p, constants)

    # solve the system
    x = x0 - ImplicitAD.implicit_linear(jacob, resid)

    # set the convergence flag
    converged[] = true

    return x
end

# nonlinear analysis function
function nlsolve!(p, constants, residual!, jacobian!)

    # unpack pre-allocated storage and the convergence flag
    @unpack x, resid, jacob, converged = constants

    # unpack nonlinear solver parameters
    @unpack show_trace, method, linesearch, ftol, iterations = constants

    # wrap the residual
    f!(resid, x) = residual!(resid, x, p, constants)

    # wrap the jacobian
    j!(jacob, x) = jacobian!(jacob, x, p, constants)

    # construct temporary storage
    df = NLsolve.OnceDifferentiable(f!, j!, x, resid, jacob)

    # solve the system
    result = NLsolve.nlsolve(df, x,
        show_trace=show_trace,
        linsolve=(x, A, b) -> x .= ImplicitAD.implicit_linear(A, b),
        method=method,
        linesearch=linesearch,
        ftol=ftol,
        iterations=iterations)

    # update the state, residual, jacobian, and convergence flag
    x .= result.zero
    resid .= df.F
    jacob .= df.DF
    converged[] = result.f_converged

    # return the result
    return result.zero
end

# mass matrix function
function mass_matrix!(jacob, x, p, constants)

    # unpack (default) parameters, parameter function, and current time
    @unpack assembly, prescribed_conditions, point_masses, pfunc, t = constants

    # unpack indices, control flags, and jacobian storage
    @unpack indices, two_dimensional, force_scaling = constants

    # overwrite default assembly and parameters (if applicable)
    parameters = pfunc(p, t)
    assembly = get(parameters, :assembly, assembly)
    prescribed_conditions = get(parameters, :prescribed_conditions, prescribed_conditions)
    point_masses = get(parameters, :point_masses, point_masses)

    # get parameters corresponding to this time step
    pcond = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(t)
    pmass = typeof(point_masses) <: AbstractDict ? point_masses : point_masses(t)

    # compute and return the jacobian
    return system_mass_matrix!(jacob, x, indices, two_dimensional, force_scaling,
        assembly, pcond, pmass)
end

# mass matrix function for an expanded system
function expanded_mass_matrix!(jacob, p, constants)

    # unpack (default) parameters, parameter function, and current time
    @unpack assembly, prescribed_conditions, point_masses, pfunc, t = constants

    # unpack indices, control flags, and jacobian storage
    @unpack indices, two_dimensional, force_scaling = constants

    # overwrite default assembly and parameters (if applicable)
    parameters = pfunc(p, t)
    assembly = get(parameters, :assembly, assembly)
    prescribed_conditions = get(parameters, :prescribed_conditions, prescribed_conditions)
    point_masses = get(parameters, :point_masses, point_masses)

    # get parameters corresponding to this time step
    pcond = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(t)
    pmass = typeof(point_masses) <: AbstractDict ? point_masses : point_masses(t)

    # compute and return the jacobian
    return expanded_system_mass_matrix!(jacob, indices, two_dimensional, force_scaling,
        assembly, pcond, pmass)
end

# copy the entire array
dual_safe_copy!(dst::AbstractArray{T,N}, src::AbstractArray{T,N}) where {T,N} = copy!(dst, src)

# copy only the primal portion of the original array
function dual_safe_copy!(dst, src)

    for i in eachindex(dst, src)
        dst[i] = nondual_value(src[i])
    end

    return dst
end

# return the primal value
nondual_value(x::ForwardDiff.Dual) = ForwardDiff.value(x)
nondual_value(x::ReverseDiff.TrackedReal) = ReverseDiff.value(x)