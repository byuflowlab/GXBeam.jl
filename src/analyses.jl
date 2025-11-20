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
 - `xpfunc = (x, p, t) -> (;)`: Similar to `pfunc`, except that parameters can also be
        defined as a function of GXBeam's state variables.  Using this function forces
        the system jacobian to be computed using automatic differentiation and switches
        the nonlinear solver to a Newton-Krylov solver (with linesearch).
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
    xpfunc = nothing,
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
        assembly, indices, two_dimensional, force_scaling, xpfunc, pfunc, t=first(time),
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
                if isnothing(xpfunc)
                    x = static_lsolve!(x, p, constants)
                else
                    x = static_matrixfree_lsolve!(x, p, constants)
                end
            else
                # use the initial state vector for the linearization
                if isnothing(xpfunc)
                    x = static_lsolve!(x0, p, constants)
                else
                    x = static_matrixfree_lsolve!(x0, p, constants)
                end
            end
        else
            if isnothing(xpfunc)
                x = ImplicitAD.implicit(static_nlsolve!, static_residual!, p, constants; drdy=static_drdy)
            else
                x = ImplicitAD.implicit(static_matrixfree_nlsolve!, static_residual!, p, constants; drdy=matrixfree_jacobian)
            end
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
function static_parameters(x, p, constants)

    # unpack default parameters, parameter function, and current time
    @unpack assembly, prescribed_conditions, distributed_loads, point_masses, gravity,
        xpfunc, pfunc, t = constants

    # overwrite default assembly and parameters (if applicable)
    parameters = isnothing(xpfunc) ? pfunc(p, t) : xpfunc(x, p, t)
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
            x0 = similar(system.x, promote_type(eltype(system), eltype(state))) .= system.x
        else
            x0 = similar(system.x, promote_type(eltype(system), eltype(state), eltype(p))) .= system.x
        end
        # combine constants and parameters
        assembly, pcond, dload, pmass, gvec = static_parameters(x0, p, constants)
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
    assembly, pcond, dload, pmass, gvec = static_parameters(x, p, constants)

    # compute and return the residual
    return static_system_residual!(resid, x, indices, two_dimensional, force_scaling,
        assembly, pcond, dload, pmass, gvec)
end

# jacobian function for a static analysis
function static_jacobian!(jacob, x, p, constants)

    # unpack indices, control flags, and jacobian storage
    @unpack indices, two_dimensional, force_scaling = constants

    # combine constants and parameters
    assembly, pcond, dload, pmass, gvec = static_parameters(x, p, constants)

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

# matrix-free linear solve for a static analysis (for use with ImplicitAD)
static_matrixfree_lsolve!(x0, p, constants) = matrixfree_lsolve!(x0, p, constants, static_residual!, static_jacobian!)

# matrix-free nonlinear solve for a static analysis (for use with ImplicitAD)
static_matrixfree_nlsolve!(p, constants) = matrixfree_nlsolve!(p, constants, static_residual!, static_jacobian!)

# returns post-processed state and rate variable vectors
function static_output!(system, x, p, constants)

    # get new assembly and parameters
    assembly, pcond, dload, pmass, gvec = static_parameters(x, p, constants)

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
 - `xpfunc = (x, p, t) -> (;)`: Similar to `pfunc`, except that parameters can also be
        defined as a function of GXBeam's state variables.  Using this function forces
        the system jacobian to be computed using automatic differentiation and switches
        the nonlinear solver to a Newton-Krylov solver (with linesearch).
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
    xpfunc = nothing,
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
        assembly, indices, structural_damping, two_dimensional, force_scaling, xpfunc, pfunc, t=first(time),
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
                    if isnothing(xpfunc)
                        x = expanded_steady_lsolve!(x, p, constants)
                    else
                        x = expanded_steady_matrixfree_lsolve!(x, p, constants)
                    end
                else
                    if isnothing(xpfunc)
                        x = expanded_steady_lsolve!(x0, p, constants)
                    else
                        x = expanded_steady_matrixfree_lsolve!(x0, p, constants)
                    end
                end
            else
                if update_linearization
                    if isnothing(xpfunc)
                        x = steady_lsolve!(x, p, constants)
                    else
                        x = steady_matrixfree_lsolve!(x, p, constants)
                    end
                else
                    if isnothing(xpfunc)
                        x = steady_lsolve!(x0, p, constants)
                    else
                        x = steady_matrixfree_lsolve!(x0, p, constants)
                    end
                end
            end
        else
            if constant_mass_matrix
                if isnothing(xpfunc)
                    x = ImplicitAD.implicit(expanded_steady_nlsolve!, expanded_steady_residual!, p, constants; drdy=expanded_steady_drdy)
                else
                    x = ImplicitAD.implicit(expanded_steady_matrixfree_nlsolve!, expanded_steady_residual!, p, constants; drdy=matrixfree_jacobian)
                end
            else
                if isnothing(xpfunc)
                    x = ImplicitAD.implicit(steady_nlsolve!, steady_residual!, p, constants; drdy=steady_drdy)
                else
                    x = ImplicitAD.implicit(steady_matrixfree_nlsolve!, steady_residual!, p, constants; drdy=matrixfree_jacobian)
                end
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
function steady_parameters(x, p, constants)

    # unpack default parameters, parameter function, and current time
    @unpack assembly, prescribed_conditions, distributed_loads, point_masses, gravity,
        linear_velocity, angular_velocity, linear_acceleration, angular_acceleration,
        xpfunc, pfunc, t = constants

    # overwrite default assembly and parameters (if applicable)
    parameters = isnothing(xpfunc) ? pfunc(p, t) : xpfunc(x, p, t)
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
            x0 = similar(system.x, promote_type(eltype(system), eltype(state))) .= system.x
        else
            x0 = similar(system.x, promote_type(eltype(system), eltype(state), eltype(p))) .= system.x
        end
        # combine constants and parameters
        assembly, pcond, dload, pmass, gvec = static_parameters(x0, p, constants)
        # set initial state variables in `x`
        set_state!(x0, system, assembly, state; prescribed_conditions=pcond)
    end

    return x0
end

# residual function for a steady analysis (in format expected by ImplicitAD)
function steady_residual!(resid, x, p, constants)
    # @show x

    # unpack indices and control flags
    @unpack indices, structural_damping, two_dimensional, force_scaling = constants

    # combine constants and parameters
    assembly, pcond, dload, pmass, gvec, vb_p, ωb_p, ab_p, αb_p = steady_parameters(x, p, constants)

    # @show dload[1]
    # @show pcond[1]
    # @show typeof(x)

    # compute and return the residual
    return steady_system_residual!(resid, x, indices, two_dimensional, force_scaling,
        structural_damping, assembly, pcond, dload, pmass, gvec, vb_p, ωb_p, ab_p, αb_p)
end

# jacobian function for a steady analysis
function steady_jacobian!(jacob, x, p, constants)

    # unpack indices, control flags, and jacobian storage
    @unpack indices, structural_damping, two_dimensional, force_scaling = constants

    # combine constants and parameters
    assembly, pcond, dload, pmass, gvec, vb_p, ωb_p, ab_p, αb_p = steady_parameters(x, p, constants)

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

# defines the linear solve (for use with ImplicitAD)
steady_matrixfree_lsolve!(x0, p, constants) = matrixfree_lsolve!(x0, p, constants, steady_residual!, steady_jacobian!)

# defines the nonlinear solver (for use with ImplicitAD)
steady_matrixfree_nlsolve!(p, constants) = matrixfree_nlsolve!(p, constants, steady_residual!, steady_jacobian!)

# returns post-processed state and rate variable vectors
function steady_output!(system, x, p, constants)

    # unpack indices
    @unpack indices = constants

    # combine constants and parameters
    assembly, pcond, dload, pmass, gvec, vb_p, ωb_p, ab_p, αb_p = steady_parameters(x, p, constants)

    # initialize state rate vector (if necessary)
    dx = typeof(system.dx) <: typeof(x) ? system.dx .= 0 : zeros(eltype(x), length(x))

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
    assembly, pcond, dload, pmass, gvec, vb_p, ωb_p, ab_p, αb_p = steady_parameters(x, p, constants)

    # compute and return the residual
    return expanded_steady_system_residual!(resid, x, indices, two_dimensional, force_scaling,
        structural_damping, assembly, pcond, dload, pmass, gvec, vb_p, ωb_p, ab_p, αb_p)
end

# jacobian function for a steady analysis
function expanded_steady_jacobian!(jacob, x, p, constants)

    # unpack indices, control flags, and jacobian storage
    @unpack indices, structural_damping, two_dimensional, force_scaling = constants

    # combine constants and parameters
    assembly, pcond, dload, pmass, gvec, vb_p, ωb_p, ab_p, αb_p = steady_parameters(x, p, constants)

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

# defines the linear solve (for use with ImplicitAD)
expanded_steady_matrixfree_lsolve!(x0, p, constants) = matrixfree_lsolve!(x0, p, constants, expanded_steady_residual!, expanded_steady_jacobian!)

# defines the nonlinear solver (for use with ImplicitAD)
expanded_steady_matrixfree_nlsolve!(p, constants) = matrixfree_nlsolve!(p, constants, expanded_steady_residual!, expanded_steady_jacobian!)

# returns post-processed state and rate variable vectors
function expanded_steady_output!(system, x, p, constants)

    # unpack indices
    @unpack indices = constants

    # combine constants and parameters
    assembly, pcond, dload, pmass, gvec, vb_p, ωb_p, ab_p, αb_p = steady_parameters(x, p, constants)

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
  - `xpfunc = (x, p, t) -> (;)`: Similar to `pfunc`, except that parameters can also be
        defined as a function of GXBeam's state variables.  Using this function forces
        the system jacobian to be computed using automatic differentiation and switches
        the nonlinear solver to a Newton-Krylov solver (with linesearch).
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
    xpfunc = nothing,
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
        assembly, indices, structural_damping, two_dimensional, force_scaling, xpfunc, pfunc, t=first(time),
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

    if isnothing(xpfunc)
        # solve for the system stiffness and mass matrices
        if constant_mass_matrix
            expanded_steady_jacobian!(K, x, p, constants)
            expanded_mass_matrix!(M, p, constants)
        else
            steady_jacobian!(K, x, p, constants)
            mass_matrix!(M, x, p, constants)
        end
    else
        if constant_mass_matrix
            autodiff_jacobian!(K, expanded_steady_residual!, x, p, constants)
            expanded_mass_matrix!(M, p, constants)
        else
            autodiff_jacobian!(K, steady_residual!, x, p, constants)
            mass_matrix!(M, x, p, constants)
        end
    end

    # update the jacobians in `system`
    dual_safe_copy!(system.K, K)
    dual_safe_copy!(system.M, M)

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
    f! = (b, x) -> ldiv!(b, Kfact, M * x)
    fc! = (b, x) -> mul!(b, M', Kfact' \ x)
    A = LinearMap{T}(f!, fc!, nx, nx; ismutating=true)

    # compute eigenvalues and eigenvectors
    λ, V = partialeigen(partialschur(A; nev=min(nx, nev), which=ArnoldiMethod.LM(), tol=1e-9)[1])

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
        KmλMfact = factorize(K' + (λ[iλ]+1e-9)'*M')

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
        KmλMfact = factorize(K' + (λ[iλ]+1e-9)'*M')

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
 - `xpfunc = (x, p, t) -> (;)`: Similar to `pfunc`, except that parameters can also be
        defined as a function of GXBeam's state variables.  Using this function forces
        the system jacobian to be computed using automatic differentiation and switches
        the nonlinear solver to a Newton-Krylov solver (with linesearch).
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
    xpfunc = nothing,
    pfunc = (p, t) -> (;),
    p = nothing,
    # eigenvalue analysis keyword arguments
    nev = 6,
    eigenvector_sensitivities=true,
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
            xpfunc = xpfunc,
            pfunc = pfunc,
            p = p,
            )
    else
        # use current state variables
        converged = true
    end

    # unpack force scaling parameter and system pointers
    @unpack force_scaling, indices = system

    # package up keyword arguments corresponding to this analysis
    constants = (;
        # assembly, indices, control flags, parameter function, and current time
        assembly, indices, structural_damping, two_dimensional, force_scaling, xpfunc, pfunc, t=first(time),
        # default parameters
        prescribed_conditions, distributed_loads, point_masses, gravity,
        linear_velocity, angular_velocity, linear_acceleration, angular_acceleration,
    )

    # initialize state vector, jacobian, and mass matrix
    x = steady_state_vector(system, initial_state, p, constants)
    K = spzeros(eltype(x), indices.nstates, indices.nstates)
    M = spzeros(eltype(x), indices.nstates, indices.nstates)

    # compute the jacobian and mass matrix
    if isnothing(xpfunc)
        if constant_mass_matrix
            expanded_steady_jacobian!(K, x, p, constants)
            expanded_mass_matrix!(M, p, constants)
        else
            steady_jacobian!(K, x, p, constants)
            mass_matrix!(M, x, p, constants)
        end
    else
        if constant_mass_matrix
            autodiff_jacobian!(K, expanded_steady_residual!, x, p, constants)
            expanded_mass_matrix!(M, p, constants)
        else
            autodiff_jacobian!(K, steady_residual!, x, p, constants)
            mass_matrix!(M, x, p, constants)
        end
    end

    # update the system storage
    xv = dual_safe_copy!(system.x, x)
    Kv = dual_safe_copy!(system.K, K)
    Mv = dual_safe_copy!(system.M, M)

    # augment parameters with state variables
    px = x
    pp = p
    p = isnothing(p) ? x : vcat(x, p)

    # compute eigenvalues and eigenvectors (without sensitivities)
    λv, Vv = solve_eigensystem(xv, Kv, Mv, nev)

    # resolve eigenvector ambiguities
    imax = zeros(Int, length(λv))
    for iλ = 1:length(λv)
        # extract eigenvector
        v = view(Vv, :, iλ)
        # apply normalization
        Vv[:,iλ] ./= norm(v)
        # resolve phase ambiguity
        imax[iλ] = findfirst(x->abs(x)>1e-4, v)
        Vv[:,iλ] .*= conj(v[imax[iλ]])/abs(v[imax[iλ]])
    end

    # correlate eigenmodes
    if !isnothing(Uprev)
        # construct correlation matrix
        C = Uprev*M*Vv
        # find correct permutation of eigenmodes
        perm, corruption = correlate_eigenmodes(C)
        # rearrange eigenmodes
        λv = λv[perm]
        Vv = Vv[:,perm]
        # also rearrange the index for resolving the phase ambiguity
        imax = imax[perm]
    end

    # compute desired sensitivities
    if eigenvector_sensitivities
        # compute eigenvalue and eigenvector sensitivities

        # store precomputed values
        constants = (; constants..., K_value=Kv, M_value=Mv, K_dual=K, M_dual=M)

        # process sensitivities for each eigenvalue individually
        λ = similar(λv, complex(eltype(x))) # λ (with sensitivities)
        V = similar(Vv, complex(eltype(x))) # V (with sensitivities)
        for iλ = 1:length(λ)
            # update stored eigenvalue and eigenvector
            constants = (; constants..., λ=λv[iλ], v=view(Vv, :, iλ), imax=imax[iλ])
            # compute sensitivities
            eigenstate = ImplicitAD.implicit(eigenproblem_solve!, eigenproblem_residual!, p, constants; drdy=eigenproblem_drdy)
            # extract eigenvalues and eigenvectors (with sensitivities) from the output
            λ[iλ] = complex(eigenstate[2*length(x)+1], eigenstate[2*length(x)+2])
            for ix = 1:length(x)
                V[ix,iλ] = complex(eigenstate[2*ix-1], eigenstate[2*ix])
            end
        end

        # compute corresponding left eigenvectors
        if left
            # find the left eigenvector corresponding to each right eigenvector
            U = left_eigenvectors(K, M, λ, V)
            # return left and right eigenvectors
            return system, λ, U, V, converged
        else
            # return only right eigenvectors
            return system, λ, V, converged
        end

    else
        # compute only eigenvalue sensitivities

        # find the left eigenvector corresponding to each right eigenvector
        Uv = left_eigenvectors(Kv, Mv, λv, Vv)

        # propagate partial derivatives: λdot = -transpose(ui)*(Kdot + Mdot*λi)*vi
        λ = similar(λv, complex(eltype(x)))
        for iλ = 1:length(λ)
            λi = λv[iλ]
            ui = transpose(view(Uv, iλ, :))
            vi = view(Vv, :, iλ)
            # NOTE: Since `ui*(K + M*λi)*vi = 0` we can augment the primal values `λi`
            # with the analytically derivatived sensitivities by adding `ui*(K + M*λi)*vi`
            λ[iλ] = λi - ui*(K + M*λi)*vi
        end

        # ignore sensitivities associated with V
        V = Vv

        # return requested output
        if left
            # ignore sensitivities associated with U
            U = Uv
            # return left and right eigenvectors
            return system, λ, U, V, converged
        else
            # return only right eigenvectors
            return system, λ, V, converged
        end

    end

end

# combines constant and variable parameters for an eigenvalue analysis
function eigenvalue_parameters(p, constants)
    # extract state vector
    nx = constants.indices.nstates
    x = view(p, 1:nx)
    # extract parameters
    p = view(p, nx+1:length(p))
    parameters = steady_parameters(x, p, constants)
    # return state vector and parameters
    return x, parameters...
end

# residual function for an eigenvalue analysis (in format expected by ImplicitAD)
function eigenproblem_residual!(resid, eigenstate, p, constants)
    # unpack indices and control flags
    @unpack indices, two_dimensional, force_scaling, structural_damping, imax = constants
    # get jacobian and mass matrix
    if eltype(resid) <: ForwardDiff.Dual
        # unpack precomputed jacobian and mass matrix
        K, M = constants.K_dual, constants.M_dual
    else
        # combine constants and parameters
        x, assembly, pcond, dload, pmass, gvec, vb_p, ωb_p, ab_p, αb_p = eigenvalue_parameters(p, constants)
        # initialize jacobian and mass matrix
        K = similar(resid, length(x), length(x))
        M = similar(resid, length(x), length(x))
        # populate jacobian and mass matrix
        steady_system_jacobian!(K, x, indices, two_dimensional, force_scaling,
            structural_damping, assembly, pcond, dload, pmass, gvec, vb_p, ωb_p, ab_p, αb_p)
        system_mass_matrix!(M, x, indices, two_dimensional, force_scaling,
            assembly, pcond, pmass)
    end
    # number of state variables
    nx = div(length(eigenstate) - 2, 2)
    # extract real and imaginary parts of the eigenvalue
    λr = eigenstate[2*nx+1]
    λi = eigenstate[2*nx+2]
    # extract real and imaginary parts of the eigenvector
    vr = view(eigenstate, 1:2:2*nx)
    vi = view(eigenstate, 2:2:2*nx)
    # eigenproblem residual: (K + M*λ)*v = 0
    # real((K + M*λ)*v) = K*real(v) + M*real(λ)*real(v) - M*imag(λ)*imag(v) = 0
    rr = view(resid, 1:2:2*nx)
    mul!(rr, K, vr)
    mul!(rr, M, vr, λr, 1)
    mul!(rr, M, vi, -λi, 1)
    # imag((K + M*λ)*v) = K*imag(v) + M*imag(λ)*real(v) + M*real(λ)*imag(v) = 0
    ri = view(resid, 2:2:2*nx)
    mul!(ri, K, vi)
    mul!(ri, M, vr, λi, 1)
    mul!(ri, M, vi, λr, 1)
    # normalization residual: v'*v - 1 = 0
    resid[2*nx+1] = vr'*vr + vi'*vi - 1
    # phase shift residual: argmax(abs, v) = 0
    resid[2*nx+2] = vi[imax]
    # return result
    return resid
end

# jacobian function for an eigenvalue analysis
function eigenproblem_jacobian!(jacob, eigenstate, p, constants)
    # unpack indices and control flags
    @unpack indices, two_dimensional, force_scaling, structural_damping, imax = constants
    # unpack precomputed jacobian and mass matrix
    K, M = constants.K_value, constants.M_value
    # number of state variables
    nx = div(length(eigenstate) - 2, 2)
    # extract real and imaginary parts of eigenvalue
    λr = eigenstate[2*nx+1]
    λi = eigenstate[2*nx+2]
    # extract real and imaginary parts of eigenvector
    vr = view(eigenstate, 1:2:2*nx)
    vi = view(eigenstate, 2:2:2*nx)
    # eigenproblem residual: (K + M*λ)*v = 0
    # residual for real part
    copy!(view(jacob, 1:2:2*nx, 1:2:2*nx), K) # K*dvr
    mul!(view(jacob, 1:2:2*nx, 1:2:2*nx), M, λr, 1, 1) # M*real(λ)*dvr
    mul!(view(jacob, 1:2:2*nx, 2:2:2*nx), M, λi, -1, 1) # -M*imag(λ)*dvi
    mul!(view(jacob, 1:2:2*nx, 2*nx+1), M, vr, 1, 1) # M*real(v)*dλr
    mul!(view(jacob, 1:2:2*nx, 2*nx+2), M, vi, -1, 1) # -M*imag(v)*dλi
    # residual for imaginary part
    copy!(view(jacob, 2:2:2*nx, 2:2:2*nx), K) # K*dvi
    mul!(view(jacob, 2:2:2*nx, 1:2:2*nx), M, λi, 1, 1) # M*imag(λ)*dvr
    mul!(view(jacob, 2:2:2*nx, 2:2:2*nx), M, λr, 1, 1) # M*real(λ)*dvi
    mul!(view(jacob, 2:2:2*nx, 2*nx+1), M, vi) # M*imag(v)*dλr
    mul!(view(jacob, 2:2:2*nx, 2*nx+2), M, vr) # M*real(v)*dλi
    # normalization residual: v'*v - 1 = 0
    # real(v'*v - 1) = real(v)'*real(v) + imag(v)'*imag(v)
    mul!(view(jacob, 2*nx+1, 1:2:2*nx), vr, 2) # real(v)'*real(v)
    mul!(view(jacob, 2*nx+1, 2:2:2*nx), vi, 2) # imag(v)'*imag(v)
    # phase shift residual: argmax(abs, vi) = 0
    # imag(v[i]) = imag(v[i])
    jacob[2*nx+2, 2*imax] = 1
    # return result
    return jacob
end

# jacobian function (in format expected by ImplicitAD)
eigenproblem_drdy(residual, x, p, constants) = eigenproblem_jacobian!(spzeros(eltype(x), length(x), length(x)), x, p, constants)

# dummy function (for use with ImplicitAD)
function eigenproblem_solve!(p, constants)
    # extract current eigenvalue and eigenvector
    @unpack λ, v = constants
    # combine eigenvalue and eigenvector into one real-valued vector
    eigenstate = zeros(real(eltype(λ)), 2*length(v)+2)
    for iv = 1:length(v)
        eigenstate[2*iv-1] = real(v[iv])
        eigenstate[2*iv] = imag(v[iv])
    end
    eigenstate[2*length(v)+1] = real(λ)
    eigenstate[2*length(v)+2] = imag(λ)
    # return result
    return eigenstate
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
 - `xpfunc = (x, p, t) -> (;)`: Similar to `pfunc`, except that parameters can also be
        defined as a function of GXBeam's state variables.  Using this function forces
        the system jacobian to be computed using automatic differentiation and switches
        the nonlinear solver to a Newton-Krylov solver (with linesearch).
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
    xpfunc = nothing,
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
            xpfunc = xpfunc,
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
        assembly, indices, structural_damping, two_dimensional, force_scaling, xpfunc, pfunc, t=first(t0),
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
        Vdot0, Omegadot0 = initial_parameters(FillArrays.Zeros(system.x), p, constants)

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
            if isnothing(xpfunc)
                x = initial_lsolve!(x, p, constants)
            else
                x = initial_matrixfree_lsolve!(x, p, constants)
            end
        else
            if isnothing(xpfunc)
                x = initial_lsolve!(x0, p, constants)
            else
                x = initial_matrixfree_lsolve!(x0, p, constants)
            end
        end
    else
        if isnothing(xpfunc)
            x = ImplicitAD.implicit(initial_nlsolve!, initial_residual!, p, constants; drdy=initial_drdy)
        else
            x = ImplicitAD.implicit(initial_matrixfree_nlsolve!, initial_residual!, p, constants; drdy=matrixfree_jacobian)
        end
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
        xx = similar(original_system.x, eltype(state)) 
        dxx = similar(original_system.dx, eltype(state)) 
        

        # construct state and rate vectors
        set_state!(xx, original_system, assembly, state; prescribed_conditions=pcond) 
        set_rate!(dxx, original_system, assembly, state; prescribed_conditions=pcond) 
    

        # copy only the primal portion of the variables . 
        dual_safe_copy!(system.x, xx)  
        dual_safe_copy!(system.dx, dxx)
         
    end

    # return the system, state, and the (de-referenced) convergence flag
    return original_system, state, converged[]
end


function initial_state_analysis!(system, assembly, t0;
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
    xpfunc = nothing,
    pfunc = (p, t) -> (;),
    p = nothing,
    )

    # perform steady state analysis instead (if requested)
    if steady_state
        error("initial_state_analysis!(): steady state not set up. ")
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
            xpfunc = xpfunc,
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
        assembly, indices, structural_damping, two_dimensional, force_scaling, xpfunc, pfunc, t=first(t0),
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
        Vdot0, Omegadot0 = initial_parameters(FillArrays.Zeros(system.x), p, constants)

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
            if isnothing(xpfunc)
                x = initial_lsolve!(x, p, constants)
            else
                x = initial_matrixfree_lsolve!(x, p, constants)
            end
        else
            if isnothing(xpfunc)
                x = initial_lsolve!(x0, p, constants)
            else
                x = initial_matrixfree_lsolve!(x0, p, constants)
            end
        end
    else
        if isnothing(xpfunc)
            x = ImplicitAD.implicit(initial_nlsolve!, initial_residual!, p, constants; drdy=initial_drdy)
            
        else
            x = ImplicitAD.implicit(initial_matrixfree_nlsolve!, initial_residual!, p, constants; drdy=matrixfree_jacobian)
        end
    end

    # NOTE: `x`, `r`, `K`, and `converged` are updated in lsolve!/nlsolve!

    # postprocess the dual number state
    state = initial_output(system, x, p, constants)

    # copy the state and rate variables to the original system
    # if eltype(original_system) <: eltype(state)
    #     # set the values of the original system directly
    #     set_state!(original_system, assembly, state; prescribed_conditions=pcond)
    #     set_rate!(original_system, assembly, state; prescribed_conditions=pcond)

    # else
    #     # initialize storage
        xx = similar(original_system.x, eltype(state)) 
        dxx = similar(original_system.dx, eltype(state)) 
        

        # construct state and rate vectors
        set_state!(xx, original_system, assembly, state; prescribed_conditions=pcond) 
        set_rate!(dxx, original_system, assembly, state; prescribed_conditions=pcond) 
    

        # copy only the primal portion of the variables . 
        dual_safe_copy!(system.x, xx)  
        dual_safe_copy!(system.dx, dxx)
         
    # end

    # return the system, state, and the (de-referenced) convergence flag
    return xx, dxx, state, converged[]
end

function initial_parameters(x, p, constants)

    # unpack default parameters, parameter function, and current time
    @unpack assembly, prescribed_conditions, distributed_loads, point_masses, gravity,
        linear_velocity, angular_velocity, linear_acceleration, angular_acceleration,
        xpfunc, pfunc, t = constants

    # also unpack initial conditions
    @unpack u0, theta0, V0, Omega0, Vdot0, Omegadot0 = constants

    # overwrite default assembly and parameters (if applicable)
    parameters = isnothing(xpfunc) ? pfunc(p, t) : xpfunc(x, p, t)
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
            x0 = similar(system.x, promote_type(eltype(system), eltype(state))) .= system.x
        else
            x0 = similar(system.x, promote_type(eltype(system), eltype(state), eltype(p))) .= system.x
        end
        # combine constants and parameters
        assembly, pcond, dload, pmass, gvec = static_parameters(x0, p, constants)
        # set initial state variables in `x`
        set_initial_state!(x0, system, rate_vars, state, pcond)
    end

    return x0
end

# defines the initial residual (for use with ImplicitAD)
function initial_residual!(resid, x, p, constants)

    # unpack indices, control flags, and parameter function
    @unpack indices, structural_damping, two_dimensional, force_scaling = constants

    @unpack rate_vars1, rate_vars2 = constants

    # combine constants and parameters
    assembly, pcond, dload, pmass, gvec, vb_p, ωb_p, ab_p, αb_p, u0, theta0, V0, Omega0,
        Vdot0, Omegadot0 = initial_parameters(x, p, constants)

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
    @unpack indices, structural_damping, two_dimensional, force_scaling = constants

    @unpack rate_vars1, rate_vars2 = constants

    # combine constants and parameters
    assembly, pcond, dload, pmass, gvec, vb_p, ωb_p, ab_p, αb_p, u0, theta0, V0, Omega0,
        Vdot0, Omegadot0 = initial_parameters(x, p, constants)

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

# defines the linear solve (for use with ImplicitAD)
initial_matrixfree_lsolve!(x0, p, constants) = matrixfree_lsolve!(x0, p, constants, initial_residual!, initial_jacobian!)

# defines the nonlinear solver (for use with ImplicitAD)
initial_matrixfree_nlsolve!(p, constants) = matrixfree_nlsolve!(p, constants, initial_residual!, initial_jacobian!)

# returns post-processed state and rate variable vectors
function initial_output(system, x, p, constants)

    # unpack indices
    @unpack indices, rate_vars1, rate_vars2, force_scaling = constants

    # combine constants and parameters
    assembly, pcond, dload, pmass, gvec, vb_p, ωb_p, ab_p, αb_p, u0, theta0, V0, Omega0,
        Vdot0, Omegadot0 = initial_parameters(x, p, constants)

    # initialize state rate vector (if necessary)
    xx = typeof(system.x) <: typeof(x) ? system.x .= x : similar(x) .= x
    dx = typeof(system.dx) <: typeof(x) ? system.dx : similar(x) .= 0

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
        set_linear_displacement!(xx, system, pcond, u, ipoint)
        set_angular_displacement!(xx, system, pcond, θ, ipoint)
        set_linear_velocity!(xx, system, V, ipoint)
        set_angular_velocity!(xx, system, Ω, ipoint)
        set_external_forces!(xx, system, pcond, F, ipoint)
        set_external_moments!(xx, system, pcond, M, ipoint)

        # insert result into the rate vector
        icol = indices.icol_point[ipoint]
        udot_udot, θdot_θdot = point_displacement_jacobians(ipoint, pcond)
        dx[icol:icol+2] = udot_udot*udot
        dx[icol+3:icol+5] = θdot_θdot*θdot
        dx[icol+6:icol+8] = Vdot
        dx[icol+9:icol+11] = Ωdot
    end


    # update the state and rate variables in `system`
    dual_safe_copy!(system.x, xx)
    dual_safe_copy!(system.dx, dx)


    # return result
    state = AssemblyState(dx, xx, system, assembly; prescribed_conditions=pcond)

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
 - `xpfunc = (x, p, t) -> (;)`: Similar to `pfunc`, except that parameters can also be
        defined as a function of GXBeam's state variables.  Using this function forces
        the system jacobian to be computed using automatic differentiation and switches
        the nonlinear solver to a Newton-Krylov solver (with linesearch).
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

function show_tracked_array(arr)
    if isa(arr[1], ReverseDiff.TrackedReal)
        newarr = [arr[i].value for i in eachindex(arr)]
        println(newarr)

    elseif isa(arr[1], ForwardDiff.Dual)
        newarr = [arr[i].value for i in eachindex(arr)]
        println(newarr)

    else
        println(arr)
    end
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
    xpfunc = nothing,
    pfunc = (p, t) -> (;),
    p = nothing,
    )

    # --- Find the Initial States and Rates --- #

    if isnothing(initial_state)
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
            xpfunc=xpfunc,
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

    # initialize convergence flag
    converged = Ref(converged)


    # package up keyword arguments corresponding to this analysis
    constants = (;
        # assembly, indices, control flags, parameter function, current time, and time step
        assembly, indices, two_dimensional, structural_damping, force_scaling, xpfunc, pfunc, t, dt,
        # pointers to the pre-allocated storage and the convergence flag
        x=system.x, resid=system.r, jacob=system.K, converged=converged,
        # default parameters
        prescribed_conditions, distributed_loads, point_masses, gravity,
        linear_velocity, angular_velocity,
        # nonlinear analysis keyword arguments
        show_trace, method, linesearch, ftol, iterations,
    )

    # --- Process the Initial States --- #
    # construct initial state and rate vector
    x0, dx0 = newmark_state_vector(system, initial_state, p, constants)


    # current state and rate vector is the initial state and rate vector
    x = x0
    dx = dx0
    # println("x0")
    # show_tracked_array(x0)
    # println("dx0")
    # show_tracked_array(dx0)

    # also update the state and rate variables in system.x
    dual_safe_copy!(system.x, x)
    dual_safe_copy!(system.dx, dx)


    # --- Initialize Time-Domain Solution --- #

    # initialize storage for each time step
    history = Vector{AssemblyState{eltype(x), Vector{PointState{eltype(x)}}, Vector{ElementState{eltype(x)}}}}(undef, length(save))
    isave = 1

    # add initial state to the solution history
    if isave in save
        history[isave] = initial_state
        isave += 1
    end

    # augment parameter vector with space for the initial states
    if isnothing(p)
        # parameter vector is non-existant
        paug = zeros(eltype(x), 12*length(assembly.points))
    else
        # parameter vector exists
        paug = zeros(eltype(x), 12*length(assembly.points) + length(p))
        # copy parameter values to the augmented parameter vector
        paug[12*length(assembly.points) + 1 : 12*length(assembly.points) + length(p)] .= p
    end

    # --- Begin Time Stepping --- #

    for it = 2:length(tvec)

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

        # get updated prescribed conditions
        parameters = isnothing(xpfunc) ? pfunc(p, t) : xpfunc(x, p, t)
        pcond = get(parameters, :prescribed_conditions, prescribed_conditions)
        pcond = typeof(pcond) <: AbstractDict ? pcond : pcond(t)

        # update parameter vector with new initialization terms
        for ipoint = 1:length(assembly.points)
            # index into parameter vector
            irate = 12 * (ipoint - 1)

            # extract state variables
            u, θ = point_displacement(x, ipoint, indices.icol_point, pcond)
            V, Ω = point_velocities(x, indices.icol_point[ipoint])

            # extract rate variables
            if it <= 2

                # println("TDA: Got here. $it")
                # if ipoint==length(assembly.points)
                #     @show u, θ
                # end

                udot = initial_state.points[ipoint].udot
                θdot = initial_state.points[ipoint].thetadot
                Vdot = initial_state.points[ipoint].Vdot
                Ωdot = initial_state.points[ipoint].Omegadot

                # if ipoint==length(assembly.points)
                #     @show udot, θdot, Vdot, Ωdot
                # end

            else
                dtprev = tvec[it-1] - tvec[it-2]
                
                udot = 2/dtprev*u - SVector(paug[irate+1], paug[irate+2], paug[irate+3])
                θdot = 2/dtprev*θ - SVector(paug[irate+4], paug[irate+5], paug[irate+6])
                Vdot = 2/dtprev*V - SVector(paug[irate+7], paug[irate+8], paug[irate+9])
                Ωdot = 2/dtprev*Ω - SVector(paug[irate+10], paug[irate+11], paug[irate+12])
            end

            # store initialization terms in the parameter vector
            paug[irate+1:irate+3] = 2/dt*u + udot
            paug[irate+4:irate+6] = 2/dt*θ + θdot
            paug[irate+7:irate+9] = 2/dt*V + Vdot
            paug[irate+10:irate+12] = 2/dt*Ω + Ωdot

            # if ipoint==length(assembly.points)&&it==2
            #     # println("TDA: $it")
            #     @show typeof(paug)
            #     println("")
            # end
        end

        # solve for the new set of state variables
        if linear
            if update_linearization
                if isnothing(xpfunc)
                    x = newmark_lsolve!(x, paug, constants)
                else
                    x = newmark_matrixfree_lsolve!(x, paug, constants)
                end

            else
                if isnothing(xpfunc)
                    x = newmark_lsolve!(x0, paug, constants)
                else
                    x = newmark_matrixfree_lsolve!(x0, paug, constants)
                end
            end

        else
            if isnothing(xpfunc)
                # @show t
                # println("this statement")

                x = ImplicitAD.implicit(newmark_nlsolve!, newmark_residual!, paug, constants; drdy=newmark_drdy)

                # println("x")
                # @show_tracked_array x

            else
                x = ImplicitAD.implicit(newmark_matrixfree_nlsolve!, newmark_residual!, paug, constants; drdy=matrixfree_jacobian)
            end
        end

        
        # add state to history
        if it in save
            # @show paug
            history[isave] = newmark_output(system, x, paug, constants)
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

    end

    return system, history, converged[]
end

"""
    initialize_system!(system, assembly, tvec; kwargs...)

Pre-allocate the history for a stepped time-domain analysis and conduct the initial condition analysis. 
"""
function initialize_system!(system::DynamicSystem, assembly, tvec;
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
    xpfunc = nothing,
    pfunc = (p, t) -> (;),
    p = nothing,
    )
    # --- Find the Initial States and Rates --- #

    if isnothing(initial_state)
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
            xpfunc=xpfunc,
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

    # initialize convergence flag
    converged = Ref(converged)


    # package up keyword arguments corresponding to this analysis
    constants = (;
        # assembly, indices, control flags, parameter function, current time, and time step
        assembly, indices, two_dimensional, structural_damping, force_scaling, xpfunc, pfunc, t, dt,
        # pointers to the pre-allocated storage and the convergence flag
        x=system.x, resid=system.r, jacob=system.K, converged=converged,
        # default parameters
        prescribed_conditions, distributed_loads, point_masses, gravity,
        linear_velocity, angular_velocity,
        # nonlinear analysis keyword arguments
        show_trace, method, linesearch, ftol, iterations,
    )

    # --- Process the Initial States --- #
    # construct initial state and rate vector
    x0, dx0 = newmark_state_vector(system, initial_state, p, constants)


    # current state and rate vector is the initial state and rate vector
    x = x0
    dx = dx0
    # println("x0")
    # show_tracked_array(x0)
    # println("dx0")
    # show_tracked_array(dx0)

    # also update the state and rate variables in system.x
    dual_safe_copy!(system.x, x)
    dual_safe_copy!(system.dx, dx)


    # --- Initialize Time-Domain Solution --- #

    # initialize storage for each time step #Todo: I should check the types at some point. 
    history = Vector{AssemblyState{eltype(x), Vector{PointState{eltype(x)}}, Vector{ElementState{eltype(x)}}}}(undef, length(save))
    isave = 1

    # add initial state to the solution history
    if isave in save
        history[isave] = initial_state
        isave += 1
    end

    # augment parameter vector with space for the initial states
    if isnothing(p)
        # parameter vector is non-existant
        # println("No parameter vector")
        paug = zeros(eltype(x), 12*length(assembly.points))
    else
        # parameter vector exists
        # println("Parameter vector exists")
        paug = zeros(eltype(x), 12*length(assembly.points) + length(p))
        # copy parameter values to the augmented parameter vector
        paug[12*length(assembly.points) + 1 : 12*length(assembly.points) + length(p)] .= p
    end


    return system, history, constants, paug, x, converged[]
end


"""
    step_system!(history, system, assembly, tvec; kwargs...)

A similar function to time_domain_analysis, but instead history is already allocated and it only takes one step with the system. 

**Inputs**
- system::DynamicSystem : The dynamic system to be stepped.
- paug: The augmented parameter vector.
- x: The current state vector.
- constants: The constants tuple.
- initial_state: The initial state of the system.
- assembly: The assembly of the system.
- tvec: The time vector.
- i: The current time step index.
- prescribed_conditions: A dictionary of prescribed conditions.
- distributed_loads: A dictionary of distributed loads.
- point_masses: A dictionary of point masses.
- linear_velocity: The linear velocity of the body frame.
- angular_velocity: The angular velocity of the body frame.
- linear_acceleration: The linear acceleration of the body frame.
- angular_acceleration: The angular acceleration of the body frame.
- gravity: The gravity vector in the body frame.
- structural_damping: A flag indicating whether to enable structural damping.
- linear: A flag indicating whether a linear analysis should be performed.
- show_trace: A flag indicating whether to display the solution progress.
- update_linearization: A flag indicating whether to update the linearization state variables for a linear analysis with the instantaneous state variables.
- xpfunc: A function that returns parameters as a function of the state variables.
- pfunc: A function that returns parameters as a function of time.
- p: Sensitivity parameters.
- converged: A reference to a boolean indicating whether the solution converged.
"""
function step_system!(system::DynamicSystem, paug, x, constants, initial_state, assembly, tvec, i;
    # general keyword arguments
    prescribed_conditions=Dict{Int,PrescribedConditions{Float64}}(),
    distributed_loads=Dict{Int,DistributedLoads{Float64}}(),
    point_masses=Dict{Int,PointMass{Float64}}(),
    linear_velocity=(@SVector zeros(3)),
    angular_velocity=(@SVector zeros(3)),
    linear_acceleration=(@SVector zeros(3)),
    angular_acceleration=(@SVector zeros(3)),
    gravity=(@SVector zeros(3)),
    ## control flag keyword arguments
    structural_damping=true,
    linear=false,
    show_trace=false,
    ## linear analysis keyword arguments
    update_linearization=false,
    ## nonlinear analysis keyword arguments
    xpfunc = nothing,
    pfunc = (p, t) -> (;),
    p = nothing,
    converged=true)


    @unpack force_scaling, indices = system

    constants = (; constants..., prescribed_conditions, distributed_loads, point_masses, linear_velocity, angular_velocity, linear_acceleration, angular_acceleration, gravity, structural_damping)


    # current time
    t = tvec[i]

    # current time step size
    if i == 1
        dt = 0 #Todo: Should this case ever happen? If so, what safeguards can I put up to protect against it? 
    else
        dt = tvec[i] - tvec[i-1]
    end


    ### print the current time
    if show_trace
        println("Solving for t=$t")
    end

    ### update the system
    system.t = t

    ### update the stored time and time step
    constants = (; constants..., t, dt)

    ### get updated prescribed conditions
    parameters = isnothing(xpfunc) ? pfunc(p, t) : xpfunc(x, p, t)
    pcond = get(parameters, :prescribed_conditions, prescribed_conditions)
    pcond = typeof(pcond) <: AbstractDict ? pcond : pcond(t)

    ### update parameter vector with new initialization terms
    for ipoint = 1:length(assembly.points)
        ### index into parameter vector
        irate = 12 * (ipoint - 1)

        # extract state variables
        u, θ = point_displacement(x, ipoint, indices.icol_point, pcond)
        V, Ω = point_velocities(x, indices.icol_point[ipoint])

        # extract rate variables
        if i <= 2 
            udot = initial_state.points[ipoint].udot 
            θdot = initial_state.points[ipoint].thetadot
            Vdot = initial_state.points[ipoint].Vdot
            Ωdot = initial_state.points[ipoint].Omegadot

        else
            dtprev = tvec[i-1] - tvec[i-2]
                
            udot = 2/dtprev*u - SVector(paug[irate+1], paug[irate+2], paug[irate+3])
            θdot = 2/dtprev*θ - SVector(paug[irate+4], paug[irate+5], paug[irate+6])
            Vdot = 2/dtprev*V - SVector(paug[irate+7], paug[irate+8], paug[irate+9])
            Ωdot = 2/dtprev*Ω - SVector(paug[irate+10], paug[irate+11], paug[irate+12])
        end

        # store initialization terms in the parameter vector
        paug[irate+1:irate+3] = 2/dt*u + udot
        paug[irate+4:irate+6] = 2/dt*θ + θdot
        paug[irate+7:irate+9] = 2/dt*V + Vdot
        paug[irate+10:irate+12] = 2/dt*Ω + Ωdot
    end

    # solve for the new set of state variables
    if linear
        if update_linearization
            if isnothing(xpfunc)
                x = newmark_lsolve!(x, paug, constants)
            else
                x = newmark_matrixfree_lsolve!(x, paug, constants)
            end

        else
            if isnothing(xpfunc)
                x = newmark_lsolve!(x0, paug, constants)
            else
                x = newmark_matrixfree_lsolve!(x0, paug, constants)
            end
        end

    else
        if isnothing(xpfunc)
            x = ImplicitAD.implicit(newmark_nlsolve!, newmark_residual!, paug, constants; drdy=newmark_drdy)
        else
            x = ImplicitAD.implicit(newmark_matrixfree_nlsolve!, newmark_residual!, paug, constants; drdy=matrixfree_jacobian)
        end
    end

    historyi = newmark_output(system, x, paug, constants)

    # stop early if unconverged
    if !converged[]
        # print error message
        if show_trace
            println("Solution failed to converge")
        end
    end

    return system, historyi, constants, paug, x, converged[]
end
    



function implicit_euler_residual!(resid, rn, inputs, p) 

    t, dt, n, nd,
    assembly, prescribed_conditions, distributed_loads, point_masses,
    gravity, linear_velocity, angular_velocity,
    xpfunc, pfunc,
    indices, force_scaling, structural_damping, two_dimensional = p

    ### Extract the states
    x = view(inputs, 1:n) #previous states


    ### Update the inputs based on the design variables. 
    if nd>0
        xd = view(inputs, 2n+1:2n+nd) #Design variables
    else
        xd = nothing
    end

    parameters = isnothing(xpfunc) ? pfunc(xd, t) : xpfunc(x, xd, t) #Todo: This will unfortunately get called at every time step. 
    assmb = get(parameters, :assembly, assembly)
    pcs = get(parameters, :prescribed_conditions, prescribed_conditions)
    dls = get(parameters, :distributed_loads, distributed_loads)
    pms = get(parameters, :point_masses, point_masses)
    grv = get(parameters, :gravity, gravity)
    lv = get(parameters, :linear_velocity, linear_velocity)
    av = get(parameters, :angular_velocity, angular_velocity)

    pcond = typeof(pcs) <: AbstractDict ? pcs : pcs(t)
    dload = typeof(dls) <: AbstractDict ? dls : dls(t)
    pmass = typeof(pms) <: AbstractDict ? pms : pms(t)
    gvec = typeof(grv) <: AbstractVector ? SVector{3}(grv) : SVector{3}(grv(t))
    vvec = typeof(lv) <: AbstractVector ? SVector{3}(lv) : SVector{3}(lv(t))
    omegavec = typeof(av) <: AbstractVector ? SVector{3}(av) : SVector{3}(av(t))


    #Extract the outputs. 
    xn = view(rn, 1:n) #New state
    dxn = view(rn, n+1:2n) #New state rates

    constraint = view(resid, 1:n) #The algebraic equations of the DAE. 
    integrand = view(resid, n+1:2n) #The differential equations of the DAE. 


    #Constraining the states and state rates to be consistent. #Todo: Something is sneaking in as a trackedReal. I suspect it's the assembly (from pfunc)
    dynamic_system_residual!(constraint, dxn, xn, indices, two_dimensional, force_scaling,
        structural_damping, assmb, pcond, dload, pmass,
        gvec, vvec, omegavec)

    #Constraining the state rates to be an implicit Euler step. 
    integrand .= @. x + dxn*dt - xn 
end

function implicit_euler_solve(inputs, p)

    dt = p[2]
    n = p[3]
    
    x = view(inputs, 1:n)
    dx = view(inputs, n+1:2*n)
    #TODO: I could make a smaller inputs for the residual. 
    
    residual! = (resid, rn) -> implicit_euler_residual!(resid, rn, inputs, p)

    #TODO: Is there a way I can provide the jacobian to speed up calcs? or minimally decrease the length of the tape. 
    result = nlsolve(residual!, vcat(x + dx.*dt, dx)) #Adding autodiff=:forward didn't decrease the length of the tape. 

    return result.zero
end

"""
    take_step(x, dx, system, assembly, t, tprev, 
    prescribed_conditions, distributed_loads, point_masses,
    gravity, linear_velocity, angular_velocity,
    xpfunc, pfunc, p,
    structural_damping, two_dimensional)

Take a single implicit Euler time step for the given system.
"""
function take_step(x, dx, system, assembly, t, tprev, 
    prescribed_conditions, distributed_loads, point_masses,
    gravity, linear_velocity, angular_velocity,
    xpfunc, pfunc, p,
    structural_damping, two_dimensional)
    
    #TODO: Make an inplace version that doesn't allocate.
    @unpack force_scaling, indices = system

    # @show x
    # @show dx

    ### Prepare the vector of design vars (which includes the states and state rates)
    if isnothing(p)
        nd = 0
        inputs = vcat(x, dx)
    else
        nd = length(p)
        inputs = vcat(x, dx, p)
    end

    ### Fixed parameters (i.e. not going to have derivatives)
    dt = t-tprev
    n = length(x)

    params = (t, dt, n, nd,
    assembly, prescribed_conditions, distributed_loads, point_masses,
    gravity, linear_velocity, angular_velocity,
    xpfunc, pfunc,
    indices, force_scaling, structural_damping, two_dimensional)

    
    x_out = implicit(implicit_euler_solve, implicit_euler_residual!, inputs, params)

    # , result.x_converged #TODO: get the convergence. 
    return x_out[1:n], x_out[n+1:end]
end



"""
    simulate(assembly, tvec; prescribed_conditions, distributed_loads, point_masses, gravity, linear_velocity, angular_velocity,
    xpfunc, pfunc, p, structural_damping, two_dimensional, verbose)

Use the take_step function to run a time_domain_analysis. 

WARNING: The implicit Euler method is not as stable as the Newmark-beta method. You may experience instability for stiff systems or large time steps.
"""
function simulate(assembly, tvec; 
    prescribed_conditions=Dict{Int,PrescribedConditions{Float64}}(),
    distributed_loads=Dict{Int,DistributedLoads{Float64}}(),
    point_masses=Dict{Int,PointMass{Float64}}(),
    gravity=(@SVector zeros(3)), linear_velocity=(@SVector zeros(3)), angular_velocity=(@SVector zeros(3)),
    xpfunc = nothing, pfunc = (p, t) -> (;), p = nothing,
    structural_damping::Bool=true, two_dimensional::Bool=false, verbose::Bool=false)

    ### Initialize the solution
    if verbose
        println("Initializing...")
    end

    system = DynamicSystem(assembly)

    nt = length(tvec)
    t0 = first(tvec)


    x0, dx0, initial_state, converged0 = initial_state_analysis!(system, assembly, t0;
    prescribed_conditions, distributed_loads, point_masses, gravity, linear_velocity, angular_velocity,
    structural_damping, two_dimensional)

    # x0 = system.x
    # dx0 = system.dx

    nx = length(x0)

    inittype = eltype(x0)
    history = Vector{AssemblyState{inittype, Vector{PointState{inittype}}, Vector{ElementState{inittype}}}}(undef, nt)
    x = zeros(nt, nx)
    dx = zeros(nt, nx)
    converged = Vector{Bool}(undef, nt)

    x[1,:] .= x0
    dx[1,:] .= dx0
    converged[1] = converged0
    history[1] = initial_state

    
    ### Loop through the solution
    if verbose
        println("Solving...")
    end
    for i = 2:nt
        if verbose
            println("Step $i")
        end
        #Time step information
        t = tvec[i]
        tprev = tvec[i-1]
        dt = t-tprev

        #TODO: I probably will need to update the prescribed_conditions, distributed_loads, gravity, linear_velocity, and anything else that I'm supposed to 
        parameters = isnothing(xpfunc) ? pfunc(p, t) : xpfunc(x[i,:], p, t)
        pcond = get(parameters, :prescribed_conditions, prescribed_conditions)
        pcond = typeof(pcond) <: AbstractDict ? pcond : pcond(t)


        x[i,:], dx[i,:] = take_step(x[i-1,:], dx[i-1,:], system, assembly, t, tprev, 
                            prescribed_conditions, distributed_loads, point_masses,
                            gravity, linear_velocity, angular_velocity,
                            xpfunc, pfunc, p,
                            structural_damping, two_dimensional)

        #calculate output state
        history[i] = AssemblyState(dx[i,:], x[i,:], system, assembly; prescribed_conditions=pcond)
    end


    return history, converged
end



# combines constant and variable parameters for a newmark-scheme time marching analysis
function newmark_parameters(p, constants)

    # unpack default parameters, parameter function, current time, and step size
    @unpack assembly, prescribed_conditions, distributed_loads, point_masses, gravity,
        linear_velocity, angular_velocity, xpfunc, pfunc, t, dt = constants

    # also unpack indices
    @unpack indices = constants

    # number of state variables and points
    np = length(assembly.points)

    # separate initialization terms from the parameter vector
    p_i = view(p, 1:12*np) # initialization terms
    p_p = view(p, 12*np+1:length(p)) # provided parameter vector

    # define state rate initialization terms
    udot_init = [SVector{3}(p_i[12*(ip-1)+1], p_i[12*(ip-1)+2], p_i[12*(ip-1)+3]) for ip = 1:length(assembly.points)]
    θdot_init = [SVector{3}(p_i[12*(ip-1)+4], p_i[12*(ip-1)+5], p_i[12*(ip-1)+6]) for ip = 1:length(assembly.points)]
    Vdot_init = [SVector{3}(p_i[12*(ip-1)+7], p_i[12*(ip-1)+8], p_i[12*(ip-1)+9]) for ip = 1:length(assembly.points)]
    Ωdot_init = [SVector{3}(p_i[12*(ip-1)+10], p_i[12*(ip-1)+11], p_i[12*(ip-1)+12]) for ip = 1:length(assembly.points)]

    # overwrite default assembly and parameters (if applicable)
    parameters = isnothing(xpfunc) ? pfunc(p_p, t) : xpfunc(x, p_p, t)
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
            x0 = similar(system.x, promote_type(eltype(system), eltype(state))) .= system.x
            dx0 = similar(system.x, promote_type(eltype(system), eltype(state))) .= system.dx
        else
            x0 = similar(system.x, promote_type(eltype(system), eltype(state), eltype(p))) .= system.x
            dx0 = similar(system.x, promote_type(eltype(system), eltype(state), eltype(p))) .= system.dx
        end
        # combine constants and parameters
        assembly, pcond, dload, pmass, gvec = static_parameters(x0, p, constants)
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

# defines the linear solve (for use with ImplicitAD)
newmark_matrixfree_lsolve!(x0, p, constants) = matrixfree_lsolve!(x0, p, constants, newmark_residual!, newmark_jacobian!)

# defines the nonlinear solver (for use with ImplicitAD)
newmark_matrixfree_nlsolve!(p, constants) = matrixfree_nlsolve!(p, constants, newmark_residual!, newmark_jacobian!)

# returns post-processed state and rate variable vectors
function newmark_output(system, x, p, constants)

    # unpack indices
    @unpack indices = constants

    # combine constants and parameters
    assembly, pcond, dload, pmass, gvec, vb_p, ωb_p, udot_init, θdot_init, Vdot_init, Ωdot_init, dt = newmark_parameters(p, constants)

    # initialize state rate vector (if necessary)
    dx = typeof(system.dx) <: typeof(x) ? system.dx .= 0 : similar(x) .= 0

    # populate state rate vector with differentiable variables
    for ipoint in eachindex(assembly.points)
        # index into parameter vector
        irate = 12 * (ipoint - 1)

        # state variables
        u, θ = point_displacement(x, ipoint, indices.icol_point, pcond)
        V, Ω = point_velocities(x, ipoint, indices.icol_point)

        # rate variables
        udot = 2/dt*u - SVector(p[irate+1], p[irate+2], p[irate+3])
        θdot = 2/dt*θ - SVector(p[irate+4], p[irate+5], p[irate+6])
        Vdot = 2/dt*V - SVector(p[irate+7], p[irate+8], p[irate+9])
        Ωdot = 2/dt*Ω - SVector(p[irate+10], p[irate+11], p[irate+12])

        # save state rates
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

    # @show x
    # @show dx

    # return result
    return AssemblyState(dx, x, system, assembly; prescribed_conditions=pcond)
end

# linear analysis function
function lsolve!(x0, p, constants, residual!, jacobian!)

    # unpack pre-allocated storage and the convergence flag
    @unpack x, resid, jacob, converged = constants

    # update the residual
    residual!(resid, x0, p, constants)

    # update the jacobian
    jacobian!(jacob, x0, p, constants)

    # update the state
    x .= x0 .- ImplicitAD.implicit_linear(jacob, resid)

    # update the convergence flag
    converged[] = true

    return x
end

# matrix free linear analysis function
function matrixfree_lsolve!(x0, p, constants, residual!, jacobian!; coupled_jacobian=matrixfree_jacobian)

    # augment the parameter vector
    paug = isnothing(p) ? x0 : vcat(x0, p)

    # contruct linear solver
    f_lsolve = matrixfree_lsolve_lsolve(paug, constants, jacobian!)

    # define solve function
    f_solve = (paug, constants) -> matrixfree_lsolve_solve!(paug, constants, residual!, f_lsolve)

    # define residual function
    f_residual = (resid, dx, paug, constants) -> matrixfree_lsolve_residual!(resid, dx, paug, constants, residual!)

    # define drdy
    f_drdy = (residual!, dx, paug, constants) = matrixfree_lsolve_drdy(residual!, dx, paug, constants)

    # update the state
    x = x0 .- ImplicitAD(f_solve, f_residual, paug, constants; drdy = f_drdy, lsolve=f_lsolve)

    return x
end

# function mylinsolve!(A, b)
#     prob = LinearProblem(A, b)
#     return solve(prob).u
# end

# nonlinear analysis function
function nlsolve!(p, constants, residual!, jacobian!)

    # unpack pre-allocated storage and the convergence flag
    @unpack x, resid, jacob, converged = constants

    # unpack nonlinear solver parameters
    @unpack show_trace, method, linesearch, ftol, iterations = constants

    # @show p

    # wrap the residual
    f!(resid, x) = residual!(resid, x, p, constants)

    # wrap the jacobian
    j!(jacob, x) = jacobian!(jacob, x, p, constants)

    # @show x[1] #always zero
    # @show x

    # construct temporary storage
    df = NLsolve.OnceDifferentiable(f!, j!, x, resid, jacob)

    # @show df.F[1] #Always zero. 
    # @show df.DF[1,1] #Always zero.


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

    # @show x[1] #Very close but slightly different. 

    # if !result.f_converged
    #     # @show fieldnames(typeof(result))

    #     #(:method, :initial_x, :zero, :residual_norm, :iterations, :x_converged, :xtol, :f_converged, :ftol, :trace, :f_calls, :g_calls)
    #     # @show result.zero
    #     # @show result.f_converged
    #     # @show result.x_converged #-> Implies that the norm difference between two successive iterations is zero. (NLsolve says this is convergence. )
    #     # @show result.residual_norm, result.xtol
    # end

    # return the result
    return result.zero
end

function matrixfree_lsolve_lsolve(paug, constants, jacobian!)

    # unpack pre-allocated storage and the convergence flag
    @unpack x, jacob, converged = constants

    # unpack nonlinear solver parameters
    @unpack show_trace, ftol, iterations = constants

    # extract the primal values of the parameter vector
    paug = nondual_value.(paug)

    # split the parameter vector
    x0 = view(paug, 1:length(x)) # initialization terms
    p = view(p, length(x)+1:length(paug)) # provided parameter vector

    # update the structural jacobian
    jacobian!(jacob, x0, p, constants)

    # use structural jacobian factorization as a preconditioner
    Pl = lu(jacob)

    return (A, b) -> IterativeSolvers.gmres(A, b; Pl=Pl,
        abstol=ftol/10, reltol=ftol/10, maxiter=iterations, verbose=show_trace)
end

function matrixfree_lsolve_solve!(paug, constants, residual!, f_lsolve)

    # unpack pre-allocated storage and the convergence flag
    @unpack x, resid, converged = constants

    # split the parameter vector
    x0 = view(paug, 1:length(x)) # initialization terms
    p = view(p, length(x)+1:length(paug)) # provided parameter vector

    # update the residual
    residual!(resid, x0, p, constants)

    # construct the coupled jacobian
    coupled_jacob = coupled_jacobian(residual!, x0, p, constants)

    # solve the system
    dx = f_lsolve(coupled_jacob, resid)

    # update the state
    x .= x0 - dx

    # update the convergence flag
    converged[] = true

    # return the result
    return dx
end

function matrixfree_lsolve_residual!(resid, dx, paug, constants, residual!)
    x0 = view(paug, 1:length(dx)) # initialization terms
    p = view(p, length(dx)+1:length(paug)) # provided parameter vector
    residual!(resid, x0, p, constants) # right hand side
    jacob = coupled_jacobian(residual!, x0, p, constants) # jacobian vector product operator
    mul!(resid, jacob, dx, 1, -1) # linear solution residual A*x - b = 0
end

function matrixfree_lsolve_drdy(residual!, dx, paug, constants)
    x0 = view(paug, 1:length(dx)) # initialization terms
    p = view(p, length(dx)+1:length(paug)) # provided parameter vector
    return coupled_jacobian(residual!, x0, p, constants) # jacobian vector product operator
end

# nonlinear analysis function
function matrixfree_nlsolve!(p, constants, residual!, jacobian!; coupled_jacobian=matrixfree_jacobian)

    # unpack pre-allocated storage and the convergence flag
    @unpack x, resid, jacob, converged = constants

    # unpack nonlinear solver parameters
    @unpack show_trace, linesearch, ftol, iterations = constants

    if show_trace
        println("Iter     f(x) inf-norm")
        println("------   --------------")
    end

    # initialize temporary storage
    xtmp = copy(x)
    dx = copy(x)

    # construct coupled jacobian
    coupled_jacob = coupled_jacobian(residual!, x, p, constants) # compute coupled_jacobian

    # perform newton-raphson iteration
    for iter = 1:iterations

        # check residual convergence
        residual!(resid, x, p, constants)
        rnorm = norm(resid, Inf)
        show_trace && @printf("%6d%14e\n", iter, rnorm)
        converged[] = rnorm < ftol
        converged[] && break # exit if converged

        # use structural jacobian factorization as a preconditioner
        jacobian!(jacob, x, p, constants)
        Pl = lu(jacob)

        # get proposed step (Newton's Method)
        dx .= 0
        IterativeSolvers.gmres!(dx, coupled_jacob, resid; initially_zero=true, Pl=Pl,
            abstol=ftol/10, reltol=ftol/10, maxiter=iterations)
        rmul!(dx, -1)

        # initial line search objective and derivative
        ϕ0 = dot(resid, resid) / 2
        dϕ0 = dot(resid, mul!(xtmp, jacob, dx))

        # line search objective function
        function ϕ(α)
            xtmp .= x .+ α .* dx # proposed state variables
            residual!(resid, xtmp, p, constants) # associated residual
            dot(resid, resid) / 2 # objective function
        end

        # line search directional derivative
        function dϕ(α)
            xtmp .= x .+ α .* dx # proposed state variables
            residual!(resid, xtmp, p, constants) # associated residual
            jacob = coupled_jacobian(residual!, xtmp, p, constants) # associated jacobian
            dot(resid, mul!(xtmp, jacob, dx)) # directional derivative
        end

        # line search objective and directional derivative
        function ϕdϕ(α)
            xtmp .= x .+ α .* dx # proposed state variables
            residual!(resid, xtmp, p, constants) # associated residual
            jacob = coupled_jacobian(residual!, xtmp, p, constants) # associated jacobian
            dot(resid, resid) / 2, dot(resid, mul!(xtmp, jacob, dx))  # objective and gradient
        end

        # perform line search
        α, ϕα = linesearch(ϕ, dϕ, ϕdϕ, 1.0, ϕ0, dϕ0)

        # apply step
        x .+= α .* dx

    end

    # return the result
    return x
end

# mass matrix function
function mass_matrix!(jacob, x, p, constants)

    # unpack (default) parameters, parameter function, and current time
    @unpack assembly, prescribed_conditions, point_masses, xpfunc, pfunc, t = constants

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
    @unpack assembly, prescribed_conditions, point_masses, xpfunc, pfunc, t = constants

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

# sparsity pattern detection
function jacobian_colors(residual!, x, p, constants)
    resid = similar(x)
    config = ForwardDiff.JacobianConfig(residual!, resid, x)
    J1 = ForwardDiff.jacobian(residual!, resid, x1, config)
    J2 = ForwardDiff.jacobian(residual!, resid, x2, config)
    J3 = ForwardDiff.jacobian(residual!, resid, x3, config)
    @. jacob = abs(J1) + abs(J2) + abs(J3)
    return SparseDiffTools.matrix_colors(jacob)
end

# automatic differentiation jacobian construction
function autodiff_jacobian!(jacob, residual!, x, p, constants; colors=1:length(x))

    f = (r, x) -> residual!(r, x, p, constants)

    return SparseDiffTools.forwarddiff_color_jacobian!(jacob, f, x, colorvec = colors)
end

# matrix-free jacobian construction
function matrixfree_jacobian(residual!, x, p, constants)
    return SparseDiffTools.JacVec((resid, x)->residual!(resid, x, p, constants), x)
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
nondual_value(x) = x
nondual_value(x::ForwardDiff.Dual) = ForwardDiff.value(x)
nondual_value(x::ReverseDiff.TrackedReal) = ReverseDiff.value(x)