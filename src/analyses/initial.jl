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
 - `linear_displacement = zeros(3)`: Initial linear displacement of the body frame.  
 - `angular_displacement = zeros(3)`: Initial angular displacement of the body frame. 
        (using Wiener Milenkovic parameters).
 - `linear_velocity = zeros(3)`: Initial linear velocity of the body frame. 
 - `angular_velocity = zeros(3)`: Initial angular velocity of the body frame. 
 - `linear_acceleration = zeros(3)`: Initial linear acceleration of the body frame.
 - `angular_acceleration = zeros(3)`: Initial angular acceleration of the body frame.
 - `gravity = [0,0,0]`: Gravity vector in the inertial frame.
     
# Initial Condition Keyword Arguments
 - `u0 = fill(zeros(3), length(assembly.points))`: Initial linear displacement of 
        each point relative to the body frame
 - `theta0 = fill(zeros(3), length(assembly.points))`: Initial angular displacement of 
        each point relative to the body frame (using Wiener-Milenkovic Parameters) 
 - `V0 = fill(zeros(3), length(assembly.points))`: Initial linear velocity of 
        each point relative to the body frame
 - `Omega0 = fill(zeros(3), length(assembly.points))`: Initial angular velocity of 
        each point relative to the body frame
 - `Vdot0 = fill(zeros(3), length(assembly.points))`: Initial linear acceleration of 
        each point relative to the body frame
 - `Omegadot0 = fill(zeros(3), length(assembly.points))`: Initial angular acceleration of 
        each point relative to the body frame

# Control Flag Keyword Arguments
 - `structural_damping = false`: Indicates whether to enable structural damping
 - `constant_mass_matrix = false`: Indicates whether to use a constant mass matrix system
 - `reset_state = true`: Flag indicating whether the system state variables should be 
        set to zero prior to performing this analysis.
 - `steady_state=false`: Flag indicating whether to initialize by performing a steady state 
        analysis.
 - `linear = false`: Flag indicating whether a linear analysis should be performed.
 - `show_trace = false`: Flag indicating whether to display the solution progress.

# Nonlinear Solver Keyword Arguments
 - `method = :newton`: Method (as defined in NLsolve) to solve nonlinear system of equations
 - `linesearch = LineSearches.BackTracking(maxstep=1e6)`: Line search used to solve the
        nonlinear system of equations
 - `ftol = 1e-9`: tolerance for solving the nonlinear system of equations
 - `iterations = 1000`: maximum iterations for solving the nonlinear system of equations

# Linear Solver Keyword Arguments
 - `linearization_state`: Linearization state variables.  Defaults to zeros.
 - `update_linearization`: Flag indicating whether to update the linearization state 
        variables for a linear analysis with the instantaneous state variables.
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
    linear_displacement=(@SVector zeros(3)),
    angular_displacement=(@SVector zeros(3)),
    linear_velocity=(@SVector zeros(3)),
    angular_velocity=(@SVector zeros(3)),
    linear_acceleration=(@SVector zeros(3)),
    angular_acceleration=(@SVector zeros(3)),
    gravity=(@SVector zeros(3)),
    # initial condition keyword arguments
    u0=fill((@SVector zeros(3)), length(assembly.points)),
    theta0=fill((@SVector zeros(3)), length(assembly.points)),
    V0=fill((@SVector zeros(3)), length(assembly.points)),
    Omega0=fill((@SVector zeros(3)), length(assembly.points)),
    Vdot0=fill((@SVector zeros(3)), length(assembly.points)),
    Omegadot0=fill((@SVector zeros(3)), length(assembly.points)),
    # control flag keyword arguments
    structural_damping=true,
    constant_mass_matrix=typeof(system)<:ExpandedSystem,
    reset_state=true,
    steady_state=false,
    linear=false,
    show_trace=false,
    # nonlinear solver keyword arguments
    method=:newton,
    linesearch=LineSearches.BackTracking(maxstep=1e6),
    ftol=1e-9,
    iterations=1000,
    # linear solver keyword arguments
    linearization_state=nothing,
    update_linearization=false,
    )

    # perform steady state analysis (if requested)
    if steady_state
        return steady_state_analysis!(system, assembly;
            # general keyword arguments
            prescribed_conditions=prescribed_conditions,
            distributed_loads=distributed_loads,
            point_masses=point_masses,
            linear_displacement=linear_displacement,
            angular_displacement=angular_displacement,
            linear_velocity=linear_velocity,
            angular_velocity=angular_velocity,
            linear_acceleration=linear_acceleration,
            angular_acceleration=angular_acceleration,
            gravity=gravity,
            time=t0,
            # control flag keyword arguments
            structural_damping=structural_damping,
            constant_mass_matrix=constant_mass_matrix,
            reset_state=reset_state,
            linear=linear,
            show_trace=show_trace,
            # nonlinear solver keyword arguments
            method=method,
            linesearch=linesearch,
            ftol=ftol,
            iterations=iterations,
            # linear solver keyword arguments
            linearization_state=linearization_state,
            update_linearization=update_linearization,
            )
    end

    # check if provided system is consistent with provided keyword arguments
    constant_mass_matrix && @assert typeof(system) <: ExpandedSystem

    # reset state, if specified
    if reset_state
        reset_state!(system)
    end

    # save provided system
    original_system = system

    # convert constant mass matrix system to a dynamic system
    if constant_mass_matrix
        # construct new system
        system = DynamicSystem(assembly; force_scaling = system.force_scaling)
        # copy state variables from the original system to the new system
        copy_state!(system, original_system, assembly;     
            prescribed_conditions=prescribed_conditions,
            distributed_loads=distributed_loads,
            point_masses=point_masses,
            linear_displacement=linear_displacement,
            angular_displacement=angular_displacement,
            linear_velocity=linear_velocity,
            angular_velocity=angular_velocity,
            linear_acceleration=linear_acceleration,
            angular_acceleration=angular_acceleration,
            gravity=gravity,
            time=t0)
    end

    # unpack pre-allocated storage and pointers
    @unpack x, r, K, M, force_scaling, indices, udot, θdot, Vdot, Ωdot = system

    # get the current time
    t0 = first(t0)

    # print the current time
    if show_trace
        println("Solving for t=$(t0)")
    end
    
    # update stored time
    system.t = t0

    # current parameters
    pcond = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(t0)
    dload = typeof(distributed_loads) <: AbstractDict ? distributed_loads : distributed_loads(t0)
    pmass = typeof(point_masses) <: AbstractDict ? point_masses : point_masses(t0)
    gvec = typeof(gravity) <: AbstractVector ? SVector{3}(gravity) : SVector{3}(gravity(t0))
    ub_p = typeof(linear_displacement) <: AbstractVector ? SVector{3}(linear_displacement) : SVector{3}(linear_displacement(t0))
    θb_p = typeof(angular_displacement) <: AbstractVector ? SVector{3}(angular_displacement) : SVector{3}(angular_displacement(t0))
    vb_p = typeof(linear_velocity) <: AbstractVector ? SVector{3}(linear_velocity) : SVector{3}(linear_velocity(t0))
    ωb_p = typeof(angular_velocity) <: AbstractVector ? SVector{3}(angular_velocity) : SVector{3}(angular_velocity(t0))
    ab_p = typeof(linear_acceleration) <: AbstractVector ? SVector{3}(linear_acceleration) : SVector{3}(linear_acceleration(t0))
    αb_p = typeof(angular_acceleration) <: AbstractVector ? SVector{3}(angular_acceleration) : SVector{3}(angular_acceleration(t0))

    # # --- Determine whether all rigid body modes are constrained --- #

    # # save original prescribed conditions
    # original_pcond = pcond

    # # find unconstrained rigid body modes
    # not_constrained = (reduce(.&, [.!p.pd for (i,p) in pcond]))

    # # constrain unconstrained rigid body modes
    # if any(not_constrained)
    #     # create copy of the original prescribed conditions
    #     pcond = copy(pcond)

    #     # unpack prescribed condition fields
    #     if haskey(pcond, 1)
    #         @unpack pd, pl, u, theta, F, M, Ff, Mf = pcond[1]
    #     else
    #         pd = @SVector zeros(Bool, 6)
    #         pl = @SVector ones(Bool, 6)
    #         u = @SVector zeros(Float64, 3)
    #         theta = @SVector zeros(Float64, 3)
    #         F = @SVector zeros(Float64, 3)
    #         M = @SVector zeros(Float64, 3)
    #         Ff = @SVector zeros(Float64, 3)
    #         Mf = @SVector zeros(Float64, 3)
    #     end

    #     # add constraints for unconstrained rigid body modes
    #     for i = 1:3
    #         # prescribe linear displacement, introduce linear acceleration state variable
    #         if not_constrained[i]
    #             pd = setindex(pd, true, i)
    #             u = setindex(u, u0[1][i], i)
    #         end
    #         # prescribe angular displacement, introduce angular acceleration state variable
    #         if not_constrained[3+i]
    #             pd = setindex(pd, true, 3+i)
    #             theta = setindex(theta, theta0[1][i], i)
    #         end
    #     end

    #     # overwrite original prescribed condition
    #     pcond[1] = PrescribedConditions(pd, pl, u, theta, F, M, Ff, Mf)
    # end

    # not_constrained = (reduce(.&, [.!p.pd for (i,p) in pcond]))
    # println(not_constrained)

    # indices corresponding to rigid body acceleration state variables
    icol_accel = body_frame_acceleration_indices(system, pcond)

    # --- Determine whether Vdot and Ωdot may be found using the equilibrium equations --- #

    # NOTE: Our system of equations cannot be solved for Vdot and Ωdot if the corresponding 
    # rows and columns of the mass matrix are zero.

    # Fortunately, though we cannot solve for these variables, their values are not used.

    # construct `rate_vars` vector
    dynamic_system_mass_matrix!(system.M, x, indices, force_scaling, assembly, pcond, pmass)
    rate_vars1 = .!(iszero.(sum(system.M, dims=1)))

    # --- Determine whether Fi and Mi may be found using the compatability equations --- #

    # NOTE: The compatability equations cannot be solved for Fi and Mi if the corresponding 
    # rows and columns of the compliance matrix are zero.
    
    # Values for these variables must be derived from the equilibrium equations.  Since the 
    # equilibrium equations cannot then be used to find Vdot and Ωdot, compatible values 
    # for Vdot and Ωdot must be provided.
    
    # replace mass matrix with compliance matrix
    elements = copy(assembly.elements)
    for (i, e) in enumerate(assembly.elements)
        assembly.elements[i] = Element(e.L, e.x, e.compliance, e.compliance, e.Cab, e.mu)
    end

    # temporarily ignore the mass of all point masses
    pmass2 = Dict{Int,PointMass{Float64}}()

    # construct `rate_vars2` vector to test if the second case applies
    dynamic_system_mass_matrix!(system.M, x, indices, force_scaling, assembly, pcond, pmass2)
    rate_vars2 = .!(iszero.(sum(system.M, dims=1)))

    # restore original element mass matrices
    for i in eachindex(assembly.elements)
        assembly.elements[i] = elements[i]
    end

    # --- Define the residual and jacobian functions --- #

    f! = (resid, x) -> initial_condition_system_residual!(resid, x, indices, rate_vars1, rate_vars2, 
        icol_accel, force_scaling, structural_damping, assembly, pcond, dload, pmass, 
        gvec, ub_p, θb_p, vb_p, ωb_p, ab_p, αb_p, u0, theta0, V0, Omega0, Vdot0, Omegadot0)

    j! = (jacob, x) -> initial_condition_system_jacobian!(jacob, x, indices, rate_vars1, rate_vars2,
        icol_accel, force_scaling, structural_damping, assembly, pcond, dload, pmass, 
        gvec, ub_p, θb_p, vb_p, ωb_p, ab_p, αb_p, u0, theta0, V0, Omega0, Vdot0, Omegadot0)

    # --- Solve for the corresponding state variables --- #

    if linear
        # perform a linear analysis
        if isnothing(linearization_state)
            x .= 0
        else
            x .= linearization_state
        end
        f!(r, x)
        j!(K, x)

        # update the solution vector            
        x .-= safe_lu(K) \ r
    
        # set the convergence flag
        converged = true
    else
        # perform a nonlinear analysis
        df = OnceDifferentiable(f!, j!, x, r, collect(K))

        result = NLsolve.nlsolve(#f!, x, autodiff=:forward,
            df, x,
            show_trace=show_trace,
            linsolve=(x, A, b) -> ldiv!(x, safe_lu(A), b),
            method=method,
            linesearch=linesearch,
            ftol=ftol,
            iterations=iterations)

        # update the solution vector and jacobian
        x .= result.zero
        K .= df.DF

        # set the convergence flag
        converged = result.f_converged
    end

    # --- Save state and rate variables associated with each point --- #

    for ipoint in eachindex(assembly.points)
        # external loads
        Fe, Me = point_loads(x, ipoint, indices.icol_point, force_scaling, pcond)
        # linear and angular displacement
        u, θ = initial_point_displacement(x, ipoint, indices.icol_point, pcond, u0, theta0, 
            rate_vars2)
        # linear and angular velocity
        V, Ω = V0[ipoint], Omega0[ipoint]
        # linear and angular displacement rates
        udot[ipoint], θdot[ipoint] = point_velocities(x, ipoint, indices.icol_point)
        # linear and angular velocity rates
        Vdot[ipoint], Ωdot[ipoint] = initial_point_velocity_rates(x, ipoint, 
            indices.icol_point, pcond, Vdot0, Omegadot0, rate_vars2)
        # add rigid body accelerations to node accelerations
        Δx = assembly.points[ipoint]
        a, α = body_frame_acceleration(x)
        Vdot[ipoint] += a + cross(α, Δx) + cross(α, u)
        Ωdot[ipoint] += α
        # set new state variables
        set_external_forces!(system, prescribed_conditions, Fe, ipoint)
        set_external_moments!(system, prescribed_conditions, Me, ipoint)
        set_linear_displacement!(system, prescribed_conditions, u, ipoint)
        set_angular_displacement!(system, prescribed_conditions, θ, ipoint)
        set_linear_velocity!(system, V, ipoint)
        set_angular_velocity!(system, Ω, ipoint)
    end

    # --- Restore constant mass matrix system (if applicable) --- #

    if !(typeof(original_system) <: DynamicSystem)
        # copy state variables from the dynamic system to the original system
        copy_state!(original_system, system, assembly;     
            prescribed_conditions=prescribed_conditions,
            distributed_loads=distributed_loads,
            point_masses=point_masses,
            linear_displacement=linear_displacement,
            angular_displacement=angular_displacement,
            linear_velocity=linear_velocity,
            angular_velocity=angular_velocity,
            linear_acceleration=linear_acceleration,
            angular_acceleration=angular_acceleration,
            gravity=gravity,
            time=t0)
    end

    return original_system, converged
end
