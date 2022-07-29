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
 - `linear_displacement = zeros(3)`: Initial linear displacement of the body frame.  
 - `angular_displacement = zeros(3)`: Initial angular displacement of the body frame. 
        (using Wiener Milenkovic parameters).
 - `linear_velocity = zeros(3)`: Initial linear velocity of the body frame. 
 - `angular_velocity = zeros(3)`: Initial angular velocity of the body frame. 
 - `linear_acceleration = zeros(3)`: Prescribed linear acceleration of the body frame.
 - `angular_acceleration = zeros(3)`: Prescribed angular acceleration of the body frame.
 - `gravity = [0,0,0]`: Gravity vector in the inertial frame.
 - `save = 1:length(time)`: Steps at which to save the time history

 # Initial Condition Keyword Arguments
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

# Control Flag Keyword Arguments
 - `structural_damping = false`: Indicates whether to enable structural damping
 - `reset_state = true`: Flag indicating whether the system state variables should be 
        set to zero prior to performing this analysis.
 - `initialize = true`: Flag indicating whether a consistent set of initial
        conditions should be found using [`initial_condition_analysis`](@ref).
 - `steady_state=false`: Flag indicating whether to initialize by performing a steady state 
        analysis.
 - `linear = false`: Flag indicating whether a linear analysis should be performed.
 - `show_trace = false`: Flag indicating whether to display the solution progress.

 # Nonlinear Solver Keyword Arguments
 - `method = :newton`: Method (as defined in NLsolve) to solve nonlinear system of equations
 - `linesearch = LineSearches.BackTracking(maxstep=1e6)`: Line search used to solve 
        nonlinear systems of equations
 - `ftol = 1e-9`: tolerance for solving the nonlinear system of equations
 - `iterations = 1000`: maximum iterations for solving the nonlinear system of equations

# Linear Solver Keyword Arguments
 - `linearization_state`: Linearization state variables.  Defaults to zeros.
 - `update_linearization`: Flag indicating whether to update the linearization state 
        variables for a linear analysis with the instantaneous state variables.
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
    linear_displacement=(@SVector zeros(3)),
    angular_displacement=(@SVector zeros(3)),
    linear_velocity=(@SVector zeros(3)),
    angular_velocity=(@SVector zeros(3)),
    linear_acceleration=(@SVector zeros(3)),
    angular_acceleration=(@SVector zeros(3)),
    gravity=(@SVector zeros(3)), 
    save=1:length(tvec),
    # initial condition keyword arguments
    u0=fill((@SVector zeros(3)), length(assembly.points)),
    theta0=fill((@SVector zeros(3)), length(assembly.points)),
    V0=fill((@SVector zeros(3)), length(assembly.points)),
    Omega0=fill((@SVector zeros(3)), length(assembly.points)),
    Vdot0=fill((@SVector zeros(3)), length(assembly.points)),
    Omegadot0=fill((@SVector zeros(3)), length(assembly.points)),
    # control flag keyword arguments
    structural_damping=true,
    reset_state=true,
    initialize=true,
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

    # reset state, if specified
    if reset_state
        reset_state!(system)
    end

    # perform initial condition analysis
    if initialize
        # perform initialization procedure
        system, converged = initial_condition_analysis!(system, assembly, tvec[1];
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
            u0=u0,
            theta0=theta0,
            V0=V0,
            Omega0=Omega0,
            Vdot0=Vdot0,
            Omegadot0=Omegadot0,
            structural_damping=structural_damping,
            reset_state=reset_state,
            linear=linear,
            steady_state=steady_state,
            show_trace=show_trace,
            method=method,
            linesearch=linesearch,
            ftol=ftol,
            iterations=iterations,
            linearization_state=linearization_state,
            )
    else
        # converged by default
        converged = true
    end

    # unpack pre-allocated storage and pointers
    @unpack x, r, K, force_scaling, indices, udot, θdot, Vdot, Ωdot = system

    # initialize storage for each time step
    history = Vector{AssemblyState{eltype(system)}}(undef, length(save))
    isave = 1

    # add initial state to the solution history
    if isave in save
        pcond = typeof(prescribed_conditions) <: AbstractDict ?
            prescribed_conditions : prescribed_conditions(tvec[1])
        history[isave] = AssemblyState(system, assembly, prescribed_conditions=pcond)
        isave += 1
    end

    # begin time stepping
    for it = 2:length(tvec)

        # print the current time
        if show_trace
            println("Solving for t=$(tvec[it])")
        end

        # update stored time
        system.t = tvec[it]

        # current time step size
        dt = tvec[it] - tvec[it-1]

        # current parameters
        pcond = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(tvec[it])
        dload = typeof(distributed_loads) <: AbstractDict ? distributed_loads : distributed_loads(tvec[it])
        pmass = typeof(point_masses) <: AbstractDict ? point_masses : point_masses(tvec[it])
        gvec = typeof(gravity) <: AbstractVector ? SVector{3}(gravity) : SVector{3}(gravity(tvec[it]))
        ab_p = typeof(linear_acceleration) <: AbstractVector ? SVector{3}(linear_acceleration) : SVector{3}(linear_acceleration(tvec[it]))
        αb_p = typeof(angular_acceleration) <: AbstractVector ? SVector{3}(angular_acceleration) : SVector{3}(angular_acceleration(tvec[it]))

        # indices corresponding to rigid body acceleration state variables
        icol_accel = body_frame_acceleration_indices(system, pcond)

        # set current state rate initialization parameters

        # state variables for the body frame motion
        ub, θb = body_frame_displacement(x)
        vb, ωb = body_frame_velocity(x)
        ab, αb = body_frame_acceleration(x)
        # initialization terms for the body frame motion
        ubdot = 2/dt*ub + vb
        θbdot = 2/dt*θb + get_Qinv(θb)*get_C(θb)*ωb
        vbdot = 2/dt*vb + ab
        ωbdot = 2/dt*ωb + αb

        # initialization terms for each point
        for ipoint = 1:length(assembly.points)
            # state variables for this point
            u, θ = point_displacement(x, ipoint, indices.icol_point, pcond)
            V, Ω = point_velocities(x, ipoint, indices.icol_point)
            # initialization terms for this point, use storage for state rates
            udot[ipoint] = 2/dt*u + udot[ipoint]
            θdot[ipoint] = 2/dt*θ + θdot[ipoint]
            Vdot[ipoint] = 2/dt*V + Vdot[ipoint]
            Ωdot[ipoint] = 2/dt*Ω + Ωdot[ipoint]
        end

        # define the residual and jacobian functions
        f! = (r, x) -> newmark_system_residual!(r, x, indices, icol_accel, force_scaling, 
            structural_damping, assembly, pcond, dload, pmass, gvec, ab_p, αb_p, 
            ubdot, θbdot, vbdot, ωbdot, udot, θdot, Vdot, Ωdot, dt)

        j! = (K, x) -> newmark_system_jacobian!(K, x, indices, icol_accel, force_scaling, 
            structural_damping, assembly, pcond, dload, pmass, gvec, ab_p, αb_p, 
            ubdot, θbdot, vbdot, ωbdot, udot, θdot, Vdot, Ωdot, dt)

        # solve for the new set of state variables
        if linear
            # perform a linear analysis
            if !update_linearization
                if isnothing(linearization_state)
                    x .= 0
                else
                    x .= linearization_state
                end
            end
            f!(r, x)
            j!(K, x)

            # update the solution vector            
            x .-= safe_lu(K) \ r
        else
            # perform a nonlinear analysis
            df = OnceDifferentiable(f!, j!, x, r, K)

            result = NLsolve.nlsolve(df, x,
                show_trace=show_trace,
                linsolve=(x, A, b) -> ldiv!(x, safe_lu(A), b),
                method=method,
                linesearch=linesearch,
                ftol=ftol,
                iterations=iterations)

            # update the solution vector and jacobian
            x .= result.zero
            K .= df.DF

            # update the convergence flag
            converged = result.f_converged
        end

        # set new state rates
        for ipoint = 1:length(assembly.points)
            u, θ = point_displacement(x, ipoint, indices.icol_point, pcond)
            V, Ω = point_velocities(x, ipoint, indices.icol_point)
            udot[ipoint] = 2/dt*u - udot[ipoint]
            θdot[ipoint] = 2/dt*θ - θdot[ipoint]
            Vdot[ipoint] = 2/dt*V - Vdot[ipoint]
            Ωdot[ipoint] = 2/dt*Ω - Ωdot[ipoint]
        end

        # add state to history
        if it in save
            history[isave] = AssemblyState(system, assembly, prescribed_conditions=pcond)
            isave += 1
        end

        # stop early if unconverged
        if !converged
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

    return system, history, converged
end