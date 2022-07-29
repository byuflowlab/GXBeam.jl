
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
 - `linear_displacement = zeros(3)`: Prescribed linear displacement of the body frame.  
        If time varying, this input may be provided as a function of time.
 - `angular_displacement = zeros(3)`: Prescribed angular displacement of the body frame. 
        (using Wiener Milenkovic parameters). If time varying, this input may be provided 
        as a function of time.
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
 - `structural_damping = false`: Indicates whether to enable structural damping
 - `constant_mass_matrix = false`: Indicates whether to use a constant mass matrix system
 - `reset_state = true`: Flag indicating whether the system state variables should be 
        set to zero prior to performing this analysis.
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
    linear_displacement=(@SVector zeros(3)),
    angular_displacement=(@SVector zeros(3)),
    linear_velocity=(@SVector zeros(3)),
    angular_velocity=(@SVector zeros(3)),
    linear_acceleration=(@SVector zeros(3)),
    angular_acceleration=(@SVector zeros(3)),
    gravity=(@SVector zeros(3)),
    time=0.0,
    # control flag keyword arguments
    structural_damping=false,
    constant_mass_matrix=typeof(system)<:ExpandedSystem,
    reset_state=true,
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

    # check if provided system is consistent with provided keyword arguments
    constant_mass_matrix && @assert typeof(system) <: ExpandedSystem

    # reset state, if specified
    if reset_state
        reset_state!(system)
    end

    # unpack pre-allocated storage and pointers
    @unpack x, r, K, force_scaling, indices = system

    # assume converged until proven otherwise
    converged = true

    # begin time stepping
    for t in time

        # print the current time
        if show_trace
            println("Solving for t=$t")
        end

        # update stored time
        system.t = t

        # current parameters
        pcond = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(t)
        dload = typeof(distributed_loads) <: AbstractDict ? distributed_loads : distributed_loads(t)
        pmass = typeof(point_masses) <: AbstractDict ? point_masses : point_masses(t)
        gvec = typeof(gravity) <: AbstractVector ? SVector{3}(gravity) : SVector{3}(gravity(t))
        ub_p = typeof(linear_displacement) <: AbstractVector ? SVector{3}(linear_displacement) : SVector{3}(linear_displacement(t))
        θb_p = typeof(angular_displacement) <: AbstractVector ? SVector{3}(angular_displacement) : SVector{3}(angular_displacement(t))
        vb_p = typeof(linear_velocity) <: AbstractVector ? SVector{3}(linear_velocity) : SVector{3}(linear_velocity(t))
        ωb_p = typeof(angular_velocity) <: AbstractVector ? SVector{3}(angular_velocity) : SVector{3}(angular_velocity(t))
        ab_p = typeof(linear_acceleration) <: AbstractVector ? SVector{3}(linear_acceleration) : SVector{3}(linear_acceleration(t))
        αb_p = typeof(angular_acceleration) <: AbstractVector ? SVector{3}(angular_acceleration) : SVector{3}(angular_acceleration(t))

        # indices corresponding to rigid body acceleration state variables
        icol_accel = body_frame_acceleration_indices(system, pcond)

        # define the residual and jacobian functions
        if constant_mass_matrix
            f! = (resid, x) -> expanded_steady_system_residual!(resid, x, indices, icol_accel, 
                force_scaling, structural_damping, assembly, pcond, dload, pmass, gvec, 
                ub_p, θb_p, vb_p, ωb_p, ab_p, αb_p)
            j! = (jacob, x) -> expanded_steady_system_jacobian!(jacob, x, indices, icol_accel, 
                force_scaling, structural_damping, assembly, pcond, dload, pmass, gvec, 
                ub_p, θb_p, vb_p, ωb_p, ab_p, αb_p)
        else
            f! = (resid, x) -> steady_state_system_residual!(resid, x, indices, icol_accel, 
                force_scaling, structural_damping, assembly, pcond, dload, pmass, gvec, 
                ub_p, θb_p, vb_p, ωb_p, ab_p, αb_p)
            j! = (jacob, x) -> steady_state_system_jacobian!(jacob, x, indices, icol_accel, 
                force_scaling, structural_damping, assembly, pcond, dload, pmass, gvec, 
                ub_p, θb_p, vb_p, ωb_p, ab_p, αb_p)
        end

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
            df = NLsolve.OnceDifferentiable(f!, j!, x, r, K)

            # result = NLsolve.nlsolve(df, x,
            result = NLsolve.nlsolve(f!, x, autodiff=:forward,
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

        # update state variable rates
        if !constant_mass_matrix
            @unpack udot, θdot, Vdot, Ωdot = system
            ab, αb = body_frame_acceleration(x)
            for ipoint = 1:length(assembly.points)
                Δx = assembly.points[ipoint]
                u, _ = point_displacement(x, ipoint, indices.icol_point, pcond)
                udot[ipoint] = @SVector zeros(3)
                θdot[ipoint] = @SVector zeros(3)
                Vdot[ipoint] = ab + cross(αb, Δx) + cross(αb, u)
                Ωdot[ipoint] = αb
            end 
        end
    end

    @assert system.x == x

    return system, converged
end
