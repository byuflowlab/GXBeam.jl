"""
    static_analysis(assembly; kwargs...)

Perform a static analysis of the system of nonlinear beams contained in
`assembly`. Return the resulting system and a flag indicating whether the
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
 - `gravity = [0,0,0]`: Gravity vector.  If time varying, this input may be provided as a 
        function of time.
 - `time = 0.0`: Current time or vector of times corresponding to each step. May be used 
        in conjunction with time varying prescribed conditions, distributed loads, and 
        body frame motion to gradually increase displacements and loads.     

# Control Flag Keyword Arguments
 - `reset_state = true`: Flag indicating whether the system state variables should be 
        set to zero prior to performing this analysis.
 - `linear = false`: Flag indicating whether a linear analysis should be performed.
 - `show_trace = false`: Flag indicating whether to display the solution progress.

# Nonlinear Solver Keyword Arguments (see [NLsolve.jl](https://github.com/JuliaNLSolvers/NLsolve.jl))
 - `method = :newton`: Solution method for nonlinear systems of equations
 - `linesearch = LineSearches.BackTracking(maxstep=1e6)`: Line search for solving nonlinear
        systems of equations
 - `ftol = 1e-9`: Tolerance for solving the nonlinear system of equations
 - `iterations = 1000`: Iteration limit when solving the nonlinear systems of equations

# Linear Solver Keyword Arguments
 - `linearization_state`: Linearization state variables.  Defaults to zeros.
 - `update_linearization`: Flag indicating whether to update the linearization state 
        variables for a linear analysis with the instantaneous state variables.
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

        # update the stored time
        system.t = t

        # current parameters
        pcond = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(t)
        dload = typeof(distributed_loads) <: AbstractDict ? distributed_loads : distributed_loads(t)
        pmass = typeof(point_masses) <: AbstractDict ? point_masses : point_masses(t)
        gvec = typeof(gravity) <: AbstractVector ? SVector{3}(gravity) : SVector{3}(gravity(t))

        # define the residual and jacobian functions
        f! = (resid, x) -> static_system_residual!(resid, x, indices, force_scaling, 
            assembly, pcond, dload, pmass, gvec)

        j! = (jacob, x) -> static_system_jacobian!(jacob, x, indices, force_scaling, 
            assembly, pcond, dload, pmass, gvec)

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
    end

    return system, converged
end
