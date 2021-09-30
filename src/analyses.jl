"""
    static_analysis(assembly; kwargs...)

Perform a static analysis of the system of nonlinear beams contained in
`assembly`. Return the resulting system and a flag indicating whether the
iteration procedure converged.

# Keyword Arguments
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
 - `point_masses = Dict{Int,Vector{PointMass{Float64}}}()`: A dictionary with keys 
        corresponding to the points at which point masses are attached and values 
        containing vectors of objects of type [`PointMass`](@ref) which describe 
        the point masses attached at those points.  If time varying, this input may
        be provided as a function of time.
 - `gravity = [0,0,0]`: Gravity vector.  If time varying, this input may be provided as a 
        function of time.
 - `linear = false`: Set to `true` for a linear analysis
 - `method = :newton`: Method (as defined in NLsolve) to solve nonlinear system of equations
 - `linesearch = LineSearches.BackTracking(maxstep=1e6)`: Line search used to solve nonlinear 
        system of equations
 - `ftol = 1e-9`: tolerance for solving nonlinear system of equations
 - `iterations = 1000`: maximum iterations for solving the nonlinear system of equations
 - `tvec = 0`: Time vector/value. May be used in conjunction with time varying
        prescribed conditions and distributed loads to gradually increase
        displacements/loads.
 - `reset_state = true`: Flag indicating whether the state variables should be
        reset prior to performing the analysis.  This keyword argument is only valid
        for the pre-allocated version of this function.
"""
function static_analysis(assembly;
    prescribed_conditions=Dict{Int,PrescribedConditions{Float64}}(),
    distributed_loads=Dict{Int,DistributedLoads{Float64}}(),
    point_masses=Dict{Int,Vector{PointMass{Float64}}}(),
    gravity=SVector(0,0,0),
    linear=false,
    method=:newton,
    linesearch=LineSearches.BackTracking(maxstep=1e6),
    ftol=1e-9,
    iterations=1000,
    tvec=0.0,
    )

    static = true

    pc = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(tvec[1])

    system = System(assembly, static; prescribed_points=keys(pc))

    return static_analysis!(system, assembly;
        prescribed_conditions=prescribed_conditions,
        distributed_loads=distributed_loads,
        point_masses=point_masses,
        gravity=gravity,
        linear=linear,
        method=method,
        linesearch=linesearch,
        ftol=ftol,
        iterations=iterations,
        tvec=tvec,
        reset_state=false)
end

"""
    static_analysis!(system, assembly; kwargs...)

Pre-allocated version of `static_analysis`.
"""
function static_analysis!(system, assembly;
    prescribed_conditions=Dict{Int,PrescribedConditions{Float64}}(),
    distributed_loads=Dict{Int,DistributedLoads{Float64}}(),
    point_masses=Dict{Int,Vector{PointMass{Float64}}}(),
    gravity=SVector(0,0,0),
    linear=false,
    method=:newton,
    linesearch=LineSearches.BackTracking(maxstep=1e6),
    ftol=1e-9,
    iterations=1000,
    tvec=0.0,
    reset_state=true)

    # check to make sure system is static
    @assert system.static == true

    # reset state, if specified
    if reset_state
        reset_state!(system)
    end

    # unpack pre-allocated storage and pointers
    x = system.x
    F = system.r
    J = system.K
    force_scaling = system.force_scaling
    mass_scaling = system.mass_scaling
    irow_point = system.irow_point
    irow_elem = system.irow_elem
    irow_elem1 = system.irow_elem1
    irow_elem2 = system.irow_elem2
    icol_point = system.icol_point
    icol_elem = system.icol_elem

    converged = true
    for t in tvec

        # update stored time
        system.t = t

        # current parameters
        pcond = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(t)
        dload = typeof(distributed_loads) <: AbstractDict ? distributed_loads : distributed_loads(t)
        pmass = typeof(point_masses) <: AbstractDict ? point_masses : point_masses(t)
        gvec = typeof(gravity) <: AbstractVector ? SVector{3}(gravity) : SVector{3}(gravity(tvec[it]))

        # solve the system of equations
        f! = (F, x) -> system_residual!(F, x, assembly, pcond, dload, pmass, gvec, force_scaling,
            mass_scaling, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem)

        j! = (J, x) -> system_jacobian!(J, x, assembly, pcond, dload, pmass, gvec, force_scaling,
            mass_scaling, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem)

        if linear
            # linear analysis
            x .= 0.0
            f!(F, x)
            j!(J, x)
            x .-= safe_lu(J) \ F
        else
            # nonlinear analysis
            df = NLsolve.OnceDifferentiable(f!, j!, x, F, J)

            result = NLsolve.nlsolve(df, x,
                linsolve=(x, A, b) -> ldiv!(x, safe_lu(A), b),
                method=method,
                linesearch=linesearch,
                ftol=ftol,
                iterations=iterations)

            # update the solution
            x .= result.zero
            J .= df.DF

            # update convergence flag
            converged = result.f_converged
        end
    end

    return system, converged
end

"""
    steady_state_analysis(assembly; kwargs...)

Perform a steady-state analysis for the system of nonlinear beams contained in
`assembly`.  Return the resulting system and a flag indicating whether the
iteration procedure converged.

# Keyword Arguments
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
 - `point_masses = Dict{Int,Vector{PointMass{Float64}}}()`: A dictionary with keys 
        corresponding to the points at which point masses are attached and values 
        containing vectors of objects of type [`PointMass`](@ref) which describe 
        the point masses attached at those points.  If time varying, this input may
        be provided as a function of time.
 - `gravity = [0,0,0]`: Gravity vector.  If time varying, this input may be provided as a 
        function of time.            
 - `linear = false`: Set to `true` for a linear analysis
 - `method = :newton`: Method (as defined in NLsolve) to solve nonlinear system of equations
 - `linesearch = LineSearches.LineSearches.BackTracking(maxstep=1e6)`: Line search used to solve nonlinear system of equations
 - `ftol = 1e-9`: tolerance for solving nonlinear system of equations
 - `iterations = 1000`: maximum iterations for solving the nonlinear system of equations
 - `origin = zeros(3)`: Global frame origin vector. If time varying, this input
        may be provided as a function of time.
 - `linear_velocity = zeros(3)`: Global frame linear velocity vector. If time
        varying, this vector may be provided as a function of time.
 - `angular_velocity = zeros(3)`: Global frame angular velocity vector. If time
        varying, this vector may be provided as a function of time.
 - `tvec = 0.0`: Time vector/value. May be used in conjunction with time varying
    prescribed conditions, distributed loads, and global motion to gradually
    increase displacements/loads.
 - `reset_state = true`: Flag indicating whether the state variables should be
    reset prior to performing the analysis.  This keyword argument is only valid
    for the pre-allocated version of this function.
"""
function steady_state_analysis(assembly;
    prescribed_conditions=Dict{Int,PrescribedConditions{Float64}}(),
    distributed_loads=Dict{Int,DistributedLoads{Float64}}(),
    point_masses=Dict{Int,Vector{PointMass{Float64}}}(),
    gravity=SVector(0,0,0),
    linear=false,
    method=:newton,
    linesearch=LineSearches.BackTracking(maxstep=1e6),
    ftol=1e-9,
    iterations=1000,
    origin=(@SVector zeros(3)),
    linear_velocity=(@SVector zeros(3)),
    angular_velocity=(@SVector zeros(3)),
    tvec=0.0,
    )

    static = false

    pc = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(tvec[1])

    system = System(assembly, static; prescribed_points=keys(pc))

    return steady_state_analysis!(system, assembly;
        prescribed_conditions=prescribed_conditions,
        distributed_loads=distributed_loads,
        point_masses=point_masses,
        gravity=gravity,
        linear=linear,
        method=method,
        linesearch=linesearch,
        ftol=ftol,
        iterations=iterations,
        origin=origin,
        linear_velocity=linear_velocity,
        angular_velocity=angular_velocity,
        tvec=tvec,
        reset_state=false,
        )
end

"""
    steady_state_analysis!(system, assembly; kwargs...)

Pre-allocated version of `steady_state_analysis`.
"""
function steady_state_analysis!(system, assembly;
    prescribed_conditions=Dict{Int,PrescribedConditions{Float64}}(),
    distributed_loads=Dict{Int,DistributedLoads{Float64}}(),
    point_masses=Dict{Int,Vector{PointMass{Float64}}}(),
    gravity=SVector(0,0,0),
    linear=false,
    method=:newton,
    linesearch=LineSearches.BackTracking(maxstep=1e6),
    ftol=1e-9,
    iterations=1000,
    origin=(@SVector zeros(3)),
    linear_velocity=(@SVector zeros(3)),
    angular_velocity=(@SVector zeros(3)),
    tvec=0.0,
    reset_state=true,
    )

    # check to make sure the simulation is dynamic
    @assert system.static == false

    # reset state, if specified
    if reset_state
        reset_state!(system)
    end

    # unpack pointers to pre-allocated storage
    x = system.x
    F = system.r
    J = system.K
    force_scaling = system.force_scaling
    mass_scaling = system.mass_scaling
    irow_point = system.irow_point
    irow_elem = system.irow_elem
    irow_elem1 = system.irow_elem1
    irow_elem2 = system.irow_elem2
    icol_point = system.icol_point
    icol_elem = system.icol_elem

    converged = true
    for t in tvec

        # update stored time
        system.t = t

        # current parameters
        pcond = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(t)
        dload = typeof(distributed_loads) <: AbstractDict ? distributed_loads : distributed_loads(t)
        pmass = typeof(point_masses) <: AbstractDict ? point_masses : point_masses(t)
        gvec = typeof(gravity) <: AbstractVector ? SVector{3}(gravity) : SVector{3}(gravity(tvec[it]))
        x0 = typeof(origin) <: AbstractVector ? SVector{3}(origin) : SVector{3}(origin(t))
        v0 = typeof(linear_velocity) <: AbstractVector ? SVector{3}(linear_velocity) : SVector{3}(linear_velocity(t))
        ω0 = typeof(angular_velocity) <: AbstractVector ? SVector{3}(angular_velocity) : SVector{3}(angular_velocity(t))

        f! = (F, x) -> system_residual!(F, x, assembly, pcond, dload, pmass, gvec, force_scaling, 
            mass_scaling, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, 
            icol_elem, x0, v0, ω0)

        j! = (J, x) -> system_jacobian!(J, x, assembly, pcond, dload, pmass, gvec, force_scaling, 
            mass_scaling, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, 
            icol_elem, x0, v0, ω0)

        # solve the system of equations
        if linear
            # linear analysis
            x .= 0.0
            f!(F, x)
            j!(J, x)
            x .-= safe_lu(J) \ F
        else
            # nonlinear analysis
            df = NLsolve.OnceDifferentiable(f!, j!, x, F, J)

            result = NLsolve.nlsolve(df, x,
                linsolve=(x, A, b) -> ldiv!(x, safe_lu(A), b),
                method=method,
                linesearch=linesearch,
                ftol=ftol,
                iterations=iterations)

            # update the solution
            x .= result.zero
            J .= df.DF

            # update the convergence flag
            convergence = result.f_converged
        end
    end

    return system, converged
end

"""
    eigenvalue_analysis(assembly; kwargs...)

Compute the eigenvalues and eigenvectors of the system of nonlinear beams
contained in `assembly`.  Return the modified system, eigenvalues, eigenvectors,
and a convergence flag indicating whether the corresponding steady-state analysis
converged.

# Keyword Arguments
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
 - `point_masses = Dict{Int,Vector{PointMass{Float64}}}()`: A dictionary with keys 
        corresponding to the points at which point masses are attached and values 
        containing vectors of objects of type [`PointMass`](@ref) which describe 
        the point masses attached at those points.  If time varying, this input may
        be provided as a function of time.
 - `gravity = [0,0,0]`: Gravity vector.  If time varying, this input may be provided as a 
        function of time.            
 - `linear = false`: Set to `true` for a linear analysis
 - `method = :newton`: Method (as defined in NLsolve) to solve nonlinear system of equations
 - `linesearch = LineSearches.LineSearches.BackTracking(maxstep=1e6)`: Line search used to solve nonlinear system of equations
 - `ftol = 1e-9`: tolerance for solving nonlinear system of equations
 - `iterations = 1000`: maximum iterations for solving the nonlinear system of equations
 - `reset_state = true`: Flag indicating whether the state variables should be
    reset prior to performing the steady-state analysis.  This keyword argument
    is only valid for the pre-allocated version of this function.
 - `find_steady_state = reset_state && !linear`: Flag indicating whether the
    steady state solution should be found prior to performing the eigenvalue analysis.
 - `origin = zeros(3)`: Global frame origin.
    If time varying, this vector may be provided as a function of time.
 - `linear_velocity = zeros(3)`: Global frame linear velocity vector.
    If time varying, this vector may be provided as a function of time.
    May be provided either as a constant or as a function of time.
 - `angular_velocity = zeros(3)`: Global frame angular velocity vector.
    If time varying, this vector may be provided as a function of time.
 - `tvec`: Time vector. May be used in conjunction with time varying
    prescribed conditions, distributed loads, and global motion to gradually
    increase displacements/loads during the steady-state analysis.
 - `nev = 6`: Number of eigenvalues to compute
"""
function eigenvalue_analysis(assembly;
    prescribed_conditions=Dict{Int,PrescribedConditions{Float64}}(),
    distributed_loads=Dict{Int,DistributedLoads{Float64}}(),
    point_masses=Dict{Int,Vector{PointMass{Float64}}}(),
    gravity=SVector(0,0,0),
    method=:newton,
    linear=false,
    linesearch=LineSearches.BackTracking(maxstep=1e6),
    ftol=1e-9,
    iterations=1000,
    find_steady_state=!linear,
    origin=(@SVector zeros(3)),
    linear_velocity=(@SVector zeros(3)),
    angular_velocity=(@SVector zeros(3)),
    tvec=0.0,
    nev=6
    )

    static = false

    pc = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(tvec[1])

    system = System(assembly, static; prescribed_points=keys(pc))

    return eigenvalue_analysis!(system, assembly;
        prescribed_conditions=prescribed_conditions,
        distributed_loads=distributed_loads,
        point_masses=point_masses,
        gravity=gravity,
        linear=linear,
        method=method,
        linesearch=linesearch,
        ftol=ftol,
        iterations=iterations,
        reset_state=false,
        find_steady_state=find_steady_state,
        origin=origin,
        linear_velocity=linear_velocity,
        angular_velocity=angular_velocity,
        tvec=tvec,
        nev=nev,
        )
end

"""
    eigenvalue_analysis!(system, assembly; kwargs...)

Pre-allocated version of `eigenvalue_analysis`.  Uses the state variables stored in
`system` as an initial guess for iterating to find the steady state solution.
"""
function eigenvalue_analysis!(system, assembly;
    prescribed_conditions=Dict{Int,PrescribedConditions{Float64}}(),
    distributed_loads=Dict{Int,DistributedLoads{Float64}}(),
    point_masses=Dict{Int,Vector{PointMass{Float64}}}(),
    gravity=SVector(0,0,0),
    linear=false,
    method=:newton,
    linesearch=LineSearches.BackTracking(maxstep=1e6),
    ftol=1e-9,
    iterations=1000,
    reset_state=true,
    find_steady_state=!linear && reset_state,
    origin=(@SVector zeros(3)),
    linear_velocity=(@SVector zeros(3)),
    angular_velocity=(@SVector zeros(3)),
    tvec=0.0,
    nev=6,
    )

    if reset_state
        reset_state!(system)
    end

    # perform steady state analysis (if nonlinear)
    if find_steady_state
        system, converged = steady_state_analysis!(system, assembly;
            prescribed_conditions=prescribed_conditions,
            distributed_loads=distributed_loads,
            point_masses=point_masses,
            gravity=gravity,
            linear=linear,
            method=method,
            linesearch=linesearch,
            ftol=ftol,
            iterations=iterations,
            origin=origin,
            linear_velocity=linear_velocity,
            angular_velocity=angular_velocity,
            tvec=tvec,
            )
    else
        # converged by default
        converged = true
    end

    # unpack state vector, stiffness, and mass matrices
    x = system.x # populated during steady state solution
    K = system.K # needs to be updated
    M = system.M # still needs to be populated

    # unpack scaling parameter
    force_scaling = system.force_scaling
    mass_scaling = system.mass_scaling

    # also unpack system indices
    irow_point = system.irow_point
    irow_elem = system.irow_elem
    irow_elem1 = system.irow_elem1
    irow_elem2 = system.irow_elem2
    icol_point = system.icol_point
    icol_elem = system.icol_elem

    # current time
    t = system.t

    # current parameters
    pcond = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(t)
    dload = typeof(distributed_loads) <: AbstractDict ? distributed_loads : distributed_loads(t)
    pmass = typeof(point_masses) <: AbstractDict ? point_masses : point_masses(t)
    gvec = typeof(gravity) <: AbstractVector ? SVector{3}(gravity) : SVector{3}(gravity(tvec[it]))
    x0 = typeof(origin) <: AbstractVector ? SVector{3}(origin) : SVector{3}(origin(t))
    v0 = typeof(linear_velocity) <: AbstractVector ? SVector{3}(linear_velocity) : SVector{3}(linear_velocity(t))
    ω0 = typeof(angular_velocity) <: AbstractVector ? SVector{3}(angular_velocity) : SVector{3}(angular_velocity(t))

    # solve for the system stiffness matrix
    K = system_jacobian!(K, x, assembly, pcond, dload, pmass, gvec, force_scaling, mass_scaling,
        irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem, x0, v0, ω0)

    # solve for the system mass matrix
    M = system_mass_matrix!(M, x, assembly, force_scaling, mass_scaling, irow_point, irow_elem,
        irow_elem1, irow_elem2, icol_point, icol_elem)

    # construct linear map
    T = eltype(system)
    nx = length(x)
    Kfact = safe_lu(K)
    f! = (b, x) -> ldiv!(b, Kfact, M * x)
    fc! = (b, x) -> mul!(b, M', Kfact' \ x)
    A = LinearMap{T}(f!, fc!, nx, nx; ismutating=true)

    # compute eigenvalues and eigenvectors
    λ, V = partialeigen(partialschur(A; nev=min(nx, nev), which=LM())[1])

    # sort eigenvalues by magnitude
    perm = sortperm(λ, by=(λ) -> (abs(λ), imag(λ)), rev=true)
    λ .= λ[perm]
    V .= V[:,perm]

    # eigenvalues are actually -1/λ, no modification necessary for eigenvectors
    λ .= -1 ./ λ

    return system, λ, V, converged
end

"""
    initial_condition_analysis(assembly, t0; kwargs...)

Perform an analysis to obtain a consistent set of initial conditions.  Return the
final system with the new initial conditions.

# Keyword Arguments
 - `prescribed_conditions: A dictionary with keys corresponding to the points at
        which prescribed conditions are applied and values of type
        [`PrescribedConditions`](@ref) which describe the prescribed conditions
        at those points.  If time varying, this input may be provided as a
        function of time.
 - `distributed_loads: A dictionary with keys corresponding to the elements to
        which distributed loads are applied and values of type
        [`DistributedLoads`](@ref) which describe the distributed loads at those
        points.  If time varying, this input may be provided as a function of
        time.
 - `point_masses = Dict{Int,Vector{PointMass{Float64}}}()`: A dictionary with keys 
        corresponding to the points at which point masses are attached and values 
        containing vectors of objects of type [`PointMass`](@ref) which describe 
        the point masses attached at those points.  If time varying, this input may
        be provided as a function of time.
 - `gravity = [0,0,0]`: Gravity vector.  If time varying, this input may be provided as a 
        function of time.
 - `linear = false`: Set to `true` for a linear analysis
 - `method = :newton`: Method (as defined in NLsolve) to solve nonlinear system of equations
 - `linesearch = LineSearches.LineSearches.BackTracking(maxstep=1e6)`: Line search used to solve nonlinear system of equations
 - `ftol = 1e-9`: tolerance for solving nonlinear system of equations
 - `iterations = 1000`: maximum iterations for solving the nonlinear system of equations
 - `reset_state = true`: Flag indicating whether the state variables should be
    reset prior to performing the analysis.  This keyword argument is only valid
    for the pre-allocated version of this function.
 - `origin = zeros(3)`: Global frame origin.
    If time varying, this vector may be provided as a function of time.
 - `linear_velocity = zeros(3)`: Global frame linear velocity vector.
    If time varying, this vector may be provided as a function of time.
    May be provided either as a constant or as a function of time.
 - `angular_velocity = zeros(3)`: Global frame angular velocity vector.
    If time varying, this vector may be provided as a function of time.
 - `u0=fill(zeros(3), length(assembly.elements))`: Initial displacment of each beam element,
 - `theta0=fill(zeros(3), length(assembly.elements))`: Initial angular displacement of each beam element,
 - `udot0=fill(zeros(3), length(assembly.elements))`: Initial time derivative with respect to `u`
 - `thetadot0=fill(zeros(3), length(assembly.elements))`: Initial time derivative with respect to `theta`
 - `save=1:length(tvec)`: Steps at which to save the time history
"""
function initial_condition_analysis(assembly, t0;
    prescribed_conditions=Dict{Int,PrescribedConditions{Float64}}(),
    distributed_loads=Dict{Int,DistributedLoads{Float64}}(),
    point_masses=Dict{Int,Vector{PointMass{Float64}}}(),
    gravity=SVector(0,0,0),
    linear=false,
    method=:newton,
    linesearch=LineSearches.BackTracking(maxstep=1e6),
    ftol=1e-9,
    iterations=1000,
    origin=(@SVector zeros(3)),
    linear_velocity=(@SVector zeros(3)),
    angular_velocity=(@SVector zeros(3)),
    u0=fill((@SVector zeros(3)), length(assembly.elements)),
    theta0=fill((@SVector zeros(3)), length(assembly.elements)),
    udot0=fill((@SVector zeros(3)), length(assembly.elements)),
    thetadot0=fill((@SVector zeros(3)), length(assembly.elements)),
    )

    static = false

    pc = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(t0)

    system = System(assembly, static; prescribed_points=keys(pc))

    return initial_condition_analysis!(system, assembly, t0;
        prescribed_conditions=prescribed_conditions,
        distributed_loads=distributed_loads,
        point_masses=point_masses,
        gravity=gravity,
        linear=linear,
        method=method,
        linesearch=linesearch,
        ftol=ftol,
        iterations=iterations,
        reset_state=false,
        origin=origin,
        linear_velocity=linear_velocity,
        angular_velocity=angular_velocity,
        u0=u0,
        theta0=theta0,
        udot0=udot0,
        thetadot0=thetadot0,
        )
end

"""
    initial_condition_analysis!(system, assembly, t0; kwargs...)

Pre-allocated version of `initial_condition_analysis`.
"""
function initial_condition_analysis!(system, assembly, t0;
    prescribed_conditions=Dict{Int,PrescribedConditions{Float64}}(),
    distributed_loads=Dict{Int,DistributedLoads{Float64}}(),
    point_masses=Dict{Int,Vector{PointMass{Float64}}}(),
    gravity=SVector(0,0,0),
    linear=false,
    method=:newton,
    linesearch=LineSearches.BackTracking(maxstep=1e6),
    ftol=1e-9,
    iterations=1000,
    reset_state=true,
    origin=(@SVector zeros(3)),
    linear_velocity=(@SVector zeros(3)),
    angular_velocity=(@SVector zeros(3)),
    u0=fill((@SVector zeros(3)), length(assembly.elements)),
    theta0=fill((@SVector zeros(3)), length(assembly.elements)),
    udot0=fill((@SVector zeros(3)), length(assembly.elements)),
    thetadot0=fill((@SVector zeros(3)), length(assembly.elements)),
    )

    # check to make sure the simulation is dynamic
    @assert system.static == false

    if reset_state
        reset_state!(system)
    end

    # unpack pre-allocated storage and pointers for system
    x = system.x
    F = system.r
    J = system.K
    force_scaling = system.force_scaling
    mass_scaling = system.mass_scaling
    irow_point = system.irow_point
    irow_elem = system.irow_elem
    irow_elem1 = system.irow_elem1
    irow_elem2 = system.irow_elem2
    icol_point = system.icol_point
    icol_elem = system.icol_elem
    udot = system.udot
    θdot = system.θdot
    Pdot = system.Pdot
    Hdot = system.Hdot

    nelem = length(assembly.elements)

    # set current time step
    system.t = t0

    # set current parameters
    pcond = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(t0)
    dload = typeof(distributed_loads) <: AbstractDict ? distributed_loads : distributed_loads(t0)
    pmass = typeof(point_masses) <: AbstractDict ? point_masses : point_masses(t)
    gvec = typeof(gravity) <: AbstractVector ? SVector{3}(gravity) : SVector{3}(gravity(tvec[it]))
    x0 = typeof(origin) <: AbstractVector ? SVector{3}(origin) : SVector{3}(origin(t0))
    v0 = typeof(linear_velocity) <: AbstractVector ? SVector{3}(linear_velocity) : SVector{3}(linear_velocity(t0))
    ω0 = typeof(angular_velocity) <: AbstractVector ? SVector{3}(angular_velocity) : SVector{3}(angular_velocity(t0))

    # construct residual and jacobian functions
    f! = (F, x) -> system_residual!(F, x, assembly, pcond, dload, pmass, gvec, force_scaling,
        mass_scaling, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem,
        x0, v0, ω0, u0, theta0, udot0, thetadot0)

    j! = (J, x) -> system_jacobian!(J, x, assembly, pcond, dload, pmass, gvec, force_scaling,
        mass_scaling, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem,
        x0, v0, ω0, u0, theta0, udot0, thetadot0)

    # solve system of equations
    if linear
        # linear analysis
        x .= 0.0
        f!(F, x)
        j!(J, x)
        x .-= safe_lu(J) \ F
    else
        # nonlinear analysis
        df = OnceDifferentiable(f!, j!, x, F, J)

        result = NLsolve.nlsolve(df, x,
            linsolve=(x, A, b) -> ldiv!(x, safe_lu(A), b),
            method=method,
            linesearch=linesearch,
            ftol=ftol,
            iterations=iterations)

        x .= result.zero
        J .= df.DF
    end

    # get convergence flag
    converged = result.f_converged

    # save states and state rates
    for ielem = 1:nelem
        icol = icol_elem[ielem]
        # extract rotation parameters for this beam element
        C = get_C(theta0[ielem])
        Cab = assembly.elements[ielem].Cab
        CtCab = C' * Cab
        # save states and state rates
        udot[ielem] = udot0[ielem]
        θdot[ielem] = thetadot0[ielem]
        Pdot[ielem] = CtCab' * SVector(x[icol], x[icol + 1], x[icol + 2]) .* mass_scaling
        Hdot[ielem] = CtCab' * SVector(x[icol + 3], x[icol + 4], x[icol + 5]) .* mass_scaling
# restore original state vector
        x[icol:icol + 2] .= u0[ielem]
        x[icol + 3:icol + 5] .= theta0[ielem]
    end

    return system, converged
end

"""
    time_domain_analysis(assembly, tvec; kwargs...)

Perform a time-domain analysis for the system of nonlinear beams contained in
`assembly` using the time vector `tvec`.  Return the final system, a post-processed
solution history, and a convergence flag indicating whether the iterations
converged for each time step.

# Keyword Arguments
 - `prescribed_conditions: A dictionary with keys corresponding to the points at
        which prescribed conditions are applied and values of type
        [`PrescribedConditions`](@ref) which describe the prescribed conditions
        at those points.  If time varying, this input may be provided as a
        function of time.
 - `distributed_loads: A dictionary with keys corresponding to the elements to
        which distributed loads are applied and values of type
        [`DistributedLoads`](@ref) which describe the distributed loads at those
        points.  If time varying, this input may be provided as a function of
        time.
 - `point_masses = Dict{Int,Vector{PointMass{Float64}}}()`: A dictionary with keys 
        corresponding to the points at which point masses are attached and values 
        containing vectors of objects of type [`PointMass`](@ref) which describe 
        the point masses attached at those points.  If time varying, this input may
        be provided as a function of time.
 - `gravity = [0,0,0]`: Gravity vector.  If time varying, this input may be provided as a 
        function of time.
 - `linear = false`: Set to `true` for a linear analysis
 - `method = :newton`: Method (as defined in NLsolve) to solve nonlinear system of equations
 - `linesearch = LineSearches.LineSearches.BackTracking(maxstep=1e6)`: Line search used to solve nonlinear system of equations
 - `ftol = 1e-9`: tolerance for solving nonlinear system of equations
 - `iterations = 1000`: maximum iterations for solving the nonlinear system of equations
 - `reset_state = true`: Flag indicating whether the state variables should be
    reset prior to performing the analysis.  This keyword argument is only valid
    for the pre-allocated version of this function.
 - `initialize = true`: Flag indicating whether a consistent set of initial
    conditions should be found using [`initial_condition_analysis`](@ref). If
    `false`, the keyword arguments `u0`, `theta0`, `udot0` and `thetadot0` will
    be ignored and the system state vector will be used as the initial state
    variables.
 - `origin`: Global frame origin vector. If time varying, this input
        may be provided as a function of time.
 - `linear_velocity`: Global frame linear velocity vector. If time
        varying, this vector may be provided as a function of time.
 - `angular_velocity`: Global frame angular velocity vector. If time
        varying, this vector may be provided as a function of time.
 - `u0=fill(zeros(3), length(assembly.elements))`: Initial displacment of each beam element,
 - `theta0=fill(zeros(3), length(assembly.elements))`: Initial angular displacement of each beam element,
 - `udot0=fill(zeros(3), length(assembly.elements))`: Initial time derivative with respect to `u`
 - `thetadot0=fill(zeros(3), length(assembly.elements))`: Initial time derivative with respect to `theta`
 - `save=1:length(tvec)`: Steps at which to save the time history
"""
function time_domain_analysis(assembly, tvec;
    prescribed_conditions=Dict{Int,PrescribedConditions{Float64}}(),
    distributed_loads=Dict{Int,DistributedLoads{Float64}}(),
    point_masses=Dict{Int,Vector{PointMass{Float64}}}(),
    gravity=SVector(0,0,0),
    linear=false,
    method=:newton,
    linesearch=LineSearches.BackTracking(maxstep=1e6),
    ftol=1e-9,
    iterations=1000,
    initialize=true,
    origin=(@SVector zeros(3)),
    linear_velocity=(@SVector zeros(3)),
    angular_velocity=(@SVector zeros(3)),
    u0=fill((@SVector zeros(3)), length(assembly.elements)),
    theta0=fill((@SVector zeros(3)), length(assembly.elements)),
    udot0=fill((@SVector zeros(3)), length(assembly.elements)),
    thetadot0=fill((@SVector zeros(3)), length(assembly.elements)),
    save=1:length(tvec)
    )

    static = false

    pc = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(tvec[1])

    system = System(assembly, static; prescribed_points=keys(pc))

    return time_domain_analysis!(system, assembly, tvec;
        prescribed_conditions=prescribed_conditions,
        distributed_loads=distributed_loads,
        point_masses=point_masses,
        gravity=gravity,
        linear=linear,
        method=method,
        linesearch=linesearch,
        ftol=ftol,
        iterations=iterations,
        reset_state=false,
        initialize=initialize,
        origin=origin,
        linear_velocity=linear_velocity,
        angular_velocity=angular_velocity,
        u0=u0,
        theta0=theta0,
        udot0=udot0,
        thetadot0=thetadot0,
        save=save,
        )
end

    """
    time_domain_analysis!(system, assembly, tvec; kwargs...)

Pre-allocated version of `time_domain_analysis`.
"""
function time_domain_analysis!(system, assembly, tvec;
    prescribed_conditions=Dict{Int,PrescribedConditions{Float64}}(),
    distributed_loads=Dict{Int,DistributedLoads{Float64}}(),
    point_masses=Dict{Int,Vector{PointMass{Float64}}}(),
    gravity=SVector(0,0,0),
    linear=false,
    method=:newton,
    linesearch=LineSearches.BackTracking(maxstep=1e6),
    ftol=1e-9,
    iterations=1000,
    reset_state=true,
    initialize=true,
    origin=(@SVector zeros(3)),
    linear_velocity=(@SVector zeros(3)),
    angular_velocity=(@SVector zeros(3)),
    u0=fill((@SVector zeros(3)), length(assembly.elements)),
    theta0=fill((@SVector zeros(3)), length(assembly.elements)),
    udot0=fill((@SVector zeros(3)), length(assembly.elements)),
    thetadot0=fill((@SVector zeros(3)), length(assembly.elements)),
    save=1:length(tvec)
    )

    # check to make sure the simulation is dynamic
    @assert system.static == false

    if reset_state
        reset_state!(system)
    end

    # perform initial condition analysis
    if initialize
        system, converged = initial_condition_analysis!(system, assembly, tvec[1];
            prescribed_conditions=prescribed_conditions,
            distributed_loads=distributed_loads,
            point_masses=point_masses,
            gravity=gravity,
            linear=linear,
            method=method,
            linesearch=linesearch,
            ftol=ftol,
            iterations=iterations,
            reset_state=false,
            origin=origin,
            linear_velocity=linear_velocity,
            angular_velocity=angular_velocity,
            u0=u0,
            theta0=theta0,
            udot0=udot0,
            thetadot0=thetadot0,
            )
    else
        # converged by default
        converged = true
    end

    # unpack pre-allocated storage and pointers for system
    x = system.x
    F = system.r
    J = system.K
    force_scaling = system.force_scaling
    mass_scaling = system.mass_scaling
    irow_point = system.irow_point
    irow_elem = system.irow_elem
    irow_elem1 = system.irow_elem1
    irow_elem2 = system.irow_elem2
    icol_point = system.icol_point
    icol_elem = system.icol_elem
    udot = system.udot
    θdot = system.θdot
    Pdot = system.Pdot
    Hdot = system.Hdot

    # number of beam elements
    nelem = length(assembly.elements)

    # initialize storage for each time step
    isave = 1
    history = Vector{AssemblyState{eltype(system)}}(undef, length(save))

    # add initial state to the solution history
    if isave in save
        pcond = typeof(prescribed_conditions) <: AbstractDict ?
            prescribed_conditions : prescribed_conditions(tvec[1])
        history[isave] = AssemblyState(system, assembly, prescribed_conditions=pcond)
        isave += 1
    end

    # --- Begin Time Domain Simulation --- #

    for it = 2:length(tvec)

        # update current time
        system.t = tvec[it]

        # current time step size
        dt = tvec[it] - tvec[it - 1]

        # current parameters
        pcond = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(tvec[it])
        dload = typeof(distributed_loads) <: AbstractDict ? distributed_loads : distributed_loads(tvec[it])
        pmass = typeof(point_masses) <: AbstractDict ? point_masses : point_masses(t)
        gvec = typeof(gravity) <: AbstractVector ? SVector{3}(gravity) : SVector{3}(gravity(tvec[it]))
        x0 = typeof(origin) <: AbstractVector ? SVector{3}(origin) : SVector{3}(origin(tvec[it]))
        v0 = typeof(linear_velocity) <: AbstractVector ? SVector{3}(linear_velocity) : SVector{3}(linear_velocity(tvec[it]))
        ω0 = typeof(angular_velocity) <: AbstractVector ? SVector{3}(angular_velocity) : SVector{3}(angular_velocity(tvec[it]))

        # current initialization parameters
        for ielem = 1:nelem
            icol = icol_elem[ielem]
            # get beam element states
            u = SVector(x[icol], x[icol + 1], x[icol + 2])
            θ = SVector(x[icol + 3], x[icol + 4], x[icol + 5])
            P = SVector(x[icol + 12], x[icol + 13], x[icol + 14]) .* mass_scaling
            H = SVector(x[icol + 15], x[icol + 16], x[icol + 17]) .* mass_scaling
            # extract rotation parameters
            C = get_C(θ)
            Cab = assembly.elements[ielem].Cab
            CtCab = C' * Cab
            # store `udot_init` in `udot`
            udot[ielem] = 2 / dt * u + udot[ielem]
            # store `θdot_init` in `θdot`
            θdot[ielem] = 2 / dt * θ + θdot[ielem]
            # store `CtCabPdot_init` in `Pdot`
            Pdot[ielem] = 2 / dt * CtCab * P + CtCab * Pdot[ielem]
            # store `CtCabHdot_init` in `Hdot`
            Hdot[ielem] = 2 / dt * CtCab * H + CtCab * Hdot[ielem]
        end

        # solve for the state variables at the next time step
        f! = (F, x) -> system_residual!(F, x, assembly, pcond, dload, pmass, gvec, force_scaling,
            mass_scaling, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem,
            x0, v0, ω0, udot, θdot, Pdot, Hdot, dt)

        j! = (J, x) -> system_jacobian!(J, x, assembly, pcond, dload, pmass, gvec, force_scaling,
            mass_scaling, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem,
            x0, v0, ω0, udot, θdot, Pdot, Hdot, dt)

        # solve system of equations
        if linear
            # linear analysis
            x .= 0.0
            f!(F, x)
            j!(J, x)
            x .-= safe_lu(J) \ F
        else
            df = OnceDifferentiable(f!, j!, x, F, J)

            result = NLsolve.nlsolve(df, x,
                linsolve=(x, A, b) -> ldiv!(x, safe_lu(A), b),
                method=method,
                linesearch=linesearch,
                ftol=ftol,
                iterations=iterations)

            x .= result.zero
            J .= df.DF
        end

        # set new state rates
        for ielem = 1:nelem
            icol = icol_elem[ielem]
            # get beam element states
            u = SVector(x[icol], x[icol + 1], x[icol + 2])
            θ = SVector(x[icol + 3], x[icol + 4], x[icol + 5])
            P = SVector(x[icol + 12], x[icol + 13], x[icol + 14]) .* mass_scaling
            H = SVector(x[icol + 15], x[icol + 16], x[icol + 17]) .* mass_scaling
            # extract rotation parameters
            C = get_C(θ)
            Cab = assembly.elements[ielem].Cab
            CtCab = C' * Cab
            # save state rates
            udot[ielem] = 2 / dt * u - udot[ielem]
            θdot[ielem] = 2 / dt * θ - θdot[ielem]
            Pdot[ielem] = 2 / dt * P - CtCab' * Pdot[ielem]
            Hdot[ielem] = 2 / dt * H - CtCab' * Hdot[ielem]
        end

        # add state to history
        if it in save
            history[isave] = AssemblyState(system, assembly, prescribed_conditions=prescribed_conditions)
            isave += 1
        end

        # stop early if unconverged
        if !linear && !result.f_converged
            converged = false
            break
        end

    end

    return system, history, converged
end
