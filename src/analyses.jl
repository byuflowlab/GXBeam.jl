"""
    static_analysis(assembly; kwargs...)

Perform a static analysis of the system of nonlinear beams contained in
`assembly`. Return the resulting system and a flag indicating whether the
iteration procedure converged.

# Keyword Arguments
 - `prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}()`: Dictionary
     holding PrescribedConditions composite types for the points in `keys(prescribed_conditions)`
 - `distributed_loads = Dict{Int,DistributedLoads{Float64}}()`: Dictionary holding
     DistributedLoads composite types for the beam elements in `keys(distributed_loads)`
 - `linear = false`: Set to `true` for a linear analysis
 - `method = :newton`: Method (as defined in NLsolve) to solve nonlinear system of equations
 - `linesearch = LineSearches.BackTracking()`: Line search used to solve nonlinear system of equations
 - `ftol = 1e-9`: tolerance for solving nonlinear system of equations
 - `iterations = 1000`: maximum iterations for solving the nonlinear system of equations
 - `nstep = 1`: Number of time steps. May be used in conjunction with time varying
     prescribed conditions and distributed loads to gradually increase
    displacements/loads.
"""
function static_analysis(assembly;
    prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}(),
    distributed_loads = Dict{Int,DistributedLoads{Float64}}(),
    linear = false,
    method = :newton,
    linesearch = BackTracking(),
    ftol = 1e-9,
    iterations = 1000,
    nstep = 1,
    )

    static = true

    system = System(assembly, keys(prescribed_conditions), static)

    return static_analysis!(system, assembly;
        prescribed_conditions = prescribed_conditions,
        distributed_loads = distributed_loads,
        linear = linear,
        method = method,
        linesearch = linesearch,
        ftol = ftol,
        iterations = iterations,
        nstep = nstep)
end

"""
    static_analysis!(system, assembly; kwargs...)

Pre-allocated version of `static_analysis`.  Uses the state variables stored in
`system` as an initial guess for iterating.
"""
function static_analysis!(system, assembly;
    prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}(),
    distributed_loads = Dict{Int,DistributedLoads{Float64}}(),
    linear = false,
    method = :newton,
    linesearch = BackTracking(),
    ftol = 1e-9,
    iterations = 1000,
    nstep = 1)

    # check to make sure system is static
    @assert system.static == true

    # unpack pre-allocated storage and pointers
    x = system.x
    F = system.r
    J = system.K
    irow_pt = system.irow_pt
    irow_beam = system.irow_beam
    irow_beam1 = system.irow_beam1
    irow_beam2 = system.irow_beam2
    icol_pt = system.icol_pt
    icol_beam = system.icol_beam
    current_step_ref = system.current_step

    converged = true
    for istep = 1:nstep

        # update current time step
        current_step_ref[] = istep

        # solve the system of equations
        f! = (F, x) -> system_residual!(F, x, assembly, prescribed_conditions,
            distributed_loads, istep, irow_pt, irow_beam, irow_beam1,
            irow_beam2, icol_pt, icol_beam)

        j! = (J, x) -> system_jacobian!(J, x, assembly, prescribed_conditions,
            distributed_loads, istep, irow_pt, irow_beam, irow_beam1, irow_beam2,
            icol_pt, icol_beam)

        if linear
            # linear analysis
            x .= 0.0
            f!(F, x)
            j!(J, x)
            x .-= J\F
        else
            # nonlinear analysis
            df = NLsolve.OnceDifferentiable(f!, j!, x, F, J)

            result = NLsolve.nlsolve(df, x,
                linsolve=(x, A, b) -> ldiv!(x, lu(A), b),
                method=method,
                linesearch=linesearch,
                ftol=ftol,
                iterations=iterations)

            # update the solution
            x .= result.zero

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
 - `prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}()`: Dictionary
     holding PrescribedConditions composite types for the points in `keys(prescribed_conditions)`
 - `distributed_loads = Dict{Int,DistributedLoads{Float64}}()`: Dictionary holding
     DistributedLoads composite types for the beam elements in `keys(distributed_loads)`
 - `linear = false`: Set to `true` for a linear analysis
 - `method = :newton`: Method (as defined in NLsolve) to solve nonlinear system of equations
 - `linesearch = LineSearches.BackTracking()`: Line search used to solve nonlinear system of equations
 - `ftol = 1e-9`: tolerance for solving nonlinear system of equations
 - `iterations = 1000`: maximum iterations for solving the nonlinear system of equations
 - `nstep = 1`: Number of time steps. May be used in conjunction with time varying
     prescribed conditions, distributed loads, and global motion to gradually
     increase displacements/loads.
 - `origin = zeros(3)`: Global frame origin
 - `linear_velocity = fill(zeros(3), nstep)`: Global frame linear velocity for each time step.
 - `angular_velocity = fill(zeros(3), nstep)`: Global frame angular velocity for each time step
"""
function steady_state_analysis(assembly;
    prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}(),
    distributed_loads = Dict{Int,DistributedLoads{Float64}}(),
    linear = false,
    method = :newton,
    linesearch = BackTracking(),
    ftol = 1e-9,
    iterations = 1000,
    nstep = 1,
    origin = (@SVector zeros(3)),
    linear_velocity = fill((@SVector zeros(3)), nstep),
    angular_velocity = fill((@SVector zeros(3)), nstep),
    )

    static = false

    system = System(assembly, keys(prescribed_conditions), static)

    return steady_state_analysis!(system, assembly;
        prescribed_conditions = prescribed_conditions,
        distributed_loads = distributed_loads,
        linear = linear,
        method = method,
        linesearch = linesearch,
        ftol = ftol,
        iterations = iterations,
        nstep = nstep,
        origin = origin,
        linear_velocity = linear_velocity,
        angular_velocity = angular_velocity
        )
end

"""
    steady_state_analysis!(system, assembly; kwargs...)

Pre-allocated version of `steady_state_analysis`.  Uses the state variables stored in
`system` as an initial guess for iterating.
"""
function steady_state_analysis!(system, assembly;
    prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}(),
    distributed_loads = Dict{Int,DistributedLoads{Float64}}(),
    linear = false,
    method = :newton,
    linesearch = BackTracking(),
    ftol = 1e-9,
    iterations = 1000,
    nstep = 1,
    origin = (@SVector zeros(3)),
    linear_velocity = fill((@SVector zeros(3)), nstep),
    angular_velocity = fill((@SVector zeros(3)), nstep),
    )

    # check to make sure the simulation is dynamic
    @assert system.static == false

    # unpack pointers to pre-allocated storage
    x = system.x
    F = system.r
    J = system.K
    irow_pt = system.irow_pt
    irow_beam = system.irow_beam
    irow_beam1 = system.irow_beam1
    irow_beam2 = system.irow_beam2
    icol_pt = system.icol_pt
    icol_beam = system.icol_beam
    current_step_ref = system.current_step

    # process global frame motion inputs (in case only one value was provided)
    x0 = SVector{3}(origin)
    if eltype(linear_velocity) <: Number
        v0 = fill(SVector{3}(linear_velocity), nstep)
    else
        v0 = SVector{3}.(linear_velocity)
    end
    if eltype(angular_velocity) <: Number
        ω0 = fill(SVector{3}(angular_velocity), nstep)
    else
        ω0 = SVector{3}.(angular_velocity)
    end

    converged = true
    for istep = 1:nstep

        # update current time step
        current_step_ref[] = istep

        f! = (F, x) -> system_residual!(F, x, assembly, prescribed_conditions,
            distributed_loads, istep, irow_pt, irow_beam, irow_beam1, irow_beam2,
            icol_pt, icol_beam, x0, v0[istep], ω0[istep])

        j! = (J, x) -> system_jacobian!(J, x, assembly, prescribed_conditions,
            distributed_loads, istep, irow_pt, irow_beam, irow_beam1, irow_beam2,
            icol_pt, icol_beam, x0, v0[istep], ω0[istep])

        # solve the system of equations
        if linear
            # linear analysis
            x .= 0.0
            f!(F, x)
            j!(J, x)
            x .-= J\F
        else
            # nonlinear analysis
            df = NLsolve.OnceDifferentiable(f!, j!, x, F, J)

            result = NLsolve.nlsolve(df, x,
                linsolve=(x, A, b) -> ldiv!(x, lu(A), b),
                method=method,
                linesearch=linesearch,
                ftol=ftol,
                iterations=iterations)

            # update the solution
            x .= result.zero

            # update the convergence flag
            convergence = result.f_converged
        end
    end

    return system, converged
end

"""
    eigenvalue_analysis(assembly; kwargs...)

Compute the eigenvalues and eigenvectors of the system of nonlinear beams
contained in `assembly` by calling ARPACK.  Return the modified system,
eigenvalues, eigenvectors, and a convergence flag indicating whether the corresponding steady-state analysis
converged.

# Keyword Arguments
 - `prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}()`: Dictionary
     holding PrescribedConditions composite types for the points in `keys(prescribed_conditions)`
 - `distributed_loads = Dict{Int,DistributedLoads{Float64}}()`: Dictionary holding
     DistributedLoads composite types for the beam elements in `keys(distributed_loads)`
 - `linear = false`: Set to `true` for a linear analysis
 - `method = :newton`: Method (as defined in NLsolve) to solve nonlinear system of equations
 - `linesearch = LineSearches.BackTracking()`: Line search used to solve nonlinear system of equations
 - `ftol = 1e-9`: tolerance for solving nonlinear system of equations
 - `iterations = 1000`: maximum iterations for solving the nonlinear system of equations
 - `nstep = 1`: Number of time steps. May be used in conjunction with time varying
     prescribed conditions, distributed loads, and global motion to gradually
     increase displacements/loads.
 - `origin = zeros(3)`: Global frame origin
 - `linear_velocity = fill(zeros(3), nstep)`: Global frame linear velocity for each time step.
 - `angular_velocity = fill(zeros(3), nstep)`: Global frame angular velocity for each time step
 - `nev = 6`: Number of eigenvalues to compute
"""
function eigenvalue_analysis(assembly;
    prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}(),
    distributed_loads = Dict{Int,DistributedLoads{Float64}}(),
    method = :newton,
    linear = false,
    linesearch = BackTracking(),
    ftol = 1e-9,
    iterations = 1000,
    nstep = 1,
    origin = (@SVector zeros(3)),
    linear_velocity = fill((@SVector zeros(3)), nstep),
    angular_velocity = fill((@SVector zeros(3)), nstep),
    nev = 6
    )

    static = false

    system = System(assembly, keys(prescribed_conditions), static)

    return eigenvalue_analysis!(system, assembly;
        prescribed_conditions = prescribed_conditions,
        distributed_loads = distributed_loads,
        linear = linear,
        method = method,
        linesearch = linesearch,
        ftol = ftol,
        iterations = iterations,
        nstep = nstep,
        origin = origin,
        linear_velocity = linear_velocity,
        angular_velocity = angular_velocity,
        nev = nev
        )
end

"""
    eigenvalue_analysis!(system, assembly; kwargs...)

Pre-allocated version of `eigenvalue_analysis`.  Uses the state variables stored in
`system` as an initial guess for iterating to find the steady state solution.
"""
function eigenvalue_analysis!(system, assembly;
    prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}(),
    distributed_loads = Dict{Int,DistributedLoads{Float64}}(),
    linear = false,
    method = :newton,
    linesearch = BackTracking(),
    ftol = 1e-9,
    iterations = 1000,
    nstep = 1,
    origin = (@SVector zeros(3)),
    linear_velocity = fill((@SVector zeros(3)), nstep),
    angular_velocity = fill((@SVector zeros(3)), nstep),
    nev = 6,
    )

    # perform steady state analysis (if nonlinear)
    if linear
        system.x .= 0.0
        converged = true
    else
        system, converged = steady_state_analysis!(system, assembly;
            prescribed_conditions = prescribed_conditions,
            distributed_loads = distributed_loads,
            linear = linear,
            method = method,
            linesearch = linesearch,
            ftol = ftol,
            iterations = iterations,
            nstep = nstep,
            origin = origin,
            linear_velocity = linear_velocity,
            angular_velocity = angular_velocity
            )
    end

    # unpack state vector, stiffness, and mass matrices
    x = system.x # populated during steady state solution
    K = system.K # needs to be updated
    M = system.M # still needs to be populated

    # also unpack system indices
    irow_pt = system.irow_pt
    irow_beam = system.irow_beam
    irow_beam1 = system.irow_beam1
    irow_beam2 = system.irow_beam2
    icol_pt = system.icol_pt
    icol_beam = system.icol_beam

    # get current time step
    istep = system.current_step[]

    # process global frame motion inputs (in case only one value was provided)
    x0 = SVector{3}(origin)
    if eltype(linear_velocity) <: Number
        v0 = SVector{3}(linear_velocity)
    else
        v0 = SVector{3}(linear_velocity[istep])
    end
    if eltype(angular_velocity) <: Number
        ω0 = SVector{3}(angular_velocity)
    else
        ω0 = SVector{3}(angular_velocity[istep])
    end

    # solve for the system stiffness matrix
    K = system_jacobian!(K, x, assembly, prescribed_conditions,
        distributed_loads, istep, irow_pt, irow_beam,
        irow_beam1, irow_beam2, icol_pt, icol_beam, x0, v0, ω0)

    # solve for the system mass matrix
    M = system_mass_matrix!(M, x, assembly, irow_pt, irow_beam, irow_beam1,
        irow_beam2, icol_pt, icol_beam)

    # construct linear map
    T = eltype(system)
    nx = length(x)
    Kfact = lu(K)
    f! = (b, x) -> ldiv!(b, Kfact, M*x)
    fc! = (b, x) -> mul!(b, M', Kfact'\x)
    A = LinearMap{T}(f!, fc!, nx, nx; ismutating=true)

    # compute eigenvalues and eigenvectors
    λ, V, _ = Arpack.eigs(A; nev=min(nx,nev), which=:LM)

    # eigenvalues are actually 1/λ, no modification necessary for eigenvectors
    λ .= 1 ./ λ

    return system, λ, V, converged
end

"""
    time_domain_analysis(assembly, dt; kwargs...)

Perform a time-domain analysis for the system of nonlinear beams contained in
`assembly`.  Return the final system, a post-processed solution history, and a
convergence flag indicating whether the iterations converged for each time step.

# Keyword Arguments
 - `prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}()`: Dictionary
     holding PrescribedConditions composite types for the points in `keys(prescribed_conditions)`
 - `distributed_loads = Dict{Int,DistributedLoads{Float64}}()`: Dictionary holding
     DistributedLoads composite types for the beam elements in `keys(distributed_loads)`
 - `linear = false`: Set to `true` for a linear analysis
 - `method = :newton`: Method (as defined in NLsolve) to solve nonlinear system of equations
 - `linesearch = LineSearches.BackTracking()`: Line search used to solve nonlinear system of equations
 - `ftol = 1e-9`: tolerance for solving nonlinear system of equations
 - `iterations = 1000`: maximum iterations for solving the nonlinear system of equations
 - `nstep = 1`: The total length of the time vector
 - `origin = zeros(3)`: Global frame origin
 - `linear_velocity = fill(zeros(3), nstep)`: Global frame linear velocity for each time step.
 - `angular_velocity = fill(zeros(3), nstep)`: Global frame angular velocity for each time step.
 - `u0=fill(zeros(3), length(assembly.elements))`: initial displacment of each beam element,
 - `theta0=fill(zeros(3), length(assembly.elements))`: initial angular displacement of each beam element,
 - `udot0=fill(zeros(3), length(assembly.elements))`: initial time derivative with respect to `u`
 - `thetadot0=fill(zeros(3), length(assembly.elements))`: initial time derivative with respect to `theta`
 - `save=1:nstep`: Steps at which to save the time history
"""
function time_domain_analysis(assembly, dt;
    prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}(),
    distributed_loads = Dict{Int,DistributedLoads{Float64}}(),
    linear = false,
    method = :newton,
    linesearch = BackTracking(),
    ftol = 1e-9,
    iterations = 1000,
    nstep = 1,
    origin = (@SVector zeros(3)),
    linear_velocity = fill((@SVector zeros(3)), nstep),
    angular_velocity = fill((@SVector zeros(3)), nstep),
    u0 = fill((@SVector zeros(3)), length(assembly.elements)),
    theta0 = fill((@SVector zeros(3)), length(assembly.elements)),
    udot0 = fill((@SVector zeros(3)), length(assembly.elements)),
    thetadot0 = fill((@SVector zeros(3)), length(assembly.elements)),
    save = 1:nstep
    )

    static = false

    system = System(assembly, keys(prescribed_conditions), static)

    return time_domain_analysis!(system, assembly, dt;
        prescribed_conditions = prescribed_conditions,
        distributed_loads = distributed_loads,
        linear = linear,
        method = method,
        linesearch = linesearch,
        ftol = ftol,
        iterations = iterations,
        nstep = nstep,
        origin = origin,
        linear_velocity = linear_velocity,
        angular_velocity = angular_velocity,
        u0 = u0,
        theta0 = theta0,
        udot0 = udot0,
        thetadot0 = thetadot0,
        save = 1:nstep
        )
end

"""
    time_domain_analysis!(system, assembly, dt; kwargs...)

Pre-allocated version of `time_domain_analysis`.  Uses the state variables stored
in `system` as an initial guess for iterating.
"""
function time_domain_analysis!(system, assembly, dt;
    prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}(),
    distributed_loads = Dict{Int,DistributedLoads{Float64}}(),
    linear = false,
    method = :newton,
    linesearch = BackTracking(),
    ftol = 1e-9,
    iterations = 1000,
    nstep = 1,
    origin = (@SVector zeros(3)),
    linear_velocity = fill((@SVector zeros(3)), nstep),
    angular_velocity = fill((@SVector zeros(3)), nstep),
    u0 = fill((@SVector zeros(3)), length(assembly.elements)),
    theta0 = fill((@SVector zeros(3)), length(assembly.elements)),
    udot0 = fill((@SVector zeros(3)), length(assembly.elements)),
    thetadot0 = fill((@SVector zeros(3)), length(assembly.elements)),
    save = 1:nstep
    )

    # check to make sure the simulation is dynamic
    @assert system.static == false

    # unpack pre-allocated storage and pointers for system
    x = system.x
    F = system.r
    J = system.K
    irow_pt = system.irow_pt
    irow_beam = system.irow_beam
    irow_beam1 = system.irow_beam1
    irow_beam2 = system.irow_beam2
    icol_pt = system.icol_pt
    icol_beam = system.icol_beam
    udot_init = system.udot_init
    θdot_init = system.θdot_init
    CtCabPdot_init = system.CtCabPdot_init
    CtCabHdot_init = system.CtCabHdot_init
    current_step_ref = system.current_step

    nbeam = length(assembly.elements)

    # process global frame motion inputs (in case only one value was provided)
    x0 = SVector{3}(origin)
    if eltype(linear_velocity) <: Number
        v0 = fill(SVector{3}(linear_velocity), nstep)
    else
        v0 = SVector{3}.(linear_velocity)
    end
    if eltype(angular_velocity) <: Number
        ω0 = fill(SVector{3}(angular_velocity), nstep)
    else
        ω0 = SVector{3}.(angular_velocity)
    end

    # --- Initial Condition Run --- #

    # update current time step
    current_step_ref[] = 1

    # construct residual and jacobian functions
    f! = (F, x) -> system_residual!(F, x, assembly, prescribed_conditions,
        distributed_loads, 1, irow_pt, irow_beam, irow_beam1, irow_beam2,
        icol_pt, icol_beam, x0, v0[1], ω0[1], u0, theta0, udot0, thetadot0)

    j! = (J, x) -> system_jacobian!(J, x, assembly, prescribed_conditions,
        distributed_loads, 1, irow_pt, irow_beam, irow_beam1, irow_beam2,
        icol_pt, icol_beam, x0, v0[1], ω0[1], u0, theta0, udot0, thetadot0)

    # solve system of equations
    if linear
        # linear analysis
        x .= 0.0
        f!(F, x)
        j!(J, x)
        x .-= J\F
    else
        # nonlinear analysis
        df = OnceDifferentiable(f!, j!, x, F, J)

        result = NLsolve.nlsolve(df, x,
            linsolve=(x, A, b) -> ldiv!(x, lu(A), b),
            method=method,
            linesearch=linesearch,
            ftol=ftol,
            iterations=iterations)

        x .= result.zero
    end

    converged = result.f_converged

    # --- End Initial Condition Run --- #

    # now set up for the time-domain run
    for ibeam = 1:nbeam
        icol = icol_beam[ibeam]
        # calculate udot_init
        udot_init[ibeam] = 2/dt*u0[ibeam] + udot0[ibeam]
        # calculate θdot_init
        θdot_init[ibeam] = 2/dt*theta0[ibeam] + thetadot0[ibeam]
        # extract rotation parameters
        C = get_C(theta0[ibeam])
        Cab = assembly.elements[ibeam].Cab
        CtCab = C'*Cab
        # calculate CtCabPdot_init
        CtCabP = CtCab*SVector{3}(x[icol+12], x[icol+13], x[icol+14])
        CtCabPdot = SVector{3}(x[icol], x[icol+1], x[icol+2])
        CtCabPdot_init[ibeam] = 2/dt*CtCabP + CtCabPdot
        # calculate CtCabHdot_init
        CtCabH = CtCab*SVector{3}(x[icol+15], x[icol+16], x[icol+17])
        CtCabHdot = SVector{3}(x[icol+3], x[icol+4], x[icol+5])
        CtCabHdot_init[ibeam] = 2/dt*CtCabH + CtCabHdot
        # insert initial conditions for time-domain analysis
        x[icol:icol+2] .= u0[ibeam]
        x[icol+3:icol+5] .= theta0[ibeam]
    end

    # initialize storage for each time step
    isave = 1
    history = Vector{AssemblyState{eltype(system)}}(undef, length(save))

    # add initial state to the solution history (if it should be saved)
    if 1 in save
        history[isave] = AssemblyState(system, assembly, prescribed_conditions = prescribed_conditions)
        isave += 1
    end

    # --- Begin Time Domain Simulation --- #

    for istep = 2:nstep

        # update current time step
        current_step_ref[] = istep

        # solve for the state variables at the next time step
        f! = (F, x) -> system_residual!(F, x, assembly, prescribed_conditions,
            distributed_loads, istep, irow_pt, irow_beam, irow_beam1, irow_beam2,
            icol_pt, icol_beam, x0, v0[istep], ω0[istep], udot_init, θdot_init, CtCabPdot_init,
            CtCabHdot_init, dt)

        j! = (J, x) -> system_jacobian!(J, x, assembly, prescribed_conditions,
            distributed_loads, istep, irow_pt, irow_beam, irow_beam1, irow_beam2,
            icol_pt, icol_beam, x0, v0[istep], ω0[istep], udot_init, θdot_init, CtCabPdot_init,
            CtCabHdot_init, dt)

        # solve system of equations
        if linear
            # linear analysis
            x .= 0.0
            f!(F, x)
            j!(J, x)
            x .-= J\F
        else
            df = OnceDifferentiable(f!, j!, x, F, J)

            result = NLsolve.nlsolve(df, x,
                linsolve=(x, A, b) -> ldiv!(x, lu(A), b),
                method=method,
                linesearch=linesearch,
                ftol=ftol,
                iterations=iterations)

            x .= result.zero
        end

        # add state to history
        if istep in save
            history[isave] = AssemblyState(system, assembly, prescribed_conditions = prescribed_conditions)
            isave += 1
        end

        # stop early if unconverged
        if !result.f_converged
            converged = false
            break
        end

        # exit loop if done iterating
        if istep == nstep
            break
        end

        # update time derivative terms for the next time step
        for ibeam = 1:nbeam
            icol = icol_beam[ibeam]
            # calculate udot for next time step
            u = SVector(x[icol], x[icol+1], x[icol+2])
            udot = 2/dt*u - udot_init[ibeam]
            udot_init[ibeam] = 2/dt*u + udot
            # calculate θdot for next time step
            θ = SVector(x[icol+3], x[icol+4], x[icol+5])
            θdot = 2/dt*θ - θdot_init[ibeam]
            θdot_init[ibeam] = 2/dt*θ + θdot
            # extract rotation parameters
            C = get_C(θ)
            Cab = assembly.elements[ibeam].Cab
            CtCab = C'*Cab
            # calculate CtCabPdot for next time step
            CtCabP = CtCab*SVector(x[icol+12], x[icol+13], x[icol+14])
            CtCabPdot = 2/dt*CtCabP - CtCabPdot_init[ibeam]
            CtCabPdot_init[ibeam] = 2/dt*CtCabP + CtCabPdot
            # calculate CtCabHdot for next time step
            CtCabH = CtCab*SVector(x[icol+15], x[icol+16], x[icol+17])
            CtCabHdot = 2/dt*CtCabH - CtCabHdot_init[ibeam]
            CtCabHdot_init[ibeam] = 2/dt*CtCabH + CtCabHdot
        end
    end

    # --- End Time Domain Simulation --- #

    return system, history, converged
end

# this function is currently unused
function linsolve!(x, A, b)

    # get scaling vectors
    r, s = linf_norm_scaling(A)

    # scale matrix and vector
    A .= A .* r .* s'
    b .= b .* r

    # do linear solve
    ldiv!(x, lu(A), b)

    # remove scaling from result
    x .= x .* s

    # remove scaling from matrix and vector
    A .= A ./ r ./ s'
    b .= b ./ r

    return x
end
