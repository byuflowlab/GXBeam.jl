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
 - `t = 0`: Current time or time vector. May be used in conjunction with time varying
        prescribed conditions and distributed loads to gradually increase displacements and 
        loads.

# Control Flag Keyword Arguments
 - `reset_state = true`: Flag indicating whether the system state variables should be 
        set to zero prior to performing this analysis.
 - `linear = false`: Flag indicating whether a linear analysis should be performed.
 - `show_trace = false`: Flag indicating whether to display the solution progress.

# Nonlinear Solver Keyword Arguments (see [NLsolve.jl](https://github.com/JuliaNLSolvers/NLsolve.jl))
 - `method = :newton`: Solution method for nonlinear systems of equations
 - `linesearch = LineSearches.BackTracking(maxstep=1e6)`: Line search for solving nonlinear
        systems of equations
 - `ftol = 1e-9`: Tolerance for solving nonlinear systems of equations
 - `iterations = 1000`: Iteration limit when solving nonlinear systems of equations

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
    t=0.0,
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
    for ti in t

        # print the current time
        if show_trace
            println("Solving for t=$t")
        end

        # update the stored time
        system.t = ti

        # current parameters
        pcond = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(ti)
        dload = typeof(distributed_loads) <: AbstractDict ? distributed_loads : distributed_loads(ti)
        pmass = typeof(point_masses) <: AbstractDict ? point_masses : point_masses(ti)
        gvec = typeof(gravity) <: AbstractVector ? SVector{3}(gravity) : SVector{3}(gravity(ti))

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

    return system, converged
end

"""
    linearize!(system, assembly; kwargs...)

Return the state variables, jacobian matrix, and mass matrix of a linearized system using
the current system state vector.  Note that the returned vectors and matrices are aliased 
with variables in `system` so a copy should be made prior to modifying them.

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
 - `show_trace = false`: Flag indicating whether to display the solution progress.

"""
function linearize!(system, assembly;
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
    constant_mass_matrix=typeof(system) <: ExpandedSystem,
    show_trace=false,
    )

    # check if provided system is consistent with provided keyword arguments
    constant_mass_matrix && @assert typeof(system) <: ExpandedSystem

    # unpack state vector, stiffness, and mass matrices
    @unpack x, K, M, force_scaling, indices = system

    # current parameters
    pcond = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(time)
    dload = typeof(distributed_loads) <: AbstractDict ? distributed_loads : distributed_loads(time)
    pmass = typeof(point_masses) <: AbstractDict ? point_masses : point_masses(time)
    gvec = typeof(gravity) <: AbstractVector ? SVector{3}(gravity) : SVector{3}(gravity(time))
    ub_p = typeof(linear_displacement) <: AbstractVector ? SVector{3}(linear_displacement) : SVector{3}(linear_displacement(time))
    θb_p = typeof(angular_displacement) <: AbstractVector ? SVector{3}(angular_displacement) : SVector{3}(angular_displacement(time))
    vb_p = typeof(linear_velocity) <: AbstractVector ? SVector{3}(linear_velocity) : SVector{3}(linear_velocity(time))
    ωb_p = typeof(angular_velocity) <: AbstractVector ? SVector{3}(angular_velocity) : SVector{3}(angular_velocity(time))
    ab_p = typeof(linear_acceleration) <: AbstractVector ? SVector{3}(linear_acceleration) : SVector{3}(linear_acceleration(time))
    αb_p = typeof(angular_acceleration) <: AbstractVector ? SVector{3}(angular_acceleration) : SVector{3}(angular_acceleration(time))

    # indices corresponding to rigid body acceleration state variables
    icol_accel = body_frame_acceleration_indices(system, pcond)

    if constant_mass_matrix

        # solve for the system stiffness matrix
        expanded_system_jacobian!(K, x, indices, icol_accel, force_scaling, 
            structural_damping, assembly, pcond, dload, pmass, gvec, ub_p, θb_p, 
            vb_p, ωb_p, ab_p, αb_p)

        # solve for the system mass matrix
        expanded_system_mass_matrix!(M, indices, force_scaling, assembly, 
            pcond, pmass)

    else

        # solve for the system stiffness matrix
        # dynamic_system_jacobian!(K, FillArrays.Zeros(x), x, indices, icol_accel, force_scaling, 
        #     structural_damping, assembly, pcond, dload, pmass, gvec, ab_p, αb_p)
        steady_state_system_jacobian!(K, x, indices, icol_accel, force_scaling, 
            structural_damping, assembly, pcond, dload, pmass, gvec, ub_p, θb_p, 
            vb_p, ωb_p, ab_p, αb_p)

        # solve for the system mass matrix
        system_mass_matrix!(M, x, indices, force_scaling, assembly, pcond, pmass)

    end

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
            for k = 1:length(valM)
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
 - `linesearch = LineSearches.LineSearches.BackTracking(maxstep=1e6)`: Line search used to 
        solve nonlinear system of equations
 - `ftol = 1e-9`: tolerance for solving the nonlinear system of equations
 - `iterations = 1000`: maximum iterations for solving the nonlinear system of equations

# Linear Solver Keyword Arguments
 - `linearization_state`: Linearization state variables.  Defaults to zeros.
 - `update_linearization`: Flag indicating whether to update the linearization state 
        variables for a linear analysis with the instantaneous state variables.

# Eigenvalue Solution Keyword Arguments
 - `nev = 6`: Number of eigenvalues to compute
 - `steady_state = reset_state && !linear`: Flag indicating whether the steady state 
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

Pre-allocated version of `eigenvalue_analysis`.  Uses the state variables stored in
`system` as an initial guess for iterating to find the steady state solution.
"""
function eigenvalue_analysis!(system, assembly; 
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
    constant_mass_matrix=false,
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
    # eigenvalue solution keyword arguments
    nev = 6,
    steady_state=true,
    left=false,
    Uprev=nothing,
    )

    # check if provided system is consistent with provided keyword arguments
    constant_mass_matrix && @assert typeof(system) <: ExpandedSystem

    # reset state, if specified
    if reset_state
        reset_state!(system)
    end

    # set linearization state variables
    if steady_state

        if show_trace
            println("Finding a Steady-State Solution")
        end

        # perform a steady state analysis to find the linearization state variables
        system, converged = steady_state_analysis!(system, assembly;
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
            time=time,
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
    else
        # use specified linearization state variables
        if linear && !update_linearization
            if isnothing(linearization_state)
                system.x .= 0
            else
                system.x .= linearization_state
            end
        end

        # converged by default
        converged = true
    end

    # linearize the resulting system of equations
    x, K, M = linearize!(system, assembly; 
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
        time=time,
        structural_damping=structural_damping,
        constant_mass_matrix=constant_mass_matrix,
        show_trace=show_trace,
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
            time=time,
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

    # convert constant mass matrix system to a dynamic system
    if constant_mass_matrix
        # save provided system
        original_system = system
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

    # --- Determine whether all rigid body modes are constrained --- #

    # save original prescribed conditions
    original_pcond = pcond

    # find unconstrained rigid body modes
    not_constrained = reduce(.&, [p.isforce for (i,p) in pcond])

    # constrain unconstrained rigid body modes
    if any(not_constrained)
        # create copy of the original prescribed conditions
        pcond = copy(pcond)

        # unpack prescribed condition fields
        if haskey(pcond, 1)
            @unpack pd, pl, u, theta, F, M, Ff, Mf = pcond[1]
        else
            pd = @SVector zeros(Bool, 6)
            pl = @SVector ones(Bool, 6)
            u = @SVector zeros(Float64, 3)
            theta = @SVector zeros(Float64, 3)
            F = @SVector zeros(Float64, 3)
            M = @SVector zeros(Float64, 3)
            Ff = @SVector zeros(Float64, 3)
            Mf = @SVector zeros(Float64, 3)
        end

        # add constraints for unconstrained rigid body modes
        for i = 1:3
            # prescribe linear displacement, introduce linear acceleration state variable
            if not_constrained[i]
                pd = setindex(pd, true, i)
                u = setindex(u, u0[1][i], i)
            end
            # prescribe angular displacement, introduce angular acceleration state variable
            if not_constrained[3+i]
                pd = setindex(pd, true, 3+i)
                theta = setindex(theta, theta0[1][i], 3+i)
            end
        end

        # overwrite original prescribed condition
        pcond[1] = PrescribedConditions(pd, pl, u, theta, F, M, Ff, Mf)
    end

    # indices corresponding to rigid body acceleration state variables
    icol_accel = body_frame_acceleration_indices(system, pcond)

    # --- Determine whether Vdot and Ωdot are state variables for each point --- #

    # NOTE: Our system of equations cannot be solved for Vdot and Ωdot if 
    # 1. These terms are not used because the corresponding rows and columns of the 
    #    mass matrix are zero
    # 2. The compatibility equations do not impose any restrictions on the internal loads
    #    because the corresponding rows and columns of the compliance matrix are zero
    # In the former case, Vdot and Ωdot are not used.  Instead internal loads are 
    # solved using the equilibrium equations and used to define compatabile displacements.
    # In the latter case, the compatibility equations are solved to find compatable 
    # displacements and the equilibrium equations are used to find internal loads.
    # The equilibrium equations cannot therefore be used to solve for Vdot and Ωdot.
    # To determine whether either of these cases apply, we construct the array `rate_vars`
    # using a modified system mass matrix.  We then rely on the user to provide values for 
    # Vdot and and Ωdot which are compatible if either case applies.
    
    # save original beam elements
    original_elements = copy(assembly.elements)
    
    # combine the mass and compliance matrices of each beam element
    for i = 1:length(assembly.elements)
    
        # element properties
        L = assembly.elements[i].L
        xe = assembly.elements[i].x
        mass = assembly.elements[i].mass
        compliance = assembly.elements[i].compliance
        Cab = assembly.elements[i].Cab
        mu = assembly.elements[i].mu

        # modified element mass matrix
        mass = similar(mass) .= mass
        for j = 1:6
            if iszero(compliance[j,:])
                mass[j,:] .= 0.0
                mass[:,j] .= 0.0
            end
        end
        mass = SMatrix{6,6}(mass)

        # replace original element mass matrix
        assembly.elements[i] = Element(L, xe, compliance, mass, Cab, mu)
    end

    # construct `rate_vars` vector
    system_mass_matrix!(M, x, indices, force_scaling, assembly, pcond, pmass)
    rate_vars = .!(iszero.(sum(M, dims=1)))

    # restore original element mass matrices
    for i = 1:length(assembly.elements)
        assembly.elements[i] = original_elements[i]
    end

    # --- Define the residual and jacobian functions --- #

    f! = (resid, x) -> initial_condition_system_residual!(resid, x, indices, rate_vars, 
        icol_accel, force_scaling, structural_damping, assembly, pcond, dload, pmass, 
        gvec, ub_p, θb_p, vb_p, ωb_p, ab_p, αb_p, u0, theta0, V0, Omega0, Vdot0, Omegadot0)

    j! = (jacob, x) -> initial_condition_system_jacobian!(jacob, x, indices, rate_vars, 
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

        # set the convergence flag
        converged = result.f_converged
    end

    # --- Save state and rate variables associated with each point --- #

    for ipoint in eachindex(assembly.points)
        # save linear and angular displacement rates
        udot[ipoint], θdot[ipoint] = point_velocities(x, ipoint, indices.icol_point)
        # save linear and angular velocity rates
        tmp1, tmp2 = point_displacement_rates(x, ipoint, indices.icol_point, pcond)
        icol = indices.icol_point[ipoint]
        tmp1 = SVector{3}(
            rate_vars[icol+6] ? tmp1[1] : Vdot0[ipoint][1], 
            rate_vars[icol+7] ? tmp1[2] : Vdot0[ipoint][2],
            rate_vars[icol+8] ? tmp1[3] : Vdot0[ipoint][3],
        )
        tmp2 = SVector{3}(
            rate_vars[icol+9] ? tmp2[1] : Omegadot0[ipoint][1],
            rate_vars[icol+10] ? tmp2[2] : Omegadot0[ipoint][2],
            rate_vars[icol+11] ? tmp2[3] : Omegadot0[ipoint][3], 
        )
        Vdot[ipoint], Ωdot[ipoint] = tmp1, tmp2
        # save linear and angular displacement
        icol = indices.icol_point[ipoint]
        u, θ = point_displacement(x, ipoint, indices.icol_point, pcond)
        u = SVector{3}(
            rate_vars[icol+6] ? u0[ipoint][1] : u[1], 
            rate_vars[icol+7] ? u0[ipoint][2] : u[2],
            rate_vars[icol+8] ? u0[ipoint][3] : u[3],
        )
        θ = SVector{3}(
            rate_vars[icol+9] ? theta0[ipoint][1] : θ[1],
            rate_vars[icol+10] ? theta0[ipoint][2] : θ[2],
            rate_vars[icol+11] ? theta0[ipoint][3] : θ[3], 
        )
        set_linear_displacement!(system, pcond, u, ipoint)
        set_angular_displacement!(system, pcond, θ, ipoint)
        # save linear and angular velocity
        set_linear_velocity!(system, V0[ipoint], ipoint)
        set_angular_velocity!(system, Omega0[ipoint], ipoint)
    end

    # --- Restore original prescribed conditions --- #

    for ipoint in eachindex(assembly.points)
        # add rigid body accelerations to node accelerations
        Δx = assembly.points[ipoint]
        u, θ = point_displacement(x, ipoint, indices.icol_point, pcond)
        v, ω = body_frame_velocity(x)
        a, α = body_frame_acceleration(x, icol_accel, zero(ab_p), zero(αb_p))
        Vdot[ipoint] += a + cross(α, Δx) + cross(α, Δu)
        Ωdot[ipoint] += α
        # restore original state vector
        F, M = point_loads(x, ipoint, indices.icol_point, force_scaling)
        set_linear_displacement!(system, pcond, u, ipoint)
        set_angular_displacement!(system, pcond, θ, ipoint)
        set_external_forces!(system, pcond, F, ipoint)
        set_external_moments!(system, pcond, M, ipoint)
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
        ub_p = typeof(linear_displacement) <: AbstractVector ? SVector{3}(linear_displacement) : SVector{3}(linear_displacement(tvec[it]))
        θb_p = typeof(angular_displacement) <: AbstractVector ? SVector{3}(angular_displacement) : SVector{3}(angular_displacement(tvec[it]))    
        vb_p = typeof(linear_velocity) <: AbstractVector ? SVector{3}(linear_velocity) : SVector{3}(linear_velocity(tvec[it]))
        ωb_p = typeof(angular_velocity) <: AbstractVector ? SVector{3}(angular_velocity) : SVector{3}(angular_velocity(tvec[it]))
        ab_p = typeof(linear_acceleration) <: AbstractVector ? SVector{3}(linear_acceleration) : SVector{3}(linear_acceleration(tvec[it]))
        αb_p = typeof(angular_acceleration) <: AbstractVector ? SVector{3}(angular_acceleration) : SVector{3}(angular_acceleration(tvec[it]))

        # indices corresponding to rigid body acceleration state variables
        icol_accel = body_frame_acceleration_indices(system, pcond)

        # set current state rate initialization parameters
        ubdot = 2/dt*u + vb
        θbdot = 2/dt*θb + get_Qinv(θb)*get_C(θb)*ωb
        vbdot = 2/dt*v + ab
        ωbdot = 2/dt*ω + αb
        for ipoint = 1:length(assembly.points)
            # extract beam element state variables
            u, θ = point_displacement(x, ipoint, indices.icol_point, pcond)
            V, Ω = point_velocities(x, ipoint, indices.icol_point)
            # calculate state rate initialization terms, use storage for state rates
            udot[ipoint] = 2/dt*u + udot[ipoint]
            θdot[ipoint] = 2/dt*θ + θdot[ipoint]
            Vdot[ipoint] = 2/dt*V + Vdot[ipoint]
            Ωdot[ipoint] = 2/dt*Ω + Ωdot[ipoint]
        end

        # define the residual and jacobian functions
        f! = (r, x) -> newmark_system_residual!(r, x, indices, icol_accel, force_scaling, 
            structural_damping, assembly, pcond, dload, pmass, gvec, ub_p, θb_p, 
            vb_p, ωb_p, ab_p, αb_p, ubdot, θbdot, vbdot, ωbdot, udot, θdot, Vdot, Ωdot, dt)

        j! = (K, x) -> newmark_system_jacobian!(K, x, indices, icol_accel, force_scaling, 
            structural_damping, assembly, pcond, dload, pmass, gvec, ub_p, θb_p, 
            vb_p, ωb_p, ab_p, αb_p, ubdot, θbdot, vbdot, ωbdot, udot, θdot, Vdot, Ωdot, dt)

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