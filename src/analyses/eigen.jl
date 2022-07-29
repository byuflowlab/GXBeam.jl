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
        expanded_steady_system_jacobian!(K, x, indices, icol_accel, force_scaling, 
            structural_damping, assembly, pcond, dload, pmass, gvec, ub_p, θb_p, 
            vb_p, ωb_p, ab_p, αb_p)

        # solve for the system mass matrix
        expanded_steady_system_mass_matrix!(M, indices, force_scaling, assembly, 
            pcond, pmass)

    else

        # solve for the system stiffness matrix
        steady_state_system_jacobian!(K, x, indices, icol_accel, force_scaling, 
            structural_damping, assembly, pcond, dload, pmass, gvec, ub_p, θb_p, 
            vb_p, ωb_p, ab_p, αb_p)

        # solve for the system mass matrix
        steady_system_mass_matrix!(M, x, indices, force_scaling, assembly, pcond, pmass)

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
        solve the nonlinear system of equations
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

Pre-allocated version of `eigenvalue_analysis`.
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
