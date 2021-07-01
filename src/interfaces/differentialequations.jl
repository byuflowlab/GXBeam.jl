"""
    ODEProblem(system::GXBeam.System, assembly, tspan; kwargs...)

Construct a `ODEProblem` for the system of nonlinear beams
contained in `assembly` which may be used with the DifferentialEquations package.

Keyword Arguments:
 - `prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}()`:
        A dictionary with keys corresponding to the points at
        which prescribed conditions are applied and elements of type
        [`PrescribedConditions`](@ref) which describe the prescribed conditions
        at those points.  If time varying, this input may be provided as a
        function of time.
 - `distributed_loads = Dict{Int,DistributedLoads{Float64}}()`: A dictionary
        with keys corresponding to the elements to which distributed loads are
        applied and elements of type [`DistributedLoads`](@ref) which describe
        the distributed loads at those points.  If time varying, this input may
        be provided as a function of time.
 - `origin = zeros(3)`: Global frame origin vector. If time varying, this input
        may be provided as a function of time.
 - `linear_velocity = zeros(3)`: Global frame linear velocity vector. If time
        varying, this vector may be provided as a function of time.
 - `angular_velocity = zeros(3)`: Global frame angular velocity vector. If time
        varying, this vector may be provided as a function of time.
"""
function ODEProblem(system::System, assembly, tspan;
    prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}(),
    distributed_loads = Dict{Int,DistributedLoads{Float64}}(),
    origin = (@SVector zeros(3)),
    linear_velocity = (@SVector zeros(3)),
    angular_velocity = (@SVector zeros(3)),
    )

    # create ODEFunction
    func = ODEFunction(system, assembly)

    # use initial state from `system`
    u0 = copy(system.x)

    # set parameters
    p = (prescribed_conditions, distributed_loads, origin, linear_velocity, angular_velocity)

    return ODEProblem{true}(func, u0, tspan, p)
end

"""
    ODEFunction(system::GXBeam.System, assembly)

Construct a `ODEFunction` for the system of nonlinear beams
contained in `assembly` which may be used with the DifferentialEquations package.

The parameters associated with the resulting ODEFunction are defined by the tuple
`(prescribed_conditions, distributed_loads, origin, linear_velocity, angular_velocity)`
where each parameter is defined as follows:
 - `prescribed_conditions`: A dictionary with keys corresponding to the points at
        which prescribed conditions are applied and elements of type
        [`PrescribedConditions`](@ref) which describe the prescribed conditions
        at those points.  If time varying, this input may be provided as a
        function of time.
 - `distributed_loads`: A dictionary with keys corresponding to the elements to
        which distributed loads are applied and elements of type
        [`DistributedLoads`](@ref) which describe the distributed loads at those
        points.  If time varying, this input may be provided as a function of
        time.
 - `origin`: Global frame origin vector. If time varying, this input
        may be provided as a function of time.
 - `linear_velocity`: Global frame linear velocity vector. If time
        varying, this vector may be provided as a function of time.
 - `angular_velocity`: Global frame angular velocity vector. If time
        varying, this vector may be provided as a function of time.
"""
function ODEFunction(system::System, assembly)

    # check to make sure the system isn't static
    @assert !system.static

    # unpack system pointers
    irow_pt = system.irow_pt
    irow_beam = system.irow_beam
    irow_beam1 = system.irow_beam1
    irow_beam2 = system.irow_beam2
    icol_pt = system.icol_pt
    icol_beam = system.icol_beam

    # unpack scaling parameters
    force_scaling = system.force_scaling
    mass_scaling = system.mass_scaling

    # DAE function
    f = function(resid, u, p, t)

        # get current parameters
        prescribed_conditions = typeof(p[1]) <: AbstractDict ? p[1] : p[1](t)
        distributed_loads = typeof(p[2]) <: AbstractDict ? p[2] : p[2](t)
        x0 = typeof(p[3]) <: AbstractVector ? SVector{3}(p[3]) : SVector{3}(p[3](t))
        v0 = typeof(p[4]) <: AbstractVector ? SVector{3}(p[4]) : SVector{3}(p[4](t))
        ω0 = typeof(p[5]) <: AbstractVector ? SVector{3}(p[5]) : SVector{3}(p[5](t))

        # calculate residual
        steady_state_system_residual!(resid, u, assembly, prescribed_conditions,
            distributed_loads, force_scaling, mass_scaling, irow_pt, irow_beam, irow_beam1,
            irow_beam2, icol_pt, icol_beam, x0, v0, ω0)

        return resid
    end

    update_mass_matrix! = function(M, u, p, t)

        # zero out all mass matrix entries
        M .= 0.0

        # get current parameters
        prescribed_conditions = typeof(p[1]) <: AbstractDict ? p[1] : p[1](t)
        distributed_loads = typeof(p[2]) <: AbstractDict ? p[2] : p[2](t)
        x0 = typeof(p[3]) <: AbstractVector ? SVector{3}(p[3]) : SVector{3}(p[3](t))
        v0 = typeof(p[4]) <: AbstractVector ? SVector{3}(p[4]) : SVector{3}(p[4](t))
        ω0 = typeof(p[5]) <: AbstractVector ? SVector{3}(p[5]) : SVector{3}(p[5](t))

        # calculate mass matrix
        system_mass_matrix!(M, u, assembly, force_scaling, mass_scaling,
            irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam)

        return M
    end

    mass_matrix = DiffEqArrayOperator(copy(system.M), update_func = update_mass_matrix!)

    # jacobian function with respect to states/state rates
    jac = function(J, u, p, t)

        # zero out all jacobian entries
        J .= 0.0

        # get current parameters
        prescribed_conditions = typeof(p[1]) <: AbstractDict ? p[1] : p[1](t)
        distributed_loads = typeof(p[2]) <: AbstractDict ? p[2] : p[2](t)
        x0 = typeof(p[3]) <: AbstractVector ? SVector{3}(p[3]) : SVector{3}(p[3](t))
        v0 = typeof(p[4]) <: AbstractVector ? SVector{3}(p[4]) : SVector{3}(p[4](t))
        ω0 = typeof(p[5]) <: AbstractVector ? SVector{3}(p[5]) : SVector{3}(p[5](t))

        # calculate jacobian
        steady_state_system_jacobian!(J, u, assembly, prescribed_conditions,
            distributed_loads, force_scaling, mass_scaling, irow_pt, irow_beam, irow_beam1,
            irow_beam2, icol_pt, icol_beam, x0, v0, ω0)

        return J
    end

    # sparsity structure
    sparsity = get_sparsity(system, assembly)

    # jacobian prototype (use dense since sparse isn't working)
    jac_prototype = collect(system.K)

    # TODO: figure out how to use a sparse matrix here.
    # It's failing with a singular exception during the LU factorization.

    return ODEFunction{true,true}(f; mass_matrix, jac, jac_prototype, sparsity)
end

"""
    DAEProblem(system::GXBeam.System, assembly, tspan; kwargs...)

Construct a `DAEProblem` for the system of nonlinear beams
contained in `assembly` which may be used with the DifferentialEquations package.

A consistent set of initial conditions may be obtained prior to constructing the
DAEProblem using [`initial_condition_analysis!`](@ref) or by constructing a
DAEProblem after a time domain analysis.

Keyword Arguments:
 - `prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}()`:
        A dictionary with keys corresponding to the points at
        which prescribed conditions are applied and elements of type
        [`PrescribedConditions`](@ref) which describe the prescribed conditions
        at those points.  If time varying, this input may be provided as a
        function of time.
 - `distributed_loads = Dict{Int,DistributedLoads{Float64}}()`: A dictionary
        with keys corresponding to the elements to which distributed loads are
        applied and elements of type [`DistributedLoads`](@ref) which describe
        the distributed loads at those points.  If time varying, this input may
        be provided as a function of time.
 - `origin = zeros(3)`: Global frame origin vector. If time varying, this input
        may be provided as a function of time.
 - `linear_velocity = zeros(3)`: Global frame linear velocity vector. If time
        varying, this vector may be provided as a function of time.
 - `angular_velocity = zeros(3)`: Global frame angular velocity vector. If time
        varying, this vector may be provided as a function of time.
"""
function DAEProblem(system::System, assembly, tspan;
    prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}(),
    distributed_loads = Dict{Int,DistributedLoads{Float64}}(),
    origin = (@SVector zeros(3)),
    linear_velocity = (@SVector zeros(3)),
    angular_velocity = (@SVector zeros(3)),
    )

    # create DAEFunction
    func = DAEFunction(system, assembly)

    # use initial state from `system`
    u0 = copy(system.x)

    # use initial state rates from `system`
    du0 = zero(u0)
    for (ibeam, icol) in enumerate(system.icol_beam)
        du0[icol:icol+2] = system.udot[ibeam]
        du0[icol+3:icol+5] = system.θdot[ibeam]
        du0[icol+12:icol+14] = system.Pdot[ibeam]
        du0[icol+15:icol+17] = system.Hdot[ibeam]
    end

    # set parameters
    p = (prescribed_conditions, distributed_loads, origin, linear_velocity, angular_velocity)

    # get differential variables
    differential_vars = get_differential_vars(system)

    return DAEProblem{true}(func, du0, u0, tspan, p; differential_vars)
end

"""
    DAEFunction(system::GXBeam.System, assembly)

Construct a `DAEFunction` for the system of nonlinear beams
contained in `assembly` which may be used with the DifferentialEquations package.

The parameters associated with the resulting DAEFunction are defined by the tuple
`(prescribed_conditions, distributed_loads, origin, linear_velocity, angular_velocity)`
where each parameter is defined as follows:
 - `prescribed_conditions`: A dictionary with keys corresponding to the points at
        which prescribed conditions are applied and elements of type
        [`PrescribedConditions`](@ref) which describe the prescribed conditions
        at those points.  If time varying, this input may be provided as a
        function of time.
 - `distributed_loads`: A dictionary with keys corresponding to the elements to
        which distributed loads are applied and elements of type
        [`DistributedLoads`](@ref) which describe the distributed loads at those
        points.  If time varying, this input may be provided as a function of
        time.
 - `origin`: Global frame origin vector. If time varying, this input
        may be provided as a function of time.
 - `linear_velocity`: Global frame linear velocity vector. If time
        varying, this vector may be provided as a function of time.
 - `angular_velocity`: Global frame angular velocity vector. If time
        varying, this vector may be provided as a function of time.
"""
function DAEFunction(system::System, assembly)

    # check to make sure the system isn't static
    @assert !system.static

    # unpack system pointers
    irow_pt = system.irow_pt
    irow_beam = system.irow_beam
    irow_beam1 = system.irow_beam1
    irow_beam2 = system.irow_beam2
    icol_pt = system.icol_pt
    icol_beam = system.icol_beam

    # unpack scaling parameters
    force_scaling = system.force_scaling
    mass_scaling = system.mass_scaling

    # DAE function
    f = function(resid, du, u, p, t)

        # get current parameters
        prescribed_conditions = typeof(p[1]) <: AbstractDict ? p[1] : p[1](t)
        distributed_loads = typeof(p[2]) <: AbstractDict ? p[2] : p[2](t)
        x0 = typeof(p[3]) <: AbstractVector ? SVector{3}(p[3]) : SVector{3}(p[3](t))
        v0 = typeof(p[4]) <: AbstractVector ? SVector{3}(p[4]) : SVector{3}(p[4](t))
        ω0 = typeof(p[5]) <: AbstractVector ? SVector{3}(p[5]) : SVector{3}(p[5](t))

        # calculate residual
        dynamic_system_residual!(resid, u, du, assembly, prescribed_conditions,
            distributed_loads, force_scaling, mass_scaling, irow_pt, irow_beam, irow_beam1,
            irow_beam2, icol_pt, icol_beam, x0, v0, ω0)

        return resid
    end

    # jacobian function with respect to states/state rates
    jac = function(J, du, u, p, gamma, t)

        # zero out all jacobian entries
        J .= 0.0

        # get current parameters
        prescribed_conditions = typeof(p[1]) <: AbstractDict ? p[1] : p[1](t)
        distributed_loads = typeof(p[2]) <: AbstractDict ? p[2] : p[2](t)
        x0 = typeof(p[3]) <: AbstractVector ? SVector{3}(p[3]) : SVector{3}(p[3](t))
        v0 = typeof(p[4]) <: AbstractVector ? SVector{3}(p[4]) : SVector{3}(p[4](t))
        ω0 = typeof(p[5]) <: AbstractVector ? SVector{3}(p[5]) : SVector{3}(p[5](t))

        # calculate jacobian
        dynamic_system_jacobian!(J, u, du, assembly, prescribed_conditions,
            distributed_loads, force_scaling, mass_scaling, irow_pt, irow_beam, irow_beam1,
            irow_beam2, icol_pt, icol_beam, x0, v0, ω0)

        # add gamma multiplied by the mass matrix
        system_mass_matrix!(J, gamma, u, assembly, force_scaling, mass_scaling,
            irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam)

        return J
    end

    # sparsity structure
    sparsity = get_sparsity(system, assembly)

    # jacobian prototype (use dense since sparse isn't working)
    # jac_prototype = collect(system.K)

    # TODO: figure out how to use a sparse matrix here.
    # It's failing with a singular exception during the LU factorization.
    # Using `jac_prototype` also causes errors

    return DAEFunction{true,true}(f; jac, sparsity)
end

function get_differential_vars(system::System)
    differential_vars = fill(false, length(system.x))
    for icol in system.icol_beam
        differential_vars[icol:icol+2] .= true # u (for the beam)
        differential_vars[icol+3:icol+5] .= true # θ (for the beam)
        differential_vars[icol+12:icol+14] .= true # P (for the beam)
        differential_vars[icol+15:icol+17] .= true # H (for the beam)
    end
    return differential_vars
end
