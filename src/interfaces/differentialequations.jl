"""
    ODEProblem(system::GXBeam.System, assembly, tspan; kwargs...)

Construct a `ODEProblem` for the system of nonlinear beams contained in `assembly` which 
may be used with the DifferentialEquations package.
"""
function SciMLBase.ODEProblem(system::System, assembly, tspan; kwargs...)

    # create ODEFunction
    func = SciMLBase.ODEFunction(system, assembly; kwargs...)

    # use initial state from `system`
    u0 = copy(system.x)

    return SciMLBase.ODEProblem{true}(func, u0, tspan)
end

"""
    ODEFunction(system::GXBeam.System, assembly; structural_damping=true)

Construct a `ODEFunction` for the system of nonlinear beams
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
        the distributed loads on those elements.  If time varying, this input may
        be provided as a function of time.
 - `point_masses = Dict{Int,Vector{PointMass{Float64}}}()`: A dictionary with keys 
        corresponding to the points at which point masses are attached and values 
        containing vectors of objects of type [`PointMass`](@ref) which describe 
        the point masses attached at those points.  If time varying, this input may
        be provided as a function of time.
 - `structural_damping = false`: Flag indicating whether structural damping should be enabled
 - `gravity`: Gravity vector. If time varying, this input may be provided as a 
        function of time.
 - `origin = zeros(3)`: Global frame origin vector. If time varying, this input
        may be provided as a function of time.
 - `linear_velocity = zeros(3)`: Global frame linear velocity vector. If time
        varying, this vector may be provided as a function of time.
 - `angular_velocity = zeros(3)`: Global frame angular velocity vector. If time
        varying, this vector may be provided as a function of time.
 - `linear_acceleration = zeros(3)`: Global frame linear acceleration vector. If time
       varying, this vector may be provided as a function of time.
 - `angular_acceleration = zeros(3)`: Global frame angular acceleration vector. If time
       varying, this vector may be provided as a function of time.
"""
function SciMLBase.ODEFunction(system::System, assembly; 
    prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}(),
    distributed_loads = Dict{Int,DistributedLoads{Float64}}(),
    point_masses = Dict{Int,Vector{PointMass{Float64}}}(),
    structural_damping=false,
    gravity = (@SVector zeros(3)),
    origin = (@SVector zeros(3)),
    linear_velocity = (@SVector zeros(3)),
    angular_velocity = (@SVector zeros(3)),
    linear_acceleration = (@SVector zeros(3)),
    angular_acceleration = (@SVector zeros(3)),
    )

    @unpack dynamic_indices, force_scaling = system

    # DAE function
    f = function(resid, u, p, t)

        # get current parameters
        pcond = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(t)
        dload = typeof(distributed_loads) <: AbstractDict ? distributed_loads : distributed_loads(t)
        pmass = typeof(point_masses) <: AbstractDict ? point_masses : point_masses(t)
        gvec = typeof(gravity) <: AbstractVector ? SVector{3}(gravity) : SVector{3}(gravity(t))
        x0 = typeof(origin) <: AbstractVector ? SVector{3}(origin) : SVector{3}(origin(t))
        v0 = typeof(linear_velocity) <: AbstractVector ? SVector{3}(linear_velocity) : SVector{3}(linear_velocity(t))
        ω0 = typeof(angular_velocity) <: AbstractVector ? SVector{3}(angular_velocity) : SVector{3}(angular_velocity(t))
        a0 = typeof(linear_acceleration) <: AbstractVector ? SVector{3}(linear_acceleration) : SVector{3}(linear_acceleration(t))
        α0 = typeof(angular_acceleration) <: AbstractVector ? SVector{3}(angular_acceleration) : SVector{3}(angular_acceleration(t))

        # calculate residual
        steady_state_system_residual!(resid, u, dynamic_indices, force_scaling, 
            assembly, pcond, dload, pmass, gvec, x0, v0, ω0, a0, α0)

        return resid
    end

    update_mass_matrix! = function(M, u, p, t)

        # get current parameters
        pcond = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(t)
        pmass = typeof(point_masses) <: AbstractDict ? point_masses : point_masses(t)

        # zero out all mass matrix entries
        M .= 0.0

        # calculate mass matrix
        system_mass_matrix!(M, u, dynamic_indices, force_scaling, structural_damping, 
                assembly, pcond, pmass)

        return M
    end

    update_jacobian! = function(J, u, p, t)

        # get current parameters
        pcond = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(t)
        dload = typeof(distributed_loads) <: AbstractDict ? distributed_loads : distributed_loads(t)
        pmass = typeof(point_masses) <: AbstractDict ? point_masses : point_masses(t)
        gvec = typeof(gravity) <: AbstractVector ? SVector{3}(gravity) : SVector{3}(gravity(t))
        x0 = typeof(origin) <: AbstractVector ? SVector{3}(origin) : SVector{3}(origin(t))
        v0 = typeof(linear_velocity) <: AbstractVector ? SVector{3}(linear_velocity) : SVector{3}(linear_velocity(t))
        ω0 = typeof(angular_velocity) <: AbstractVector ? SVector{3}(angular_velocity) : SVector{3}(angular_velocity(t))
        a0 = typeof(linear_acceleration) <: AbstractVector ? SVector{3}(linear_acceleration) : SVector{3}(linear_acceleration(t))
        α0 = typeof(angular_acceleration) <: AbstractVector ? SVector{3}(angular_acceleration) : SVector{3}(angular_acceleration(t))

        # zero out all jacobian entries
        J .= 0.0

        # calculate jacobian
        steady_state_system_jacobian!(J, u, dynamic_indices, force_scaling, 
              assembly, pcond, dload, pmass, gvec, x0, v0, ω0, a0, α0)

        return J
    end

    return SciMLBase.ODEFunction{true,true}(f; 
        mass_matrix = SciMLBase.DiffEqArrayOperator(system.M, update_func = update_mass_matrix!),
        jac = update_jacobian!,
        jac_prototype = typeof(system.K)
        )
end

"""
    DAEProblem(system::GXBeam.System, assembly, tspan; kwargs...)

Construct a `DAEProblem` for the system of nonlinear beams contained in `assembly` which 
may be used with the DifferentialEquations package.
"""
function SciMLBase.DAEProblem(system::System, assembly, tspan; 
    prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}(),
    kwargs...)

    # create SciMLBase.DAEFunction
    func = SciMLBase.DAEFunction(system, assembly; prescribed_conditions, kwargs...)

    # use initial state from `system`
    u0 = copy(system.x)

    # use initial state rates from `system`
    du0 = zero(u0)
    for (ipoint, icol) in enumerate(system.dynamic_indices.icol_point)
        du0[icol:icol+2] = system.udot[ipoint]
        du0[icol+3:icol+5] = system.θdot[ipoint]
        du0[icol+6:icol+8] = system.Vdot[ipoint]
        du0[icol+9:icol+11] = system.Ωdot[ipoint]
    end

    # get differential variables
    differential_vars = get_differential_vars(system, assembly, prescribed_conditions)

    return SciMLBase.DAEProblem{true}(func, du0, u0, tspan; differential_vars)
end

"""
    DAEFunction(system::GXBeam.System, assembly; kwargs...)

Construct a `DAEFunction` for the system of nonlinear beams
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
        the distributed loads on those elements.  If time varying, this input may
        be provided as a function of time.
 - `point_masses = Dict{Int,Vector{PointMass{Float64}}}()`: A dictionary with keys 
        corresponding to the points at which point masses are attached and values 
        containing vectors of objects of type [`PointMass`](@ref) which describe 
        the point masses attached at those points.  If time varying, this input may
        be provided as a function of time.
 - `structural_damping = false`: Flag indicating whether structural damping should be enabled
 - `gravity`: Gravity vector. If time varying, this input may be provided as a 
        function of time.
 - `origin = zeros(3)`: Global frame origin vector. If time varying, this input
        may be provided as a function of time.
 - `linear_velocity = zeros(3)`: Global frame linear velocity vector. If time
        varying, this vector may be provided as a function of time.
 - `angular_velocity = zeros(3)`: Global frame angular velocity vector. If time
        varying, this vector may be provided as a function of time.
 - `linear_acceleration = zeros(3)`: Global frame linear acceleration vector. If time
       varying, this vector may be provided as a function of time.
 - `angular_acceleration = zeros(3)`: Global frame angular acceleration vector. If time
       varying, this vector may be provided as a function of time.
"""
function SciMLBase.DAEFunction(system::System, assembly; 
    prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}(),
    distributed_loads = Dict{Int,DistributedLoads{Float64}}(),
    point_masses = Dict{Int,Vector{PointMass{Float64}}}(),
    structural_damping = true,
    gravity = (@SVector zeros(3)),
    origin = (@SVector zeros(3)),
    linear_velocity = (@SVector zeros(3)),
    angular_velocity = (@SVector zeros(3)),
    linear_acceleration = (@SVector zeros(3)),
    angular_acceleration = (@SVector zeros(3)),
    )

    # unpack system pointers
    @unpack dynamic_indices, force_scaling = system

    # DAE function
    f = function(resid, du, u, p, t)

        # get current parameters
        pcond = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(t)
        dload = typeof(distributed_loads) <: AbstractDict ? distributed_loads : distributed_loads(t)
        pmass = typeof(point_masses) <: AbstractDict ? point_masses : point_masses(t)
        gvec = typeof(gravity) <: AbstractVector ? SVector{3}(gravity) : SVector{3}(gravity(t))
        x0 = typeof(origin) <: AbstractVector ? SVector{3}(origin) : SVector{3}(origin(t))
        v0 = typeof(linear_velocity) <: AbstractVector ? SVector{3}(linear_velocity) : SVector{3}(linear_velocity(t))
        ω0 = typeof(angular_velocity) <: AbstractVector ? SVector{3}(angular_velocity) : SVector{3}(angular_velocity(t))
        a0 = typeof(linear_acceleration) <: AbstractVector ? SVector{3}(linear_acceleration) : SVector{3}(linear_acceleration(t))
        α0 = typeof(angular_acceleration) <: AbstractVector ? SVector{3}(angular_acceleration) : SVector{3}(angular_acceleration(t))

        # calculate residual
        dynamic_system_residual!(resid, du, u, dynamic_indices, force_scaling, structural_damping, 
            assembly, pcond, dload, pmass, gvec, x0, v0, ω0, a0, α0)

        return resid
    end

    # jacobian function with respect to states/state rates
    update_jacobian! = function(J, du, u, p, gamma, t)

        # zero out all jacobian entries
        J .= 0.0

        # get current parameters
        pcond = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(t)
        dload = typeof(distributed_loads) <: AbstractDict ? distributed_loads : distributed_loads(t)
        pmass = typeof(point_masses) <: AbstractDict ? point_masses : point_masses(t)
        gvec = typeof(gravity) <: AbstractVector ? SVector{3}(gravity) : SVector{3}(gravity(t))
        x0 = typeof(origin) <: AbstractVector ? SVector{3}(origin) : SVector{3}(origin(t))
        v0 = typeof(linear_velocity) <: AbstractVector ? SVector{3}(linear_velocity) : SVector{3}(linear_velocity(t))
        ω0 = typeof(angular_velocity) <: AbstractVector ? SVector{3}(angular_velocity) : SVector{3}(angular_velocity(t))
        a0 = typeof(linear_acceleration) <: AbstractVector ? SVector{3}(linear_acceleration) : SVector{3}(linear_acceleration(t))
        α0 = typeof(angular_acceleration) <: AbstractVector ? SVector{3}(angular_acceleration) : SVector{3}(angular_acceleration(t))

        # calculate jacobian
        dynamic_system_jacobian!(J, du, u, dynamic_indices, force_scaling, structural_damping, 
            assembly, pcond, dload, pmass, gvec, x0, v0, ω0, a0, α0)

        # add gamma multiplied by the mass matrix
        system_mass_matrix!(J, gamma, u, dynamic_indices, force_scaling, structural_damping, 
            assembly, pcond, pmass)

        return J
    end

    return SciMLBase.DAEFunction{true,true}(f) # TODO: re-add jacobian here once supported
end

function get_differential_vars(system::System, assembly, prescribed_conditions)

    pcond = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(0.0)

    differential_vars = fill(false, length(system.x))

    for (ipoint, icol) in enumerate(system.dynamic_indices.icol_point)

        if haskey(pcond, ipoint)
            # displacements are differential variables, forces and moments are not
            differential_vars[icol:icol+5] .= pcond[ipoint].isforce
            # velocities are differential variables
            differential_vars[icol+6:icol+11] .= true
        else
            # displacements and velocities are differentiable variables
            differential_vars[icol:icol+11] .= true
        end
    end

    return differential_vars
end
