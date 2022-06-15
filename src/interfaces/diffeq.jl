"""
    ODEProblem(system::GXBeam.AbstractSystem, assembly, tspan; kwargs...)

Construct a `ODEProblem` for the system of nonlinear beams contained in `assembly` which 
may be used with the DifferentialEquations package.
"""
function SciMLBase.ODEProblem(system::AbstractSystem, assembly, tspan; 
    prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}(),
    distributed_loads = Dict{Int,DistributedLoads{Float64}}(),
    point_masses = Dict{Int,Vector{PointMass{Float64}}}(),
    gravity = (@SVector zeros(3)),
    origin = (@SVector zeros(3)),
    linear_velocity = (@SVector zeros(3)),
    angular_velocity = (@SVector zeros(3)),
    linear_acceleration = (@SVector zeros(3)),
    angular_acceleration = (@SVector zeros(3)),
    constant_mass_matrix = typeof(system) <: ExpandedSystem,
    structural_damping = true,
    )

    # original_system = system

    # if constant_mass_matrix
    #     # check if provided system is a constant mass matrix system
    #     if typeof(original_system) <: ExpandedSystem
    #         # use provided constant mass matrix system for the analysis
    #         system = original_system
    #     else
    #         # construct a constant mass matrix system for the analysis
    #         system = ExpandedSystem(assembly; force_scaling = system.force_scaling)
    #         # copy state variables from the original system to the constant mass matrix system
    #         copy_state!(system, original_system, assembly;     
    #             prescribed_conditions=prescribed_conditions,
    #             distributed_loads=distributed_loads,
    #             point_masses=point_masses,
    #             origin=origin,
    #             linear_velocity=linear_velocity,
    #             angular_velocity=angular_velocity,
    #             linear_acceleration=linear_acceleration,
    #             angular_acceleration=angular_acceleration,
    #             gravity=gravity,
    #             time=tspan[1])
    #     end
    # else
    #     # check if provided system is a dynamic system
    #     if typeof(original_system) <: DynamicSystem
    #         # use provided dynamic system for the analysis
    #         system = original_system
    #     else
    #         # construct a dynamic system for the analysis
    #         system = DynamicSystem(assembly; force_scaling = system.force_scaling)
    #         # copy state variables from the original system to the dynamic system
    #         copy_state!(system, original_system, assembly;     
    #             prescribed_conditions=prescribed_conditions,
    #             distributed_loads=distributed_loads,
    #             point_masses=point_masses,
    #             origin=origin,
    #             linear_velocity=linear_velocity,
    #             angular_velocity=angular_velocity,
    #             linear_acceleration=linear_acceleration,
    #             angular_acceleration=angular_acceleration,
    #             gravity=gravity,
    #             time=tspan[1])
    #     end
    # end

    # set initial state variables
    u0 = copy(system.x)

    # construct ODEFunction
    func = SciMLBase.ODEFunction(system, assembly; 
        prescribed_conditions = prescribed_conditions,
        distributed_loads = distributed_loads,
        point_masses = point_masses,
        gravity = gravity,
        origin = origin,
        linear_velocity = linear_velocity,
        angular_velocity = angular_velocity,
        linear_acceleration = linear_acceleration,
        angular_acceleration = angular_acceleration,
        constant_mass_matrix = constant_mass_matrix,
        structural_damping = structural_damping)

    return SciMLBase.ODEProblem{true}(func, u0, tspan)
end

"""
    ODEFunction(system::GXBeam.AbstractSystem, assembly; kwargs...)

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
 - `structural_damping = true`: Flag indicating whether structural damping should be enabled
 - `constant_mass_matrix = true`: Flag indicating whether to use a constant mass matrix.  
"""
function SciMLBase.ODEFunction(system::AbstractSystem, assembly; 
    prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}(),
    distributed_loads = Dict{Int,DistributedLoads{Float64}}(),
    point_masses = Dict{Int,Vector{PointMass{Float64}}}(),
    gravity = (@SVector zeros(3)),
    origin = (@SVector zeros(3)),
    linear_velocity = (@SVector zeros(3)),
    angular_velocity = (@SVector zeros(3)),
    linear_acceleration = (@SVector zeros(3)),
    angular_acceleration = (@SVector zeros(3)),
    constant_mass_matrix = typeof(system) <: ExpandedSystem,
    structural_damping = true,
    )

    force_scaling = system.force_scaling

    if constant_mass_matrix

        indices = SystemIndices(assembly.start, assembly.stop, static=false, expanded=true)

        # residual function (constant mass matrix system)
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
            expanded_system_residual!(resid, u, indices, force_scaling, structural_damping,
                assembly, pcond, dload, pmass, gvec, x0, v0, ω0, a0, α0)

            return resid
        end

        # mass matrix (constant mass matrix system)
        TF = eltype(system)
        nx = indices.nstates
        mass_matrix = spzeros(TF, nx, nx)
        pcond = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(0)
        pmass = typeof(point_masses) <: AbstractDict ? point_masses : point_masses(0)
        expanded_system_mass_matrix!(mass_matrix, -1, indices, force_scaling, assembly, pcond, pmass) 
    
        # jacobian
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
            expanded_system_jacobian!(J, u, indices, force_scaling, structural_damping, 
                assembly, pcond, dload, pmass, gvec, x0, v0, ω0, a0, α0)

            return J
        end

    else

        indices = SystemIndices(assembly.start, assembly.stop, static=false, expanded=false)

        # residual function
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
            steady_state_system_residual!(resid, u, indices, force_scaling, structural_damping,
                assembly, pcond, dload, pmass, gvec, x0, v0, ω0, a0, α0)

            return resid
        end

        # mass matrix
        update_mass_matrix! = function(M, u, p, t)

            # get current parameters
            pcond = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(t)
            pmass = typeof(point_masses) <: AbstractDict ? point_masses : point_masses(t)
    
            # zero out all mass matrix entries
            M .= 0.0
    
            # calculate mass matrix
            system_mass_matrix!(M, u, indices, force_scaling, assembly, pcond, pmass)
    
            M .*= -1
    
            return M
        end

        mass_matrix = SciMLBase.DiffEqArrayOperator(system.M, update_func = update_mass_matrix!)

        # jacobian
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
            steady_state_system_jacobian!(J, u, indices, force_scaling, structural_damping, 
                assembly, pcond, dload, pmass, gvec, x0, v0, ω0, a0, α0)

            return J
        end

    end

    return SciMLBase.ODEFunction{true,true}(f; 
        mass_matrix = mass_matrix,
        jac = update_jacobian!
        )
end

"""
    DAEProblem(system::GXBeam.DynamicSystem, assembly, tspan; kwargs...)

Construct a `DAEProblem` for the system of nonlinear beams contained in `assembly` which 
may be used with the DifferentialEquations package.
"""
function SciMLBase.DAEProblem(system::DynamicSystem, assembly, tspan; 
    prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}(),
    point_masses = Dict{Int,PointMass{Float64}}(),
    kwargs...)

    # create SciMLBase.DAEFunction
    func = SciMLBase.DAEFunction(system, assembly; 
        prescribed_conditions = prescribed_conditions,
        point_masses = point_masses,
        kwargs...)

    # use initial state from `system`
    u0 = copy(system.x)

    # use initial state rates from `system`
    du0 = zero(u0)
    for (ipoint, icol) in enumerate(system.indices.icol_point)
        du0[icol:icol+2] = system.udot[ipoint]
        du0[icol+3:icol+5] = system.θdot[ipoint]
        du0[icol+6:icol+8] = system.Vdot[ipoint]
        du0[icol+9:icol+11] = system.Ωdot[ipoint]
    end

    # get differential variables
    differential_vars = get_differential_vars(system, assembly, prescribed_conditions, point_masses)

    return SciMLBase.DAEProblem{true}(func, du0, u0, tspan; differential_vars)
end

"""
    DAEFunction(system::GXBeam.DynamicSystem, assembly; kwargs...)

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
 - `structural_damping = true`: Flag indicating whether structural damping should be enabled
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
function SciMLBase.DAEFunction(system::DynamicSystem, assembly; 
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
    @unpack force_scaling, indices = system

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
        dynamic_system_residual!(resid, du, u, indices, force_scaling, structural_damping, 
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
        dynamic_system_jacobian!(J, du, u, indices, force_scaling, structural_damping, 
            assembly, pcond, dload, pmass, gvec, x0, v0, ω0, a0, α0)

        # add gamma multiplied by the mass matrix
        system_mass_matrix!(J, gamma, u, indices, force_scaling, 
            assembly, pcond, pmass)

        return J
    end

    return SciMLBase.DAEFunction{true,true}(f) # TODO: re-add jacobian here once supported
end

function get_differential_vars(system::DynamicSystem, assembly, prescribed_conditions, point_masses)

    # unpack pre-allocated storage and pointers for system
    @unpack x, M, force_scaling, indices, t = system

    # current parameters
    pcond = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(t)
    pmass = typeof(point_masses) <: AbstractDict ? point_masses : point_masses(t)

    # solve for the system mass matrix
    system_mass_matrix!(M, x, indices, force_scaling, assembly, pcond, pmass)

    # identify differential variables
    differential_vars = dropdims(.!(iszero.(sum(M, dims=1))), dims=1)

    return differential_vars
end
