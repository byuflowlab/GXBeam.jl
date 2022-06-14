"""
    ODEProblem(system::GXBeam.AbstractSystem, assembly, tspan; kwargs...)

Construct a `ODEProblem` for the system of nonlinear beams contained in `assembly` which 
may be used with the DifferentialEquations package.
"""
function SciMLBase.ODEProblem(system::AbstractSystem, assembly, tspan; 
    prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}(),
    constant_mass_matrix = true, 
    kwargs...)

    # create ODEFunction
    func = SciMLBase.ODEFunction(system, assembly; 
        prescribed_conditions = prescribed_conditions,
        constant_mass_matrix = constant_mass_matrix, 
        kwargs...)

    if constant_mass_matrix
        u0 = get_expanded_state(system, assembly, prescribed_conditions)
    else
        u0 = copy(system.x)        
    end

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
    structural_damping=true,
    constant_mass_matrix=true,
    )

    @unpack indices, force_scaling = system

    if constant_mass_matrix

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
            expanded_system_residual!(resid, u, expanded_indices, force_scaling, structural_damping,
                assembly, pcond, dload, pmass, gvec, x0, v0, ω0, a0, α0)

            return resid
        end

        # mass matrix (constant mass matrix system)
        TF = eltype(system)
        nx = expanded_indices.nstates
        mass_matrix = spzeros(TF, nx, nx)
        pcond = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(0)
        pmass = typeof(point_masses) <: AbstractDict ? point_masses : point_masses(0)
        expanded_system_mass_matrix!(mass_matrix, -1, expanded_indices, force_scaling, assembly, pcond, pmass) 
    
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
            expanded_system_jacobian!(J, u, expanded_indices, force_scaling, structural_damping, 
                assembly, pcond, dload, pmass, gvec, x0, v0, ω0, a0, α0)

            return J
        end

    else

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
            steady_state_system_residual!(resid, u, dynamic_indices, force_scaling, structural_damping,
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
            system_mass_matrix!(M, u, dynamic_indices, force_scaling, assembly, pcond, pmass)
    
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
            steady_state_system_jacobian!(J, u, dynamic_indices, force_scaling, structural_damping, 
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
    DAEProblem(system::GXBeam.AbstractSystem, assembly, tspan; kwargs...)

Construct a `DAEProblem` for the system of nonlinear beams contained in `assembly` which 
may be used with the DifferentialEquations package.
"""
function SciMLBase.DAEProblem(system::AbstractSystem, assembly, tspan; 
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
    for (ipoint, icol) in enumerate(system.dynamic_indices.icol_point)
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
    DAEFunction(system::GXBeam.AbstractSystem, assembly; kwargs...)

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
function SciMLBase.DAEFunction(system::AbstractSystem, assembly; 
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
        system_mass_matrix!(J, gamma, u, dynamic_indices, force_scaling, 
            assembly, pcond, pmass)

        return J
    end

    return SciMLBase.DAEFunction{true,true}(f) # TODO: re-add jacobian here once supported
end

function get_differential_vars(system::AbstractSystem, assembly, prescribed_conditions, point_masses)

    # NOTE: If a point and the elements connected to it are massless, then Vdot and Ωdot for
    # the point do not appear in the system of equations and thus not differentiable variables.
    
    differential_vars = zeros(Bool, length(system.x))

    pcond = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(0.0)
    pmass = typeof(point_masses) <: AbstractDict ? point_masses : pmass(0.0)

    for ipoint = 1:length(assembly.points)

        icol = system.dynamic_indices.icol_point[ipoint]

        # check if Vdot and Ωdot are used
        hasmass = haskey(pmass, ipoint)
        for ielem = 1:length(assembly.elements)
            if ipoint == assembly.start[ielem] || ipoint == assembly.stop[ielem]
                hasmass = hasmass || (!iszero(assembly.elements[ielem].L) && !iszero(assembly.elements[ielem].mass))
            end
        end

        if haskey(pcond, ipoint)
            # displacements are differential variables, forces and moments are not
            differential_vars[icol:icol+5] .= pcond[ipoint].isforce
        else
            # displacements and velocities are differentiable variables
            differential_vars[icol:icol+11] .= true
        end

        # velocities are differentiable variables if the point or adjacent elements have mass
        differential_vars[icol+6:icol+11] .= hasmass

    end

    return differential_vars
end
