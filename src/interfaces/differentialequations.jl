"""
    ODEProblem(system::GXBeam.System, assembly, tspan; kwargs...)

Construct a `ODEProblem` for the system of nonlinear beams
contained in `assembly` which may be used with the DifferentialEquations package.

**Note that this function defines a non-constant mass matrix, which is not directly 
supported by DifferentialEquations.jl.  This function is therefore experimental**

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
function SciMLBase.ODEProblem(system::System, assembly, tspan;
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

    # create ODEFunction
    func = SciMLBase.ODEFunction(system, assembly; structural_damping)

    # use initial state from `system`
    u0 = copy(system.x)

    # set parameters
    p = (prescribed_conditions, distributed_loads, point_masses, gravity, origin, 
        linear_velocity, angular_velocity, linear_acceleration, angular_acceleration)

    return SciMLBase.ODEProblem{true}(func, u0, tspan, p)
end

"""
    ODEFunction(system::GXBeam.System, assembly; structural_damping=true)

Construct a `ODEFunction` for the system of nonlinear beams
contained in `assembly` which may be used with the DifferentialEquations package.

**Note that this function defines a non-constant mass matrix, which is not directly 
supported by DifferentialEquations.jl.  This function is therefore experimental**

The parameters associated with the resulting ODEFunction are defined by the tuple
`(prescribed_conditions, distributed_loads, point_masses, origin, linear_velocity, 
angular_velocity, linear_acceleration, angular_acceleration)` where each parameter is defined as follows:
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
 - `point_masses = Dict{Int,Vector{PointMass{Float64}}}()`: A dictionary with keys 
        corresponding to the points at which point masses are attached and values 
        containing vectors of objects of type [`PointMass`](@ref) which describe 
        the point masses attached at those points.  If time varying, this input may
        be provided as a function of time.
 - `gravity`: Gravity vector. If time varying, this input may be provided as a 
        function of time.
 - `origin`: Global frame origin vector. If time varying, this input
        may be provided as a function of time.
 - `linear_velocity`: Global frame linear velocity vector. If time
        varying, this vector may be provided as a function of time.
 - `angular_velocity`: Global frame angular velocity vector. If time
        varying, this vector may be provided as a function of time.
 - `linear_acceleration = zeros(3)`: Global frame linear acceleration vector. If time
       varying, this vector may be provided as a function of time.
 - `angular_acceleration = zeros(3)`: Global frame angular acceleration vector. If time
       varying, this vector may be provided as a function of time.
"""
function SciMLBase.ODEFunction(system::System, assembly; structural_damping=false)

    # check to make sure the system isn't static
    @assert !system.static

    # unpack system pointers
    irow_point = system.irow_point
    irow_elem = system.irow_elem
    irow_elem1 = system.irow_elem1
    irow_elem2 = system.irow_elem2
    icol_point = system.icol_point
    icol_elem = system.icol_elem

    # unpack scaling parameters
    force_scaling = system.force_scaling

    # DAE function
    f = function(resid, u, p, t)

        # get current parameters
        prescribed_conditions = typeof(p[1]) <: AbstractDict ? p[1] : p[1](t)
        distributed_loads = typeof(p[2]) <: AbstractDict ? p[2] : p[2](t)
        point_masses = typeof(p[3]) <: AbstractDict ? p[3] : p[3](t)
        gvec = typeof(p[4]) <: AbstractVector ? SVector{3}(p[4]) : SVector{3}(p[4](t))
        x0 = typeof(p[5]) <: AbstractVector ? SVector{3}(p[5]) : SVector{3}(p[5](t))
        v0 = typeof(p[6]) <: AbstractVector ? SVector{3}(p[6]) : SVector{3}(p[6](t))
        ω0 = typeof(p[7]) <: AbstractVector ? SVector{3}(p[7]) : SVector{3}(p[7](t))
        a0 = typeof(p[8]) <: AbstractVector ? SVector{3}(p[8]) : SVector{3}(p[8](t))
        α0 = typeof(p[9]) <: AbstractVector ? SVector{3}(p[9]) : SVector{3}(p[9](t))

        # calculate residual
        steady_state_system_residual!(resid, u, assembly, prescribed_conditions,
            distributed_loads, point_masses, gvec, force_scaling, irow_point, irow_elem, irow_elem1,
            irow_elem2, icol_point, icol_elem, x0, v0, ω0, a0, α0)

        return resid
    end

    update_mass_matrix! = function(M, u, p, t)

       point_masses = typeof(p[3]) <: AbstractDict ? p[3] : p[3](t)

        # zero out all mass matrix entries
        M .= 0.0

        # calculate mass matrix
        system_mass_matrix!(M, u, assembly, point_masses, structural_damping, force_scaling,
            irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem)

        return M
    end

    mass_matrix = SciMLBase.DiffEqArrayOperator(copy(system.M), update_func = update_mass_matrix!)

    # jacobian function with respect to states/state rates
    jac = function(J, u, p, t)

        # zero out all jacobian entries
        J .= 0.0

        # get current parameters
        prescribed_conditions = typeof(p[1]) <: AbstractDict ? p[1] : p[1](t)
        distributed_loads = typeof(p[2]) <: AbstractDict ? p[2] : p[2](t)
        point_masses = typeof(p[3]) <: AbstractDict ? p[3] : p[3](t)
        gvec = typeof(p[4]) <: AbstractVector ? SVector{3}(p[4]) : SVector{3}(p[4](t))
        x0 = typeof(p[5]) <: AbstractVector ? SVector{3}(p[5]) : SVector{3}(p[5](t))
        v0 = typeof(p[6]) <: AbstractVector ? SVector{3}(p[6]) : SVector{3}(p[6](t))
        ω0 = typeof(p[7]) <: AbstractVector ? SVector{3}(p[7]) : SVector{3}(p[7](t))
        a0 = typeof(p[8]) <: AbstractVector ? SVector{3}(p[8]) : SVector{3}(p[8](t))
        α0 = typeof(p[9]) <: AbstractVector ? SVector{3}(p[9]) : SVector{3}(p[9](t))

        # calculate jacobian
        steady_state_system_jacobian!(J, u, assembly, prescribed_conditions,
            distributed_loads, point_masses, gvec, force_scaling, irow_point, irow_elem, irow_elem1,
            irow_elem2, icol_point, icol_elem, x0, v0, ω0, a0, α0)

        return J
    end

    # sparsity structure
    sparsity = get_sparsity(system, assembly)

    # jacobian prototype (use dense since sparse isn't working)
    jac_prototype = collect(system.K)

    # TODO: figure out how to use a sparse matrix here.
    # It's failing with a singular exception during the LU factorization.

    return SciMLBase.ODEFunction{true,true}(f; mass_matrix, jac)
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
        the distributed loads on those elements.  If time varying, this input may
        be provided as a function of time.
 - `point_masses = Dict{Int,Vector{PointMass{Float64}}}()`: A dictionary with keys 
        corresponding to the points at which point masses are attached and values 
        containing vectors of objects of type [`PointMass`](@ref) which describe 
        the point masses attached at those points.  If time varying, this input may
        be provided as a function of time.
 - `structural_damping = true`: Flag indicating whether structural damping should be enabled
 - `gravity = zeros(3)`: Gravity vector. If time varying, this input may be provided as a 
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
function SciMLBase.DAEProblem(system::System, assembly, tspan;
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

    # create SciMLBase.DAEFunction
    func = SciMLBase.DAEFunction(system, assembly; structural_damping)

    # use initial state from `system`
    u0 = copy(system.x)

    # use initial state rates from `system`
    du0 = zero(u0)
    for (ielem, icol) in enumerate(system.icol_elem)
        du0[icol:icol+2] = system.udot[ielem]
        du0[icol+3:icol+5] = system.θdot[ielem]
        du0[icol+12:icol+14] = system.Vdot[ielem]
        du0[icol+15:icol+17] = system.Ωdot[ielem]
    end

    # set parameters
    p = (prescribed_conditions, distributed_loads, point_masses, gravity, origin, 
       linear_velocity, angular_velocity, linear_acceleration, angular_acceleration)

    # get differential variables
    differential_vars = get_differential_vars(system, assembly, structural_damping)

    return SciMLBase.DAEProblem{true}(func, du0, u0, tspan, p; differential_vars)
end

"""
    DAEFunction(system::GXBeam.System, assembly; structural_damping=true)

Construct a `DAEFunction` for the system of nonlinear beams
contained in `assembly` which may be used with the DifferentialEquations package.

The parameters associated with the resulting SciMLBase.DAEFunction are defined by the tuple
`(prescribed_conditions, distributed_loads, point_masses, origin, linear_velocity, 
angular_velocity, linear_acceleration, angular_acceleration)`
where each parameter is defined as follows:
 - `prescribed_conditions`: A dictionary with keys corresponding to the points at
        which prescribed conditions are applied and elements of type
        [`PrescribedConditions`](@ref) which describe the prescribed conditions
        at those points.  If time varying, this input may be provided as a
        function of time.
 - `distributed_loads`: A dictionary with keys corresponding to the elements to
        which distributed loads are applied and elements of type [`DistributedLoads`](@ref) 
        which describe the distributed loads on those elements.  If time varying, this 
        input may be provided as a function of time.
 - `point_masses = Dict{Int,Vector{PointMass{Float64}}}()`: A dictionary with keys 
        corresponding to the points at which point masses are attached and values 
        containing vectors of objects of type [`PointMass`](@ref) which describe 
        the point masses attached at those points.  If time varying, this input may
        be provided as a function of time.
 - `gravity`: Gravity vector. If time varying, this input may be provided as a 
        function of time.
 - `origin`: Global frame origin vector. If time varying, this input
        may be provided as a function of time.
 - `linear_velocity`: Global frame linear velocity vector. If time
        varying, this vector may be provided as a function of time.
 - `angular_velocity`: Global frame angular velocity vector. If time
        varying, this vector may be provided as a function of time.
 - `linear_acceleration = zeros(3)`: Global frame linear acceleration vector. If time
       varying, this vector may be provided as a function of time.
 - `angular_acceleration = zeros(3)`: Global frame angular acceleration vector. If time
       varying, this vector may be provided as a function of time.
"""
function SciMLBase.DAEFunction(system::System, assembly; structural_damping=true)

    # check to make sure the system isn't static
    @assert !system.static

    # unpack system pointers
    irow_point = system.irow_point
    irow_elem = system.irow_elem
    irow_elem1 = system.irow_elem1
    irow_elem2 = system.irow_elem2
    icol_point = system.icol_point
    icol_elem = system.icol_elem

    # unpack scaling parameters
    force_scaling = system.force_scaling

    # DAE function
    f = function(resid, du, u, p, t)

        # get current parameters
        prescribed_conditions = typeof(p[1]) <: AbstractDict ? p[1] : p[1](t)
        distributed_loads = typeof(p[2]) <: AbstractDict ? p[2] : p[2](t)
        point_masses = typeof(p[3]) <: AbstractDict ? p[3] : p[3](t)
        gvec = typeof(p[4]) <: AbstractVector ? SVector{3}(p[4]) : SVector{3}(p[4](t))
        x0 = typeof(p[5]) <: AbstractVector ? SVector{3}(p[5]) : SVector{3}(p[5](t))
        v0 = typeof(p[6]) <: AbstractVector ? SVector{3}(p[6]) : SVector{3}(p[6](t))
        ω0 = typeof(p[7]) <: AbstractVector ? SVector{3}(p[7]) : SVector{3}(p[7](t))
        a0 = typeof(p[8]) <: AbstractVector ? SVector{3}(p[8]) : SVector{3}(p[8](t))
        α0 = typeof(p[9]) <: AbstractVector ? SVector{3}(p[9]) : SVector{3}(p[9](t))

        # calculate residual
        dynamic_system_residual!(resid, du, u, assembly, prescribed_conditions,
            distributed_loads, point_masses, structural_damping, gvec, force_scaling, 
            irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem, 
            x0, v0, ω0, a0, α0)

        return resid
    end

    # jacobian function with respect to states/state rates
    jac = function(J, du, u, p, gamma, t)

        # zero out all jacobian entries
        J .= 0.0

        # get current parameters
        prescribed_conditions = typeof(p[1]) <: AbstractDict ? p[1] : p[1](t)
        distributed_loads = typeof(p[2]) <: AbstractDict ? p[2] : p[2](t)
        point_masses = typeof(p[3]) <: AbstractDict ? p[3] : p[3](t)
        gvec = typeof(p[4]) <: AbstractVector ? SVector{3}(p[4]) : SVector{3}(p[4](t))
        x0 = typeof(p[5]) <: AbstractVector ? SVector{3}(p[5]) : SVector{3}(p[5](t))
        v0 = typeof(p[6]) <: AbstractVector ? SVector{3}(p[6]) : SVector{3}(p[6](t))
        ω0 = typeof(p[7]) <: AbstractVector ? SVector{3}(p[7]) : SVector{3}(p[7](t))
        a0 = typeof(p[8]) <: AbstractVector ? SVector{3}(p[8]) : SVector{3}(p[8](t))
        α0 = typeof(p[9]) <: AbstractVector ? SVector{3}(p[9]) : SVector{3}(p[9](t))

        # calculate jacobian
        dynamic_system_jacobian!(J, du, u, assembly, prescribed_conditions,
            distributed_loads, point_masses, structural_damping, gvec, force_scaling, 
            irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem, 
            x0, v0, ω0, a0, α0)

        # add gamma multiplied by the mass matrix
        system_mass_matrix!(J, gamma, u, assembly, point_masses, structural_damping, 
        force_scaling, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem)

        return J
    end

    # sparsity structure
    sparsity = get_sparsity(system, assembly)

    # jacobian prototype (use dense since sparse isn't working)
    # jac_prototype = collect(system.K)

    # TODO: figure out how to use a sparse matrix here.
    # It's failing with a singular exception during the LU factorization.
    # Using `jac_prototype` also causes errors

    return SciMLBase.DAEFunction{true,true}(f) # TODO: re-add jacobian here once supported
end

function get_differential_vars(system::System, assembly::Assembly, structural_damping)
    differential_vars = fill(false, length(system.x))
    for (ielem, icol) in enumerate(system.icol_elem)
       icol = system.icol_elem[ielem]
       differential_vars[icol:icol+2] .= true # u (for the beam element)
       differential_vars[icol+3:icol+5] .= true # θ (for the beam element)
       for i = 1:6
           if structural_damping && !iszero(assembly.elements[ielem].mu[i])
              differential_vars[icol+5+i] = true
           end
       end
       differential_vars[icol+12:icol+14] .= true # V (for the beam element)
       differential_vars[icol+15:icol+17] .= true # Ω (for the beam element)
    end
    return differential_vars
end
