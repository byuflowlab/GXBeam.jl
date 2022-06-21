"""
    ODEProblem(system::GXBeam.AbstractSystem, assembly, tspan, p=(;); kwargs...)

Construct an `ODEProblem` for the system of nonlinear beams contained in `assembly` which 
may be used with the DifferentialEquations package.

# Arguments
 - `system`:  Object of type `GXBeam.AbstractSystem` which holds indices for accessing the 
    state variables and equations associated with each point and beam element in a system.
 - `assembly`: Object of type `GXBeam.Assembly` which defines an assembly of connected 
    nonlinear beam elements.
 - `tspan`: Time span over which to solve the ODE problem
 - `p`: Parameters, as defined in conjunction with the keyword argument `pfunc`.  
    Defaults to an empty named tuple.   

# Keyword Arguments
 - `pfunc = (p, t) -> p`: Function which returns a named tuple with parameters as 
    described in [`ODEFunction`](@ref).
 - `structural_damping = true`: Flag indicating whether structural damping should be enabled
 - `constant_mass_matrix = true`: Flag indicating whether to use a constant mass matrix.  
"""
function SciMLBase.ODEProblem(system::AbstractSystem, assembly, tspan, p=(;); 
    pfunc = (p, t) -> p,
    structural_damping = true,
    constant_mass_matrix = typeof(system) <: ExpandedSystem)

    # extract parameters from the parameter vector using `pfunc`
    parameters = pfunc(p, tspan[1])

    # unpack parameters
    pcond = get(parameters, :prescribed_conditions, Dict{Int,PrescribedConditions{Float64}}())
    pmass = get(parameters, :point_masses, Dict{Int,PointMass{Float64}}())

    # get parameters for the initial time
    pcond = typeof(pcond) <: AbstractDict ? pcond : pcond(tspan[1])
    pmass = typeof(pmass) <: AbstractDict ? pmass : pmass(tspan[1])

    # set initial state variables
    u0 = copy(system.x)

    # construct ODEFunction
    func = SciMLBase.ODEFunction(system, assembly, pfunc; 
        structural_damping = structural_damping,
        constant_mass_matrix = constant_mass_matrix,
        prescribed_conditions = pcond,
        point_masses = pmass)

    return SciMLBase.ODEProblem{true}(func, u0, tspan, p)
end

"""
    ODEFunction(system::GXBeam.AbstractSystem, assembly, pfunc = (p, t) -> p; kwargs...)

Construct a `ODEFunction` for the system of nonlinear beams contained in `assembly` 
which may be used with the DifferentialEquations package.

# Arguments
 - `system`:  Object of type `GXBeam.AbstractSystem` which holds indices for accessing the 
    state variables and equations associated with each point and beam element in a system.
 - `assembly`: Object of type `GXBeam.Assembly` which defines an assembly of connected 
    nonlinear beam elements.
 - `pfunc = (p, t) -> p`: Function which returns a named tuple with the fields
   - `prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}()`:
      A dictionary with keys corresponding to the points at which prescribed conditions 
      are applied and elements of type [`PrescribedConditions`](@ref) which describe the 
      prescribed conditions at those points.
   - `distributed_loads = Dict{Int,DistributedLoads{Float64}}()`: A dictionary
      with keys corresponding to the elements to which distributed loads are
      applied and elements of type [`DistributedLoads`](@ref) which describe
      the distributed loads on those elements.
   - `point_masses = Dict{Int,Vector{PointMass{Float64}}}()`: A dictionary with keys 
      corresponding to the points at which point masses are attached and values 
      containing vectors of objects of type [`PointMass`](@ref) which describe 
      the point masses attached at those points.
   - `gravity`: Gravity vector.
   - `origin = zeros(3)`: Global frame origin vector.
   - `linear_velocity = zeros(3)`: Global frame linear velocity vector.
   - `angular_velocity = zeros(3)`: Global frame angular velocity vector.
   - `linear_acceleration = zeros(3)`: Global frame linear acceleration vector.
   - `angular_acceleration = zeros(3)`: Global frame angular acceleration vector.

# Keyword Arguments
 - `structural_damping = true`: Flag indicating whether structural damping should be enabled
 - `constant_mass_matrix = true`: Flag indicating whether to use a constant mass matrix.  
 - `prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}()`: Initial prescribed 
    conditions (only used for constant mass matrix systems).  Note that the type of 
    prescribed condition for each point must remain constant for a constant mass matrix 
    system, though the magnitude of prescribed displacements and/or loads may change. 
 - `point_masses = Dict{Int,PointMass{Float64}}()`: Point masses (only used for constant 
    mass matrix systems).  Point mass properties cannot be changed when using a constant 
    mass matrix system.     
"""
function SciMLBase.ODEFunction(system::AbstractSystem, assembly, pfunc = (p, t) -> p; 
    structural_damping = true,
    constant_mass_matrix = typeof(system) <: ExpandedSystem,
    prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}(),
    point_masses = Dict{Int,PointMass{Float64}}())

    force_scaling = system.force_scaling

    if constant_mass_matrix

        indices = SystemIndices(assembly.start, assembly.stop, static=false, expanded=true)

        # residual function (constant mass matrix system)
        f = function(resid, u, p, t)

            # extract parameters from the parameter vector using `pfunc`
            parameters = pfunc(p, t)

            # unpack parameters
            pcond = get(parameters, :prescribed_conditions, Dict{Int,PrescribedConditions{Float64}}())
            dload = get(parameters, :distributed_loads, Dict{Int,DistributedLoads{Float64}}())
            pmass = get(parameters, :point_masses, Dict{Int,PointMass{Float64}}())
            gvec = get(parameters, :gravity, (@SVector zeros(3)))
            x0 = get(parameters, :origin, (@SVector zeros(3)))
            v0 = get(parameters, :linear_velocity, (@SVector zeros(3)))
            ω0 = get(parameters, :angular_velocity, (@SVector zeros(3)))
            a0 = get(parameters, :linear_acceleration, (@SVector zeros(3)))
            α0 = get(parameters, :angular_acceleration, (@SVector zeros(3)))

            # get parameters for this time step
            pcond = typeof(pcond) <: AbstractDict ? pcond : pcond(t)
            dload = typeof(dload) <: AbstractDict ? dload : dload(t)
            pmass = typeof(pmass) <: AbstractDict ? pmass : pmass(t)
            gvec = typeof(gvec) <: AbstractVector ? SVector{3}(gvec) : SVector{3}(gvec(t))
            x0 = typeof(x0) <: AbstractVector ? SVector{3}(x0) : SVector{3}(x0(t))
            v0 = typeof(v0) <: AbstractVector ? SVector{3}(v0) : SVector{3}(v0(t))
            ω0 = typeof(ω0) <: AbstractVector ? SVector{3}(ω0) : SVector{3}(ω0(t))
            a0 = typeof(a0) <: AbstractVector ? SVector{3}(a0) : SVector{3}(a0(t))
            α0 = typeof(α0) <: AbstractVector ? SVector{3}(α0) : SVector{3}(α0(t))

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

            # extract parameters from the parameter vector using `pfunc`
            parameters = pfunc(p, t)

            # unpack parameters
            pcond = get(parameters, :prescribed_conditions, Dict{Int,PrescribedConditions{Float64}}())
            dload = get(parameters, :distributed_loads, Dict{Int,DistributedLoads{Float64}}())
            pmass = get(parameters, :point_masses, Dict{Int,PointMass{Float64}}())
            gvec = get(parameters, :gravity, (@SVector zeros(3)))
            x0 = get(parameters, :origin, (@SVector zeros(3)))
            v0 = get(parameters, :linear_velocity, (@SVector zeros(3)))
            ω0 = get(parameters, :angular_velocity, (@SVector zeros(3)))
            a0 = get(parameters, :linear_acceleration, (@SVector zeros(3)))
            α0 = get(parameters, :angular_acceleration, (@SVector zeros(3)))

            # get parameters for this time step
            pcond = typeof(pcond) <: AbstractDict ? pcond : pcond(t)
            dload = typeof(dload) <: AbstractDict ? dload : dload(t)
            pmass = typeof(pmass) <: AbstractDict ? pmass : pmass(t)
            gvec = typeof(gvec) <: AbstractVector ? SVector{3}(gvec) : SVector{3}(gvec(t))
            x0 = typeof(x0) <: AbstractVector ? SVector{3}(x0) : SVector{3}(x0(t))
            v0 = typeof(v0) <: AbstractVector ? SVector{3}(v0) : SVector{3}(v0(t))
            ω0 = typeof(ω0) <: AbstractVector ? SVector{3}(ω0) : SVector{3}(ω0(t))
            a0 = typeof(a0) <: AbstractVector ? SVector{3}(a0) : SVector{3}(a0(t))
            α0 = typeof(α0) <: AbstractVector ? SVector{3}(α0) : SVector{3}(α0(t))

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

            # extract parameters from the parameter vector using `pfunc`
            parameters = pfunc(p, t)

            # unpack parameters
            pcond = get(parameters, :prescribed_conditions, Dict{Int,PrescribedConditions{Float64}}())
            dload = get(parameters, :distributed_loads, Dict{Int,DistributedLoads{Float64}}())
            pmass = get(parameters, :point_masses, Dict{Int,PointMass{Float64}}())
            gvec = get(parameters, :gravity, (@SVector zeros(3)))
            x0 = get(parameters, :origin, (@SVector zeros(3)))
            v0 = get(parameters, :linear_velocity, (@SVector zeros(3)))
            ω0 = get(parameters, :angular_velocity, (@SVector zeros(3)))
            a0 = get(parameters, :linear_acceleration, (@SVector zeros(3)))
            α0 = get(parameters, :angular_acceleration, (@SVector zeros(3)))

            # get parameters for this time step
            pcond = typeof(pcond) <: AbstractDict ? pcond : pcond(t)
            dload = typeof(dload) <: AbstractDict ? dload : dload(t)
            pmass = typeof(pmass) <: AbstractDict ? pmass : pmass(t)
            gvec = typeof(gvec) <: AbstractVector ? SVector{3}(gvec) : SVector{3}(gvec(t))
            x0 = typeof(x0) <: AbstractVector ? SVector{3}(x0) : SVector{3}(x0(t))
            v0 = typeof(v0) <: AbstractVector ? SVector{3}(v0) : SVector{3}(v0(t))
            ω0 = typeof(ω0) <: AbstractVector ? SVector{3}(ω0) : SVector{3}(ω0(t))
            a0 = typeof(a0) <: AbstractVector ? SVector{3}(a0) : SVector{3}(a0(t))
            α0 = typeof(α0) <: AbstractVector ? SVector{3}(α0) : SVector{3}(α0(t))

            # calculate residual
            steady_state_system_residual!(resid, u, indices, force_scaling, structural_damping,
                assembly, pcond, dload, pmass, gvec, x0, v0, ω0, a0, α0)

            return resid
        end

        # mass matrix
        update_mass_matrix! = function(M, u, p, t)

            # extract parameters from the parameter vector using `pfunc`
            parameters = pfunc(p, t)

            # unpack parameters
            pcond = get(parameters, :prescribed_conditions, Dict{Int,PrescribedConditions{Float64}}())
            pmass = get(parameters, :point_masses, Dict{Int,PointMass{Float64}}())

            # get parameters for this time step
            pcond = typeof(pcond) <: AbstractDict ? pcond : pcond(t)
            pmass = typeof(pmass) <: AbstractDict ? pmass : pmass(t)

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

            # extract parameters from the parameter vector using `pfunc`
            parameters = pfunc(p, t)

            # unpack parameters
            pcond = get(parameters, :prescribed_conditions, Dict{Int,PrescribedConditions{Float64}}())
            dload = get(parameters, :distributed_loads, Dict{Int,DistributedLoads{Float64}}())
            pmass = get(parameters, :point_masses, Dict{Int,PointMass{Float64}}())
            gvec = get(parameters, :gravity, (@SVector zeros(3)))
            x0 = get(parameters, :origin, (@SVector zeros(3)))
            v0 = get(parameters, :linear_velocity, (@SVector zeros(3)))
            ω0 = get(parameters, :angular_velocity, (@SVector zeros(3)))
            a0 = get(parameters, :linear_acceleration, (@SVector zeros(3)))
            α0 = get(parameters, :angular_acceleration, (@SVector zeros(3)))

            # get parameters for this time step
            pcond = typeof(pcond) <: AbstractDict ? pcond : pcond(t)
            dload = typeof(dload) <: AbstractDict ? dload : dload(t)
            pmass = typeof(pmass) <: AbstractDict ? pmass : pmass(t)
            gvec = typeof(gvec) <: AbstractVector ? SVector{3}(gvec) : SVector{3}(gvec(t))
            x0 = typeof(x0) <: AbstractVector ? SVector{3}(x0) : SVector{3}(x0(t))
            v0 = typeof(v0) <: AbstractVector ? SVector{3}(v0) : SVector{3}(v0(t))
            ω0 = typeof(ω0) <: AbstractVector ? SVector{3}(ω0) : SVector{3}(ω0(t))
            a0 = typeof(a0) <: AbstractVector ? SVector{3}(a0) : SVector{3}(a0(t))
            α0 = typeof(α0) <: AbstractVector ? SVector{3}(α0) : SVector{3}(α0(t))

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
    DAEProblem(system::GXBeam.DynamicSystem, assembly, tspan, p=(;); kwargs...)

Construct a `DAEProblem` for the system of nonlinear beams contained in `assembly` which 
may be used with the DifferentialEquations package.

# Arguments
 - `system`:  Object of type `GXBeam.AbstractSystem` which holds indices for accessing the 
    state variables and equations associated with each point and beam element in a system.
 - `assembly`: Object of type `GXBeam.Assembly` which defines an assembly of connected 
    nonlinear beam elements.
 - `tspan`: Time span over which to solve the ODE problem    
 - `p`: Parameters, as defined in conjunction with the keyword argument `pfunc`.  
    Defaults to an empty named tuple.

# Keyword Arguments
 - `pfunc = (p, t) -> p`: Function which returns a named tuple with parameters as 
    described in [`DAEFunction`](@ref).
 - `structural_damping = true`: Flag indicating whether structural damping should be enabled
"""

function SciMLBase.DAEProblem(system::AbstractSystem, assembly, tspan, p=(;); 
    pfunc = (p, t) -> p,
    structural_damping = true)

    # extract parameters from the parameter vector using `pfunc`
    parameters = pfunc(p, tspan[1])

    # unpack parameters
    pcond = get(parameters, :prescribed_conditions, Dict{Int,PrescribedConditions{Float64}}())
    pmass = get(parameters, :point_masses, Dict{Int,PointMass{Float64}}())

    # get parameters for the initial time
    pcond = typeof(pcond) <: AbstractDict ? pcond : pcond(tspan[1])
    pmass = typeof(pmass) <: AbstractDict ? pmass : pmass(tspan[1])

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
    differential_vars = get_differential_vars(system, assembly, pcond, pmass)

    # create SciMLBase.DAEFunction
    func = SciMLBase.DAEFunction(system, assembly, pfunc; structural_damping)

    return SciMLBase.DAEProblem{true}(func, du0, u0, tspan, p; differential_vars)
end

"""
    DAEFunction(system::GXBeam.DynamicSystem, assembly; kwargs...)

Construct a `DAEFunction` for the system of nonlinear beams contained in `assembly` 
which may be used with the DifferentialEquations package.

# Arguments
 - `system`:  Object of type `GXBeam.AbstractSystem` which holds indices for accessing the 
    state variables and equations associated with each point and beam element in a system.
 - `assembly`: Object of type `GXBeam.Assembly` which defines an assembly of connected 
    nonlinear beam elements.
 - `pfunc = (p, t) -> p`: Function which returns a named tuple with the fields
   - `prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}()`:
      A dictionary with keys corresponding to the points at which prescribed conditions 
      are applied and elements of type [`PrescribedConditions`](@ref) which describe the 
      prescribed conditions at those points.
   - `distributed_loads = Dict{Int,DistributedLoads{Float64}}()`: A dictionary
      with keys corresponding to the elements to which distributed loads are
      applied and elements of type [`DistributedLoads`](@ref) which describe
      the distributed loads on those elements.
   - `point_masses = Dict{Int,Vector{PointMass{Float64}}}()`: A dictionary with keys 
      corresponding to the points at which point masses are attached and values 
      containing vectors of objects of type [`PointMass`](@ref) which describe 
      the point masses attached at those points.
   - `gravity`: Gravity vector.
   - `origin = zeros(3)`: Global frame origin vector.
   - `linear_velocity = zeros(3)`: Global frame linear velocity vector.
   - `angular_velocity = zeros(3)`: Global frame angular velocity vector.
   - `linear_acceleration = zeros(3)`: Global frame linear acceleration vector.
   - `angular_acceleration = zeros(3)`: Global frame angular acceleration vector.

Keyword Arguments:
 - `structural_damping = true`: Flag indicating whether structural damping should be enabled
"""
function SciMLBase.DAEFunction(system::DynamicSystem, assembly, pfunc = (p, t) -> p; 
    structural_damping = true)

    # unpack system pointers
    @unpack force_scaling, indices = system

    # DAE function
    f = function(resid, du, u, p, t)

        # extract parameters from the parameter vector using `pfunc`
        parameters = pfunc(p, t)

        # unpack parameters
        pcond = get(parameters, :prescribed_conditions, Dict{Int,PrescribedConditions{Float64}}())
        dload = get(parameters, :distributed_loads, Dict{Int,DistributedLoads{Float64}}())
        pmass = get(parameters, :point_masses, Dict{Int,PointMass{Float64}}())
        gvec = get(parameters, :gravity, (@SVector zeros(3)))
        x0 = get(parameters, :origin, (@SVector zeros(3)))
        v0 = get(parameters, :linear_velocity, (@SVector zeros(3)))
        ω0 = get(parameters, :angular_velocity, (@SVector zeros(3)))
        a0 = get(parameters, :linear_acceleration, (@SVector zeros(3)))
        α0 = get(parameters, :angular_acceleration, (@SVector zeros(3)))

        # get parameters for this time step
        pcond = typeof(pcond) <: AbstractDict ? pcond : pcond(t)
        dload = typeof(dload) <: AbstractDict ? dload : dload(t)
        pmass = typeof(pmass) <: AbstractDict ? pmass : pmass(t)
        gvec = typeof(gvec) <: AbstractVector ? SVector{3}(gvec) : SVector{3}(gvec(t))
        x0 = typeof(x0) <: AbstractVector ? SVector{3}(x0) : SVector{3}(x0(t))
        v0 = typeof(v0) <: AbstractVector ? SVector{3}(v0) : SVector{3}(v0(t))
        ω0 = typeof(ω0) <: AbstractVector ? SVector{3}(ω0) : SVector{3}(ω0(t))
        a0 = typeof(a0) <: AbstractVector ? SVector{3}(a0) : SVector{3}(a0(t))
        α0 = typeof(α0) <: AbstractVector ? SVector{3}(α0) : SVector{3}(α0(t))

        # calculate residual
        dynamic_system_residual!(resid, du, u, indices, force_scaling, structural_damping, 
            assembly, pcond, dload, pmass, gvec, x0, v0, ω0, a0, α0)

        return resid
    end

    # jacobian function with respect to states/state rates
    update_jacobian! = function(J, du, u, p, gamma, t)

        # zero out all jacobian entries
        J .= 0.0

        # extract parameters from the parameter vector using `pfunc`
        parameters = pfunc(p, t)

        # unpack parameters
        pcond = get(parameters, :prescribed_conditions, Dict{Int,PrescribedConditions{Float64}}())
        dload = get(parameters, :distributed_loads, Dict{Int,DistributedLoads{Float64}}())
        pmass = get(parameters, :point_masses, Dict{Int,PointMass{Float64}}())
        gvec = get(parameters, :gravity, (@SVector zeros(3)))
        x0 = get(parameters, :origin, (@SVector zeros(3)))
        v0 = get(parameters, :linear_velocity, (@SVector zeros(3)))
        ω0 = get(parameters, :angular_velocity, (@SVector zeros(3)))
        a0 = get(parameters, :linear_acceleration, (@SVector zeros(3)))
        α0 = get(parameters, :angular_acceleration, (@SVector zeros(3)))

        # get parameters for this time step
        pcond = typeof(pcond) <: AbstractDict ? pcond : pcond(t)
        dload = typeof(dload) <: AbstractDict ? dload : dload(t)
        pmass = typeof(pmass) <: AbstractDict ? pmass : pmass(t)
        gvec = typeof(gvec) <: AbstractVector ? SVector{3}(gvec) : SVector{3}(gvec(t))
        x0 = typeof(x0) <: AbstractVector ? SVector{3}(x0) : SVector{3}(x0(t))
        v0 = typeof(v0) <: AbstractVector ? SVector{3}(v0) : SVector{3}(v0(t))
        ω0 = typeof(ω0) <: AbstractVector ? SVector{3}(ω0) : SVector{3}(ω0(t))
        a0 = typeof(a0) <: AbstractVector ? SVector{3}(a0) : SVector{3}(a0(t))
        α0 = typeof(α0) <: AbstractVector ? SVector{3}(α0) : SVector{3}(α0(t))

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
