"""
    ODEProblem(system::GXBeam.AbstractSystem, assembly, tspan; kwargs...)

Construct an `ODEProblem` for the system of nonlinear beams contained in `assembly` which
may be used with the DifferentialEquations package.

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
 - `linear_velocity = zeros(3)`: Prescribed linear velocity of the body frame.
 - `angular_velocity = zeros(3)`: Prescribed angular velocity of the body frame.
 - `gravity = [0,0,0]`: Gravity vector in the body frame.  If time varying, this input
        may be provided as a function of time.

# Control Flag Keyword Arguments
 - `initial_state = nothing`: Object of type `AssemblyState`, which defines the initial
        states and state rates corresponding to the analysis.  By default, this input is
        calculated using either `steady_state_analysis` or `initial_condition_analysis`.
 - `structural_damping = true`: Flag indicating whether to enable structural damping
 - `two_dimensional = false`: Flag indicating whether to constrain results to the x-y plane
 - `constant_mass_matrix = true`: Flag indicating whether to use a constant mass matrix.
 - `sparse = false`: Flag indicating whether to use a sparse jacobian.

 # Sensitivity Analysis Keyword Arguments
 - `pfunc = (p, t) -> (;)`: Function which returns a named tuple with fields corresponding
         to updated versions of the arguments `assembly`, `prescribed_conditions`,
         `distributed_loads`, `point_masses`, `linear_velocity`, `angular_velocity`, and
         `gravity`. Only fields contained in the resulting named tuple will be overwritten.
 - `p`: Sensitivity parameters, as defined in conjunction with the keyword argument `pfunc`.
         While not necessary, using `pfunc` and `p` to define the arguments to this function
         allows automatic differentiation sensitivities to be computed more efficiently

Additional keyword arguments are passed forward to ODEProblem.
"""
function SciMLBase.ODEProblem(system::AbstractSystem, assembly, tspan, p=(;);
    pfunc = (p, t) -> p,
    two_dimensional = false,
    structural_damping = true,
    constant_mass_matrix = typeof(system) <: ExpandedSystem,
    sparse = false,
    kwargs...)

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
    func = SciMLBase.ODEFunction(system, assembly, pfunc, p;
        two_dimensional, structural_damping, constant_mass_matrix, sparse)

    return SciMLBase.ODEProblem{true}(func, u0, tspan, p; kwargs...)
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
   - `linear_velocity = zeros(3)`: Global frame linear velocity vector.
   - `angular_velocity = zeros(3)`: Global frame angular velocity vector.
 - `p`: Parameters, as defined in conjunction with the keyword argument `pfunc`.
    Defaults to an empty named tuple.

# Keyword Arguments
 - `two_dimensional = false`: Flag indicating whether to constrain results to the x-y plane
 - `structural_damping = true`: Flag indicating whether structural damping should be enabled
 - `constant_mass_matrix = true`: Flag indicating whether to use a constant mass matrix.
 - `sparse = false`: Flag indicating whether to use a sparse jacobian.
"""
function SciMLBase.ODEFunction(system::AbstractSystem, assembly, pfunc=(p, t) -> p, p0=(;);
    two_dimensional = false, structural_damping = true, sparse=false,
    constant_mass_matrix = typeof(system) <: ExpandedSystem)

    for ielem = 1:length(assembly.elements)
        @assert !iszero(assembly.elements[ielem].L) "Zero length elements cannot be used "*
            "with DifferentialEquations"
        @assert !iszero(assembly.elements[ielem].compliance[1,:]) &&
            !iszero(assembly.elements[ielem].compliance[2,:]) &&
            !iszero(assembly.elements[ielem].compliance[3,:]) &&
            !iszero(assembly.elements[ielem].compliance[4,:]) &&
            !iszero(assembly.elements[ielem].compliance[5,:]) &&
            !iszero(assembly.elements[ielem].compliance[6,:]) "Compliance matrix must "*
            "be invertible when using ODEFunction."
        @assert !iszero(assembly.elements[ielem].mass[1,:]) &&
            !iszero(assembly.elements[ielem].mass[2,:]) &&
            !iszero(assembly.elements[ielem].mass[3,:]) &&
            !iszero(assembly.elements[ielem].mass[4,:]) &&
            !iszero(assembly.elements[ielem].mass[5,:]) &&
            !iszero(assembly.elements[ielem].mass[6,:]) "Mass matrix must "*
            "be invertible when using ODEFunction."
    end

    # get parameters at time t=0.0
    pcond, dload, pmass, gvec, vb_p, ωb_p = extract_parameters(pfunc, p0, 0.0)

    # extract force scaling
    force_scaling = system.force_scaling

    if constant_mass_matrix

        indices = SystemIndices(assembly.start, assembly.stop, static=false, expanded=true)
        du = zeros(indices.nstates)

        # residual function (constant mass matrix system)
        f = function(resid, u, p, t)

            # extract parameters from the parameter vector using `pfunc`
            pcond, dload, pmass, gvec, vb_p, ωb_p = extract_parameters(pfunc, p, t)

            # calculate residual
            expanded_dynamic_system_residual!(resid, du, u, indices, two_dimensional, force_scaling,
                structural_damping, assembly, pcond, dload, pmass, gvec, vb_p, ωb_p)

            return resid
        end

        # mass matrix (constant mass matrix system)
        TF = eltype(system)
        nx = indices.nstates
        mass_matrix = zeros(TF, nx, nx)
        expanded_system_mass_matrix!(mass_matrix, -1, indices, two_dimensional, force_scaling, assembly, pcond, pmass)

        # jacobian
        update_jacobian! = function(J, u, p, t)

            # zero out all jacobian entries
            J .= 0.0

            # extract parameters from the parameter vector using `pfunc`
            pcond, dload, pmass, gvec, vb_p, ωb_p = extract_parameters(pfunc, p, t)

            # calculate jacobian
            expanded_dynamic_system_jacobian!(J, du, u, indices, two_dimensional, force_scaling,
                structural_damping, assembly, pcond, dload, pmass, gvec, vb_p, ωb_p)

            return J
        end

        # jacobian prototype
        if sparse
            jac_prototype = spzeros(TF, nx, nx)
            expanded_dynamic_system_jacobian!(jac_prototype, rand(nx), rand(nx), indices, two_dimensional, force_scaling,
                structural_damping, assembly, pcond, dload, pmass, rand(3), rand(3), rand(3))
        else
            jac_prototype = nothing
        end

    else

        indices = SystemIndices(assembly.start, assembly.stop, static=false, expanded=false)
        du = zeros(indices.nstates)
        u = rand(indices.nstates)

        # residual function
        f = function(resid, u, p, t)

            # extract parameters from the parameter vector using `pfunc`
            pcond, dload, pmass, gvec, vb_p, ωb_p = extract_parameters(pfunc, p, t)

            # calculate residual
            dynamic_system_residual!(resid, du, u, indices, two_dimensional, force_scaling,
                structural_damping, assembly, pcond, dload, pmass, gvec, vb_p, ωb_p)

            return resid
        end

        # mass matrix
        update_mass_matrix! = function(M, u, p, t)

            # zero out all mass matrix entries
            M .= 0.0

            # extract parameters from the parameter vector using `pfunc`
            pcond, dload, pmass, gvec, vb_p, ωb_p = extract_parameters(pfunc, p, t)

            # calculate mass matrix
            system_mass_matrix!(M, u, indices, two_dimensional, force_scaling, assembly, pcond, pmass)

            # change sign of mass matrix
            M .*= -1

            return M
        end

        mass_matrix = SciMLBase.DiffEqArrayOperator(system.M, update_func = update_mass_matrix!)

        # jacobian
        update_jacobian! = function(J, u, p, t)

            # zero out all jacobian entries
            J .= 0.0

            # extract parameters from the parameter vector using `pfunc`
            pcond, dload, pmass, gvec, vb_p, ωb_p = extract_parameters(pfunc, p, t)

            # calculate jacobian
            dynamic_system_jacobian!(J, du, u, indices, two_dimensional, force_scaling,
                structural_damping, assembly, pcond, dload, pmass, gvec, vb_p, ωb_p)

            return J
        end

        # jacobian prototype
        if sparse
            TF = eltype(system)
            nx = indices.nstates
            jac_prototype = spzeros(TF, nx, nx)
            dynamic_system_jacobian!(jac_prototype, du, u, indices, two_dimensional, force_scaling,
                structural_damping, assembly, pcond, dload, pmass, gvec, vb_p, ωb_p)
        else
            jac_prototype = nothing
        end

    end

    return SciMLBase.ODEFunction{true,true}(f; mass_matrix = mass_matrix,
        jac = update_jacobian!, jac_prototype = jac_prototype)
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
 - `two_dimensional = false`: Flag indicating whether to constrain results to the x-y plane
 - `structural_damping = true`: Flag indicating whether structural damping should be enabled
"""
function SciMLBase.DAEProblem(system::AbstractSystem, assembly, tspan, p=(;);
    pfunc = (p, t) -> p, two_dimensional=false, structural_damping = true)

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
    differential_vars = get_differential_vars(system, assembly, pcond, pmass, two_dimensional)

    # create SciMLBase.DAEFunction
    func = SciMLBase.DAEFunction(system, assembly, pfunc; two_dimensional, structural_damping)

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
   - `linear_velocity = zeros(3)`: Global frame linear velocity vector.
   - `angular_velocity = zeros(3)`: Global frame angular velocity vector.

Keyword Arguments:
 - `two_dimensional = false`: Flag indicating whether to constrain results to the x-y plane
 - `structural_damping = true`: Flag indicating whether structural damping should be enabled
"""
function SciMLBase.DAEFunction(system::DynamicSystem, assembly, pfunc = (p, t) -> p;
    two_dimensional = false, structural_damping = true)

    # unpack system pointers
    @unpack force_scaling, indices = system

    # DAE function
    f = function(resid, du, u, p, t)

        # extract parameters from the parameter vector using `pfunc`
        pcond, dload, pmass, gvec, vb_p, ωb_p = extract_parameters(pfunc, p, t)

        # calculate residual
        dynamic_system_residual!(resid, du, u, indices, two_dimensional, force_scaling,
            structural_damping, assembly, pcond, dload, pmass, gvec, vb_p, ωb_p)

        return resid
    end

    # jacobian function with respect to states/state rates
    update_jacobian! = function(J, du, u, p, gamma, t)

        # zero out all jacobian entries
        J .= 0.0

        # extract parameters from the parameter vector using `pfunc`
        pcond, dload, pmass, gvec, vb_p, ωb_p = extract_parameters(pfunc, p, t)

        # calculate jacobian
        dynamic_system_jacobian!(J, du, u, indices, two_dimensional, force_scaling,
            structural_damping, assembly, pcond, dload, pmass, gvec, vb_p, ωb_p)

        # add gamma multiplied by the mass matrix
        system_mass_matrix!(J, gamma, u, indices, two_dimensional, force_scaling,
            assembly, pcond, pmass)

        return J
    end

    return SciMLBase.DAEFunction{true,true}(f) # TODO: re-add jacobian here once supported
end

function extract_parameters(pfunc, p, t)

    # extract parameters from the parameter vector using `pfunc`
    parameters = pfunc(p, t)

    # unpack parameters
    pcond = get(parameters, :prescribed_conditions, Dict{Int,PrescribedConditions{Float64}}())
    dload = get(parameters, :distributed_loads, Dict{Int,DistributedLoads{Float64}}())
    pmass = get(parameters, :point_masses, Dict{Int,PointMass{Float64}}())
    gvec = get(parameters, :gravity, (@SVector zeros(3)))
    vb_p = get(parameters, :linear_velocity, (@SVector zeros(3)))
    ωb_p = get(parameters, :angular_velocity, (@SVector zeros(3)))

    # get parameters for this time step
    pcond = typeof(pcond) <: AbstractDict ? pcond : pcond(t)
    dload = typeof(dload) <: AbstractDict ? dload : dload(t)
    pmass = typeof(pmass) <: AbstractDict ? pmass : pmass(t)
    gvec = typeof(gvec) <: AbstractVector ? SVector{3}(gvec) : SVector{3}(gvec(t))
    vb_p = typeof(vb_p) <: AbstractVector ? SVector{3}(vb_p) : SVector{3}(vb_p(t))
    ωb_p = typeof(ωb_p) <: AbstractVector ? SVector{3}(ωb_p) : SVector{3}(ωb_p(t))

    return (; pcond, dload, pmass, gvec, vb_p, ωb_p)
end

function get_differential_vars(system::DynamicSystem, assembly, prescribed_conditions,
    point_masses, two_dimensional)

    # unpack pre-allocated storage and pointers for system
    @unpack x, M, force_scaling, indices, t = system

    # current parameters
    pcond = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(t)
    pmass = typeof(point_masses) <: AbstractDict ? point_masses : point_masses(t)

    # solve for the system mass matrix
    system_mass_matrix!(M, x, indices, two_dimensional, force_scaling, assembly, pcond, pmass)

    # identify differential variables
    differential_vars = dropdims(.!(iszero.(sum(M, dims=1))), dims=1)

    return differential_vars
end
