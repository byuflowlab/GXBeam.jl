"""
    ODEProblem(system::GXBeam.AbstractSystem, assembly, tspan; kwargs...)

Construct an `ODEProblem` corresponding to the system of nonlinear beams contained in
`assembly`.  This problem may be solved using the `DifferentialEquations` package

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
        If time varying, this input may be provided as a function of time.
 - `angular_velocity = zeros(3)`: Prescribed angular velocity of the body frame.
        If time varying, this input may be provided as a function of time.
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
 - `xpfunc = (x, p, t) -> (;)`: Similar to `pfunc`, except that parameters can also be
        defined as a function of GXBeam's state variables.  Using this function forces
        the system jacobian to be computed using automatic differentiation.
 - `pfunc = (p, t) -> (;)`: Function which returns a named tuple with fields corresponding
         to updated versions of the arguments `assembly`, `prescribed_conditions`,
         `distributed_loads`, `point_masses`, `linear_velocity`, `angular_velocity`, and
         `gravity`. Only fields contained in the resulting named tuple will be overwritten.
 - `p`: Sensitivity parameters, as defined in conjunction with the keyword argument `pfunc`.
         While not necessary, using `pfunc` and `p` to define the arguments to this function
         allows automatic differentiation sensitivities to be computed more efficiently

Additional keyword arguments are passed on to the ODEProblem constructor.
"""
function SciMLBase.ODEProblem(system::AbstractSystem, assembly, tspan;
    # general keyword arguments
    prescribed_conditions=Dict{Int,PrescribedConditions{Float64}}(),
    distributed_loads=Dict{Int,DistributedLoads{Float64}}(),
    point_masses=Dict{Int,PointMass{Float64}}(),
    linear_velocity=(@SVector zeros(3)),
    angular_velocity=(@SVector zeros(3)),
    gravity=(@SVector zeros(3)),
    # control flag keyword arguments
    initial_state=nothing,
    structural_damping=true,
    two_dimensional=false,
    constant_mass_matrix=typeof(system) <: ExpandedSystem,
    sparse=false,
    # sensitivity analysis keyword arguments
    xpfunc = nothing,
    pfunc = (p, t) -> (;),
    p = nothing,
    # additional keyword arguments (passed to ODEProblem constructor)
    kwargs...)

    # check if provided system is consistent with the provided keyword arguments
    if constant_mass_matrix
        @assert typeof(system) <: ExpandedSystem
    else
        @assert typeof(system) <: DynamicSystem
    end

    if isnothing(initial_state)
        # use stored state vector
        u0 = system.x
    else
        # initialize new state vector
        u0 = similar(system.x, promote_type(eltype(system), eltype(initial_state))) .= system.x

        # set current time
        t = tspan[1]

        # get prescribed conditions
        parameters = isnothing(xpfunc) ? pfunc(p, t) : xpfunc(x, p, t)
        pcond = get(parameters, :prescribed_conditions, prescribed_conditions)
        pcond = typeof(pcond) <: AbstractDict ? pcond : pcond(t)

        # set state variables to provided values
        set_state!(u0, system, initial_state, pcond)
    end

    # construct ODEFunction
    func = SciMLBase.ODEFunction(system, assembly;
        # general keyword arguments
        prescribed_conditions=prescribed_conditions,
        distributed_loads=distributed_loads,
        point_masses=point_masses,
        linear_velocity=linear_velocity,
        angular_velocity=angular_velocity,
        gravity=gravity,
        # control flag keyword arguments
        structural_damping=structural_damping,
        two_dimensional=two_dimensional,
        constant_mass_matrix=constant_mass_matrix,
        sparse=sparse,
        # sensitivity analysis keyword arguments
        xpfunc=xpfunc,
        pfunc=pfunc,
        p=p)

    # return ODEProblem
    return SciMLBase.ODEProblem{true}(func, u0, tspan, p; kwargs...)
end

"""
    SciMLBase.ODEFunction(system::AbstractSystem, assembly;
    # general keyword arguments
    prescribed_conditions=Dict{Int,PrescribedConditions{Float64}}(),
    distributed_loads=Dict{Int,DistributedLoads{Float64}}(),
    point_masses=Dict{Int,PointMass{Float64}}(),
    linear_velocity=(@SVector zeros(3)),
    angular_velocity=(@SVector zeros(3)),
    gravity=(@SVector zeros(3)),
    # control flag keyword arguments
    structural_damping=true,
    two_dimensional=false,
    constant_mass_matrix=true,
    sparse=false,
    # sensitivity analysis keyword arguments
    xpfunc = nothing,
    pfunc = (p, t) -> (;),
    p = nothing,
    # additional keyword arguments (passed to ODEFunction constructor)
    kwargs...)

Construct a `ODEFunction` for the system of nonlinear beams contained in `assembly`
which may be used with the DifferentialEquations package.

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

Additional keyword arguments are passed on to the ODEFunction constructor.
"""
function SciMLBase.ODEFunction(system::AbstractSystem, assembly;
    # general keyword arguments
    prescribed_conditions=Dict{Int,PrescribedConditions{Float64}}(),
    distributed_loads=Dict{Int,DistributedLoads{Float64}}(),
    point_masses=Dict{Int,PointMass{Float64}}(),
    linear_velocity=(@SVector zeros(3)),
    angular_velocity=(@SVector zeros(3)),
    gravity=(@SVector zeros(3)),
    # control flag keyword arguments
    structural_damping=true,
    two_dimensional=false,
    constant_mass_matrix=true,
    sparse=false,
    # sensitivity analysis keyword arguments
    xpfunc = nothing,
    pfunc = (p, t) -> (;),
    p = nothing,
    # additional keyword arguments (passed to ODEFunction constructor)
    kwargs...)

    # --- Input Argument Checks --- #

    # check if provided system is consistent with the provided keyword arguments
    if constant_mass_matrix
        @assert typeof(system) <: ExpandedSystem
    else
        @assert typeof(system) <: DynamicSystem
    end

    # check that compliance and mass matrices are invertible
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

    # unpack force scaling term
    force_scaling = system.force_scaling

    # --- Define mass matrix, residual, and jacobian --- #

    if constant_mass_matrix

        # system indices
        indices = SystemIndices(assembly.start, assembly.stop, static=false, expanded=true)

        # default keyword arguments
        constants = (; assembly, indices, two_dimensional, structural_damping, force_scaling,
            xpfunc, pfunc, prescribed_conditions, distributed_loads, point_masses,
            linear_velocity, angular_velocity, gravity, t=0.0)

        # set state rate vector
        du = zeros(indices.nstates) # state rate vector must be zero
        u = rand(indices.nstates) # state vector can be anything
        t = rand() # time can be anything

        # mass matrix (constant mass matrix system)
        TF = eltype(system)
        nx = indices.nstates
        mass_matrix = zeros(TF, nx, nx)
        expanded_mass_matrix!(mass_matrix, p, constants)
        mass_matrix .*= -1

        # residual function
        f = (resid, u, p, t) -> expanded_dynamic_residual!(resid, du, u, p, (; constants..., t))

        # jacobian function
        update_jacobian! = (jacob, u, p, t) -> expanded_dynamic_jacobian!(jacob, du, u, p, (; constants..., t))

        # jacobian prototype
        if sparse
            jac_prototype = spzeros(TF, nx, nx)
            update_jacobian!(jac_prototype, u, p, t)
        else
            jac_prototype = nothing
        end

    else

        # system indices
        indices = SystemIndices(assembly.start, assembly.stop, static=false, expanded=false)

        # default keyword arguments
        constants = (; assembly, indices, two_dimensional, structural_damping, force_scaling,
            xpfunc, pfunc, prescribed_conditions, distributed_loads, point_masses,
            linear_velocity, angular_velocity, gravity, t=0.0)

        # set state rate vector
        du = zeros(indices.nstates) # state rate vector must be zero
        u = rand(indices.nstates) # state vector can be anything
        t = rand() # time can be anything

        # mass matrix
        TF = eltype(system)
        nx = indices.nstates
        
        update_mass_matrix! = (jacob, x, p, t) -> begin
            # zero out all mass matrix entries
            jacob .= 0.0 

            # compute mass matrix
            mass_matrix!(jacob, x, p, (; constants..., t))

            # change sign of mass matrix
            jacob .*= -1
        end
        
        M = zeros(TF, nx, nx)
        mass_matrix = SciMLBase.MatrixOperator(M, update_func! = update_mass_matrix!) 

        # residual function
        f = (resid, u, p, t) -> dynamic_residual!(resid, du, u, p, (; constants..., t))

        # jacobian function
        update_jacobian! = (jacob, u, p, t) -> dynamic_jacobian!(jacob, du, u, p, (; constants..., t))

        # jacobian prototype
        if sparse
            jac_prototype = spzeros(TF, nx, nx)
            update_jacobian!(jac_prototype, u, p, t)
        else
            jac_prototype = nothing
        end

    end

    if isnothing(xpfunc)
        odefunc = SciMLBase.ODEFunction{true,true}(f; mass_matrix = mass_matrix,
            jac = update_jacobian!, jac_prototype = jac_prototype)
    else
        odefunc = SciMLBase.ODEFunction{true,true}(f; mass_matrix = mass_matrix)
    end

    return odefunc
end

"""
    DAEProblem(system::GXBeam.AbstractSystem, assembly, tspan; kwargs...)

Construct an `DAEProblem` corresponding to the system of nonlinear beams contained in
`assembly`.  This problem may be solved using the `DifferentialEquations` package

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
        If time varying, this input may be provided as a function of time.
 - `angular_velocity = zeros(3)`: Prescribed angular velocity of the body frame.
        If time varying, this input may be provided as a function of time.
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
 - `xpfunc = (x, p, t) -> (;)`: Similar to `pfunc`, except that parameters can also be
        defined as a function of GXBeam's state variables.  Using this function forces
        the system jacobian to be computed using automatic differentiation.
 - `pfunc = (p, t) -> (;)`: Function which returns a named tuple with fields corresponding
         to updated versions of the arguments `assembly`, `prescribed_conditions`,
         `distributed_loads`, `point_masses`, `linear_velocity`, `angular_velocity`, and
         `gravity`. Only fields contained in the resulting named tuple will be overwritten.
 - `p`: Sensitivity parameters, as defined in conjunction with the keyword argument `pfunc`.
         While not necessary, using `pfunc` and `p` to define the arguments to this function
         allows automatic differentiation sensitivities to be computed more efficiently

Additional keyword arguments are passed on to the DAEProblem constructor.
"""
function SciMLBase.DAEProblem(system::AbstractSystem, assembly, tspan;
    # general keyword arguments
    prescribed_conditions=Dict{Int,PrescribedConditions{Float64}}(),
    distributed_loads=Dict{Int,DistributedLoads{Float64}}(),
    point_masses=Dict{Int,PointMass{Float64}}(),
    linear_velocity=(@SVector zeros(3)),
    angular_velocity=(@SVector zeros(3)),
    gravity=(@SVector zeros(3)),
    # control flag keyword arguments
    initial_state=nothing,
    structural_damping=true,
    two_dimensional=false,
    constant_mass_matrix=false,
    sparse=false,
    # sensitivity analysis keyword arguments
    xpfunc = nothing,
    pfunc = (p, t) -> (;),
    p = nothing,
    # additional keyword arguments (passed to ODEProblem constructor)
    kwargs...)

    if isnothing(initial_state)
        # use stored state and rate vector
        dx0 = system.dx
        x0 = system.x
    else
        # initialize new state and rate vector
        dx0 = similar(system.dx, promote_type(eltype(system), eltype(initial_state))) .= system.dx
        x0 = similar(system.x, promote_type(eltype(system), eltype(initial_state))) .= system.x
    end

    # set current time
    t = tspan[1]

    # get prescribed conditions and point masses
    parameters = isnothing(xpfunc) ? pfunc(p, t) : xpfunc(x0, p, t)
    pcond = get(parameters, :prescribed_conditions, prescribed_conditions)
    pmass = get(parameters, :point_masses, point_masses)
    pcond = typeof(pcond) <: AbstractDict ? pcond : pcond(t)
    pmass = typeof(point_masses) <: AbstractDict ? pmass : pmass(t)

    if !isnothing(initial_state)
        # set state and rate variables to provided values
        set_state!(x0, system, initial_state, pcond)
        set_rate!(dx0, system, initial_state, pcond)
    end

    # define differential variables
    if constant_mass_matrix
        indices = SystemIndices(assembly.start, assembly.stop, static=false, expanded=true)
        differential_vars = expanded_differential_vars(indices, two_dimensional, assembly, pcond, pmass)
    else
        indices = SystemIndices(assembly.start, assembly.stop, static=false, expanded=false)
        differential_vars = dynamic_differential_vars(x0, indices, two_dimensional, assembly, pcond, pmass)
    end

    # construct DAEFunction
    func = SciMLBase.DAEFunction(system, assembly;
        # general keyword arguments
        prescribed_conditions=prescribed_conditions,
        distributed_loads=distributed_loads,
        point_masses=point_masses,
        linear_velocity=linear_velocity,
        angular_velocity=angular_velocity,
        gravity=gravity,
        # control flag keyword arguments
        structural_damping=structural_damping,
        two_dimensional=two_dimensional,
        constant_mass_matrix=constant_mass_matrix,
        sparse=sparse,
        # sensitivity analysis keyword arguments
        xpfunc=xpfunc,
        pfunc=pfunc,
        p=p,)

    # return DAEProblem
    return SciMLBase.DAEProblem{true}(func, dx0, x0, tspan, p; differential_vars, kwargs...)
end

"""
    DAEFunction(system::GXBeam.AbstractSystem, assembly; kwargs...)

Construct a `DAEFunction` for the system of nonlinear beams contained in `assembly`
which may be used with the DifferentialEquations package.

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
 - `structural_damping = true`: Flag indicating whether to enable structural damping
 - `two_dimensional = false`: Flag indicating whether to constrain results to the x-y plane
 - `constant_mass_matrix = true`: Flag indicating whether to use a constant mass matrix.
 - `sparse = false`: Flag indicating whether to use a sparse jacobian.

# Sensitivity Analysis Keyword Arguments
 - `xpfunc = (x, p, t) -> (;)`: Similar to `pfunc`, except that parameters can also be
        defined as a function of GXBeam's state variables.  Using this function forces
        the system jacobian to be computed using automatic differentiation.
 - `pfunc = (p, t) -> (;)`: Function which returns a named tuple with fields corresponding
        to updated versions of the arguments `assembly`, `prescribed_conditions`,
        `distributed_loads`, `point_masses`, `linear_velocity`, `angular_velocity`, and
        `gravity`. Only fields contained in the resulting named tuple will be overwritten.
 - `p`: Sensitivity parameters, as defined in conjunction with the keyword argument `pfunc`.
        While not necessary, using `pfunc` and `p` to define the arguments to this function
        allows automatic differentiation sensitivities to be computed more efficiently

Additional keyword arguments are passed on to the DAEFunction constructor.
"""
function SciMLBase.DAEFunction(system::AbstractSystem, assembly;
    # general keyword arguments
    prescribed_conditions=Dict{Int,PrescribedConditions{Float64}}(),
    distributed_loads=Dict{Int,DistributedLoads{Float64}}(),
    point_masses=Dict{Int,PointMass{Float64}}(),
    linear_velocity=(@SVector zeros(3)),
    angular_velocity=(@SVector zeros(3)),
    gravity=(@SVector zeros(3)),
    # control flag keyword arguments
    structural_damping=true,
    two_dimensional=false,
    constant_mass_matrix=false,
    sparse=false,
    # sensitivity analysis keyword arguments
    xpfunc = nothing,
    pfunc = (p, t) -> (;),
    p = nothing,
    # additional keyword arguments (passed to ODEFunction constructor)
    kwargs...)

    # --- Input Argument Checks --- #

    for ielem = 1:length(assembly.elements)
        @assert !iszero(assembly.elements[ielem].L) "Zero length elements cannot be used "*
            "with DifferentialEquations"
        @assert !iszero(assembly.elements[ielem].compliance[1,:]) &&
            !iszero(assembly.elements[ielem].compliance[2,:]) &&
            !iszero(assembly.elements[ielem].compliance[3,:]) &&
            !iszero(assembly.elements[ielem].compliance[4,:]) &&
            !iszero(assembly.elements[ielem].compliance[5,:]) &&
            !iszero(assembly.elements[ielem].compliance[6,:]) "Compliance matrix must "*
            "be invertible when using DAEFunction."
        @assert !iszero(assembly.elements[ielem].mass[1,:]) &&
            !iszero(assembly.elements[ielem].mass[2,:]) &&
            !iszero(assembly.elements[ielem].mass[3,:]) &&
            !iszero(assembly.elements[ielem].mass[4,:]) &&
            !iszero(assembly.elements[ielem].mass[5,:]) &&
            !iszero(assembly.elements[ielem].mass[6,:]) "Mass matrix must "*
            "be invertible when using DAEFunction."
    end

    # check if provided system is consistent with the provided keyword arguments
    if constant_mass_matrix
        @assert typeof(system) <: ExpandedSystem
    else
        @assert typeof(system) <: DynamicSystem
    end

    # unpack force scaling term
    force_scaling = system.force_scaling

    # --- Define mass matrix, residual, and jacobian --- #

    if constant_mass_matrix

        # system indices
        indices = SystemIndices(assembly.start, assembly.stop, static=false, expanded=true)

        # default keyword arguments
        constants = (; assembly, indices, two_dimensional, structural_damping, force_scaling,
            xpfunc, pfunc, prescribed_conditions, distributed_loads, point_masses,
            linear_velocity, angular_velocity, gravity, t=0.0)

        # residual function
        f = (resid, du, u, p, t) -> expanded_dynamic_residual!(resid, du, u, p, (; constants..., t))

        # jacobian function
        update_jacobian! = function(J, du, u, p, gamma, t)

            # zero out all jacobian entries
            J .= 0.0

            # unpack indices and control flags
            @unpack indices, structural_damping, two_dimensional, force_scaling = constants

            # combine constants and parameters
            assembly, pcond, dload, pmass, gvec, vb_p, ωb_p = dynamic_parameters(x, p, constants)

            # compute and return the residual
            expanded_dynamic_system_jacobian!(resid, du, u, indices, two_dimensional, force_scaling,
                structural_damping, assembly, pcond, dload, pmass, gvec, vb_p, ωb_p)

            # add gamma multiplied by the mass matrix
            expanded_system_mass_matrix!(J, gamma, indices, two_dimensional, force_scaling,
                assembly, pcond, pmass)

            return J
        end

        # jacobian prototype
        if sparse
            jac_prototype = spzeros(eltype(system), indices.nstates, indices.nstates)
        else
            jac_prototype = nothing
        end

    else

        # system indices
        indices = SystemIndices(assembly.start, assembly.stop, static=false, expanded=false)

        # default keyword arguments
        constants = (; assembly, indices, two_dimensional, structural_damping, force_scaling,
            xpfunc, pfunc, prescribed_conditions, distributed_loads, point_masses,
            linear_velocity, angular_velocity, gravity, t=0.0)

        # residual function
        f = (resid, du, u, p, t) -> dynamic_residual!(resid, du, u, p, (; constants..., t))

        # jacobian function
        update_jacobian! = function(J, du, u, p, gamma, t)

            # zero out all jacobian entries
            J .= 0.0

            # unpack indices and control flags
            @unpack indices, structural_damping, two_dimensional, force_scaling = constants

            # combine constants and parameters
            assembly, pcond, dload, pmass, gvec, vb_p, ωb_p = dynamic_parameters(x, p, constants)

            # update acceleration state variable indices
            update_body_acceleration_indices!(indices, pcond)

            # compute and return the residual
            dynamic_system_jacobian!(resid, dx, x, indices, two_dimensional, force_scaling,
                structural_damping, assembly, pcond, dload, pmass, gvec, vb_p, ωb_p)

            # add gamma multiplied by the mass matrix
            system_mass_matrix!(J, gamma, u, indices, two_dimensional, force_scaling,
                assembly, pcond, pmass)

            return J
        end

        # jacobian prototype
        if sparse
            jac_prototype = spzeros(eltype(system), indices.nstates, indices.nstates)
        else
            jac_prototype = nothing
        end

    end

    return SciMLBase.DAEFunction{true,true}(f) # TODO: re-add jacobian here once supported
end

# combines constant and variable parameters for a dynamic analysis
function dynamic_parameters(x, p, constants)

    # unpack default parameters, parameter function, and current time
    @unpack assembly, prescribed_conditions, distributed_loads, point_masses, gravity,
        linear_velocity, angular_velocity, xpfunc, pfunc, t = constants

    # overwrite default assembly and parameters (if applicable)
    parameters = isnothing(xpfunc) ? pfunc(p, t) : xpfunc(x, p, t)
    assembly = get(parameters, :assembly, assembly)
    prescribed_conditions = get(parameters, :prescribed_conditions, prescribed_conditions)
    distributed_loads = get(parameters, :distributed_loads, distributed_loads)
    point_masses = get(parameters, :point_masses, point_masses)
    gravity = get(parameters, :gravity, gravity)
    linear_velocity = get(parameters, :linear_velocity, linear_velocity)
    angular_velocity = get(parameters, :angular_velocity, angular_velocity)

    # get parameters corresponding to this time step
    pcond = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(t)
    dload = typeof(distributed_loads) <: AbstractDict ? distributed_loads : distributed_loads(t)
    pmass = typeof(point_masses) <: AbstractDict ? point_masses : point_masses(t)
    gvec = typeof(gravity) <: AbstractVector ? SVector{3}(gravity) : SVector{3}(gravity(t))
    vb_p = typeof(linear_velocity) <: AbstractVector ? SVector{3}(linear_velocity) : SVector{3}(linear_velocity(t))
    ωb_p = typeof(angular_velocity) <: AbstractVector ? SVector{3}(angular_velocity) : SVector{3}(angular_velocity(t))

    # update body acceleration frame indices
    update_body_acceleration_indices!(constants.indices, pcond)

    return assembly, pcond, dload, pmass, gvec, vb_p, ωb_p
end

# residual function for a dynamic analysis
function dynamic_residual!(resid, dx, x, p, constants)

    # unpack indices and control flags
    @unpack indices, structural_damping, two_dimensional, force_scaling = constants

    # combine constants and parameters
    assembly, pcond, dload, pmass, gvec, vb_p, ωb_p = dynamic_parameters(x, p, constants)

    # update acceleration state variable indices
    update_body_acceleration_indices!(indices, pcond)

    # compute and return the residual
    return dynamic_system_residual!(resid, dx, x, indices, two_dimensional, force_scaling,
        structural_damping, assembly, pcond, dload, pmass, gvec, vb_p, ωb_p)
end

# jacobian function for a dynamic analysis
function dynamic_jacobian!(jacob, dx, x, p, constants)

    # zero out all jacobian entries
    jacob .= 0.0

    # unpack indices and control flags
    @unpack indices, structural_damping, two_dimensional, force_scaling = constants

    # combine constants and parameters
    assembly, pcond, dload, pmass, gvec, vb_p, ωb_p = dynamic_parameters(x, p, constants)

    # update acceleration state variable indices
    update_body_acceleration_indices!(indices, pcond)

    # compute and return the residual
    return dynamic_system_jacobian!(jacob, dx, x, indices, two_dimensional, force_scaling,
        structural_damping, assembly, pcond, dload, pmass, gvec, vb_p, ωb_p)
end

# differential variables for a constant mass matrix system
function dynamic_differential_vars(x, indices, two_dimensional, assembly, pcond, pmass)

    # get floating point type
    TF = eltype(assembly)

    # intialize temporary mass matrix
    M = spzeros(TF, indices.nstates, indices.nstates)

    # use arbitrary force scaling parameter
    force_scaling = 1.0

    # solve for the system mass matrix
    system_mass_matrix!(M, x, indices, two_dimensional, force_scaling, assembly, pcond, pmass)

    # identify differential variables
    differential_vars = dropdims(.!(iszero.(sum(M, dims=1))), dims=1)

    return differential_vars
end

# residual function for a constant mass matrix system
function expanded_dynamic_residual!(resid, dx, x, p, constants)

    # unpack indices and control flags
    @unpack indices, structural_damping, two_dimensional, force_scaling = constants

    # combine constants and parameters
    assembly, pcond, dload, pmass, gvec, vb_p, ωb_p = dynamic_parameters(x, p, constants)

    # compute and return the residual
    return expanded_dynamic_system_residual!(resid, dx, x, indices, two_dimensional, force_scaling,
        structural_damping, assembly, pcond, dload, pmass, gvec, vb_p, ωb_p)
end

# jacobian function for a constant mass matrix system
function expanded_dynamic_jacobian!(jacob, dx, x, p, constants)

    # zero out all jacobian entries
    jacob .= 0.0

    # unpack indices and control flags
    @unpack indices, structural_damping, two_dimensional, force_scaling = constants

    # combine constants and parameters
    assembly, pcond, dload, pmass, gvec, vb_p, ωb_p = dynamic_parameters(x, p, constants)

    # compute and return the residual
    result = expanded_dynamic_system_jacobian!(jacob, dx, x, indices, two_dimensional, force_scaling,
        structural_damping, assembly, pcond, dload, pmass, gvec, vb_p, ωb_p)

    return result
end

# differential variables for a constant mass matrix system
function expanded_differential_vars(indices, two_dimensional, assembly, pcond, pmass)

    # get floating point type
    TF = eltype(assembly)

    # intialize temporary mass matrix
    M = spzeros(TF, indices.nstates, indices.nstates)

    # use arbitrary force scaling parameter
    force_scaling = 1.0

    # solve for the system mass matrix
    expanded_system_mass_matrix!(M, indices, two_dimensional, force_scaling, assembly, pcond, pmass)

    # identify differential variables
    differential_vars = dropdims(.!(iszero.(sum(M, dims=1))), dims=1)

    return differential_vars
end