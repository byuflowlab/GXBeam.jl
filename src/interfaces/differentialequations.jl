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
"""
function SciMLBase.ODEProblem(system::System, assembly, tspan; kwargs...)

    # use initial state from `system`
    u0 = vcat(system.x, zeros(12*length(assembly.elements)))

    # create ODEFunction
    func = SciMLBase.ODEFunction(system, assembly; kwargs...)

    return SciMLBase.ODEProblem{true}(func, u0, tspan)
end

"""
    ODEFunction(system::GXBeam.System, assembly; kwargs...)

Construct a `ODEFunction` for the system of nonlinear beams contained in `assembly` which 
may be used with the DifferentialEquations package.

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
"""
function SciMLBase.ODEFunction(system::System, assembly;
    prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}(),
    distributed_loads = Dict{Int,DistributedLoads{Float64}}(),
    point_masses = Dict{Int,Vector{PointMass{Float64}}}(),
    gravity = (@SVector zeros(3)),
    origin = (@SVector zeros(3)),
    linear_velocity = (@SVector zeros(3)),
    angular_velocity = (@SVector zeros(3)),
    linear_acceleration = (@SVector zeros(3)),
    angular_acceleration = (@SVector zeros(3)),
    )

    # system dimensions
    N = length(system.x)
    nelem = length(assembly.elements)
    npoint = length(assembly.points)

    # system pointers
    irow_point = system.irow_point
    irow_elem = system.irow_elem
    irow_elem1 = system.irow_elem1
    irow_elem2 = system.irow_elem2
    icol_point = system.icol_point
    icol_elem = system.icol_elem

    # unpack scaling parameters
    force_scaling = system.force_scaling

    # DAE function
    f = function(dx, x, p, t)

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

        for ielem = 1:nelem
    
            # extract element
            elem = assembly.elements[ielem]

            # get pointers for element
            icol = icol_elem[ielem]
            irow_e = irow_elem[ielem]
            irow_e1 = irow_elem1[ielem]
            irow_p1 = irow_point[assembly.start[ielem]]
            irow_e2 = irow_elem2[ielem]
            irow_p2 = irow_point[assembly.stop[ielem]]

            # compute steady state constraints
            steady_state_element_residual!(dx, x, ielem, elem, dload, 
                pmass, gvec, force_scaling, icol, irow_e, irow_e1, irow_p1, irow_e2, irow_p2, 
                x0, v0, ω0, a0, α0)

            # extract state variable element resultants
            u1 = SVector(x[N+12*(ielem-1)+1], x[N+12*(ielem-1)+2], x[N+12*(ielem-1)+3])
            ψ1 = SVector(x[N+12*(ielem-1)+4], x[N+12*(ielem-1)+5], x[N+12*(ielem-1)+6])
            u2 = SVector(x[N+12*(ielem-1)+7], x[N+12*(ielem-1)+8], x[N+12*(ielem-1)+9])
            ψ2 = SVector(x[N+12*(ielem-1)+10], x[N+12*(ielem-1)+11], x[N+12*(ielem-1)+12])

            # extract calculated element resultants
            f_u1 = SVector(dx[irow_p1], dx[irow_p1+1], dx[irow_p1+2])
            f_ψ1 = SVector(dx[irow_p1+3], dx[irow_p1+4], dx[irow_p1+5])
            f_u2 = SVector(dx[irow_p2], dx[irow_p2+1], dx[irow_p2+2])
            f_ψ2 = SVector(dx[irow_p2+3], dx[irow_p2+4], dx[irow_p2+5])
            
            # correct element resultants
            θ = SVector(x[icol+3], x[icol+4], x[icol+5])
            V = SVector(x[icol+12], x[icol+13], x[icol+14])
            Ω = SVector(x[icol+15], x[icol+16], x[icol+17])
            P = element_linear_momentum(elem, V, Ω)
            H = element_angular_momentum(elem, V, Ω)
            C = get_C(θ)
            Cab = elem.Cab
            θdot = SVector(dx[irow_e+3], dx[irow_e+4], dx[irow_e+5])
            Cdot = get_C_t(C, θ, θdot)
            CbaC = C*Cab'

            tmp1 = elem.L/2*Cdot'*Cab*P/force_scaling
            tmp2 = elem.L/2*Cdot'*Cab*H/force_scaling

            f_u1 += tmp1
            f_ψ1 += tmp2
            f_u2 += tmp1
            f_ψ2 += tmp2

            # construct element equilibrium constraints
            dx[N+12*(ielem-1)+1:N+12*(ielem-1)+3] .= CbaC*(f_u1 - u1)
            dx[N+12*(ielem-1)+4:N+12*(ielem-1)+6] .= CbaC*(f_ψ1 - ψ1)
            dx[N+12*(ielem-1)+7:N+12*(ielem-1)+9] .= CbaC*(f_u2 - u2)
            dx[N+12*(ielem-1)+10:N+12*(ielem-1)+12] .= CbaC*(f_ψ2 - ψ2)

            # initialize node equilibrium constraints
            dx[irow_p1:irow_p1+2] .= u1
            dx[irow_p1+3:irow_p1+5] .= ψ1
            dx[irow_p2:irow_p2+2] .= u2
            dx[irow_p2+3:irow_p2+5] .= ψ2
            
        end
    
        # add to node equilibrium and compatability constraints
        for ipoint = 1:npoint
    
            # skip if the unknowns have been eliminated from the system of equations
            if icol_point[ipoint] <= 0
                continue
            end
    
            icol = icol_point[ipoint]
            irow_p = irow_point[ipoint]
    
            point_residual!(dx, x, ipoint, assembly, pcond, force_scaling, icol, irow_p, 
                irow_elem1, irow_elem2)
        end
    
        return dx
    end

    # mass matrix
    mass_matrix = similar(system.M, N+12*nelem, N+12*nelem)
    mass_matrix .= 0
    system_mass_matrix!(mass_matrix, zeros(N), assembly, point_masses, force_scaling, 
    irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem)
    for ielem = 1:nelem
        # get pointers for element
        icol = icol_elem[ielem]
        irow_p1 = irow_point[assembly.start[ielem]]
        irow_p2 = irow_point[assembly.stop[ielem]]
        # move mass matrix entries to new equations
        mass_matrix[N+12*(ielem-1)+1:N+12*(ielem-1)+6, icol:icol+17] .= mass_matrix[irow_p1:irow_p1+5, icol:icol+17]
        mass_matrix[N+12*(ielem-1)+7:N+12*(ielem-1)+12, icol:icol+17] .= mass_matrix[irow_p2:irow_p2+5, icol:icol+17]
        # zero out previous contribution
        mass_matrix[irow_p1:irow_p1+5, icol:icol+17] .= 0
        mass_matrix[irow_p2:irow_p2+5, icol:icol+17] .= 0
    end
    dropzeros!(mass_matrix)

    return SciMLBase.ODEFunction{true,true}(f; mass_matrix = collect(mass_matrix))
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
    gravity = (@SVector zeros(3)),
    origin = (@SVector zeros(3)),
    linear_velocity = (@SVector zeros(3)),
    angular_velocity = (@SVector zeros(3)),
    linear_acceleration = (@SVector zeros(3)),
    angular_acceleration = (@SVector zeros(3)),
    )

    # create SciMLBase.DAEFunction
    func = SciMLBase.DAEFunction(system, assembly)

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
    differential_vars = get_differential_vars(system)

    return SciMLBase.DAEProblem{true}(func, du0, u0, tspan, p; differential_vars)
end

"""
    DAEFunction(system::GXBeam.System, assembly)

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
function SciMLBase.DAEFunction(system::System, assembly)

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
            distributed_loads, point_masses, gvec, force_scaling, irow_point, irow_elem, 
            irow_elem1, irow_elem2, icol_point, icol_elem, x0, v0, ω0, a0, α0)

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
            distributed_loads, point_masses, gvec, force_scaling, irow_point, irow_elem, 
            irow_elem1, irow_elem2, icol_point, icol_elem, x0, v0, ω0, a0, α0)

        # add gamma multiplied by the mass matrix
        system_mass_matrix!(J, gamma, u, assembly, point_masses, force_scaling,
            irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem)

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

function get_differential_vars(system::System)
    differential_vars = fill(false, length(system.x))
    for icol in system.icol_elem
        differential_vars[icol:icol+2] .= true # u (for the beam element)
        differential_vars[icol+3:icol+5] .= true # θ (for the beam element)
        differential_vars[icol+12:icol+14] .= true # V (for the beam element)
        differential_vars[icol+15:icol+17] .= true # Ω (for the beam element)
    end
    return differential_vars
end
