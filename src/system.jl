"""
    System{TF, TV<:AbstractVector{TF}, TM<:AbstractMatrix{TF}, TTF<:AbstractVector{TF}}

Contains the system state, residual vector, and jacobian matrices as well as
pointers to be able to access their contents.  Also contains additional storage
needed for time domain simulations.

# Fields:
 - `static`: Flag indicating whether system matrices are only valid for static analyses
 - `x`: State vector
 - `r`: Residual vector
 - `K`: System jacobian matrix with respect to the state variables
 - `M`: System jacobian matrix with respect to the time derivative of the state variables
 - `irow_pt`: Row index of first equilibrium equation for each point
 - `irow_beam`: Row index of first equation for just this beam element
 - `irow_beam1`: Row index of first equation for the left side of each beam
 - `irow_beam2`: Row index of first equation for the right side of each beam
 - `icol_pt`: Row/Column index of first state variable for each point
 - `icol_beam`: Row/Column index of first state variable for each beam element
 - `current_step`: Current time step
 - `udot_init`: `2/dt*u + udot` for each beam element from the previous time step
 - `θdot_init`: `2/dt*θ + θdot` for each beam element from the previous time step
 - `CtCabPdot_init`: `2/dt*C'*Cab*P + C'*Cab*Pdot` for each beam element from the previous time step
 - `CtCabHdot_init`: `2/dt*C'*Cab*H + C'*Cab*Hdot` for each beam element from the previous time step
"""
struct System{TF, TV<:AbstractVector{TF}, TM<:AbstractMatrix{TF}}
    static::Bool
    x::TV
    r::TV
    K::TM
    M::TM
    irow_pt::Vector{Int}
    irow_beam::Vector{Int}
    irow_beam1::Vector{Int}
    irow_beam2::Vector{Int}
    icol_pt::Vector{Int}
    icol_beam::Vector{Int}
    current_step::Ref{Int}
    udot_init::Vector{SVector{3,TF}}
    θdot_init::Vector{SVector{3,TF}}
    CtCabPdot_init::Vector{SVector{3,TF}}
    CtCabHdot_init::Vector{SVector{3,TF}}
end
Base.eltype(::System{TF, TV, TM}) where {TF, TV, TM} = TF

"""
    System([TF=eltype(assembly),] assembly, points, static)

Initialize an object of type `System` which stores the system state, residual vector,
current time function values,and jacobian matrices as well as pointers to
be able to access their contents.

# Arguments:
 - `TF:` (optional) Used to specify floating point type used by resulting `System` object
 - `assembly`: Assembly of rigidly connected nonlinear beam elements
 - `points`: Point indices which should be preserved in the system of equations.
        All points with prescribed conditions should be included.
 - `static`: Flag indicating whether system matrices will be used for static simulations
"""
System(assembly, points, static) = System(eltype(assembly), assembly, points, static)

function System(TF::Type{<:AbstractFloat}, assembly, points, static)

    # get number of beams
    nbeam = length(assembly.elements)

    # get the number of beams connected to each point.
    n_connections = point_connections(assembly)

    # get the size and pointers into the system matrices
    n, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam =
        system_indices(assembly, points, n_connections, static)

    # initialize system matrices
    x = zeros(TF, n)
    r = zeros(TF, n)
    K = spzeros(TF, n, n)
    M = spzeros(TF, n, n)

    # initialize pointer to value of current time step
    current_step = Ref(0)

    # initialize storage for time domain simulations
    udot_init = [@SVector zeros(TF, 3) for i = 1:nbeam]
    θdot_init = [@SVector zeros(TF, 3) for i = 1:nbeam]
    CtCabPdot_init = [@SVector zeros(TF, 3) for i = 1:nbeam]
    CtCabHdot_init = [@SVector zeros(TF, 3) for i = 1:nbeam]

    return System(static, x, r, K, M, irow_pt, irow_beam, irow_beam1, irow_beam2,
        icol_pt, icol_beam, current_step, udot_init, θdot_init, CtCabPdot_init,
        CtCabHdot_init)
end

"""
   reset_state!(system)

Reset the state variables in `system` (stored in `system.x`) to zero.
"""
function reset_state!(system)
    system.x .= 0
    return system
end

"""
    point_connections(assembly)

Count the number of beams connected to each point
"""
function point_connections(assembly)

    npt = length(assembly.points)

    n_connections = Array{Int, 1}(undef, npt)

    for ipt = 1:npt
        n_connections[ipt] = count(x -> x == ipt, assembly.start) + count(x -> x == ipt, assembly.stop)
    end

    return n_connections
end

"""
    system_indices(assembly, points, n_connections, static)

Solve for the row indices of the first equilibrium or compatability equations for
each point and side of each beam element.  Also solve for the row/column index of
each point and beam state variable.

Note that this function includes the following logic which reduces the size of
the system of equations where possible (without sacrificing any accuracy):

If only two beams meet at a point, the 6 unknowns associated with that point as
well as the 6 compatability equations are eliminated from the system, except if
specified in the array `points`.  Points for which unknowns have been
eliminated are assigned a column index of -1.  Beams for which the compatability
equations have been eliminated are also assigned an index of -1

# Arguments:
 - `assembly`: Assembly of rigidly connected nonlinear beam elements
 - `points`: Point indices which should be preserved in the system matrices
 - `n_connections`: Number of connections to each point
 - `static`: flag indicating whether analysis is static

# Return Arguments:
 - `n`: total number of equations/unknowns in the system
 - `irow_pt`: Row index of first equilibrium equation for each point
 - `irow_beam`: Row index of first equation for just this beam element
 - `irow_beam1`: Row index of first equation for the left side of each beam
 - `irow_beam2`: Row index of first equation for the right side of each beam
 - `icol_pt`: Column index of first state variable for each point
 - `icol_beam`: Column index of first state variable for each beam element
"""
function system_indices(assembly, points, n_connections, static)

    npt = length(assembly.points)
    nbeam = length(assembly.elements)

    assigned = fill(false, npt)

    irow_pt = Array{Int, 1}(undef, npt)
    irow_beam = Array{Int, 1}(undef, nbeam)
    irow_beam1 = Array{Int, 1}(undef, nbeam)
    irow_beam2 = Array{Int, 1}(undef, nbeam)

    icol_pt = Array{Int, 1}(undef, npt)
    icol_beam = Array{Int, 1}(undef, nbeam)

    irow = 1
    icol = 1
    for ibeam = 1:nbeam
        ipt = assembly.start[ibeam]
        if !assigned[ipt]
            assigned[ipt] = true
            # 6 equilibrium equations + 6 compatability equations
            irow_pt[ipt] = irow
            irow_beam1[ibeam] = irow
            irow += 12
            # add unknowns for each point
            if n_connections[ipt] == 2 && !(ipt in points)
                # no additional unknowns when only two beams meet at a point
                icol_pt[ipt] = -1
            else
                # 6 additional unknowns for each point
                icol_pt[ipt] = icol
                icol += 6
            end
        else
            # add compatability equations
            if n_connections[ipt] == 2 && !(ipt in points)
                # no additional equations when only two beams meet at a point
                irow_beam1[ibeam] = -1
            else
                # 6 additional compatibility equations
                irow_beam1[ibeam] = irow
                irow += 6
            end
        end

        # 12 unknowns for each element
        icol_beam[ibeam] = icol
        icol += 12

        if static
            irow_beam[ibeam] = -1
        else
            # 6 linear and angular momentum residual equations for each element
            irow_beam[ibeam] = irow
            irow += 6
            # 6 additional unknowns for each member for unsteady analyses
            icol += 6
        end

        ipt = assembly.stop[ibeam]
        if !assigned[ipt]
            assigned[ipt] = true
            # 6 equilibrium equations + 6 compatability equations
            irow_pt[ipt] = irow
            irow_beam2[ibeam] = irow
            irow += 12
            # add unknowns for each point
            if n_connections[ipt] == 2 && !(ipt in points)
                # no additional unknowns when only two beams meet at a point
                icol_pt[ipt] = -1
            else
                # 6 additional unknowns for each point
                icol_pt[ipt] = icol
                icol += 6
            end
        else
            if n_connections[ipt] == 2 && !(ipt in points)
                # no additional compatability equations when only two beams meet at a point
                irow_beam2[ibeam] = -1
            else
                # 6 additional compatibility equations
                irow_beam2[ibeam] = irow
                irow += 6
            end
        end
    end

    # number of unknowns/equations
    n = irow - 1

    return n, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam
end

"""
    system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads,
        istep, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam)
    system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads,
        istep, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam, x0, v0, ω0)
    system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads,
        istep, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam, x0, v0, ω0,
        u, θ, udot, θdot)
    system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads,
        istep, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam, x0, v0, ω0,
        udot, θdot_init, CtCabPdot, CtCabHdot, dt)

Populate the residual vector `resid` with the results of the residual equations
for the system.

There are four implementations corresponding to the following analysis types:
 - Static
 - Dynamic - Steady State
 - Dynamic - Initial Step (for initializing time domain simulations)
 - Dynamic - Time Marching

 See "GEBT: A general-purpose nonlinear analysis tool for composite beams" by
 Wenbin Yu and "Geometrically nonlinear analysis of composite beams using
 Wiener-Milenković parameters" by Qi Wang and Wenbin Yu.

# Arguments
 - `resid`: System residual vector
 - `x`: Current state variables of the system
 - `assembly`: Assembly of rigidly connected nonlinear beam elements
 - `prescribed_conditions`: Dictionary of prescribed conditions at all time steps
 - `distributed_loads`: Dictionary of distributed loads at all time steps
 - `istep`: Current time step
 - `irow_pt`: Row index of first equilibrium equation for each point
 - `irow_beam`: Row index of first equation for just this beam element
 - `irow_beam1`: Row index of first equation for the left side of each beam
 - `irow_beam2`: Row index of first equation for the right side of each beam
 - `icol_pt`: Column index of first state variable for each point
 - `icol_beam`: Column index of first state variable for each beam element

# Additional Arguments for Dynamic Analyses
 - `x0`: Global frame origin (for the current time step)
 - `v0`: Global frame linear velocity (for the current time step)
 - `ω0`: Global frame angular velocity (for the current time step)

# Additional Arguments for Initial Step Analyses
 - `u`: deflection variables for each beam element
 - `θ`: rotation variables for each beam element
 - `udot`: time derivative of u for each beam element
 - `θdot`: time derivative of θ for each beam element

# Additional Arguments for Time Marching Analyses
 - `udot`: `-2/dt*u - udot` for each beam element from the previous time step
 - `θdot_init`: `-2/dt*θ - θdot` for each beam element from the previous time step
 - `CtCabPdot`: `-2/dt*C'*Cab*P - C'*Cab*Pdot` for each beam element from the previous time step
 - `CtCabHdot`: `-2/dt*C'*Cab*H - C'*Cab*Hdot` for each beam element from the previous time step
 - `dt`: time step size
"""
system_residual!

# static
function system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads,
    istep, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam)

    npoint = length(assembly.points)
    nbeam = length(assembly.elements)

    # add contributions to residual equations from the elements
    for ibeam = 1:nbeam

        # get pointers for element
        icol = icol_beam[ibeam]
        irow_b = irow_beam[ibeam]
        irow_b1 = irow_beam1[ibeam]
        irow_p1 = irow_pt[assembly.start[ibeam]]
        irow_b2 = irow_beam2[ibeam]
        irow_p2 = irow_pt[assembly.stop[ibeam]]

        resid = element_residual!(resid, x, ibeam, assembly.elements[ibeam],
             distributed_loads, istep, icol, irow_b, irow_b1,
            irow_p1, irow_b2, irow_p2)
    end

    # add contributions to the residual equations from the prescribed point conditions
    for ipoint = 1:npoint

        # skip if the unknowns have been eliminated from the system of equations
        if icol_pt[ipoint] <= 0
            continue
        end

        icol = icol_pt[ipoint]
        irow_p = irow_pt[ipoint]

        resid = point_residual!(resid, x, ipoint, assembly, prescribed_conditions,
            istep, icol, irow_p, irow_beam1, irow_beam2)
    end

    return resid
end

# dynamic
function system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads,
    istep, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam,
    x0, v0, ω0)

    npoint = length(assembly.points)
    nbeam = length(assembly.elements)

    # add contributions to residual equations from the elements
    for ibeam = 1:nbeam

        # get pointers for element
        icol = icol_beam[ibeam]
        irow_b = irow_beam[ibeam]
        irow_b1 = irow_beam1[ibeam]
        irow_p1 = irow_pt[assembly.start[ibeam]]
        irow_b2 = irow_beam2[ibeam]
        irow_p2 = irow_pt[assembly.stop[ibeam]]

        resid = element_residual!(resid, x, ibeam, assembly.elements[ibeam],
             distributed_loads, istep, icol, irow_b, irow_b1,
            irow_p1, irow_b2, irow_p2, x0, v0, ω0)
    end

    # add contributions to the residual equations from the prescribed point conditions
    for ipoint = 1:npoint

        # skip if the unknowns have been eliminated from the system of equations
        if icol_pt[ipoint] <= 0
            continue
        end

        icol = icol_pt[ipoint]
        irow_p = irow_pt[ipoint]

        resid = point_residual!(resid, x, ipoint, assembly, prescribed_conditions,
            istep, icol, irow_p, irow_beam1, irow_beam2)
    end

    return resid
end

# initial step
function system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads,
    istep, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam,
    x0, v0, ω0, u, θ, udot, θdot)

    npoint = length(assembly.points)
    nbeam = length(assembly.elements)

    # add contributions to residual equations from the elements
    for ibeam = 1:nbeam

        # get pointers for element
        icol = icol_beam[ibeam]
        irow_b = irow_beam[ibeam]
        irow_b1 = irow_beam1[ibeam]
        irow_p1 = irow_pt[assembly.start[ibeam]]
        irow_b2 = irow_beam2[ibeam]
        irow_p2 = irow_pt[assembly.stop[ibeam]]

        resid = element_residual!(resid, x, ibeam, assembly.elements[ibeam],
             distributed_loads, istep, icol, irow_b, irow_b1,
            irow_p1, irow_b2, irow_p2, x0, v0, ω0, u, θ, udot, θdot)
    end

    # add contributions to the residual equations from the prescribed point conditions
    for ipoint = 1:npoint

        # skip if the unknowns have been eliminated from the system of equations
        if icol_pt[ipoint] <= 0
            continue
        end

        icol = icol_pt[ipoint]
        irow_p = irow_pt[ipoint]

        resid = point_residual!(resid, x, ipoint, assembly, prescribed_conditions,
            istep, icol, irow_p, irow_beam1, irow_beam2)
    end

    return resid
end

# time marching
function system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads,
    istep, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam,
    x0, v0, ω0, udot, θdot_init, CtCabPdot, CtCabHdot, dt)

    nbeam = length(assembly.elements)
    npoint = length(assembly.points)

    # add contributions to residual equations from the beam elements
    for ibeam = 1:nbeam

        # get pointers for element
        icol = icol_beam[ibeam]
        irow_b = irow_beam[ibeam]
        irow_b1 = irow_beam1[ibeam]
        irow_p1 = irow_pt[assembly.start[ibeam]]
        irow_b2 = irow_beam2[ibeam]
        irow_p2 = irow_pt[assembly.stop[ibeam]]

        resid = element_residual!(resid, x, ibeam, assembly.elements[ibeam],
             distributed_loads, istep, icol, irow_b, irow_b1,
            irow_p1, irow_b2, irow_p2, x0, v0, ω0, udot, θdot_init, CtCabPdot,
            CtCabHdot, dt)
    end

    # add contributions to the residual equations from the prescribed point conditions
    for ipoint = 1:npoint

        # skip if the unknowns have been eliminated from the system of equations
        if icol_pt[ipoint] <= 0
            continue
        end

        icol = icol_pt[ipoint]
        irow_p = irow_pt[ipoint]

        resid = point_residual!(resid, x, ipoint, assembly, prescribed_conditions,
            istep, icol, irow_p, irow_beam1, irow_beam2)
    end

    return resid
end

"""
    system_jacobian!(jacob, x, assembly, prescribed_conditions, distributed_loads,
        istep, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam)
    system_jacobian!(jacob, x, assembly, prescribed_conditions, distributed_loads,
        istep, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam,
        x0, v0, ω0)
    system_jacobian!(jacob, x, assembly, prescribed_conditions, distributed_loads,
        istep, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam,
        x0, v0, ω0, u, θ, udot, θdot)
    system_jacobian!(jacob, x, assembly, prescribed_conditions, distributed_loads,
        istep, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam,
        x0, v0, ω0, udot_init, θdot_init, CtCabPdot_init, CtCabHdot_init, dt)

Populate the jacobian matrix `jacob` with the jacobian of the residual vector
with respect to the state variables.

There are four implementations corresponding to the following analysis types:
 - Static
 - Dynamic - Steady State
 - Dynamic - Initial Step (for initializing time domain simulations)
 - Dynamic - Time Marching

 See "GEBT: A general-purpose nonlinear analysis tool for composite beams" by
 Wenbin Yu and "Geometrically nonlinear analysis of composite beams using
 Wiener-Milenković parameters" by Qi Wang and Wenbin Yu.

# Arguments
 - `jacob`: Jacobian matrix
 - `x`: Vector containing current state variables of the system
 - `assembly`: Assembly of rigidly connected nonlinear beam elements
 - `prescribed_conditions`: Dictionary of prescribed conditions at all time steps
 - `distributed_loads`: Dictionary of distributed loads at all time steps
 - `istep`: Current time step
 - `irow_pt`: Row index of first equilibrium equation for each point
 - `irow_beam`: Row index of first equation for just this beam element
 - `irow_beam1`: Row index of first equation for the left side of each beam
 - `irow_beam2`: Row index of first equation for the right side of each beam
 - `icol_pt`: Column index of first state variable for each point
 - `icol_beam`: Column index of first state variable for each beam element

# Additional Arguments for Dynamic Analyses
 - `x0`: Global frame origin (for the current time step)
 - `v0`: Global frame linear velocity (for the current time step)
 - `ω0`: Global frame angular velocity (for the current time step)

# Additional Arguments for Initial Step Analyses
 - `u`: deflection variables for each beam element
 - `θ`: rotation variables for each beam element
 - `udot`: time derivative of u for each beam element
 - `θdot`: time derivative of θ for each beam element

# Additional Arguments for Time Marching Analyses
 - `udot`: `-2/dt*u - udot` for each beam element from the previous time step
 - `θdot_init`: `-2/dt*θ - θdot` for each beam element from the previous time step
 - `CtCabPdot`: `-2/dt*C'*Cab*P - C'*Cab*Pdot` for each beam element from the previous time step
 - `CtCabHdot`: `-2/dt*C'*Cab*H - C'*Cab*Hdot` for each beam element from the previous time step
 - `dt`: time step size
"""
system_jacobian!

# static
function system_jacobian!(jacob, x, assembly, prescribed_conditions, distributed_loads,
    istep, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam)

    npoint = length(assembly.points)
    nbeam = length(assembly.elements)

    # add contributions to residual equations from the elements
    for ibeam = 1:nbeam

        icol = icol_beam[ibeam]
        irow_b = irow_beam[ibeam]
        irow_b1 = irow_beam1[ibeam]
        irow_p1 = irow_pt[assembly.start[ibeam]]
        irow_b2 = irow_beam2[ibeam]
        irow_p2 = irow_pt[assembly.stop[ibeam]]

        jacob = element_jacobian!(jacob, x, ibeam, assembly.elements[ibeam],
            distributed_loads, istep, icol, irow_b, irow_b1,
            irow_p1, irow_b2, irow_p2)
    end

    # add contributions to the system jacobian matrix from the prescribed point conditions
    for ipoint = 1:npoint

        # skip if the unknowns have been eliminated from the system of equations
        if icol_pt[ipoint] <= 0
            continue
        end

        icol = icol_pt[ipoint]
        irow_p = irow_pt[ipoint]

        jacob = point_jacobian!(jacob, x, ipoint, assembly, prescribed_conditions,
            istep, icol, irow_p, irow_beam1, irow_beam2)

    end

    return jacob
end

# dynamic
function system_jacobian!(jacob, x, assembly, prescribed_conditions, distributed_loads,
    istep, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam,
    x0, v0, ω0)

    npoint = length(assembly.points)
    nbeam = length(assembly.elements)

    # add contributions to residual equations from the elements
    for ibeam = 1:nbeam

        icol = icol_beam[ibeam]
        irow_b = irow_beam[ibeam]
        irow_b1 = irow_beam1[ibeam]
        irow_p1 = irow_pt[assembly.start[ibeam]]
        irow_b2 = irow_beam2[ibeam]
        irow_p2 = irow_pt[assembly.stop[ibeam]]

        jacob = element_jacobian!(jacob, x, ibeam, assembly.elements[ibeam],
            distributed_loads, istep, icol, irow_b, irow_b1,
            irow_p1, irow_b2, irow_p2, x0, v0, ω0)
    end

    # add contributions to the system jacobian matrix from the prescribed point conditions
    for ipoint = 1:npoint

        # skip if the unknowns have been eliminated from the system of equations
        if icol_pt[ipoint] <= 0
            continue
        end

        icol = icol_pt[ipoint]
        irow_p = irow_pt[ipoint]

        jacob = point_jacobian!(jacob, x, ipoint, assembly, prescribed_conditions,
            istep, icol, irow_p, irow_beam1, irow_beam2)
    end

    return jacob
end

# initial step
function system_jacobian!(jacob, x, assembly, prescribed_conditions, distributed_loads,
    istep, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam,
    x0, v0, ω0, u, θ, udot, θdot)

    npoint = length(assembly.points)
    nbeam = length(assembly.elements)

    # add contributions to residual equations from the elements
    for ibeam = 1:nbeam

        icol = icol_beam[ibeam]
        irow_b = irow_beam[ibeam]
        irow_b1 = irow_beam1[ibeam]
        irow_p1 = irow_pt[assembly.start[ibeam]]
        irow_b2 = irow_beam2[ibeam]
        irow_p2 = irow_pt[assembly.stop[ibeam]]

        jacob = element_jacobian!(jacob, x, ibeam, assembly.elements[ibeam],
            distributed_loads, icol, irow_b, irow_b1,
            irow_p1, irow_b2, irow_p2, x0, v0, ω0, u, θ, udot, θdot)
    end

    # add contributions to the system jacobian matrix from the prescribed point conditions
    for ipoint = 1:npoint

        # skip if the unknowns have been eliminated from the system of equations
        if icol_pt[ipoint] <= 0
            continue
        end

        icol = icol_pt[ipoint]
        irow_p = irow_pt[ipoint]

        jacob = point_jacobian!(jacob, x, ipoint, assembly, prescribed_conditions,
            istep, icol, irow_p, irow_beam1, irow_beam2)
    end

    return jacob
end

# time marching
function system_jacobian!(jacob, x, assembly, prescribed_conditions, distributed_loads,
    istep, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam,
    x0, v0, ω0, udot_init, θdot_init, CtCabPdot_init, CtCabHdot_init, dt)

    npoint = length(assembly.points)
    nbeam = length(assembly.elements)

    # add contributions to residual equations from the elements
    for ibeam = 1:nbeam

        icol = icol_beam[ibeam]
        irow_b = irow_beam[ibeam]
        irow_b1 = irow_beam1[ibeam]
        irow_p1 = irow_pt[assembly.start[ibeam]]
        irow_b2 = irow_beam2[ibeam]
        irow_p2 = irow_pt[assembly.stop[ibeam]]

        jacob = element_jacobian!(jacob, x, ibeam, assembly.elements[ibeam],
            distributed_loads, istep, icol, irow_b, irow_b1,
            irow_p1, irow_b2, irow_p2, x0, v0, ω0, udot_init, θdot_init,
            CtCabPdot_init, CtCabHdot_init, dt)
    end

    # add contributions to the system jacobian matrix from the prescribed point conditions
    for ipoint = 1:npoint

        # skip if the unknowns have been eliminated from the system of equations
        if icol_pt[ipoint] <= 0
            continue
        end

        icol = icol_pt[ipoint]
        irow_p = irow_pt[ipoint]

        jacob = point_jacobian!(jacob, x, ipoint, assembly, prescribed_conditions,
            istep, icol, irow_p, irow_beam1, irow_beam2)
    end

    # zero out near-zero values ( < eps() )
    jacob = droptol!(jacob, eps(eltype(jacob)))

    return jacob
end

"""
    system_mass_matrix!(jacob, x, assembly, irow_pt, irow_beam, irow_beam1,
        irow_beam2, icol_pt, icol_beam)

Populate the system "mass matrix", the jacobian of the residual vector
with respect to the time derivatives of the state variables.

See "GEBT: A general-purpose nonlinear analysis tool for composite beams" by
Wenbin Yu and "Geometrically nonlinear analysis of composite beams using
Wiener-Milenković parameters" by Qi Wang and Wenbin Yu.

# Arguments
 - `jacob`: Jacobian matrix
 - `x`: Vector containing current state variables of the system
 - `assembly`: Assembly of rigidly connected nonlinear beam elements
 - `irow_pt`: Row index of first equilibrium equation for each point
 - `irow_beam`: Row index of first equation for just this beam element
 - `irow_beam1`: Row index of first equation for the left side of each beam
 - `irow_beam2`: Row index of first equation for the right side of each beam
 - `icol_pt`: Column index of first state variable for each point
 - `icol_beam`: Column index of first state variable for each beam element
"""
function system_mass_matrix!(jacob, x, assembly, irow_pt, irow_beam, irow_beam1,
    irow_beam2, icol_pt, icol_beam)

    npoint = length(assembly.points)
    nbeam = length(assembly.elements)

    # add contributions to "mass matrix" from the beam elements
    for ibeam = 1:nbeam

        icol = icol_beam[ibeam]
        irow_b = irow_beam[ibeam]
        irow_b1 = irow_beam1[ibeam]
        irow_p1 = irow_pt[assembly.start[ibeam]]
        irow_b2 = irow_beam2[ibeam]
        irow_p2 = irow_pt[assembly.stop[ibeam]]

        element_mass_matrix!(jacob, x, assembly.elements[ibeam], icol, irow_b,
            irow_b1, irow_p1, irow_b2, irow_p2)
    end

    # no contributions to "mass matrix" from point state variables

    # zero out near-zero values ( < eps() )
    jacob = droptol!(jacob, eps(eltype(jacob)))

    return jacob
end
