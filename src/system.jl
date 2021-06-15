"""
    System{TF, TV<:AbstractVector{TF}, TM<:AbstractMatrix{TF}}

Contains the system state, residual vector, and jacobian matrices as well as
pointers to be able to access their contents.  Also contains additional storage
needed for time domain simulations.

# Fields:
 - `static`: Flag indicating whether system matrices are only valid for static analyses
 - `x`: State vector
 - `r`: Residual vector
 - `K`: System jacobian matrix with respect to the state variables
 - `M`: System jacobian matrix with respect to the time derivative of the state variables
 - `force_scaling`: Scaling for state variables corresponding to forces/moments
 - `mass_scaling`: Scaling for state variables corresponding to masses/inertias
 - `irow_pt`: Row index of first equilibrium equation for each point
 - `irow_beam`: Row index of first equation for just this beam element
 - `irow_beam1`: Row index of first equation for the left side of each beam
 - `irow_beam2`: Row index of first equation for the right side of each beam
 - `icol_pt`: Row/Column index of first state variable for each point
 - `icol_beam`: Row/Column index of first state variable for each beam element
 - `udot`: `2/dt*u + udot` for each beam element from the previous time step
 - `θdot`: `2/dt*θ + θdot` for each beam element from the previous time step
 - `Pdot`: `2/dt*C'*Cab*P + C'*Cab*Pdot` for each beam element from the previous time step
 - `Hdot`: `2/dt*C'*Cab*H + C'*Cab*Hdot` for each beam element from the previous time step
 - `t`: Current system time
"""
mutable struct System{TF, TV<:AbstractVector{TF}, TM<:AbstractMatrix{TF}}
    static::Bool
    x::TV
    r::TV
    K::TM
    M::TM
    force_scaling::TF
    mass_scaling::TF
    irow_pt::Vector{Int}
    irow_beam::Vector{Int}
    irow_beam1::Vector{Int}
    irow_beam2::Vector{Int}
    icol_pt::Vector{Int}
    icol_beam::Vector{Int}
    udot::Vector{SVector{3,TF}}
    θdot::Vector{SVector{3,TF}}
    Pdot::Vector{SVector{3,TF}}
    Hdot::Vector{SVector{3,TF}}
    t::TF
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

function System(TF, assembly, points, static)

    # get number of beams
    nbeam = length(assembly.elements)

    # get the number of beams connected to each point.
    n_connections = point_connections(assembly)

    # get the size and pointers into the system matrices
    n, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam =
        system_indices(assembly, points, n_connections, static)

    # set the force scaling based on the average compliance matrix value
    compliance_entries = vcat([vcat(beam.C11..., beam.C12..., beam.C22...) for beam in assembly.elements]...)
    compliance_nonzero_indices = findall(xi -> abs(xi) > eps(TF), compliance_entries)
    if isempty(compliance_nonzero_indices)
        force_scaling = 1.0
    else
        nonzero_compliance_entries = compliance_entries[compliance_nonzero_indices]
        force_scaling = nextpow(2.0, length(nonzero_compliance_entries)/sum(nonzero_compliance_entries)/100)
    end

    # set the mass scaling based on the average inverse mass matrix value
    minv_entries = vcat([vcat(beam.minv11..., beam.minv12..., beam.minv22...) for beam in assembly.elements]...)
    minv_nonzero_indices = findall(xi -> abs(xi) > eps(TF), minv_entries)
    if isempty(minv_nonzero_indices)
        mass_scaling = 1.0
    else
        nonzero_minv_entries = minv_entries[minv_nonzero_indices]
        mass_scaling = nextpow(2.0, length(nonzero_minv_entries)/sum(nonzero_minv_entries))
    end

    # initialize system matrices
    x = zeros(TF, n)
    r = zeros(TF, n)
    K = spzeros(TF, n, n)
    M = spzeros(TF, n, n)

    # initialize storage for time domain simulations
    udot = [@SVector zeros(TF, 3) for i = 1:nbeam]
    θdot = [@SVector zeros(TF, 3) for i = 1:nbeam]
    Pdot = [@SVector zeros(TF, 3) for i = 1:nbeam]
    Hdot = [@SVector zeros(TF, 3) for i = 1:nbeam]

    # initialize current time
    t = 0.0

    # set system types
    TV = promote_type(typeof(x), typeof(r))
    TM = promote_type(typeof(K), typeof(M))

    return System{TF, TV, TM}(static, x, r, K, M, force_scaling, mass_scaling,
        irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam, udot, θdot,
        Pdot, Hdot, t)
end

"""
   system_state(system)

Return a vector containing the state variables of `system`.
"""
system_state(system) = system.x

"""
   set_state!([x,] system, prescribed_conditions; kwargs...)

Set the state variables in `system` (or in the vector `x`) to the provided values.
If values are not provided for a given keyword argument, then the state variables
corresponding to the keyword argument are not updated.

# Keyword Arguments
 - `u_b`: Vector containing initial displacements (in the global frame) for each beam element.
 - `theta_b`: Vector containing rotation variables (in the global frame) for each beam element.
 - `F_b`: Vector containing resultant forces (in the beam coordinate frame) for each beam element.
 - `M_b`: Vector containing resultant moments (in the beam coordinate frame) for each beam element.
 - `P_b`: Vector containing linear momenta (in the beam coordinate frame) for each beam element.
 - `H_b`: Vector containing angular momenta (in the beam coordinate frame) for each beam element.
 - `u_p`: Vector containing initial displacements (in the global frame) for each point.
 - `theta_p`: Vector containing rotation variables (in the global frame) for each point.
 - `F_p`: Vector containing externally applied forces (in the global frame) for each point.
 - `M_p`: Vector containing externally applied moments (in the global frame) for each point.
"""
set_state!

function set_state!(system; kwargs...)
    x = set_state!(system.x, system; kwargs...)
    return system
end

function set_state!(x, system; u_b = nothing, theta_b = nothing,
    F_b = nothing, M_b = nothing, P_b = nothing, H_b = nothing,
    u_p = nothing, theta_p = nothing, F_p = nothing, M_p = nothing)

    nelem = length(system.icol_beam)
    for ielem = 1:nelem
        if !isnothing(u_b)
            set_element_deflections!(x, system, u_b, ielem)
        end
        if !isnothing(theta_b)
            set_element_rotations!(x, system, theta_b, ielem)
        end
        if !isnothing(F_b)
            set_element_forces!(x, system, F_b, ielem)
        end
        if !isnothing(M_b)
            set_element_moments!(x, system, M_b, ielem)
        end
        if !isnothing(P_b)
            set_element_linear_momenta!(x, system, P_b, ielem)
        end
        if !isnothing(H_b)
            set_element_angular_momenta!(x, system, H_b, ielem)
        end
    end

    npoint = length(system.icol_point)
    for ipoint = 1:npoint
        if !isnothing(u_p)
            set_point_deflections!(x, system, u_p, ipoint, prescribed_conditions)
        end
        if !isnothing(theta_p)
            set_point_rotations!(x, system, theta_p, ipoint, prescribed_conditions)
        end
        if !isnothing(F_p)
            set_point_forces!(x, system, F_p, ipoint, prescribed_conditions)
        end
        if !isnothing(M_p)
            set_point_moments!(x, system, M_p, ipoint, prescribed_conditions)
        end
    end

    return x
end

"""
   set_element_deflections!([x,] system, u_b, ielem)

Set the state variables in `system` (or in the vector `x`) corresponding to the
deflections of element `ielem` to the provided values.
"""
set_element_deflections!

function set_element_deflections!(x, system, u_b, ielem)
    icol = icol_beam[ielem]
    x[icol:icol+2] .= u_b
    return x
end

function set_element_deflections!(system, u_b, ielem)
    x = set_element_deflections!(system.x, system, u_b, ielem)
    return system
end

"""
   set_element_rotations!([x,] system, θ_b, ielem)

Set the state variables in `system` (or in the vector `x`) corresponding to the
rotations of element `ielem` to the provided values.
"""
set_element_rotations!

function set_element_rotations!(x, system, theta_b, ielem)
    icol = icol_beam[ielem]
    x[icol+3:icol+5] .= theta_b
    return x
end

function set_element_rotations!(system, theta_b, ielem)
    x = set_element_rotations!(system.x, system, theta_b, ielem)
    return system
end

"""
   set_element_forces!([x,] system, F_b, ielem)

Set the state variables in `system` (or in the vector `x`) corresponding to the
resultant forces of element `ielem` to the provided values.
"""
set_element_forces!

function set_element_forces!(x, system, F_b, ielem)
    icol = icol_beam[ielem]
    x[icol+6:icol+8] .= F_b ./ system.force_scaling
    return x
end

function set_element_forces!(system, F_b, ielem)
    x = set_element_forces!(system.x, system, F_b, ielem)
    return system
end

"""
   set_element_moments!([x,] system, u_b, ielem)

Set the state variables in `system` (or in the vector `x`) corresponding to the
resultant moments of element `ielem` to the provided values.
"""
set_element_moments!

function set_element_moments!(x, system, M_b, ielem)
    icol = icol_beam[ielem]
    x[icol+9:icol+11] .= M_b ./ system.force_scaling
    return x
end

function set_element_moments!(system, M_b, ielem)
    x = set_element_moments!(system.x, system, M_b, ielem)
    return system
end

"""
   set_element_linear_momenta!([x,] system, u_b, ielem)

Set the state variables in `system` (or in the vector `x`) corresponding to the
linear_momenta of element `ielem` to the provided values.
"""
set_element_linear_momenta!

function set_element_linear_momenta!(x, system, P_b, ielem)
    @assert !system.static
    icol = icol_beam[ielem]
    x[icol+12:icol+14] .= P_b ./ system.mass_scaling
    return x
end

function set_element_linear_momenta!(system, P_b, ielem)
    x = set_element_linear_momenta!(system.x, system, P_b, ielem)
    return system
end

"""
   set_element_angular_momenta!([x,] system, u_b, ielem)

Set the state variables in `system` (or in the vector `x`) corresponding to the
angular_momenta of element `ielem` to the provided values.
"""
set_element_angular_momenta!

function set_element_angular_momenta!(x, system, H_b, ielem)
    @assert !system.static
    icol = icol_beam[ielem]
    x[icol+15:icol+17] .= H_b ./ system.mass_scaling
    return x
end

function set_element_angular_momenta!(system, H_b, ielem)
    x = set_element_angular_momenta!(system.x, system, H_b, ielem)
    return system
end

"""
   set_point_deflections!([x,] system, u_b, ipoint, prescribed_conditions)

Set the state variables in `system` (or in the vector `x`) corresponding to the
deflections (or externally applied forces) at point `ipoint` to the provided values.
"""
set_point_deflections!

function set_point_deflections!(x, system, u_p, ipoint, prescribed_conditions)
    icol = icol_point[ipoint]
    if icol_pt > 0
        # point variables are state variables
        for k = 1:3
            if ipoint in keys(prescribed_conditions)
                # prescribed conditions exist for this point
                if prescribed_conditions[ipoint].force_dof[k]
                    # forces are prescribed
                    x[icol+k-1] .= u_p[k]
                end
            else
                # no prescribed condition for this point
                x[icol+k-1] .= u_p[k]
            end
        end
    end
    return x
end

function set_point_deflections!(system, u_p, ipoint, prescribed_conditions)
    x = set_point_deflections!(system.x, system, u_p, ipoint, prescribed_conditions)
    return system
end

"""
   set_point_rotations!([x,] system, u_b, ipoint, prescribed_conditions)

Set the state variables in `system` (or in the vector `x`) corresponding to the
rotations (or externally applied moments) at point `ipoint` to the provided values.
"""
set_point_rotations!

function set_point_rotations!(x, system, theta_p, ipoint, prescribed_conditions)
    icol = icol_point[ipoint]
    if icol_pt > 0
        # point variables are state variables
        for k = 1:3
            if ipoint in keys(prescribed_conditions)
                # prescribed conditions exist for this point
                if prescribed_conditions[ipoint].force_dof[3+k]
                    # forces are prescribed
                    x[icol+k+2] .= theta_p[k]
                end
            else
                # no prescribed condition for this point
                x[icol+k+2] .= theta_p[k]
            end
        end
    end
    return x
end

function set_point_rotations!(system, theta_p, ipoint, prescribed_conditions)
    x = set_point_rotations!(system.x, system, theta_p, ipoint, prescribed_conditions)
    return system
end

"""
   set_point_forces!([x,] system, F_p, ipoint, prescribed_conditions)

Set the state variables in `system` (or in the vector `x`) corresponding to the
external forces applied on point `ipoint` to the provided values.
"""
set_point_forces!

function set_point_forces!(x, system, F_p, ipoint, prescribed_conditions)
    icol = icol_point[ipoint]
    if icol_pt > 0
        # point variables are state variables
        for k = 1:3
            if ipoint in keys(prescribed_conditions)
                # prescribed conditions exist for this point
                if !prescribed_conditions[ipoint].force_dof[k]
                    # forces are prescribed
                    x[icol+k-1] .= F_p[k] ./ system.force_scaling
                end
            end
        end
    end
    return x
end

function set_point_forces!(system, F_p, ipoint, prescribed_conditions)
    x = set_point_forces!(system.x, system, F_p, ipoint, prescribed_conditions)
    return system
end

"""
   set_point_moments!([x,] system, M_p, ipoint, prescribed_conditions)

Set the state variables in `system` (or in the vector `x`) corresponding to the
external moments applied on point `ipoint` to the provided values.
"""
set_point_moments!

function set_point_moments!(x, system, M_p, ipoint, prescribed_conditions)
    icol = icol_point[ipoint]
    if icol_pt > 0
        # point variables are state variables
        for k = 1:3
            if ipoint in keys(prescribed_conditions)
                # prescribed conditions exist for this point
                if !prescribed_conditions[ipoint].force_dof[3+k]
                    # forces are prescribed
                    x[icol+k+2] .= M_p[k] ./ system.force_scaling
                end
            end
        end
    end
    return x
end

function set_point_moments!(system, M_p, ipoint, prescribed_conditions)
    x = set_point_moments!(system.x, system, M_p, ipoint, prescribed_conditions)
    return system
end

"""
   reset_state!(system)

Reset the state variables in `system` (stored in `system.x`).
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
    get_sparsity(system, assembly)

Return a matrix indicating the sparsity structure of the jacobian matrix.
"""
function get_sparsity(system, assembly)

    n = length(system.x)
    irow_pt = system.irow_pt
    irow_beam = system.irow_beam
    irow_beam1 = system.irow_beam1
    irow_beam2 = system.irow_beam2
    icol_pt = system.icol_pt
    icol_beam = system.icol_beam

    sparsity = spzeros(Bool, n, n)

    # --- add constributions from beam element state variables --- #

    nbeam = length(icol_beam)

    for ibeam = 1:nbeam
        # indices for this beam
        icol = icol_beam[ibeam]
        irow_b = irow_beam[ibeam]
        irow_b1 = irow_beam1[ibeam]
        irow_p1 = irow_pt[assembly.start[ibeam]]
        irow_b2 = irow_beam2[ibeam]
        irow_p2 = irow_pt[assembly.stop[ibeam]]

        # --- static, beam element jacobian contributions --- #

        # left endpoint equilibrium equations
        sparsity[irow_p1:irow_p1+2, icol:icol+8] .= true
        sparsity[irow_p1+3:irow_p1+5, icol+3:icol+11] .= true
        # left point compatability equations
        irow = ifelse(irow_b1 == irow_p1 || irow_b1 <= 0, irow_p1+6, irow_b1)
        sparsity[irow:irow+2, icol:icol+11] .= true
        sparsity[irow+3:irow+5, icol+3:icol+11] .= true
        # right endpoint equilibrium equations
        sparsity[irow_p2:irow_p2+2, icol:icol+8] .= true
        sparsity[irow_p2+3:irow_p2+5, icol+3:icol+11] .= true
        # right point compatability equations
        irow = ifelse(irow_b2 == irow_p2 || irow_b2 <= 0, irow_p2 + 6, irow_b2)
        sparsity[irow:irow+2, icol:icol+11] .= true
        sparsity[irow+3:irow+5, icol+3:icol+11] .= true

        # --- dynamic, beam element jacobian contributions --- #
        if !system.static
            # left endpoint equilibrium equations
            sparsity[irow_p1:irow_p1+2, icol+12:icol+14] .= true
            sparsity[irow_p1+3:irow_p1+5, icol+12:icol+17] .= true
            # right endpoint equilibrium equations
            sparsity[irow_p2:irow_p2+2, icol+12:icol+14] .= true
            sparsity[irow_p2+3:irow_p2+5, icol+12:icol+17] .= true
            # element residual equations
            sparsity[irow_b:irow_b+2, icol:icol+5] .= true
            sparsity[irow_b:irow_b+2, icol+12:icol+17] .= true
            sparsity[irow_b+3:irow_b+5, icol+3:icol+5] .= true
            sparsity[irow_b+3:irow_b+5, icol+12:icol+17] .= true
        end
    end

    # --- add constributions from point state variables --- #

    npoint = length(icol_pt)

    for ipoint = 1:npoint

        # skip if the unknowns have been eliminated from the system of equations
        if icol_pt[ipoint] <= 0
            continue
        end

        icol = icol_pt[ipoint]
        irow_p = irow_pt[ipoint]

        for ibeam = 1:nbeam
            # check left side of beam
            if ipoint == assembly.start[ibeam]
                # add jacobian entries for the beam endpoint
                irow_b = irow_beam1[ibeam]
                if irow_b == irow_p
                    for i = 1:6
                        for j = 4:6
                            sparsity[irow_b+i-1, icol+j-1] = true
                        end
                        sparsity[irow_b+i+5, icol+i-1] = true
                        sparsity[irow_b+i-1, icol+i-1] = true
                    end
                else
                    for i = 1:6
                        sparsity[irow_b+i-1, icol+i-1] = true
                    end
                end
            end
            # check right side of beam
            if ipoint == assembly.stop[ibeam]
                # add jacobian entries for the beam endpoint
                irow_b = irow_beam2[ibeam]
                if irow_b == irow_p
                    for i = 1:6
                        for j = 4:6
                            sparsity[irow_b+i-1, icol+j-1] = true
                        end
                        sparsity[irow_b+i+5, icol+i-1] = true
                        sparsity[irow_b+i-1, icol+i-1] = true
                    end
                else
                    for i = 1:6
                        sparsity[irow_b+i-1, icol+i-1] = true
                    end
                end
            end
        end
    end

    return sparsity
end

"""
    system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads,
        force_scaling, mass_scaling, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam)
    system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads,
        force_scaling, mass_scaling, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam,
        x0, v0, ω0)
    system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads,
        force_scaling, mass_scaling, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam,
        x0, v0, ω0, u, θ, udot, θdot)
    system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads,
        force_scaling, mass_scaling, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam,
        x0, v0, ω0, udot_init, θdot_init, CtCabPdot_init, CtCabHdot_init, dt)

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
 - `force_scaling`: Scaling parameter for forces/moments
 - `mass_scaling`: Scaling parameter for masses/inertias
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
 - `udot_init`: `-2/dt*u - udot` for each beam element from the previous time step
 - `θdot_init`: `-2/dt*θ - θdot` for each beam element from the previous time step
 - `CtCabPdot_init`: `-2/dt*C'*Cab*P - C'*Cab*Pdot` for each beam element from the previous time step
 - `CtCabHdot_init`: `-2/dt*C'*Cab*H - C'*Cab*Hdot` for each beam element from the previous time step
 - `dt`: Time step size. If set to `nothing`, `udot_init`,
"""
system_residual!

# static
function system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads,
    force_scaling, mass_scaling, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam)

    return static_system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads,
        force_scaling, mass_scaling, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam)
end

# dynamic - steady state
function system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads,
    force_scaling, mass_scaling, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam,
    x0, v0, ω0)

    return steady_state_system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads,
        force_scaling, mass_scaling, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam,
        x0, v0, ω0)
end

# dynamic - initial step
function system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads,
    force_scaling, mass_scaling, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam,
    x0, v0, ω0, u, θ, udot, θdot)

    return initial_step_system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads,
    force_scaling, mass_scaling, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam,
    x0, v0, ω0, u, θ, udot, θdot)
end

# dynamic - newmark scheme time marching
function system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads,
    force_scaling, mass_scaling, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam,
    x0, v0, ω0, udot_init, θdot_init, CtCabPdot_init, CtCabHdot_init, dt)

    # dynamic - newmark scheme time marching
    return newmark_system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads,
        force_scaling, mass_scaling, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam,
        x0, v0, ω0, udot_init, θdot_init, CtCabPdot_init, CtCabHdot_init, dt)
end

# static
function static_system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads,
    force_scaling, mass_scaling, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam)

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

        resid = static_element_residual!(resid, x, ibeam, assembly.elements[ibeam],
             distributed_loads, force_scaling, mass_scaling, icol, irow_b, irow_b1,
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
            force_scaling, icol, irow_p, irow_beam1, irow_beam2)
    end

    return resid
end

# dynamic - steady state
function steady_state_system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads,
    force_scaling, mass_scaling, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam,
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

        resid = steady_state_element_residual!(resid, x, ibeam, assembly.elements[ibeam],
             distributed_loads, force_scaling, mass_scaling, icol, irow_b, irow_b1,
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
            force_scaling, icol, irow_p, irow_beam1, irow_beam2)
    end

    return resid
end

# dynamic - initial step
function initial_step_system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads,
    force_scaling, mass_scaling, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam,
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

        resid = initial_step_element_residual!(resid, x, ibeam, assembly.elements[ibeam],
             distributed_loads, force_scaling, mass_scaling, icol, irow_b, irow_b1,
            irow_p1, irow_b2, irow_p2, x0, v0, ω0,
            u[ibeam], θ[ibeam], udot[ibeam], θdot[ibeam])
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
            force_scaling, icol, irow_p, irow_beam1, irow_beam2)
    end

    return resid
end

# dynamic - newmark scheme time marching
function newmark_system_residual!(resid, x, assembly, prescribed_conditions, distributed_loads,
    force_scaling, mass_scaling, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam,
    x0, v0, ω0, udot_init, θdot_init, CtCabPdot_init, CtCabHdot_init, dt)

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

        resid = newmark_element_residual!(resid, x, ibeam, assembly.elements[ibeam],
             distributed_loads, force_scaling, mass_scaling, icol, irow_b, irow_b1, irow_p1, irow_b2, irow_p2,
             x0, v0, ω0,
             udot_init[ibeam], θdot_init[ibeam],
             CtCabPdot_init[ibeam], CtCabHdot_init[ibeam], dt)
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
            force_scaling, icol, irow_p, irow_beam1, irow_beam2)
    end

    return resid
end

# dynamic - general
function dynamic_system_residual!(resid, x, dx, assembly, prescribed_conditions, distributed_loads,
    force_scaling, mass_scaling, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam,
    x0, v0, ω0)

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

        # set state rates for element
        udot = SVector(dx[icol], dx[icol+1], dx[icol+2])
        θdot = SVector(dx[icol+3], dx[icol+4], dx[icol+5])
        Pdot = SVector(dx[icol+12], dx[icol+13], dx[icol+14]) .* mass_scaling
        Hdot = SVector(dx[icol+15], dx[icol+16], dx[icol+17]) .* mass_scaling

        resid = dynamic_element_residual!(resid, x, ibeam, assembly.elements[ibeam],
             distributed_loads, force_scaling, mass_scaling, icol, irow_b, irow_b1, irow_p1, irow_b2, irow_p2,
             x0, v0, ω0, udot, θdot, Pdot, Hdot)

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
            force_scaling, icol, irow_p, irow_beam1, irow_beam2)
    end

    return resid
end

"""
    system_jacobian!(jacob, x, assembly, prescribed_conditions, distributed_loads,
        force_scaling, mass_scaling, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam)
    system_jacobian!(jacob, x, assembly, prescribed_conditions, distributed_loads,
        force_scaling, mass_scaling, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam,
        x0, v0, ω0)
    system_jacobian!(jacob, x, assembly, prescribed_conditions, distributed_loads,
        force_scaling, mass_scaling, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam,
        x0, v0, ω0, u, θ, udot, θdot)
    system_jacobian!(jacob, x, assembly, prescribed_conditions, distributed_loads,
        force_scaling, mass_scaling, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam,
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
 - `force_scaling`: Scaling parameter for forces/moments
 - `mass_scaling`: Scaling parameter for masses/inertia
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
@inline function system_jacobian!(jacob, x, assembly, prescribed_conditions, distributed_loads,
    force_scaling, mass_scaling, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam)

    return static_system_jacobian!(jacob, x, assembly, prescribed_conditions, distributed_loads,
        force_scaling, mass_scaling, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam)
end

# dynamic - steady state
@inline function system_jacobian!(jacob, x, assembly, prescribed_conditions, distributed_loads,
    force_scaling, mass_scaling, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam,
    x0, v0, ω0)

    return steady_state_system_jacobian!(jacob, x, assembly, prescribed_conditions, distributed_loads,
        force_scaling, mass_scaling, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam,
        x0, v0, ω0)
end

# dynamic - initial step
@inline function system_jacobian!(jacob, x, assembly, prescribed_conditions, distributed_loads,
    force_scaling, mass_scaling, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam,
    x0, v0, ω0, u, θ, udot, θdot)

    return initial_step_system_jacobian!(jacob, x, assembly, prescribed_conditions, distributed_loads,
        force_scaling, mass_scaling, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam,
        x0, v0, ω0, u, θ, udot, θdot)
end

# dynamic - newmark scheme time marching
@inline function system_jacobian!(jacob, x, assembly, prescribed_conditions,
    distributed_loads, force_scaling, mass_scaling,
    irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam,
    x0, v0, ω0, udot_init, θdot_init, CtCabPdot_init, CtCabHdot_init, dt)

    return newmark_system_jacobian!(jacob, x, assembly, prescribed_conditions,
        distributed_loads, force_scaling, mass_scaling,
        irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam,
        x0, v0, ω0, udot_init, θdot_init, CtCabPdot_init, CtCabHdot_init, dt)
end

# static
@inline function static_system_jacobian!(jacob, x, assembly,
    prescribed_conditions, distributed_loads, force_scaling, mass_scaling,
    irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam)

    jacob .= 0

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
            distributed_loads, force_scaling, mass_scaling, icol, irow_b, irow_b1,
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
            force_scaling, icol, irow_p, irow_beam1, irow_beam2)

    end

    return jacob
end

# dynamic - steady state
@inline function steady_state_system_jacobian!(jacob, x, assembly,
    prescribed_conditions, distributed_loads, force_scaling, mass_scaling,
    irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam,
    x0, v0, ω0)

    jacob .= 0

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
            distributed_loads, force_scaling, mass_scaling, icol, irow_b, irow_b1,
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
            force_scaling, icol, irow_p, irow_beam1, irow_beam2)
    end

    return jacob
end

# dynamic - initial step
@inline function initial_step_system_jacobian!(jacob, x, assembly,
    prescribed_conditions, distributed_loads, force_scaling, mass_scaling,
    irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam,
    x0, v0, ω0, u, θ, udot, θdot)

    jacob .= 0

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
            distributed_loads, force_scaling, mass_scaling, icol, irow_b, irow_b1,
            irow_p1, irow_b2, irow_p2, x0, v0, ω0,
            u[ibeam], θ[ibeam], udot[ibeam], θdot[ibeam])
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
            force_scaling, icol, irow_p, irow_beam1, irow_beam2)
    end

    return jacob
end

# dynamic - newmark scheme time marching
@inline function newmark_system_jacobian!(jacob, x, assembly,
    prescribed_conditions, distributed_loads, force_scaling, mass_scaling,
    irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam,
    x0, v0, ω0, udot_init, θdot_init, CtCabPdot_init, CtCabHdot_init, dt)

    jacob .= 0

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
            distributed_loads, force_scaling, mass_scaling, icol, irow_b, irow_b1,
            irow_p1, irow_b2, irow_p2, x0, v0, ω0,
            udot_init[ibeam], θdot_init[ibeam],
            CtCabPdot_init[ibeam], CtCabHdot_init[ibeam], dt)
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
            force_scaling, icol, irow_p, irow_beam1, irow_beam2)
    end

    # # zero out near-zero values ( < eps() )
    # jacob = droptol!(jacob, eps(eltype(jacob)))

    return jacob
end

# dynamic - general
@inline function dynamic_system_jacobian!(jacob, x, dx, assembly,
    prescribed_conditions, distributed_loads, force_scaling, mass_scaling,
    irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam,
    x0, v0, ω0)

    jacob .= 0

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

        # set state rates for element
        udot = SVector(dx[icol], dx[icol+1], dx[icol+2])
        θdot = SVector(dx[icol+3], dx[icol+4], dx[icol+5])
        Pdot = SVector(dx[icol+12], dx[icol+13], dx[icol+14]) .* mass_scaling
        Hdot = SVector(dx[icol+15], dx[icol+16], dx[icol+17]) .* mass_scaling

        jacob = dynamic_element_jacobian!(jacob, x, ibeam, assembly.elements[ibeam],
            distributed_loads, force_scaling, mass_scaling, icol, irow_b, irow_b1,
            irow_p1, irow_b2, irow_p2, x0, v0, ω0,
            udot, θdot, Pdot, Hdot)
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
            force_scaling, icol, irow_p, irow_beam1, irow_beam2)
    end

    # # zero out near-zero values ( < eps() )
    # jacob = droptol!(jacob, eps(eltype(jacob)))

    return jacob
end

"""
    system_mass_matrix!(jacob, x, assembly, force_scaling, mass_scaling, irow_pt, irow_beam,
        irow_beam1, irow_beam2, icol_pt, icol_beam)

Populate the system "mass matrix", the jacobian of the residual vector
with respect to the time derivatives of the state variables.

See "GEBT: A general-purpose nonlinear analysis tool for composite beams" by
Wenbin Yu and "Geometrically nonlinear analysis of composite beams using
Wiener-Milenković parameters" by Qi Wang and Wenbin Yu.

# Arguments
 - `jacob`: Jacobian matrix
 - `x`: Vector containing current state variables of the system
 - `assembly`: Assembly of rigidly connected nonlinear beam elements
 - `force_scaling`: Scaling parameter for forces
 - `mass_scaling`: Scaling parameter for masses
 - `irow_pt`: Row index of first equilibrium equation for each point
 - `irow_beam`: Row index of first equation for just this beam element
 - `irow_beam1`: Row index of first equation for the left side of each beam
 - `irow_beam2`: Row index of first equation for the right side of each beam
 - `icol_pt`: Column index of first state variable for each point
 - `icol_beam`: Column index of first state variable for each beam element
"""
function system_mass_matrix!(jacob, x, assembly, force_scaling, mass_scaling, irow_pt,
    irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam)

    jacob .= 0

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

        element_mass_matrix!(jacob, x, assembly.elements[ibeam], force_scaling,
            mass_scaling, icol, irow_b, irow_b1, irow_p1, irow_b2, irow_p2)
    end

    # no contributions to "mass matrix" from point state variables

    # # zero out near-zero values ( < eps() )
    # jacob = droptol!(jacob, eps(eltype(jacob)))

    return jacob
end

"""
    system_mass_matrix!(jacob, gamma, x, dx, assembly, force_scaling, mass_scaling, irow_pt,
        irow_beam, irow_beam1, irow_beam2, icol_pt, icol_beam)

Add the system mass matrix to `jacob`, scaled by the scaling parameter `gamma`.
"""
function system_mass_matrix!(jacob, gamma, x, assembly, force_scaling, mass_scaling, irow_pt, irow_beam,
    irow_beam1, irow_beam2, icol_pt, icol_beam)

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

        element_mass_matrix!(jacob, gamma, x, assembly.elements[ibeam],
            force_scaling, mass_scaling, icol, irow_b, irow_b1, irow_p1, irow_b2, irow_p2)
    end

    # no contributions to "mass matrix" from point state variables

    return jacob
end
