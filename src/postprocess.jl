"""
    BodyState{TF}

Holds the state variables for the motion of the body frame

# Fields:
 - `u`: Linear deflection
 - `theta`: Angular deflection (Wiener-Milenkovic parameters)
 - `v`: Linear velocity
 - `omega`: Angular velocity
 - `a`: Linear acceleration
 - `alpha`: Angular acceleration
"""
struct BodyState{TF}
    u::SVector{3, TF}
    theta::SVector{3, TF}
    v::SVector{3, TF}
    omega::SVector{3, TF}
    a::SVector{3, TF}
    alpha::SVector{3, TF}
end

"""
    PointState

Holds the state variables for a point

# Fields:
 - `u`: Linear deflection
 - `theta`: Angular deflection (Wiener-Milenkovic parameters)
 - `V`: Linear velocity
 - `Omega`: Angular velocity
 - `F`: External forces
 - `M`: External moments
"""
struct PointState{TF}
    u::SVector{3, TF}
    theta::SVector{3, TF}
    V::SVector{3, TF}
    Omega::SVector{3, TF}
    F::SVector{3, TF}
    M::SVector{3, TF}
end

"""
    ElementState

Holds the state variables for an element

# Fields:
 - `u`: Linear deflection
 - `theta`: Angular deflection (Wiener-Milenkovic parameters)
 - `V`: Linear velocity
 - `Omega`: Angular velocity
 - `Fi`: Internal forces
 - `Mi`: Internal moments
"""
struct ElementState{TF}
    u::SVector{3, TF}
    theta::SVector{3, TF}
    V::SVector{3, TF}
    Omega::SVector{3, TF}
    Fi::SVector{3, TF}
    Mi::SVector{3, TF}
end

"""
    AssemblyState{TF, TB<:BodyState{TF}, TP<:AbstractVector{PointState{TF}},
        TE<:AbstractVector{ElementState{TF}}}

Struct for storing state variables for the points and elements in an assembly.

# Fields:
 - `body::TB`: Object of type [`BodyState`](@ref)
 - `points::TP`: Array of [`PointState`](@ref) for each point in the assembly
 - `elements::TE`: Array of [`ElementState`](@ref) for each element in the assembly
"""
struct AssemblyState{TF, TB<:BodyState{TF}, TP<:AbstractVector{PointState{TF}}, TE<:AbstractVector{ElementState{TF}}}
    body::TB   
    points::TP
    elements::TE
end
Base.eltype(::AssemblyState{TF, TB, TP, TE}) where {TF, TB, TP, TE} = TF

"""
    AssemblyState(system, assembly, x = system.x; prescribed_conditions = Dict())

Post-process the system state given the solution vector `x`.  Return an object
of type `AssemblyState` that defines the state of the assembly for the time step.

If `prescribed_conditions` is not provided, all point state variables are assumed
to be displacements/rotations, rather than their actual identities as used in the
analysis.
"""
function AssemblyState(system, assembly, x = system.x; 
    prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}(),
    kwargs...)

    body = extract_body_state(system, x)

    points = extract_point_states(system, assembly, x; prescribed_conditions)

    elements = extract_element_states(system, assembly, x; prescribed_conditions)

    return AssemblyState(body, points, elements)
end

"""
    extract_body_state(system, x = system.x)

Return the state variables corresponding to the motion of the body frame (see 
[`BodyState`](@ref)) given the solution vector `x`.
"""
extract_body_state

function extract_body_state(system, x = system.x)

    # displacement and velocity state variables
    u, θ = body_frame_displacement(x)
    v, ω = body_frame_velocity(x)
    a, α = body_frame_acceleration(x)

    # convert rotation parameter to Wiener-Milenkovic parameters
    scaling = rotation_parameter_scaling(θ)
    θ *= scaling

    # promote all variables to the same type
    u, θ, v, ω, a, α = promote(u, θ, v, ω, a, α)

    return BodyState(u, θ, v, ω, a, α)
end

"""
    extract_point_state(system, assembly, ipoint, x = system.x;
        prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}())

Return the state variables corresponding to point `ipoint` (see [`PointState`](@ref))
given the solution vector `x`.

If `prescribed_conditions` is not provided, all point state variables are assumed
to be displacements/rotations, rather than their actual identities as used in the
analysis.
"""
extract_point_state

function extract_point_state(system::StaticSystem, assembly, ipoint, x = system.x;
    prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}())

    @unpack force_scaling, indices, t = system

    # check state variable length
    @assert length(x) == length(system.x) "state vector length does not match with the provided system"

    # get current prescribed conditions
    pc = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(t)

    # extract point state variables
    u, θ = point_displacement(x, ipoint, indices.icol_point, pc)
    F, M = point_loads(x, ipoint, indices.icol_point, force_scaling, pc)
    V, Ω = zero(u), zero(θ)
    
    # convert rotation parameter to Wiener-Milenkovic parameters
    scaling = rotation_parameter_scaling(θ)
    θ *= scaling

    # promote all variables to the same type
    u, θ, V, Ω, F, M = promote(u, θ, V, Ω, F, M)

    return PointState(u, θ, V, Ω, F, M)
end

function extract_point_state(system::DynamicSystem, assembly, ipoint, x = system.x;
    prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}())

    @unpack force_scaling, indices, t = system

    # check state variable length
    @assert length(x) == length(system.x) "state vector length does not match with the provided system"

    # get current prescribed conditions
    pc = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(t)

    # extract point state variables
    u, θ = point_displacement(x, ipoint, indices.icol_point, pc)
    F, M = point_loads(x, ipoint, indices.icol_point, force_scaling, pc)
    V, Ω = point_velocities(x, ipoint, indices.icol_point)

    # convert rotation parameter to Wiener-Milenkovic parameters
    scaling = rotation_parameter_scaling(θ)
    θ *= scaling

    # promote all variables to the same type
    u, θ, V, Ω, F, M = promote(u, θ, V, Ω, F, M)

    return PointState(u, θ, V, Ω, F, M)
end

function extract_point_state(system::ExpandedSystem, assembly, ipoint, x = system.x;
    prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}())

    @unpack force_scaling, indices, t = system

    # check state variable length
    @assert length(x) == length(system.x) "state vector length does not match with the provided system"

    # get current prescribed conditions
    pc = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(t)

    # extract point state variables
    u, θ = point_displacement(x, ipoint, indices.icol_point, pc)
    F, M = point_loads(x, ipoint, indices.icol_point, force_scaling, pc)
    V, Ω = point_velocities(x, ipoint, indices.icol_point)

    # rotate linear and angular velocities into the body frame
    C = get_C(θ)
    V = C'*V
    Ω = C'*Ω

    # convert rotation parameters to Wiener-Milenkovic parameters
    scaling = rotation_parameter_scaling(θ)
    θ *= scaling

    # promote all variables to the same type
    u, θ, V, Ω, F, M = promote(u, θ, V, Ω, F, M)

    return PointState(u, θ, V, Ω, F, M)
end

"""
    extract_point_states(system, assembly, x = system.x;
        prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}())

Return the state variables corresponding to each point (see [`PointState`](@ref))
given the solution vector `x`.

If `prescribed_conditions` is not provided, all point state variables are assumed
to be displacements/rotations, rather than their actual identities as used in the
analysis.
"""
function extract_point_states(system, assembly, x = system.x; kwargs...)
    TF = promote_type(eltype(system), eltype(x))
    points = Vector{PointState{TF}}(undef, length(assembly.points))
    return extract_point_states!(points, system, assembly, x; kwargs...)
end

"""
    extract_point_states!(points, system, assembly, x = system.x;
        prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}())

Pre-allocated version of [`extract_point_states`](@ref)
"""
function extract_point_states!(points, system, assembly, x = system.x;
    prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}())

    for ipoint = 1:length(points)
        points[ipoint] = extract_point_state(system, assembly, ipoint, x; prescribed_conditions)
    end

    return points
end

"""
    extract_element_state(system, assembly, ielem, x = system.x;
        prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}())

Return the state variables corresponding to element `ielem` (see [`ElementState`](@ref))
given the solution vector `x`.
"""
extract_element_state

function extract_element_state(system::StaticSystem, assembly, ielem, x = system.x; 
    prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}())

    # system variables
    @unpack force_scaling, indices, t = system

    # check state variable length
    @assert length(x) == length(system.x) "state vector length does not match with the provided system"

    # current prescribed conditions
    pc = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(t)

    # linear and angular displacement
    u1, θ1 = point_displacement(x, assembly.start[ielem], indices.icol_point, pc)
    u2, θ2 = point_displacement(x, assembly.stop[ielem], indices.icol_point, pc)
    u = (u1 + u2)/2
    θ = (θ1 + θ2)/2

    # internal forces and moments
    F, M = element_loads(x, ielem, indices.icol_elem, force_scaling)

    # linear and angular velocity
    V = zero(u)
    Ω = zero(θ)

    # convert rotation parameter to Wiener-Milenkovic parameters
    scaling = rotation_parameter_scaling(θ)
    θ *= scaling

    # promote all variables to the same type
    u, θ, V, Ω, F, M = promote(u, θ, V, Ω, F, M)

    return ElementState(u, θ, V, Ω, F, M)
end

function extract_element_state(system::DynamicSystem, assembly, ielem, x = system.x; 
    prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}())

    # system variables
    @unpack force_scaling, indices, t = system

    # check state variable length
    @assert length(x) == length(system.x) "state vector length does not match with the provided system"

    # current prescribed conditions
    pc = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(t)

    # linear and angular displacement
    u1, θ1 = point_displacement(x, assembly.start[ielem], indices.icol_point, pc)
    u2, θ2 = point_displacement(x, assembly.stop[ielem], indices.icol_point, pc)
    u = (u1 + u2)/2
    θ = (θ1 + θ2)/2

    # internal forces and moments
    F, M = element_loads(x, ielem, indices.icol_elem, force_scaling)

    # linear and angular velocity
    V1, Ω1 = point_velocities(x, assembly.start[ielem], indices.icol_point)
    V2, Ω2 = point_velocities(x, assembly.stop[ielem], indices.icol_point)
    V = (V1 + V2)/2
    Ω = (Ω1 + Ω2)/2

    # convert rotation parameter to Wiener-Milenkovic parameters
    scaling = rotation_parameter_scaling(θ)
    θ *= scaling

    # promote all variables to the same type
    u, θ, V, Ω, F, M = promote(u, θ, V, Ω, F, M)

    return ElementState(u, θ, V, Ω, F, M)
end

function extract_element_state(system::ExpandedSystem, assembly, ielem, x = system.x; 
    prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}())

    # system variables
    @unpack force_scaling, indices, t = system

    # check state variable length
    @assert length(x) == length(system.x) "state vector length does not match with the provided system"

    # current prescribed conditions
    pc = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions : prescribed_conditions(t)

    # linear and angular displacement
    u1, θ1 = point_displacement(x, assembly.start[ielem], indices.icol_point, pc)
    u2, θ2 = point_displacement(x, assembly.stop[ielem], indices.icol_point, pc)
    u = (u1 + u2)/2
    θ = (θ1 + θ2)/2

    # element forces and moments
    F1, M1, F2, M2 = expanded_element_loads(x, ielem, indices.icol_elem, force_scaling)
    F = (F1 + F2)/2
    M = (M1 + M2)/2

    # linear and angular displacement of the body frame
    ub, θb = body_frame_displacement(x)

    # linear and angular velocity of the body frame
    vb, ωb = body_frame_velocity(x)

    # linear and angular acceleration of the body frame
    ab, αb = body_frame_acceleration(x)

    # distance from the rotation center
    Δx = assembly.elements[ielem].x

    # linear and angular velocity
    v = vb + cross(ωb, Δx) + cross(ωb, u)# + udot = CtCab*V
    ω = ωb# + Cab'*Q*θdot = CtCab*Ω

    # linear and angular acceleration
    a = ab + cross(αb, Δx) + cross(αb, u)# + cross(ω, V) + uddot = CtCabdot*V + CtCab*Vdot
    α = αb# + Cab'*Qdot*θdot + Cab'*Q*θddot = CtCabdot*Ωdot + CtCab*Ωdot

    # linear and angular velocity **excluding contributions from body frame motion**
    V, Ω = expanded_element_velocities(x, ielem, indices.icol_elem)

    # rotate linear and angular velocites into the body frame
    C = get_C(θ)
    Cab = assembly.elements[ielem].Cab
    CtCab = C'*Cab
    V = CtCab*V
    Ω = CtCab*Ω

    # add contributions from body frame motion to velocities
    V += v
    Ω += ω

    # convert rotation parameter to Wiener-Milenkovic parameters
    scaling = rotation_parameter_scaling(θ)
    θ *= scaling

    # promote all variables to the same type
    u, θ, V, Ω, F, M = promote(u, θ, V, Ω, F, M)

    return ElementState(u, θ, V, Ω, F, M)
end

"""
    extract_element_states(system, assembly, x = system.x;
        prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}())

Return the state variables corresponding to each element (see [`ElementState`](@ref))
given the solution vector `x`.
"""
function extract_element_states(system, assembly, x = system.x; kwargs...)
    TF = promote_type(eltype(system), eltype(x))
    elements = Vector{ElementState{TF}}(undef, length(assembly.elements))
    return extract_element_states!(elements, system, assembly, x; kwargs...)
end

"""
    extract_element_states!(elements, system, assembly, x = system.x;
        prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}())

Pre-allocated version of [`extract_element_states`](@ref)
"""
function extract_element_states!(elements, system, assembly, x = system.x; kwargs...)
    for ielem = 1:length(elements)
        elements[ielem] = extract_element_state(system, assembly, ielem, x; kwargs...)
    end
    return elements
end

"""
    rotate(xyz, r, theta)

Rotate the vectors in `xyz` about point `r` using the Wiener-Milenkovic
parameters in `theta`.
"""
rotate(xyz, r, theta) = rotate!(copy(xyz), r, theta)

"""
    rotate!(xyz, r, theta)

Pre-allocated version of [`rotate`](@ref)
"""
function rotate!(xyz, r, theta)
    # reshape cross section points (if necessary)
    rxyz = reshape(xyz, 3, :)
    # convert inputs to static arrays
    r = SVector{3}(r)
    theta = SVector{3}(theta)
    # create rotation matrix
    Ct = wiener_milenkovic(theta)'
    # rotate each point
    for ipt = 1:size(rxyz, 2)
        p = SVector(rxyz[1,ipt], rxyz[2,ipt], rxyz[3,ipt])
        rxyz[:,ipt] .= Ct*(p - r) + r
    end
    # return the result
    return xyz
end

"""
    translate(xyz, u)

Translate the points in `xyz` by the displacements in `u`.
"""
translate(xyz, u) = translate!(copy(xyz), u)

"""
    translate!(xyz, u)

Pre-allocated version of [`translate`](@ref)
"""
function translate!(xyz, u)
    # reshape cross section points (if necessary)
    rxyz = reshape(xyz, 3, :)
    # convert inputs to static arrays
    u = SVector{3}(u)
    # translate each point
    for ipt = 1:size(rxyz, 2)
        p = SVector(rxyz[1,ipt], rxyz[2,ipt], rxyz[3,ipt])
        rxyz[:,ipt] .= p .+ u
    end
    # return the result
    return xyz
end

"""
    deform_cross_section(xyz, r, u, theta)

Rotate the points in `xyz` (of shape (3, :)) about point `r` using
the Wiener-Milenkovic parameters in `theta`, then translate the points by the
displacements in `u`.
"""
deform_cross_section(xyz, r, u, theta) = deform_cross_section!(copy(xyz), r, u, theta)

"""
    deform_cross_section!(xyz, r, u, theta)

Pre-allocated version of [`deform_cross_section`](@ref)
"""
deform_cross_section!(xyz, r, u, theta) = translate!(rotate!(xyz, r, theta), u)

"""
    write_vtk(name, assembly::Assembly; kwargs...)
    write_vtk(name, assembly::Assembly, state::AssemblyState; kwargs...)
    write_vtk(name, assembly::Assembly, history::Vector{<:AssemblyState}], dt;
        kwargs...)

Write the deformed geometry (and associated data) to a VTK file for visualization
using ParaView.

The `state` argument may be omitted to write the original geometry to a VTK file
without any associated data.

If the solution time `history` is provided, the time step must also be provided

# Keyword Arguments
 - `body_motion = true`: Flag indicating whether to visualize body frame motion.
 - `sections = nothing`: Cross section geometry corresponding to each point,
    defined in a frame aligned with the body frame but centered around the
    corresponding point. Defined as an array with shape `(3, ncross, np)` where `ncross`
    is the number of points in each cross section and `np` is the number of points.
 - `scaling=1.0`: Parameter to scale the deflections (only valid if state is provided)
 - `metadata=Dict()`: Dictionary of metadata for the file(s)
"""
function write_vtk(name, assembly; sections=nothing, metadata=Dict())

    # get problem dimensions
    npoint = length(assembly.points)
    ncross = isnothing(sections) ? 1 : size(sections, 2)
    nelem = length(assembly.elements)

    if isnothing(sections)
        # extract point locations
        points = Matrix{eltype(assembly)}(undef, 3, npoint)
        for ip = 1:npoint
            for i = 1:3
                points[i,ip] = assembly.points[ip][i]
            end
        end

        # create cells
        cells = [MeshCell(PolyData.Lines(), [assembly.start[i], assembly.stop[i]]) for i = 1:nelem]

    else

        li = LinearIndices((ncross, npoint))

        # extract cell point locations
        points = Matrix{eltype(assembly)}(undef, 3, ncross*npoint)
        if ndims(sections) > 2
            for ip = 1:npoint
                for ic = 1:ncross
                    points[:,li[ic,ip]] = assembly.points[ip] + sections[:,ic,ip]
                end
            end
        else
            for ip = 1:npoint
                for ic = 1:ncross
                    points[:,li[ic,ip]] = assembly.points[ip] + sections[:,ic]
                end
            end
        end

        # construct triangle strip for each beam element
        cells = Vector{MeshCell{VTKCellType, Vector{Int64}}}(undef, nelem)
        for ielem = 1:nelem
            # index of key point corresponding to the start of the beam element
            ipt1 = assembly.start[ielem]
            # index of key point corresponding to the end of the beam element
            ipt2 = assembly.stop[ielem]
            # triangle strip points
            connectivity = Vector{Int}(undef, ncross*2)
            for ic = 1:ncross
                connectivity[2*ic-1] = li[ic, ipt1]
                connectivity[2*ic] = li[ic, ipt2]
            end
            cells[ielem] = MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP, connectivity)
        end

    end

    # write vtk file
    vtk_grid(name, points, cells) do vtkfile

        # add metadata
        for (key, value) in pairs(metadata)
            vtkfile[string(key)] = value
        end

        # add local axis data
        axis_name = ["x-axis", "y-axis", "z-axis"]
        axis_vector = [e1, e2, e3]
        for i = 1:3
            data = Matrix{eltype(assembly)}(undef, 3, nelem)
            for ielem = 1:nelem
                data[:, ielem] .= assembly.elements[ielem].Cab*axis_vector[i]
            end
            vtkfile[axis_name[i], VTKCellData()] = data
        end
    end

    return nothing
end

function write_vtk(name, assembly, state; body_motion = true, sections = nothing, 
    scaling=1.0, metadata=Dict())

    # get problem dimensions
    npoint = length(assembly.points)
    ncross = isnothing(sections) ? 1 : size(sections, 2)
    nelem = length(assembly.elements)

    # extract body position and orientation
    if body_motion
        ub = scaling*state.body.u
        θb = scaling*state.body.theta
    else
        ub = @SVector zeros(3)
        θb = @SVector zeros(3)
    end

    # rotation matrix for body frame rotation
    Ct_b = wiener_milenkovic(θb)'

    if isnothing(sections)
        # extract point locations
        points = Matrix{eltype(assembly)}(undef, 3, npoint)
        for ip = 1:npoint
            for i = 1:3
                points[i,ip] = ub + Ct_b*(assembly.points[ip][i] + scaling*(state.points[ip].u[i]))
            end
        end

        # create cells
        cells = [MeshCell(PolyData.Lines(), [assembly.start[i], assembly.stop[i]]) for i = 1:nelem]
    else
        li = LinearIndices((ncross, npoint))

        # extract cell point locations
        points = Matrix{eltype(assembly)}(undef, 3, ncross*npoint)
        if ndims(sections) > 2
            for ip = 1:npoint
                for ic = 1:ncross
                    u = scaling*state.points[ip].u
                    c = scaling*state.points[ip].theta
                    Ct_p = wiener_milenkovic(c)'
                    points[:,li[ic,ip]] = ub + Ct_b*(assembly.points[ip] + u + Ct_p*sections[:,ic,ip])
                end
            end
        else
            for ip = 1:npoint
                for ic = 1:ncross
                    u = scaling*state.points[ip].u
                    c = scaling*state.points[ip].theta
                    Ct_p = wiener_milenkovic(c)'
                    points[:,li[ic,ip]] = ub + Ct_b*(assembly.points[ip] + u + Ct_p*sections[:,ic])
                end
            end
        end

        # construct triangle strip for each beam element
        cells = Vector{MeshCell{VTKCellType, Vector{Int64}}}(undef, nelem)
        for ielem = 1:nelem
            # index of key point corresponding to the start of the beam element
            ipt1 = assembly.start[ielem]
            # index of key point corresponding to the end of the beam element
            ipt2 = assembly.stop[ielem]
            # triangle strip points
            connectivity = Vector{Int}(undef, ncross*2)
            for ic = 1:ncross
                connectivity[2*ic-1] = li[ic, ipt1]
                connectivity[2*ic] = li[ic, ipt2]
            end
            cells[ielem] = MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP, connectivity)
        end
    end

    # write vtk file
    vtk_grid(name, points, cells) do vtkfile

        vtkfile["scaling"] = scaling

        # add metadata
        for (key, value) in pairs(metadata)
            vtkfile[string(key)] = value
        end

        # add local axis data
        axis_name = ["x-axis", "y-axis", "z-axis"]
        axis_vector = [e1, e2, e3]
        for i = 1:3
            data = Matrix{eltype(assembly)}(undef, 3, nelem)
            for ielem = 1:nelem
                CtCab = wiener_milenkovic(state.elements[ielem].theta)'*assembly.elements[ielem].Cab
                data[:, ielem] .= CtCab * axis_vector[i]
            end
            vtkfile[axis_name[i], VTKCellData()] = data
        end

        # add body frame data
        for field in fieldnames(Body)
            vtkfile[string(field)] = getproperty(state.body, field)
        end

        # add point data
        for field in fieldnames(PointState)
            data = Matrix{eltype(assembly)}(undef, 3, npoint*ncross)
            li = LinearIndices((ncross, npoint))
            for ip = 1:npoint
                data[:, li[:,ip]] .= getproperty(state.points[ip], field)
            end
            vtkfile[string(field), VTKPointData()] = data
        end

        # add cell data
        for field in fieldnames(ElementState)
            data = Matrix{eltype(assembly)}(undef, 3, nelem)
            for ielem = 1:nelem
                data[:, ielem] .= getproperty(state.elements[ielem], field)
            end
            vtkfile[string(field), VTKCellData()] = data
        end

    end

    return nothing
end

function write_vtk(name, assembly, history, t; sections=nothing, scaling=1.0,
    metadata=Dict())

    # get problem dimensions
    npoint = length(assembly.points)
    ncross = isnothing(sections) ? 1 : size(sections, 2)
    nelem = length(assembly.elements)

    paraview_collection(name) do pvd

        for (current_step, state) in enumerate(history)

            # extract body position and orientation
            if body_motion
                ub = scaling*state.body.u
                θb = scaling*state.body.theta
            else
                ub = @SVector zeros(3)
                θb = @SVector zeros(3)
            end

            # rotation matrix for body frame rotation
            Ct_b = wiener_milenkovic(θb)'

            if isnothing(sections)
                # extract point locations
                points = Matrix{eltype(assembly)}(undef, 3, npoint)
                for ip = 1:npoint
                    for i = 1:3
                        points[i,ip] = ub + Ct_b*(assembly.points[ip][i] + scaling*state.points[ip].u[i])
                    end
                end

                # create cells
                cells = [MeshCell(PolyData.Lines(), [assembly.start[i], assembly.stop[i]]) for i = 1:nelem]
            else

                li = LinearIndices((ncross, npoint))

                # extract cell point locations
                points = Matrix{eltype(assembly)}(undef, 3, ncross*npoint)
                if ndims(sections) > 2
                    for ip = 1:npoint
                        for ic = 1:ncross
                            u = scaling*state.points[ip].u
                            c = scaling*state.points[ip].theta
                            Ct_p = wiener_milenkovic(c)'
                            points[:,li[ic,ip]] = ub + Ct_b*(assembly.points[ip] + u + Ct_p*sections[:,ic,ip])
                        end
                    end
                else
                    for ip = 1:npoint
                        for ic = 1:ncross
                            u = scaling*state.points[ip].u
                            c = scaling*state.points[ip].theta
                            Ct_p = wiener_milenkovic(c)'
                            points[:,li[ic,ip]] = ub + Ct_b*(assembly.points[ip] + u + Ct_p*sections[:,ic])
                        end
                    end
                end

                # construct triangle strip for each beam element
                cells = Vector{MeshCell{VTKCellType, Vector{Int64}}}(undef, nelem)
                for ielem = 1:nelem
                    # index of key point corresponding to the start of the beam element
                    ipt1 = assembly.start[ielem]
                    # index of key point corresponding to the end of the beam element
                    ipt2 = assembly.stop[ielem]
                    # triangle strip points
                    connectivity = Vector{Int}(undef, ncross*2)
                    for ic = 1:ncross
                        connectivity[2*ic-1] = li[ic, ipt1]
                        connectivity[2*ic] = li[ic, ipt2]
                    end
                    cells[ielem] = MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP, connectivity)
                end
            end

            # write vtk file
            vtkfile = vtk_grid(name*"-step$current_step", points, cells)

            # add metadata
            vtkfile["scaling"] = scaling
            vtkfile["time"] = t[current_step]
            for (key, value) in pairs(metadata)
                vtkfile[string(key)] = value
            end

            # add local axis data
            axis_name = ["x-axis", "y-axis", "z-axis"]
            axis_vector = [e1, e2, e3]
            for i = 1:3
                data = Matrix{eltype(assembly)}(undef, 3, nelem)
                for ielem = 1:nelem
                    CtCab = wiener_milenkovic(state.elements[ielem].theta)' * assembly.elements[ielem].Cab
                    data[:, ielem] .= CtCab * axis_vector[i]
                end
                vtkfile[axis_name[i], VTKCellData()] = data
            end

            # add body frame data
            for field in fieldnames(Body)
                vtkfile[string(field)] = getproperty(state.body, field)
            end

            # add point data
            for field in fieldnames(PointState)
                data = Matrix{eltype(assembly)}(undef, 3, npoint*ncross)
                li = LinearIndices((ncross, npoint))
                for ip = 1:npoint
                    data[:, li[:,ip]] .= getproperty(state.points[ip], field)
                end
                vtkfile[string(field), VTKPointData()] = data
            end

            # add cell data
            for field in fieldnames(ElementState)
                data = Matrix{eltype(assembly)}(undef, 3, nelem)
                for ielem = 1:nelem
                    data[:, ielem] .= getproperty(state.elements[ielem], field)
                end
                vtkfile[string(field), VTKCellData()] = data
            end

            pvd[t[current_step]] = vtkfile
        end
    end

    return nothing
end

"""
    write_vtk(name, assembly::Assembly, [state::AssemblyState, ]λ::Number,
        eigenstate::AssemblyState; scaling=1.0, mode_scaling=1.0, cycles=1,
        steps=100)

Write a series of files corresponding to the elastic motion of the `assembly`
about the deformed state encoded in `state` defined by the eigenvalue `λ` and
the eigenvector encoded in `eigenstate` over the time period specified by `time`.

The steady-state deflections can be scaled with `scaling` and the eigenmode
deflections can be scaled using `mode_scaling`.

The current time is encoded in the metadata tag "time"
"""
function write_vtk(name, assembly, state, λ, eigenstate;
    sections = nothing,
    scaling=1.0,
    mode_scaling=1.0,
    cycles = 1,
    steps = 100)

    npoint = length(assembly.points)
    ncross = isnothing(sections) ? 1 : size(sections, 2)
    nelem = length(assembly.elements)

    damping = real(λ)
    period = 2*pi/imag(λ)

    if damping <= 0.0
        start = -period*cycles
        stop = 0.0
    else
        start = 0.0
        stop = period*cycles
    end

    time=range(start, stop, length=steps)

    paraview_collection(name) do pvd

        for (current_step, t) in enumerate(time)

            # extract body position and orientation
            if body_motion
                ub = scaling*state.body.u + mode_scaling*real(state.body.u)*exp(λ*t)
                θb = scaling*state.body.theta + mode_scaling*real(state.body.theta)*exp(λ*t)
            else
                ub = @SVector zeros(3)
                θb = @SVector zeros(3)
            end

            # rotation matrix for body frame rotation
            Ct_b = wiener_milenkovic(θb)'

            if isnothing(sections)

                # extract point locations
                points = Matrix{eltype(assembly)}(undef, 3, npoint)
                for ip = 1:npoint
                    for i = 1:3
                        points[i,ip] = ub + Ct_b*(assembly.points[ip][i] +
                            scaling*state.points[ip].u[i] +
                            mode_scaling*real(eigenstate.points[ip].u[i]*exp(λ*t)))
                    end
                end

                # create cells
                cells = [MeshCell(PolyData.Lines(), [assembly.start[i], assembly.stop[i]]) for i = 1:nelem]

            else

                li = LinearIndices((ncross, npoint))

                # extract cell point locations
                points = Matrix{eltype(assembly)}(undef, 3, ncross*npoint)
                if ndims(sections) > 2
                    for ip = 1:npoint
                        for ic = 1:ncross
                            u = mode_scaling*real.(eigenstate.points[ip].u .* exp(λ*t))
                            c = mode_scaling*real.(eigenstate.points[ip].theta .* exp(λ*t))
                            Ct = wiener_milenkovic(c)'
                            points[:,li[ic,ip]] = ub + Ct_b*(assembly.points[ip] + u + Ct*sections[:,ic,ip])
                        end
                    end
                else
                    for ip = 1:npoint
                        for ic = 1:ncross
                            u = mode_scaling*real.(eigenstate.points[ip].u .* exp(λ*t))
                            c = mode_scaling*real.(eigenstate.points[ip].theta .* exp(λ*t))
                            Ct = wiener_milenkovic(c)'
                            points[:,li[ic,ip]] = ub + Ct_b*(assembly.points[ip] + u + Ct*sections[:,ic])
                        end
                    end
                end

                # construct triangle strip for each beam element
                cells = Vector{MeshCell{VTKCellType, Vector{Int64}}}(undef, nelem)
                for ielem = 1:nelem
                    # index of key point corresponding to the start of the beam element
                    ipt1 = assembly.start[ielem]
                    # index of key point corresponding to the end of the beam element
                    ipt2 = assembly.stop[ielem]
                    # triangle strip points
                    connectivity = Vector{Int}(undef, ncross*2)
                    for ic = 1:ncross
                        connectivity[2*ic-1] = li[ic, ipt1]
                        connectivity[2*ic] = li[ic, ipt2]
                    end
                    cells[ielem] = MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP, connectivity)
                end

            end

            # write vtk file
            vtkfile = vtk_grid(name*"-step$(current_step)", points, cells)

            # add metadata
            vtkfile["scaling"] = scaling
            vtkfile["time"] = t
            vtkfile["phase"] = 2*pi*t/period

            # add local axis data
            axis_name = ["x-axis", "y-axis", "z-axis"]
            axis_vector = [e1, e2, e3]
            for i = 1:3
                data = Matrix{eltype(assembly)}(undef, 3, nelem)
                for ielem = 1:nelem
                    CtCab = wiener_milenkovic(real(state.elements[ielem].theta .* exp(λ*t)))' *
                        assembly.elements[ielem].Cab
                    data[:, ielem] .= CtCab * axis_vector[i]
                end
                vtkfile[axis_name[i], VTKCellData()] = data
            end

            # add body data
            for field in fieldnames(Body)
                vtkfile[string(field)] = getproperty(state.body, field) +
                    mode_scaling*real(getproperty(eigenstate.body, field) .* exp(λ*t))
            end

            # add point data
            for field in fieldnames(PointState)
                data = Matrix{eltype(assembly)}(undef, 3, ncross*npoint)
                li = LinearIndices((ncross, npoint))
                for ip = 1:npoint
                    data[:, li[:,ip]] .= getproperty(state.points[ip], field) +
                         mode_scaling*real(getproperty(eigenstate.points[ip], field) .* exp(λ*t))
                end
                vtkfile[string(field), VTKPointData()] = data
            end

            # add cell data
            for field in fieldnames(ElementState)
                data = Matrix{eltype(assembly)}(undef, 3, nelem)
                for ielem = 1:nelem
                    data[:, ielem] .= getproperty(state.elements[ielem], field) +
                         mode_scaling*real(getproperty(eigenstate.elements[ielem], field).* exp(λ*t))
                end
                vtkfile[string(field), VTKCellData()] = data
            end

            pvd[t] = vtkfile
        end
    end

    return nothing
end

