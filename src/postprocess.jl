"""
    PointState

Holds the state variables for a point

# Fields:
 - `u`: Displacement variables for the point (in the global coordinate frame)
 - `theta`: Wiener-Milenkovic rotational displacement variables for the point
 - `F`: Externally applied forces on the point (in the global coordinate frame)
 - `M`: Externally applied moments on the point (in the global coordinate frame)
"""
struct PointState{TF}
    u::SVector{3, TF}
    theta::SVector{3, TF}
    F::SVector{3, TF}
    M::SVector{3, TF}
end

"""
    ElementState

Holds the state variables for an element

# Fields:
 - `u`: Displacement variables for the element (in the global coordinate frame)
 - theta: Wiener-Milenkovic rotational displacement variables for the element
 - `F`: Resultant forces for the element (in the deformed beam coordinate frame)
 - `M`: Resultant moments for the element (in the deformed beam coordinate frame)
 - `P`: Linear momenta of the element (in the deformed beam coordinate frame)
 - `H`: Angular momenta of the element (in the deformed beam coordinate frame)
"""
struct ElementState{TF}
    u::SVector{3, TF}
    theta::SVector{3, TF}
    F::SVector{3, TF}
    M::SVector{3, TF}
    P::SVector{3, TF}
    H::SVector{3, TF}
end

"""
    AssemblyState{TF, TP<:AbstractVector{PointState{TF}}, TE<:AbstractVector{ElementState{TF}}}

Struct for storing state variables for the points and elements in an assembly.

# Fields:
 - `points::TP`: Array of `PointState`s for each point in the assembly
 - `elements::TE`: Array of `ElementState`s for each element in the assembly
"""
struct AssemblyState{TF, TP<:AbstractVector{PointState{TF}}, TE<:AbstractVector{ElementState{TF}}}
    points::TP
    elements::TE
end
Base.eltype(::AssemblyState{TF, TP, TE}) where {TF, TP, TE} = TF

"""
    AssemblyState(system, assembly, x = system.x;
        prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}())

Post-process the system state given the solution vector `x`.  Return an object
of type `AssemblyState` that defines the state of the assembly for the time step.

If `prescribed_conditions` is not provided, all point state variables are assumed
to be displacements/rotations, rather than their actual identities as used in the
analysis.
"""
function AssemblyState(system, assembly, x = system.x;
    prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}())

    points = extract_point_states(system, assembly, x; prescribed_conditions)

    elements = extract_element_states(system, x)

    return AssemblyState(points, elements)
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
function extract_point_state(system, assembly, ipoint, x = system.x;
    prescribed_conditions = Dict{Int,PrescribedConditions{Float64}}())

    pc = typeof(prescribed_conditions) <: AbstractDict ? prescribed_conditions :
        prescribed_conditions(system.t)

    force_scaling = system.force_scaling

    icol_p = system.icol_pt[ipoint]

    if icol_p <= 0
        # point variables are not state variables, solve for point variables

        # find the first beam connected to the point, checking both sides
        ibeam = findfirst(x -> x == ipoint, assembly.start)
        if !isnothing(ibeam)
            side = -1
        else
            ibeam = findfirst(x -> x == ipoint, assembly.stop)
            side = 1
        end

        beam = assembly.elements[ibeam]

        # get beam state variables
        icol_b = system.icol_beam[ibeam]
        u_b = SVector(x[icol_b   ], x[icol_b+1 ], x[icol_b+2 ])
        θ_b = SVector(x[icol_b+3 ], x[icol_b+4 ], x[icol_b+5 ])
        F_b = SVector(x[icol_b+6 ], x[icol_b+7 ], x[icol_b+8 ]) .* force_scaling
        M_b = SVector(x[icol_b+9 ], x[icol_b+10], x[icol_b+11]) .* force_scaling

        # get beam element properties
        ΔL_b = beam.L
        Ct_b = get_C(θ_b)'
        Cab_b = beam.Cab
        γ_b = element_strain(beam, F_b, M_b)
        κ_b = element_curvature(beam, F_b, M_b)

        # solve for the displacements/rotations of the point
        u = u_b + side*ΔL_b/2*(Ct_b*Cab_b*(e1 + γ_b) - Cab_b*e1)
        theta = θ_b + side*ΔL_b/2*get_Qinv(θ_b)*Cab_b*κ_b
        F = zero(u)
        M = zero(u)
    else
        # point variables are state variables, extract point state variables
        prescribed = haskey(pc, ipoint)
        if prescribed
            u, theta, F, M = point_variables(x, icol_p, pc[ipoint], force_scaling)
        else
            u, theta, F, M = point_variables(x, icol_p)
        end
    end

    # convert rotation parameter to Wiener-Milenkovic parameters
    scaling = rotation_parameter_scaling(theta)
    theta *= scaling

    # promote all variables to the same type
    u, theta, F, M = promote(u, theta, F, M)

    return PointState(u, theta, F, M)
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
    extract_element_state(system, ielem, x = system.x)

Return the state variables corresponding to element `ielem` (see [`ElementState`](@ref))
given the solution vector `x`.
"""
function extract_element_state(system, ielem, x = system.x)

    force_scaling = system.force_scaling
    mass_scaling = system.mass_scaling

    icol = system.icol_beam[ielem]
    static = system.irow_beam[ielem] <= 0
    u = SVector(x[icol], x[icol+1], x[icol+2])
    theta = SVector(x[icol+3], x[icol+4], x[icol+5])
    F = SVector(x[icol+6], x[icol+7], x[icol+8]) .* force_scaling
    M = SVector(x[icol+9], x[icol+10], x[icol+11]) .* force_scaling
    P = ifelse(static, zero(u), SVector(x[icol+12], x[icol+13], x[icol+14]) .* mass_scaling)
    H = ifelse(static, zero(u), SVector(x[icol+15], x[icol+16], x[icol+17]) .* mass_scaling)

    # convert rotation parameter to Wiener-Milenkovic parameters
    scaling = rotation_parameter_scaling(theta)
    theta *= scaling

    # promote all variables to the same type
    u, theta, F, M, P, H = promote(u, theta, F, M, P, H)

    return ElementState(u, theta, F, M, P, H)
end

"""
    extract_element_states(system, x = system.x)

Return the state variables corresponding to each element (see [`ElementState`](@ref))
given the solution vector `x`.
"""
function extract_element_states(system, x = system.x)
    TF = promote_type(eltype(system), eltype(x))
    elements = Vector{ElementState{TF}}(undef, length(system.icol_beam))
    return extract_element_states!(elements, system, x)
end

"""
    extract_element_states!(elements, system, x = system.x)

Pre-allocated version of [`extract_element_states`](@ref)
"""
function extract_element_states!(elements, system, x = system.x)
    for ielem = 1:length(elements)
        elements[ielem] = extract_element_state(system, ielem, x)
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
    cross_section_velocities(xyz, r, v, ω)

Calculate the velocities of the points in `xyz` given the linear velocity `v`,
and the angular velocity `ω` about point `r`.
"""
function cross_section_velocities(xyz, r, v, ω)
    TF = promote_type(eltype(xyz), eltype(r), eltype(v), eltype(ω))
    vxyz = similar(xyz, TF)
    return cross_section_velocities!(vxyz, xyz, r, v, ω)
end

"""
    cross_section_velocities!(vxyz, xyz, r, v, ω)

Pre-allocated version of [`cross_section_velocities`](@ref)
"""
function cross_section_velocities!(vxyz, xyz, r, v, ω)
    # reshape cross section points (if necessary)
    rxyz = reshape(xyz, 3, :)
    rvxyz = reshape(vxyz, 3, :)
    # convert inputs to static arrays
    r = SVector{3}(r)
    v = SVector{3}(v)
    ω = SVector{3}(ω)
    # calculate velocities at each point
    for ipt = 1:size(rxyz, 2)
        p = SVector(rxyz[1,ipt], rxyz[2,ipt], rxyz[3,ipt])
        rvxyz[:,ipt] .= v .+ cross(ω, p - r)
    end
    # return the result
    return vxyz
end

"""
    write_vtk(name, assembly::Assembly; kwargs...)
    write_vtk(name, assembly::Assembly, state::AssemblyState; kwargs...)
    write_vtk(name, assembly::Assembly, history::Vector{<:AssemblyState}], dt; kwargs...)

Write the deformed geometry (and associated data) to a VTK file for visualization
using ParaView.

The `state` argument may be omitted to write the original geometry to a VTK file
without any associated data.

If the solution time `history` is provided, the time step must also be provided

# Keyword Arguments
 - `sections = nothing`: Cross section geometry corresponding to each point,
    defined in a frame aligned with the global frame but centered around the
    corresponding point. Defined as an array with shape `(3, ncross, np)` where `ncross`
    is the number of points in each cross section and `np` is the number of points.
 - `scaling=1.0`: Parameter to scale the deflections (only valid if state is provided)
 - `metadata=Dict()`: Dictionary of metadata for the file(s)
"""
function write_vtk(name, assembly; sections=nothing, metadata=Dict())

    # get problem dimensions
    npoint = length(assembly.points)
    ncross = isnothing(sections) ? 1 : size(sections, 2)
    nbeam = length(assembly.elements)

    if isnothing(sections)
        # extract point locations
        points = Matrix{eltype(assembly)}(undef, 3, npoint)
        for ip = 1:npoint
            for i = 1:3
                points[i,ip] = assembly.points[ip][i]
            end
        end

        # create cells
        cells = [MeshCell(PolyData.Lines(), [assembly.start[i], assembly.stop[i]]) for i = 1:nbeam]

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
        cells = Vector{MeshCell{VTKCellType, Vector{Int64}}}(undef, nbeam)
        for ibeam = 1:nbeam
            # index of key point corresponding to the start of the beam element
            ipt1 = assembly.start[ibeam]
            # index of key point corresponding to the end of the beam element
            ipt2 = assembly.stop[ibeam]
            # triangle strip points
            connectivity = Vector{Int}(undef, ncross*2)
            for ic = 1:ncross
                connectivity[2*ic-1] = li[ic, ipt1]
                connectivity[2*ic] = li[ic, ipt2]
            end
            cells[ibeam] = MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP, connectivity)
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
            data = Matrix{eltype(assembly)}(undef, 3, nbeam)
            for ibeam = 1:nbeam
                data[:, ibeam] .= assembly.elements[ibeam].Cab*axis_vector[i]
            end
            vtkfile[axis_name[i], VTKCellData()] = data
        end
    end

    return nothing
end

function write_vtk(name, assembly, state; sections = nothing, scaling=1.0,
    metadata=Dict())

    # get problem dimensions
    npoint = length(assembly.points)
    ncross = isnothing(sections) ? 1 : size(sections, 2)
    nbeam = length(assembly.elements)

    if isnothing(sections)
        # extract point locations
        points = Matrix{eltype(assembly)}(undef, 3, npoint)
        for ip = 1:npoint
            for i = 1:3
                points[i,ip] = assembly.points[ip][i] + scaling*state.points[ip].u[i]
            end
        end

        # create cells
        cells = [MeshCell(PolyData.Lines(), [assembly.start[i], assembly.stop[i]]) for i = 1:nbeam]
    else
        li = LinearIndices((ncross, npoint))

        # extract cell point locations
        points = Matrix{eltype(assembly)}(undef, 3, ncross*npoint)
        if ndims(sections) > 2
            for ip = 1:npoint
                for ic = 1:ncross
                    u = scaling*state.points[ip].u
                    c = scaling*state.points[ip].theta
                    Ct = wiener_milenkovic(c)'
                    points[:,li[ic,ip]] = assembly.points[ip] + u + Ct*sections[:,ic,ip]
                end
            end
        else
            for ip = 1:npoint
                for ic = 1:ncross
                    u = scaling*state.points[ip].u
                    c = scaling*state.points[ip].theta
                    Ct = wiener_milenkovic(c)'
                    points[:,li[ic,ip]] = assembly.points[ip] + u + Ct*sections[:,ic]
                end
            end
        end

        # construct triangle strip for each beam element
        cells = Vector{MeshCell{VTKCellType, Vector{Int64}}}(undef, nbeam)
        for ibeam = 1:nbeam
            # index of key point corresponding to the start of the beam element
            ipt1 = assembly.start[ibeam]
            # index of key point corresponding to the end of the beam element
            ipt2 = assembly.stop[ibeam]
            # triangle strip points
            connectivity = Vector{Int}(undef, ncross*2)
            for ic = 1:ncross
                connectivity[2*ic-1] = li[ic, ipt1]
                connectivity[2*ic] = li[ic, ipt2]
            end
            cells[ibeam] = MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP, connectivity)
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
            data = Matrix{eltype(assembly)}(undef, 3, nbeam)
            for ibeam = 1:nbeam
                CtCab = wiener_milenkovic(state.elements[ibeam].theta)'*assembly.elements[ibeam].Cab
                data[:, ibeam] .= CtCab * axis_vector[i]
            end
            vtkfile[axis_name[i], VTKCellData()] = data
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
            data = Matrix{eltype(assembly)}(undef, 3, nbeam)
            for ibeam = 1:nbeam
                data[:, ibeam] .= getproperty(state.elements[ibeam], field)
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
    nbeam = length(assembly.elements)

    paraview_collection(name) do pvd

        for (current_step, state) in enumerate(history)

            if isnothing(sections)
                # extract point locations
                points = Matrix{eltype(assembly)}(undef, 3, npoint)
                for ip = 1:npoint
                    for i = 1:3
                        points[i,ip] = assembly.points[ip][i] + scaling*state.points[ip].u[i]
                    end
                end

                # create cells
                cells = [MeshCell(PolyData.Lines(), [assembly.start[i], assembly.stop[i]]) for i = 1:nbeam]
            else

                li = LinearIndices((ncross, npoint))

                # extract cell point locations
                points = Matrix{eltype(assembly)}(undef, 3, ncross*npoint)
                if ndims(sections) > 2
                    for ip = 1:npoint
                        for ic = 1:ncross
                            u = scaling*state.points[ip].u
                            c = scaling*state.points[ip].theta
                            Ct = wiener_milenkovic(c)'
                            points[:,li[ic,ip]] = assembly.points[ip] + u + Ct*sections[:,ic,ip]
                        end
                    end
                else
                    for ip = 1:npoint
                        for ic = 1:ncross
                            u = scaling*state.points[ip].u
                            c = scaling*state.points[ip].theta
                            Ct = wiener_milenkovic(c)'
                            points[:,li[ic,ip]] = assembly.points[ip] + u + Ct*sections[:,ic]
                        end
                    end
                end

                # construct triangle strip for each beam element
                cells = Vector{MeshCell{VTKCellType, Vector{Int64}}}(undef, nbeam)
                for ibeam = 1:nbeam
                    # index of key point corresponding to the start of the beam element
                    ipt1 = assembly.start[ibeam]
                    # index of key point corresponding to the end of the beam element
                    ipt2 = assembly.stop[ibeam]
                    # triangle strip points
                    connectivity = Vector{Int}(undef, ncross*2)
                    for ic = 1:ncross
                        connectivity[2*ic-1] = li[ic, ipt1]
                        connectivity[2*ic] = li[ic, ipt2]
                    end
                    cells[ibeam] = MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP, connectivity)
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
                data = Matrix{eltype(assembly)}(undef, 3, nbeam)
                for ibeam = 1:nbeam
                    CtCab = wiener_milenkovic(state.elements[ibeam].theta)' * assembly.elements[ibeam].Cab
                    data[:, ibeam] .= CtCab * axis_vector[i]
                end
                vtkfile[axis_name[i], VTKCellData()] = data
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
                data = Matrix{eltype(assembly)}(undef, 3, nbeam)
                for ibeam = 1:nbeam
                    data[:, ibeam] .= getproperty(state.elements[ibeam], field)
                end
                vtkfile[string(field), VTKCellData()] = data
            end

            pvd[t[current_step]] = vtkfile
        end
    end

    return nothing
end

"""
    write_vtk(name, assembly::Assembly, [state::AssemblyState, ]λ::Number, eigenstate::AssemblyState;
    scaling=1.0, mode_scaling=1.0, cycles=1, steps=100)

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
    nbeam = length(assembly.elements)

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

            if isnothing(sections)

                # extract point locations
                points = Matrix{eltype(assembly)}(undef, 3, npoint)
                for ip = 1:npoint
                    for i = 1:3
                        points[i,ip] = assembly.points[ip][i] +
                            scaling*state.points[ip].u[i] +
                            mode_scaling*real(eigenstate.points[ip].u[i]*exp(λ*t))
                    end
                end

                # create cells
                cells = [MeshCell(PolyData.Lines(), [assembly.start[i], assembly.stop[i]]) for i = 1:nbeam]

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
                            points[:,li[ic,ip]] = assembly.points[ip] + u + Ct*sections[:,ic,ip]
                        end
                    end
                else
                    for ip = 1:npoint
                        for ic = 1:ncross
                            u = mode_scaling*real.(eigenstate.points[ip].u .* exp(λ*t))
                            c = mode_scaling*real.(eigenstate.points[ip].theta .* exp(λ*t))
                            Ct = wiener_milenkovic(c)'
                            points[:,li[ic,ip]] = assembly.points[ip] + u + Ct*sections[:,ic]
                        end
                    end
                end

                # construct triangle strip for each beam element
                cells = Vector{MeshCell{VTKCellType, Vector{Int64}}}(undef, nbeam)
                for ibeam = 1:nbeam
                    # index of key point corresponding to the start of the beam element
                    ipt1 = assembly.start[ibeam]
                    # index of key point corresponding to the end of the beam element
                    ipt2 = assembly.stop[ibeam]
                    # triangle strip points
                    connectivity = Vector{Int}(undef, ncross*2)
                    for ic = 1:ncross
                        connectivity[2*ic-1] = li[ic, ipt1]
                        connectivity[2*ic] = li[ic, ipt2]
                    end
                    cells[ibeam] = MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP, connectivity)
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
                data = Matrix{eltype(assembly)}(undef, 3, nbeam)
                for ibeam = 1:nbeam
                    CtCab = wiener_milenkovic(real(state.elements[ibeam].theta .* exp(λ*t)))' *
                        assembly.elements[ibeam].Cab
                    data[:, ibeam] .= CtCab * axis_vector[i]
                end
                vtkfile[axis_name[i], VTKCellData()] = data
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
                data = Matrix{eltype(assembly)}(undef, 3, nbeam)
                for ibeam = 1:nbeam
                    data[:, ibeam] .= getproperty(state.elements[ibeam], field) +
                         mode_scaling*real(getproperty(eigenstate.elements[ibeam], field).* exp(λ*t))
                end
                vtkfile[string(field), VTKCellData()] = data
            end

            pvd[t] = vtkfile
        end
    end

    return nothing
end


"""
    left_eigenvectors(system, λ, V)
    left_eigenvectors(K, M, λ, V)

Compute the left eigenvector matrix `U` for the `system` using inverse power
iteration given the eigenvalues `λ` and the corresponding right eigenvector
matrix `V`.

The complex conjugate of each left eigenvector is stored in each row of the
matrix `U`

Left and right eigenvectors satisfy the following M-orthogonality condition:
 - u'*M*v = 1 if u and v correspond to the same eigenvalue
 - u'*M*v = 0 if u and v correspond to different eigenvalues
This means that U*M*V = I

This function assumes that `system` has not been modified since the eigenvalues
and right eigenvectors were computed.
"""
left_eigenvectors(system, λ, V) = left_eigenvectors(system.K, system.M, λ, V)

function left_eigenvectors(K, M, λ, V)

    # problem type and dimensions
    TC = eltype(V)
    nx = size(V,1)
    nev = size(V,2)

    # allocate storage
    U = rand(TC, nev, nx)
    u = Vector{TC}(undef, nx)
    tmp = Vector{TC}(undef, nx)

    # # get entries in M
    # iM, jM, valM = findnz(M)

    # compute eigenvectors for each eigenvalue
    for iλ = 1:nev

        # factorize (K + λ*M)'
        KmλMfact = factorize(K' + λ[iλ]'*M')

        # initialize left eigenvector
        for i = 1:nx
            u[i] = U[iλ,i]
        end

        # perform a few iterations to converge the left eigenvector
        for ipass = 1:3
            # get updated u
            mul!(tmp, M, u)
            ldiv!(u, KmλMfact, tmp)
            # normalize u
            unorm = zero(TC)
            for i in axes(M, 1), j in axes(M, 2)
                unorm += conj(u[i])*M[i,j]*V[j,iλ]
            end
            rdiv!(u, conj(unorm))
        end

        # store conjugate of final eigenvector
        for i = 1:nx
            U[iλ,i] = conj(u[i])
        end
    end

    return U
end

function left_eigenvectors(K, M::SparseMatrixCSC, λ, V)

    # problem type and dimensions
    TC = eltype(V)
    nx = size(V,1)
    nev = size(V,2)

    # allocate storage
    U = rand(TC, nev, nx)
    u = Vector{TC}(undef, nx)
    tmp = Vector{TC}(undef, nx)

    # get entries in M
    iM, jM, valM = findnz(M)

    # compute eigenvectors for each eigenvalue
    for iλ = 1:nev

        # factorize (K + λ*M)'
        KmλMfact = factorize(K' + λ[iλ]'*M')

        # initialize left eigenvector
        for i = 1:nx
            u[i] = U[iλ,i]
        end

        # perform a few iterations to converge the left eigenvector
        for ipass = 1:3
            # get updated u
            mul!(tmp, M, u)
            ldiv!(u, KmλMfact, tmp)
            # normalize u
            unorm = zero(TC)
            for k = 1:length(valM)
                unorm += conj(u[iM[k]])*valM[k]*V[jM[k],iλ]
            end
            rdiv!(u, conj(unorm))
        end

        # store conjugate of final eigenvector
        for i = 1:nx
            U[iλ,i] = conj(u[i])
        end
    end

    return U
end


"""
    correlate_eigenmodes(C)

Return the permutation and the associated corruption index vector which associates
eigenmodes from the current iteration with those of the previous iteration given
the correlation matrix `C`.

The correlation matrix can take one of the following forms (in order of preference):
 - `C = U_p*M*V`
 - `C = U*M_p*V_p`
 - `C = V_p'*V`
 - `C = V'*V_p`
where `U` is a matrix of conjugated left eigenvectors, `M` is the system mass
matrix, `V` is a matrix of right eigenvectors, and `()_p` indicates a variable
from the previous iteration.

Note that the following two forms of the correlation matrix seem to be significantly
inferior to their counterparts listed above: `C = U*M*V_p` and `C = U_p*M_p*V`.
This is likely due to the way in which the left eigenvector matrix is calculated.

The corruption index is the largest magnitude in a given row of `C`
that was not chosen divided by the magnitude of the chosen eigenmode.  It is most
meaningful when using one of the forms of the correlation matrix that uses left
eigenvectors since correct eigenmodes will have magnitudes close to 1 and
incorrect eigenmodes will have magnitudes close to 0.

If the new mode number is already assigned, the next highest unassigned mode
number is used.  In this case a corruption index higher than 1 will be returned,
otherwise the values of the corruption index will always be bounded by 0 and 1.

See "New Mode Tracking Methods in Aeroelastic Analysis" by Eldred, Vankayya, and
Anderson.
"""
function correlate_eigenmodes(C)

    # get row permutation that puts maximum values on the diagonals
    nev = size(C, 1)
    tmp = Vector{real(eltype(C))}(undef, nev)
    perm = zeros(Int, nev)
    corruption = Vector{real(eltype(C))}(undef, nev)
    for i = 1:nev

        # rank each value in this row
        for j = 1:nev
            tmp[j] = abs(C[i,j])
        end
        ranked_modes = sortperm(tmp, rev=true)

        # choose the best fit that is not yet assigned
        i1 = ranked_modes[findfirst((x) -> !(x in perm), ranked_modes)]
        i2 = i1 == ranked_modes[1] ? ranked_modes[2] : ranked_modes[1]

        # assign best eigenmode fit, create corruption index
        perm[i] = i1
        corruption[i] = tmp[i2]/tmp[i1]
    end

    return perm, corruption
end
