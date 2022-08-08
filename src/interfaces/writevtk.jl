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
            points[:,ip] = assembly.points[ip]
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

function write_vtk(name, assembly, state; sections = nothing, 
    scaling=1.0, metadata=Dict())

    # get problem dimensions
    npoint = length(assembly.points)
    ncross = isnothing(sections) ? 1 : size(sections, 2)
    nelem = length(assembly.elements)

    if isnothing(sections)
        # extract point locations
        points = Matrix{eltype(assembly)}(undef, 3, npoint)
        for ip = 1:npoint
            points[:,ip] = assembly.points[ip] + scaling*(state.points[ip].u)
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
                    points[:,li[ic,ip]] = assembly.points[ip] + u + Ct_p*sections[:,ic,ip]
                end
            end
        else
            for ip = 1:npoint
                for ic = 1:ncross
                    u = scaling*state.points[ip].u
                    c = scaling*state.points[ip].theta
                    Ct_p = wiener_milenkovic(c)'
                    points[:,li[ic,ip]] = assembly.points[ip] + u + Ct_p*sections[:,ic]
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

        vtkfile["scaling"] = float(scaling)

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

function write_vtk(name, assembly, history, t; sections=nothing, scaling=1.0, metadata=Dict())

    # get problem dimensions
    npoint = length(assembly.points)
    ncross = isnothing(sections) ? 1 : size(sections, 2)
    nelem = length(assembly.elements)

    paraview_collection(name) do pvd

        for (current_step, state) in enumerate(history)

            if isnothing(sections)
                # extract point locations
                points = Matrix{eltype(assembly)}(undef, 3, npoint)
                for ip = 1:npoint
                    points[:,ip] = assembly.points[ip] + scaling*state.points[ip].u
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
                            points[:,li[ic,ip]] = assembly.points[ip] + u + Ct_p*sections[:,ic,ip]
                        end
                    end
                else
                    for ip = 1:npoint
                        for ic = 1:ncross
                            u = scaling*state.points[ip].u
                            c = scaling*state.points[ip].theta
                            Ct_p = wiener_milenkovic(c)'
                            points[:,li[ic,ip]] = assembly.points[ip] + u + Ct_p*sections[:,ic]
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
            vtkfile["scaling"] = float(scaling)
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

            if isnothing(sections)

                # extract point locations
                points = Matrix{eltype(assembly)}(undef, 3, npoint)
                for ip = 1:npoint
                    points[:,ip] = assembly.points[ip] + scaling*state.points[ip].u +
                        mode_scaling*real.(eigenstate.points[ip].u*exp(λ*t))
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

