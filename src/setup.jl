# determines row and column indices for the assembly
function system_indexes(points, beams, steady)

    npt = length(points)
    nbeam = length(beams)

    irow_pt = Array{Int, 1}(undef, npt)
    irow_beam1 = Array{Int, 1}(undef, nbeam)
    irow_beam2 = Array{Int, 1}(undef, nbeam)

    icol_pt = Array{Int, 1}(undef, npt)
    icol_beam = Array{Int, 1}(undef, nbeam)

    irow = 1
    icol = 1
    for ibeam = 1:nbeam
        ipt = beams[ibeam].pt1
        if !isassigned(irow_pt, ipt)
            # 6 equilibrium equations + 6 compatability equations
            irow_pt[ipt] = irow
            irow_beam1[ipt] = irow
            irow += 12
            # 6 unknowns for each point
            icol_pt[ipt] = icol
            icol += 6
        else
            # 6 compatibility equations
            irow_beam1[ipt] = irow
            irow += 6
        end

        # 12 unknowns for each element
        icol_beam[ibeam] = icol
        icol += 12

        if !steady
            # 6 linear and angular momentum residual equations for each element
            irow += 6
            # 6 additional unknowns for each member for unsteady analyses
            icol += 6
        end

        ipt = beams[ibeam].pt2
        if !isassigned(irow_pt, ipt)
            # 6 equilibrium equations + 6 compatability equations
            irow_pt[ipt] = irow
            irow_beam2[ipt] = irow
            irow += 12
            # 6 unknowns for each point
            icol_pt[ipt] = icol
            icol += 6
        else
            # 6 compatibility equations
            irow_beam2[ipt] = irow
            irow += 6
        end
    end

    return irow_pt, irow_beam1, irow_beam2, icol_pt, icol_beam
end
