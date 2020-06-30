
"""
    system_indices(points, beams, static)

Solve for the row indices of the first equilibrium or compatability equations for
each point and side of each beam element.  Also solve for the row/column index of
each point and beam state variable.

Return Values:
 - `N`: Number of rows/columns in system jacobian
 - irow_pt: Row index of first equilibrium equation for each point
 - irow_beam1: Row index of first equation for the left side of each beam
 - irow_beam2: Row index of first equation for the right side of each beam
 - icol_pt: Column index of first state variable for each point
 - icol_beam: Column index of first state variable for each beam element
"""
function system_indices(points, beams, static)

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

        if !static
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

    nrow = irow-1
    ncol = icol-1
    @assert nrow == ncol

    return nrow, irow_pt, irow_beam1, irow_beam2, icol_pt, icol_beam
end
