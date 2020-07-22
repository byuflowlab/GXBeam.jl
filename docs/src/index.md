# GEBT.jl

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://flow.byu.edu/GEBT.jl/dev)
![](https://github.com/byuflowlab/GEBT.jl/workflows/Run%20tests/badge.svg)

*(Almost) Pure Julia Implementation of Geometrically Exact Beam Theory*

Author: Taylor McDonnell

**GEBT.jl** is an (almost) pure Julia implementation of Geometrically Exact Beam Theory, based on the similarly named open source [GEBT code by Wenbin Yu](https://cdmhub.org/resources/367) and its associated papers[[1]](@ref References)[[2]](@ref References).  The "almost" here refers to the fact that the Fortran library ARPACK is used for eigenvalue computations.  Otherwise the code is written with pure Julia and should work with custom types and automatic differentiation packages such as [ForwardDiff](https://github.com/JuliaDiff/ForwardDiff.jl).

## Package Features
 - Performs multiple types of analyses including:
    - Nonlinear static analyses
    - Nonlinear steady-state dynamic analyses
    - Nonlinear eigenvalue analyses (by linearizing about a steady state condition)
    - Nonlinear time-marching dynamic analyses
 - Accurately models arbitrary systems of interconnected highly flexible composite beams.
    - Captures all geometric nonlinearities due to large deflections and rotations
    - Capable of using the full 6x6 Timoshenko beam stiffness matrix
 - Models arbitrary distributed forces/moments on beam elements using:
    - Dead forces/moments (which do not rotate as the beam element rotates)
    - Follower forces/moments (which rotate as the beam element rotates)
 - Models arbitrary prescribed forces/moments and/or displacements/rotations at the connection points between beam elements using:
    - Dead forces/moments (which do not rotate as the point rotates)
    - Follower forces/moments (which rotate as the point rotates)
 - Provides time functions to define the variation of the applied loads/displacements over time for both steady and unsteady simulations.
 - Capable of using arbitrary units (as long as they are compatible)
 - Simple result visualization using [WriteVTK](https://github.com/jipolanco/WriteVTK.jl)
 - Thoroughly validated against published analytical and computational results.

## Installation

Enter the package manager by typing `]` and then run the following:

```julia
pkg> add https://github.com/byuflowlab/GEBT.jl
```

## Performance

This code has been optimized to be highly performant, primarily by maintaining type stability and minimizing allocations.  As a result the performance of this package rivals (and sometimes beats) that of the Fortran implementation of GEBT provided by Wenbin Yu.  At this point, differences in performance between the two codes can be primarily attributed to the performance of the sparse linear system solver in each.

## Usage

See the [examples](@ref Examples)

Note that while the theory is identical to the Wenbin Yu's code, some of the implementation details vary.

## References
[1] Yu, W., & Blair, M. (2012).
GEBT: A general-purpose nonlinear analysis tool for composite beams.
Composite Structures, 94(9), 2677-2689.

[2] Wang, Q., & Yu, W. (2017).
Geometrically nonlinear analysis of composite beams using Wiener-Milenković parameters.
Journal of Renewable and Sustainable Energy, 9(3), 033306.