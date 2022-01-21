# GXBeam

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://flow.byu.edu/GXBeam.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://flow.byu.edu/GXBeam.jl/dev)
![](https://github.com/byuflowlab/GXBeam.jl/workflows/Run%20tests/badge.svg)
[![status](https://joss.theoj.org/papers/13cb0c41e9834510c6acf732bdfa8c05/status.svg)](https://joss.theoj.org/papers/13cb0c41e9834510c6acf732bdfa8c05)

*Pure Julia Implementation of Geometrically Exact Beam Theory*

Author: Taylor McDonnell

**GXBeam** is a pure Julia implementation of Geometrically Exact Beam Theory, originally based on the open source code [GEBT](https://cdmhub.org/resources/367) and its associated papers[[1]](#1)[[2]](#2).

As a sample of one of the many things this package can do, here's a time domain simulation of the dynamic response of a joined wing subjected to a simulated gust, scaled up in order to visualize the deflections:
![](docs/src/assets/dynamic-joined-wing.gif)

And here's a dynamic simulation of a wind turbine subjected to a sinusoidal tip load.
![](docs/src/assets/dynamic-wind-turbine.gif)

## Package Features
 - Performs multiple types of analyses including:
    - Linear/Nonlinear static analyses
    - Linear/Nonlinear steady-state dynamic analyses
    - Linear/Nonlinear eigenvalue analyses (by linearizing about a steady state condition)
    - Linear/Nonlinear time-marching dynamic analyses
 - Accurately models arbitrary systems of interconnected highly flexible composite beams.
    - Captures all geometric nonlinearities due to large deflections and rotations
    - Capable of using the full 6x6 Timoshenko beam stiffness matrix
    - Singularity-free rotational deflections of any magnitude using only 3 rotational parameters
 - Models arbitrary time-varying distributed forces/moments on beam elements using:
    - Dead forces/moments (which do not rotate as the beam element rotates)
    - Follower forces/moments (which rotate as the beam element rotates)
    - Forces/moments due to the presence of rigidly attached point masses
    - Forces/moments due to gravitational loads
    - Forces/moments due to body frame linear/angular velocities and accelerations
 - Models arbitrary time-varying prescribed forces/moments and/or displacements/rotations at the connection points between beam elements using:
    - Dead forces/moments (which do not rotate as the point rotates)
    - Follower forces/moments (which rotate as the point rotates)
 - Capable of using arbitrary units (as long as they are compatible)
 - Simple result visualization using [WriteVTK](https://github.com/jipolanco/WriteVTK.jl)
 - Built-in [DifferentialEquations](https://github.com/SciML/DifferentialEquations.jl) interface for time domain simulations.
 - Extensively validated against published analytical and computational results.  See the [examples in the documentation](https://flow.byu.edu/GXBeam.jl/dev/examples/).

## Installation

Enter the package manager by typing `]` and then run the following:

```julia
pkg> add GXBeam
```

## Performance

This code has been optimized to be highly performant.  In our tests we found that GXBeam outperforms GEBT by a significant margin across all analysis types, as seen in the following table.  More details about the specific cases which we test may be found by inspecting the input files and scripts for these tests in the `benchmark` folder.

| Package | Steady Analysis | Eigenvalue Analysis | Time Marching Analysis |
|---- | ----| --- | --- |
| GEBT | 13.722 ms | 33.712 ms | 26.870 s |
| GXBeam | 4.716 ms | 18.478 ms | 9.019 s |

## Usage

See the [documentation](https://flow.byu.edu/GXBeam.jl/dev)

## Limitations

By using the simplest possible shape functions (constant or linear shape functions), this package avoids using numerical quadrature except when integrating applied distributed loads (which can be pre-integrated).  As a result, element properties are approximated as constant throughout each beam element and a relatively large number of beam elements may be necessary to achieve grid-independent results.  More details about the convergence of this package may be found in the [examples](https://flow.byu.edu/GXBeam.jl/dev/examples/#Nonlinear-Analysis-of-a-Cantilever-Subjected-to-a-Constant-Moment).

This package does not currently model cross section warping, and therefore should not be used to model open cross sections (such as I, C, or L-beams).  The one exception to this rule is if the beam's width is much greater than its height, in which case the beam may be considered to be strip-like (like a helicopter blade).  

This package relies on the results of linear cross-sectional analyses.  Most notably, it does not model the nonlinear component of the Trapeze effect, which is the tendency of a beam to untwist when subjected to axial tension.  This nonlinear effect is typically most important when modeling rotating structures such as helicopter blades due to the presence of large centrifugal forces.  It is also more important when modeling strip-like beams than for modeling closed cross-section beams due to their low torsional rigidity.

## References
<a id="1">[1]</a>
Yu, W., & Blair, M. (2012).
GEBT: A general-purpose nonlinear analysis tool for composite beams.
Composite Structures, 94(9), 2677-2689.

<a id="2">[2]</a>
Wang, Q., & Yu, W. (2017).
Geometrically nonlinear analysis of composite beams using Wiener-Milenković parameters.
Journal of Renewable and Sustainable Energy, 9(3), 033306.

<a id="3">[3]</a> 
Hodges, D. (2006).
Nonlinear Composite Beam Theory.
American Institute of Aeronautics and Astronautics.
