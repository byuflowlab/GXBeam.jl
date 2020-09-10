# Library

```@contents
Pages = ["library.md"]
Depth = 3
```

## Public API

### Creating an Assembly

```@docs
curve_length
discretize_beam
Assembly(points, start, stop)
```

### Defining Distributed Loads
```@docs
DistributedLoads(assembly, ibeam)
```

### Defining Prescribed Conditions

```@docs
PrescribedConditions()
```

### Pre-Initializing Memory for an Analysis

```@docs
System(assembly, points, static)
reset_state!
```

### Performing an Analysis

```@docs
static_analysis
static_analysis!
steady_state_analysis
steady_state_analysis!
eigenvalue_analysis
eigenvalue_analysis!
time_domain_analysis
time_domain_analysis!
```

### Post-Processing

```@docs
AssemblyState(system, assembly)
AssemblyState
left_eigenvectors
correlate_eigenmodes
wiener_milenkovic
write_vtk
```

## Private API

### Math

```@docs
GXBeam.tilde
GXBeam.rotation_parameter_scaling
GXBeam.get_C
GXBeam.get_C_t
GXBeam.get_C_θ
GXBeam.get_C_θdot
GXBeam.get_Q
GXBeam.get_Q_θ
GXBeam.get_Qinv
GXBeam.get_Qinv_θ
GXBeam.mul3
GXBeam.gauss_quadrature
```

### Points

```@docs
GXBeam.point_variables
GXBeam.insert_point_residual!
GXBeam.point_residual!
GXBeam.point_follower_jacobians
GXBeam.insert_point_jacobian!
GXBeam.point_jacobian!
GXBeam.PointState
```

### Elements

```@docs
GXBeam.Element
GXBeam.element_strain
GXBeam.element_curvature
GXBeam.element_linear_velocity
GXBeam.element_angular_velocity
GXBeam.element_properties
GXBeam.dynamic_element_properties
GXBeam.element_equations
GXBeam.insert_element_residual!
GXBeam.element_residual!
GXBeam.element_jacobian_equations
GXBeam.insert_element_jacobian!
GXBeam.element_jacobian!
GXBeam.element_mass_matrix_properties
GXBeam.element_mass_matrix_equations
GXBeam.insert_element_mass_matrix!
GXBeam.element_mass_matrix!
GXBeam.ElementState
```

### Loads

```@docs
GXBeam.PrescribedConditions
GXBeam.DistributedLoads
```

### System

```@docs
GXBeam.Assembly
GXBeam.curve_triad
GXBeam.curve_coordinates
GXBeam.System
GXBeam.point_connections
GXBeam.system_indices
GXBeam.system_residual!
GXBeam.system_jacobian!
GXBeam.system_mass_matrix!
```

## Index

```@index
```
