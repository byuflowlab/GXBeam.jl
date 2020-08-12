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

### Defining Prescribed Conditions

```@docs
PrescribedConditions()
DistributedLoads(assembly, ibeam)
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
GEBT.tilde
GEBT.rotation_parameter_scaling
GEBT.get_C
GEBT.get_C_t
GEBT.get_C_θ
GEBT.get_C_θdot
GEBT.get_Q
GEBT.get_Q_θ
GEBT.get_Qinv
GEBT.get_Qinv_θ
GEBT.mul3
GEBT.gauss_quadrature
```

### Points

```@docs
GEBT.point_variables
GEBT.insert_point_residual!
GEBT.point_residual!
GEBT.point_follower_jacobians
GEBT.insert_point_jacobian!
GEBT.point_jacobian!
GEBT.PointState
```

### Elements

```@docs
GEBT.Element
GEBT.element_strain
GEBT.element_curvature
GEBT.element_linear_velocity
GEBT.element_angular_velocity
GEBT.element_properties
GEBT.dynamic_element_properties
GEBT.element_equations
GEBT.insert_element_residual!
GEBT.element_residual!
GEBT.element_jacobian_equations
GEBT.insert_element_jacobian!
GEBT.element_jacobian!
GEBT.element_mass_matrix_properties
GEBT.element_mass_matrix_equations
GEBT.insert_element_mass_matrix!
GEBT.element_mass_matrix!
GEBT.ElementState
```

### Loads

```@docs
GEBT.PrescribedConditions
GEBT.DistributedLoads
```

### System

```@docs
GEBT.Assembly
GEBT.curve_triad
GEBT.curve_coordinates
GEBT.System
GEBT.point_connections
GEBT.system_indices
GEBT.system_residual!
GEBT.system_jacobian!
GEBT.system_mass_matrix!
```

## Index

```@index
```
