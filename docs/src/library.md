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
Assembly
```

### Defining Prescribed Conditions

```@docs
PrescribedConditions
DistributedLoads
TimeFunction
```

### Pre-Initializing Memory for an Analysis

```@docs
System
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
AssemblyState
left_eigenvectors
wiener_milenkovic
write_vtk
```

## Private API

### Math

```@docs
GEBT.tilde
GEBT.get_C
GEBT.get_C_t
GEBT.get_C_c
GEBT.get_C_cdot
GEBT.get_Q
GEBT.get_Q_c
GEBT.get_Qinv
GEBT.get_Qinv_c
GEBT.mul3
```

### Points

```@docs
GEBT.point_variables
GEBT.insert_point_residual!
GEBT.point_residual!
GEBT.point_follower_jacobians
GEBT.insert_point_jacobians!
GEBT.point_jacobian!
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
```

### Loads

```@docs
GEBT.quadrature_points
GEBT.integrate_element_loads
GEBT.follower_load_jacobians
```

### System

```@docs
GEBT.point_connections
GEBT.system_indices
GEBT.system_residual!
GEBT.system_jacobian!
GEBT.system_mass_matrix!
```

## Index

```@index
```
