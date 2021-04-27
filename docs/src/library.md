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
combine_loads
```

### Defining Prescribed Conditions

```@docs
PrescribedConditions()
```

### Pre-Initializing Memory for an Analysis

```@docs
System(assembly, points, static)
system_state
reset_state!
set_state!
set_element_deflections!
set_element_rotations!
set_element_forces!
set_element_moments!
set_element_linear_momenta!
set_element_angular_momenta!
set_point_deflections!
set_point_rotations!
set_point_forces!
set_point_moments!
```

### Performing an Analysis

```@docs
static_analysis
static_analysis!
steady_state_analysis
steady_state_analysis!
eigenvalue_analysis
eigenvalue_analysis!
initial_condition_analysis
initial_condition_analysis!
time_domain_analysis
time_domain_analysis!
```

### Post-Processing

```@docs
AssemblyState(system, assembly)
AssemblyState
PointState
extract_element_state
extract_element_states
extract_element_states!
extract_point_state
extract_point_states
extract_point_states!
left_eigenvectors
correlate_eigenmodes
wiener_milenkovic
angular_velocities
rotate
rotate!
translate
translate!
deform_cross_section
deform_cross_section!
cross_section_velocities
cross_section_velocities!
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
```

### Elements

```@docs
GXBeam.Element
GXBeam.element_strain
GXBeam.element_curvature
GXBeam.element_linear_velocity
GXBeam.element_angular_velocity
GXBeam.element_properties
GXBeam.element_dynamic_properties
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
GXBeam.get_sparsity
GXBeam.system_indices
GXBeam.system_residual!
GXBeam.system_jacobian!
GXBeam.system_mass_matrix!
```

## Index

```@index
```
