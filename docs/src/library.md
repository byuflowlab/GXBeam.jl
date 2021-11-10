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

### Defining Point Masses
```@docs
PointMass
PointMass(m, p, J)
combine_masses
```

### Defining Distributed Loads
```@docs
DistributedLoads
DistributedLoads(assembly, ibeam)
combine_loads
```

### Defining Prescribed Conditions

```@docs
PrescribedConditions
PrescribedConditions()
```

### Pre-Initializing Memory for an Analysis

```@docs
System(assembly, static)
system_state
reset_state!
set_state!
set_element_deflection!
set_element_rotation!
set_element_forces!
set_element_moments!
set_element_linear_velocity!
set_element_angular_velocity!
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
AssemblyState
PointState
ElementState
AssemblyState(system, assembly)
extract_element_state
extract_element_states
extract_element_states!
extract_point_state
extract_point_states
extract_point_states!
left_eigenvectors
correlate_eigenmodes
wiener_milenkovic
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
GXBeam.get_C_t_θ
GXBeam.get_C_t_θdot
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
GXBeam.element_strain
GXBeam.element_curvature
GXBeam.element_linear_momentum
GXBeam.element_angular_momentum
GXBeam.static_element_state_variables
GXBeam.dynamic_element_state_variables
GXBeam.static_element_equations
GXBeam.steady_state_element_equations
GXBeam.dynamic_element_equations
GXBeam.static_insert_element_residual!
GXBeam.dynamic_insert_element_residual!
GXBeam.static_element_residual!
GXBeam.steady_state_element_residual!
GXBeam.initial_condition_element_residual!
GXBeam.dynamic_element_residual!
GXBeam.newmark_element_residual!
GXBeam.static_element_jacobian_equations
GXBeam.initial_condition_element_jacobian_equations
GXBeam.steady_state_element_jacobian_equations
GXBeam.newmark_element_jacobian_equations
GXBeam.dynamic_element_jacobian_equations
GXBeam.static_insert_element_jacobian!
GXBeam.initial_condition_insert_element_jacobian!
GXBeam.dynamic_insert_element_jacobian!
GXBeam.static_element_jacobian!
GXBeam.steady_state_element_jacobian!
GXBeam.initial_condition_element_jacobian!
GXBeam.newmark_element_jacobian!
GXBeam.dynamic_element_jacobian!
GXBeam.element_mass_matrix_equations
GXBeam.insert_element_mass_matrix!
GXBeam.element_mass_matrix!
```

### Loads

```@docs
GXBeam.acceleration_loads
```

### System

```@docs
GXBeam.Assembly
GXBeam.curve_triad
GXBeam.curve_coordinates
GXBeam.System
GXBeam.static_system_residual!
GXBeam.initial_condition_system_residual!
GXBeam.steady_state_system_residual!
GXBeam.newmark_system_residual!
GXBeam.dynamic_system_residual!
GXBeam.static_system_jacobian!
GXBeam.steady_state_system_jacobian!
GXBeam.initial_condition_system_jacobian!
GXBeam.newmark_system_jacobian!
GXBeam.dynamic_system_jacobian!
GXBeam.get_sparsity
GXBeam.system_indices
GXBeam.system_mass_matrix!
```

## Index

```@index
```
