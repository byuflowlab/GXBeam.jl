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
Assembly(points, start, stop)
```

### Section Properties

```@docs
Material
Node
MeshElement
Layer
afmesh
initialize_cache
compliance_matrix
mass_matrix
plotmesh
strain_recovery
plotsoln
tsai_wu
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
AbstractSystem
StaticSystem
StaticSystem(assembly)
DynamicSystem
DynamicSystem(assembly)
ExpandedSystem
ExpandedSystem(assembly)
reset_state!
set_state!
set_rate!
set_linear_displacement!
set_angular_displacement!
set_external_forces!
set_external_moments!
set_linear_velocity!
set_angular_velocity!
set_internal_forces!
set_internal_moments!
set_start_forces!
set_start_moments!
set_end_forces!
set_end_moments!
set_point_linear_velocity!
set_point_angular_velocity!
set_element_linear_velocity!
set_element_angular_velocity!
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
extract_point_state
extract_point_states
linearize!
solve_eigensystem
left_eigenvectors
correlate_eigenmodes
wiener_milenkovic
rotate
rotate!
translate
translate!
deform_cross_section
deform_cross_section!
write_vtk
```

## Private API

### Math

```@docs
GXBeam.tilde
GXBeam.transform_properties
GXBeam.rotation_parameter_scaling
GXBeam.get_C
GXBeam.get_C_θ
GXBeam.get_Q
GXBeam.get_Q_θ
GXBeam.get_Qinv
GXBeam.get_Qinv_θ
GXBeam.get_ΔQ
GXBeam.get_ΔQ_θ
GXBeam.mul3
GXBeam.gauss_quadrature
```

### Body Frame

```@docs
GXBeam.update_body_acceleration_indices!
GXBeam.body_accelerations
```

### Points

```@docs
GXBeam.point_loads
GXBeam.point_load_jacobians
GXBeam.point_displacement
GXBeam.point_displacement_jacobians
GXBeam.point_displacement_rates
GXBeam.point_velocity_rates
GXBeam.point_velocities
GXBeam.initial_point_displacement
GXBeam.initial_point_velocity_rates
GXBeam.initial_point_displacement_jacobian
GXBeam.initial_point_velocity_rate_jacobian
GXBeam.static_point_properties
GXBeam.steady_point_properties
GXBeam.initial_point_properties
GXBeam.newmark_point_properties
GXBeam.dynamic_point_properties
GXBeam.expanded_steady_point_properties
GXBeam.expanded_dynamic_point_properties
GXBeam.point_velocity_residuals
GXBeam.expanded_point_velocity_residuals
GXBeam.static_point_resultants
GXBeam.dynamic_point_resultants
GXBeam.expanded_point_resultants
GXBeam.static_point_residual!
GXBeam.steady_point_residual!
GXBeam.initial_point_residual!
GXBeam.newmark_point_residual!
GXBeam.dynamic_point_residual!
GXBeam.expanded_steady_point_residual!
GXBeam.expanded_dynamic_point_residual!
GXBeam.static_point_jacobian_properties
GXBeam.steady_point_jacobian_properties
GXBeam.initial_point_jacobian_properties
GXBeam.newmark_point_jacobian_properties
GXBeam.dynamic_point_jacobian_properties
GXBeam.expanded_steady_point_jacobian_properties
GXBeam.expanded_dynamic_point_jacobian_properties
GXBeam.mass_matrix_point_jacobian_properties
GXBeam.expanded_mass_matrix_point_jacobian_properties
GXBeam.static_point_resultant_jacobians
GXBeam.steady_point_resultant_jacobians
GXBeam.initial_point_resultant_jacobians
GXBeam.dynamic_point_resultant_jacobians
GXBeam.expanded_steady_point_resultant_jacobians
GXBeam.expanded_dynamic_point_resultant_jacobians
GXBeam.mass_matrix_point_resultant_jacobians
GXBeam.point_velocity_jacobians
GXBeam.initial_point_velocity_jacobians
GXBeam.newmark_point_velocity_jacobians
GXBeam.dynamic_point_velocity_jacobians
GXBeam.expanded_point_velocity_jacobians
GXBeam.mass_matrix_point_velocity_jacobians
GXBeam.insert_static_point_jacobians!
GXBeam.insert_steady_point_jacobians!
GXBeam.insert_initial_point_jacobians!
GXBeam.insert_dynamic_point_jacobians!
GXBeam.insert_expanded_steady_point_jacobians!
GXBeam.insert_expanded_dynamic_point_jacobians!
GXBeam.insert_mass_matrix_point_jacobians!
GXBeam.insert_expanded_mass_matrix_point_jacobians!
GXBeam.static_point_jacobian!
GXBeam.steady_point_jacobian!
GXBeam.initial_point_jacobian!
GXBeam.newmark_point_jacobian!
GXBeam.dynamic_point_jacobian!
GXBeam.expanded_steady_point_jacobian!
GXBeam.expanded_dynamic_point_jacobian!
GXBeam.mass_matrix_point_jacobian!
GXBeam.expanded_mass_matrix_point_jacobian!
```

### Elements

```@docs
GXBeam.Element
GXBeam.element_loads
GXBeam.expanded_element_loads
GXBeam.expanded_element_velocities
GXBeam.static_element_properties
GXBeam.steady_element_properties
GXBeam.initial_element_properties
GXBeam.newmark_element_properties
GXBeam.dynamic_element_properties
GXBeam.expanded_steady_element_properties
GXBeam.expanded_dynamic_element_properties
GXBeam.compatibility_residuals
GXBeam.expanded_element_velocity_residuals
GXBeam.expanded_element_equilibrium_residuals
GXBeam.static_element_resultants
GXBeam.dynamic_element_resultants
GXBeam.expanded_element_resultants
GXBeam.insert_element_residuals!
GXBeam.insert_expanded_element_residuals!
GXBeam.static_element_residual!
GXBeam.steady_element_residual!
GXBeam.initial_element_residual!
GXBeam.newmark_element_residual!
GXBeam.dynamic_element_residual!
GXBeam.expanded_steady_element_residual!
GXBeam.expanded_dynamic_element_residual!
GXBeam.static_element_jacobian_properties
GXBeam.steady_element_jacobian_properties
GXBeam.initial_element_jacobian_properties
GXBeam.newmark_element_jacobian_properties
GXBeam.dynamic_element_jacobian_properties
GXBeam.expanded_steady_element_jacobian_properties
GXBeam.expanded_dynamic_element_jacobian_properties
GXBeam.mass_matrix_element_jacobian_properties
GXBeam.expanded_mass_matrix_element_jacobian_properties
GXBeam.expanded_element_velocity_jacobians
GXBeam.expanded_steady_element_equilibrium_jacobians
GXBeam.expanded_dynamic_element_equilibrium_jacobians
GXBeam.expanded_mass_matrix_element_equilibrium_jacobians
GXBeam.static_element_resultant_jacobians
GXBeam.steady_element_resultant_jacobians
GXBeam.initial_element_resultant_jacobians
GXBeam.dynamic_element_resultant_jacobians
GXBeam.expanded_element_resultant_jacobians
GXBeam.mass_matrix_element_resultant_jacobians
GXBeam.static_element_jacobian!
GXBeam.steady_element_jacobian!
GXBeam.initial_element_jacobian!
GXBeam.newmark_element_jacobian!
GXBeam.dynamic_element_jacobian!
GXBeam.expanded_steady_element_jacobian!
GXBeam.expanded_dynamic_element_jacobian!
GXBeam.mass_matrix_element_jacobian!
GXBeam.expanded_mass_matrix_element_jacobian!
```

### System

```@docs
GXBeam.SystemIndices
GXBeam.default_force_scaling
GXBeam.curve_triad
GXBeam.curve_coordinates
GXBeam.set_initial_state!
GXBeam.two_dimensional_residual!
GXBeam.static_system_residual!
GXBeam.initial_system_residual!
GXBeam.steady_system_residual!
GXBeam.newmark_system_residual!
GXBeam.dynamic_system_residual!
GXBeam.expanded_steady_system_residual!
GXBeam.expanded_dynamic_system_residual!
GXBeam.two_dimensional_jacobian!
GXBeam.static_system_jacobian!
GXBeam.steady_system_jacobian!
GXBeam.initial_system_jacobian!
GXBeam.newmark_system_jacobian!
GXBeam.dynamic_system_jacobian!
GXBeam.expanded_steady_system_jacobian!
GXBeam.expanded_dynamic_system_jacobian!
GXBeam.system_mass_matrix!
GXBeam.expanded_system_mass_matrix
GXBeam.expanded_system_mass_matrix!
```

### Section
```@docs
GXBeam.SectionCache
GXBeam.area_and_centroid_of_element
GXBeam.redistribute_thickness
GXBeam.combine_halfs
GXBeam.nodes_half
GXBeam.rotate_ply_to_element
GXBeam.insertpoint
GXBeam.tangential
GXBeam.element_orientation
GXBeam.element_submatrix
GXBeam.addwebs
GXBeam.web_intersections
GXBeam.parseairfoil
GXBeam.node2idx
GXBeam.node2idx!
GXBeam.reorder
GXBeam.find_inner_surface
GXBeam.elementQ
GXBeam.resample
GXBeam.stiffness
GXBeam.preprocess_layers
GXBeam.linearsolve
GXBeam.rotate_element_to_beam
GXBeam.element_integrand
GXBeam.te_inner_intersection
```

## Index

```@index
```
