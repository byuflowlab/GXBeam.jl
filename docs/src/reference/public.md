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
initialize_system!
step_system!
take_step
simulate
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