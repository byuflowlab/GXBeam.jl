
# Theory

## Structural Damping

When stiffness proportional structural damping is used, the relationship between stresses and strains takes the following form:

```math
\begin{bmatrix} F \\  M \end{bmatrix} = S \begin{bmatrix} \gamma \\ \kappa \end{bmatrix} + \mu S \begin{bmatrix} \dot{\gamma} \\ \dot{\kappa} \end{bmatrix}  \\
```

where ``\mu`` is a diagonal matrix with damping coefficients.

Rearranging this equations yields a new expression for the strains and curvatures when stiffness proportional structural damping is assumed.

```math
\begin{bmatrix} \gamma \\ \kappa \end{bmatrix} = S^{-1} \begin{bmatrix} F \\  M \end{bmatrix} - \mu \begin{bmatrix} \dot{\gamma} \\ \dot{\kappa} \end{bmatrix}
```

### Strain Rates

The strain rates may be expressed as a function of our state variables using the strain compatability relationship.  We first rearrange the strain compatability equation to yield the following expression for the element strain.

```math
\gamma_i = C^{ba} C \frac{\Delta u}{\Delta L} + C^{ba} C C^{ab} e_1 - e_1 
```

We then differentiate this equation with respect to time to obtain an expression for the strain rates:

```math
\dot{\gamma}_i = C^{ba} C \frac{\Delta \dot{u}}{\Delta L} + C^{ba} \dot{C} \frac{\Delta u}{\Delta L} + C^{ba} \dot{C} C^{ab} e_1 
```

The derivative of the rotation matrix ``\dot{C}`` can be defined in terms of angular velocity as ``\dot{C} = -\widetilde{(\Omega-\omega)} C``.  The linear displacement rates ``\dot{u}`` can be defined in terms of linear velocity as ``\dot{u} = V - v - \tilde{\omega} u``.  Using these two expressions allows us to write the strain rates in terms of the state variables.

```math
\dot{\gamma}_i = C^{ba} C \frac{\Delta \dot{u}}{\Delta L} - C^{ba} \widetilde{(\Omega-\omega)} C \frac{\Delta u}{\Delta L} -  C^{ba} \widetilde{(\Omega-\omega)} C C^{ab} e_1 
```

### Curvature Rates

The curvature rates may be expressed as a function of our state variables using the curvature compatability relationship.  We first rearrange the curvature compatability equation to yield the following expression for the element curvature.

```math
\kappa_i = C^{ba} Q \frac{\Delta \theta}{\Delta L}
```

Differentiating this equation with respect to time yields:

```math
\dot{\kappa}_i = C^{ba} Q \frac{\Delta \dot{\theta}}{\Delta L} + C^{ba} \dot{Q} \frac{\Delta \theta}{\Delta L}
```

We can expand ``\dot{Q}`` using its partial derivatives 

```math
\dot{Q} = Q_{\theta_1} \dot{\theta_1} + Q_{\theta_2} \dot{\theta_2} + Q_{\theta_3} \dot{\theta_3}
```

so that

```math
\dot{Q} \Delta \theta = 
(Q_{\theta_1} \Delta \theta) \dot{\theta_1} + 
(Q_{\theta_2} \Delta \theta) \dot{\theta_2} + 
(Q_{\theta_3} \Delta \theta) \dot{\theta_3} 
\equiv \Delta Q \dot{\theta}
```

Our expression for the curvature rates is then

```math
\dot{\kappa}_i = C^{ba} Q \frac{\Delta \dot{\theta}}{\Delta L} + C^{ba} \Delta Q \frac{\dot{\theta}}{\Delta L}
```

We can then use the following expression for the angular displacement rates to express the curvature rates as a function of our state variables.

```math
\dot{\theta} = Q^{-1} C (\Omega - \omega)
```