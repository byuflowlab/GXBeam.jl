# Theory

To understand the underlying theory for this package, we recommend you review the provided [references](@ref References).  This page describes some of the additional features introduced by this package.

```@contents
Pages = ["theory.md"]
Depth = 3
```

## Constant Mass Matrix Form

The governing equations associated with this package are a set of fully implicit differential algebraic equations.  While these equations may be solved as is, it is often more efficient to reformulate the governing equations into constant mass matrix form.

$$
M \dot{x} = f(x, p, t)
$$

To do so, we introduce the new algebraic variables $F_i^-$, $M_i^-$, $F_i^+$, $M_i^+$ which are defined using the following constraints

$$
f_{u_i}^- - F_i^- = 0 \quad 
f_{\psi_i}^- - M_i^- = 0 \quad 
f_{u_i}^+ - F_i^+ = 0 \quad 
f_{\psi_N}^+ - M_i^+ = 0 \quad 
$$

If beam element $i$ ends and beam elements $j$ and $k$ start at a connection point $C$, the corresponding equilibrium equations for the connection point is then
$$
F_i^+ + F_j^- + F_k^- - F_C^* = 0 \quad
M_i^+ + M_j^- + M_k^- - M_C^* = 0
$$
where $F_C^*$ and $M_C^*$ are the concentrated forces and moments applied at the connection point, respectively.  Similar equilibrium equations may be constructed for any given node.

The constraints associated with any given beam element are then

$$
f_{u_i}^{-} - F_i^- = 0 \\
f_{\psi_i}^{-} - M_i^- = 0 \\ 
f_{F_i}^{-} - \hat{u}_i = 0 \\
f_{M_i}^{-} - \hat{\theta}_i = 0 \\
f_{u_i}^{+} - F_{i}^+ = 0 \\
f_{\psi_i}^{+} - M_{i}^+ = 0 \\ 
f_{F_i}^{+} + \hat{u}_{i+1} = 0 \\ 
f_{M_i}^{+} + \hat{\theta}_{i+1} = 0 \\
f_P = 0 \\
f_H = 0 \\
$$
where 
$$
f_{u_i}^{\mp} = 
\mp C^T C^{ab} F_i
- \bar{f}_i^{\mp} 
+ \frac{\Delta L}{2} 
[ \tilde{\omega} C^T C^{ab} P_i + \dot{\overline{C^T C^{ab} P_i}}] 
\\
f_{\psi_i}^{\mp} = 
\mp C^T C^{ab} M_i 
- \bar{m}_i^{\mp} 
+ \frac{\Delta L}{2} 
[ 
  \tilde{\omega} C^T C^{ab} H_i 
  + \dot{\overline{C^T C^{ab} H_i}}
  + C^T C^{ab} (\tilde{V}_i P_i - (\tilde{e}_i + \tilde{\gamma}_i) F_i) 
] \\
f_{F_i}^{\mp} = \pm u_i - \frac{\Delta L}{2} [C^T C^{ab} (e_i + \gamma_i) - C^{ab} e_i] \\
f_{M_i}^{\mp} = \pm \theta_i - \frac{\Delta L}{2} Q_a^{-1} C^{ab} \kappa_i \\
f_{P_i} = C^T C^{ab} V_i - v_i - \tilde{\omega_a} u_i - \dot{u}_i \\
f_{H_i} = \Omega_i - C^{ba} C \omega_a - C^{ba} Q_a \dot{\theta} \\
$$
and
$$
\gamma = S^{11}_i F_i + S^{12}_i M_i \\
\kappa = S^{21}_i F_i + S^{22}_i M_i \\
P = M^{11}_i V_i + M^{12}_i \Omega_i \\
H = M^{21}_i V_i + M^{22}_i \Omega_i \\
$$

We can now perform some algebraic manipulations to be able to represent the constraints for any given beam element in constant mass matrix form.

$$
-\frac{\Delta L}{2} M^{11}_i \dot{V}_i - \frac{\Delta L}{2} M^{12}_i \dot{\Omega}_i = \mp F_i + C^{ba} C \left(
- \bar{f}_i^{\mp} 
+ \frac{\Delta L}{2} 
[ \tilde{\omega} C^T C^{ab} P_i + \dot{C}^T C^{ab} P_i] - F_i^\mp \right)
\\
-\frac{\Delta L}{2} M^{21}_i \dot{V}_i - \frac{\Delta L}{2} M^{22}_i \dot{\Omega}_i = \mp M_i + C^{ba} C \left( 
 
- \bar{m}_i^{\mp} 
+ \frac{\Delta L}{2} 
[ 
  \tilde{\omega} C^T C^{ab} H_i 
  + \dot{C}^T C^{ab} H_i
  + C^T C^{ab} (\tilde{V}_i P_i - (\tilde{e}_i + \tilde{\gamma}_i) F_i) 
] - M_i^\mp
\right) \\
0 = \pm u_i - \frac{\Delta L}{2} [C^T C^{ab} (e_i + \gamma_i) - C^{ab} e_i] \mp u_i^\mp \\
0 = \pm \theta_i - \frac{\Delta L}{2} Q_a^{-1} C^{ab} \kappa_i \mp \theta_i^\mp \\
\dot{u}_i = C^T C^{ab} V_i - v_i - \tilde{\omega_a} u_i  \\
\dot{\theta_i} = Q_a^{-1} C^{ab} \left( \Omega_i - C^{ba} C \omega_a \right) \\
$$