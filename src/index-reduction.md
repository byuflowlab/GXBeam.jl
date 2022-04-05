# Index Reduction

To simulate the behavior of a given set of beam elements and nodes, this package constructs and solves a set of differential algebraic equations.  While these equations may be solved as is, it can be sometimes be more convenient or efficient to express a given set of differential algebraic equations as ordinary differential equations.  We would therefore like to be able to represent the differential algebraic equations modeled in this package as ordinary differential equations.  

We first consider the following equations for a single cantilevered beam element.

$$
f_{u_1}^{-} - F_1^* = 0 \\
f_{\psi_1}^{-} - M_1^* = 0 \\ 
f_{V_1} = 0 \\
f_{\Omega_1} = 0 \\
\hat{u}_2 - \hat{u}_1 - \Delta u = 0 \\
\hat{\theta}_2 - \hat{\theta}_1 - \Delta \theta = 0 \\
f_{u_1}^{+} - F_{N+1}^* = 0 \\
f_{\psi_1}^{+} - M_{N+1}^* = 0 \\ 
f_{V_2} = 0 \\
f_{\Omega_2} = 0 \\
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
\Delta u = \Delta L \left(C^T C^{ab} (e_i + \gamma_i) - C^{ab} e_i\right) \\
\Delta \theta = \Delta L \left( Q_a^{-1} C^{ab} \kappa_i \right) \\
f_{V_i} = C^T C^{ab} V_i - v_i - \tilde{\omega_a} u_i - \dot{u}_i \\
f_{\Omega_i} = \Omega_i - C^{ba} C \omega_a - C^{ba} Q_a \dot{\theta} \\
\gamma = S^{11}_1 F_i + S^{12}_1 M_i \\
\kappa = S^{21}_i F_i + S^{22}_1 M_i \\
P = M^{11}_i V_i + M^{12}_i \Omega_i \\
H = M^{21}_1 V_i + M^{22}_1 \Omega_i \\
$$

The last two governing equations can be manipulated to define explicit expressions for $\dot{u}$ and $\dot{\theta}$.

$$
\dot{u}_i = C^T C^{ab} V_i - v_i - \tilde{\omega_a} u_i \\
\dot{\theta} = Q_a^{-1} C^{ab} \left(\Omega_i - C^{ba} C \omega_a \right)  \\
$$

Differentiating the last governing equation again with respect to time yields an explicit expression for $\ddot{\theta}$.

$$
\ddot{\theta} 
= Q_a^{-1} C^{ab} \left(
\dot{\Omega_1} 
- C^{ba} \dot{C} \omega_a 
- C^{ba} C \dot{\omega}_a
- C^{ba} \dot{Q_a} \dot{\theta}
\right)   
$$

We would like to introduce state rate terms for each of our algebraic variables. We can do so by differentiating the first eight governing equations.

$$
% eqn 1.
- \dot{C}^T C^{ab} F_1
- C^T C^{ab} \dot{F}_1
- \dot{\bar{f}_1^{-}} 
+ \frac{\Delta L}{2} 
\left( 
  \dot{\tilde{\omega}} C^T C^{ab} P_1
  + \tilde{\omega} \dot{C}^T C^{ab} P_1
  + \tilde{\omega} C^T C^{ab} \dot{P}_1
  + \ddot{C}^T C^{ab} P_1 
  + 2 \dot{C}^T C^{ab} \dot{P}_1 
  + C^T C^{ab} \ddot{P}_1 
\right) 
- \dot{F_1^*} 
= 0 \\
% eqn 2.
- \dot{C}^T C^{ab} M_1
- C^T C^{ab} \dot{M}_1
- \dot{\bar{m}_1^{-}}
+ \frac{\Delta L}{2} 
\left(
  \dot{\tilde{\omega}} C^T C^{ab} H_1
  + \tilde{\omega} \dot{C}^T C^{ab} H_1
  + \tilde{\omega} C^T C^{ab} \dot{H}_1
  + \ddot{C}^T C^{ab} H_1
  + 2\dot{C}^T C^{ab} \dot{H}_1
  + C^T C^{ab} \ddot{H}_1
  + \dot{C}^T C^{ab} 
  ( 
    \tilde{V}_1 P_1 
    - (\tilde{e}_1 + \tilde{\gamma}_1) F_1
  )
  + C^T C^{ab} 
  (
    \dot{\tilde{V}}_1 P_1 
    + \tilde{V}_1 \dot{P}_1 
    - \dot{\tilde{\gamma}}_1 F_1
    - (\tilde{e}_1 + \tilde{\gamma}_1) \dot{F}_1
  )
\right) 
- \dot{M_1^*} 
= 0 \\
% eqn 3.
\dot{\hat{u}_2} - \dot{\hat{u}_1} - \Delta L
\left(
  \dot{C}^T C^{ab} (e_1 + \gamma_1)
  + C^T C^{ab} \dot{\gamma}_1
\right) 
= 0 \\
% eqn 4.
\dot{\hat{\theta}_2} 
- \dot{\hat{\theta}_1} 
- \Delta L
\left(
  - Q_a^{-1} \dot{Q_a} Q_a^{-1} C^{ab} \kappa_1
  + Q_a^{-1} C^{ab} \dot{\kappa}_1
\right)
 
= 0 \\
% eqn 5.
\dot{C}^T C^{ab} F_1
+ C^T C^{ab} \dot{F}_1
- \dot{\bar{f}_1^{+}} 
+ \frac{\Delta L}{2} 
\left( 
  \dot{\tilde{\omega}} C^T C^{ab} P_1
  + \tilde{\omega} \dot{C}^T C^{ab} P_1
  + \tilde{\omega} C^T C^{ab} \dot{P}_1
  + \ddot{C}^T C^{ab} P_1 
  + 2 \dot{C}^T C^{ab} \dot{P}_1 
  + C^T C^{ab} \ddot{P}_1 
\right) 
- \dot{F_1^*} 
= 0 \\
% eqn 6.
\dot{C}^T C^{ab} M_1
C^T C^{ab} \dot{M}_1
- \dot{\bar{m}_1^{+}}
+ \frac{\Delta L}{2} 
\left(
  \dot{\tilde{\omega}} C^T C^{ab} H_1
  + \tilde{\omega} \dot{C}^T C^{ab} H_1
  + \tilde{\omega} C^T C^{ab} \dot{H}_1
  + \ddot{C}^T C^{ab} H_1
  + 2\dot{C}^T C^{ab} \dot{H}_1
  + C^T C^{ab} \ddot{H}_1
  + \dot{C}^T C^{ab} 
  ( 
    \tilde{V}_1 P_1 
    - (\tilde{e}_1 + \tilde{\gamma}_1) F_1
  )
  + C^T C^{ab} 
  (
    \dot{\tilde{V}}_1 P_1 
    + \tilde{V}_1 \dot{P}_1 
    - \dot{\tilde{\gamma}}_1 F_1
    - (\tilde{e}_1 + \tilde{\gamma}_1) \dot{F}_1
  )
\right) 
- \dot{M_1^*} 
= 0 \\
$$

At this point, we define the following state variable vector for our new system of equations.

$$
x = \begin{bmatrix} 
 \dot{F^*_1} & \dot{M^*_1} & 
 F_1 & M_1 & 
 V_1 & \Omega_1 & 
 \dot{V}_1 & \dot{\Omega}_1 &
 \dot{\hat{u}}_2 & \dot{\hat{\theta}}_2  &  
\end{bmatrix}^T
$$

The governing equations for our new system of equations are

$$
% eqn 1.
- \dot{C}^T C^{ab} F_1
- C^T C^{ab} \dot{F}_1
- \dot{\bar{f}_1^{-}} 
+ \frac{\Delta L}{2} 
\left( 
  \dot{\tilde{\omega}} C^T C^{ab} P_1
  + \tilde{\omega} \dot{C}^T C^{ab} P_1
  + \tilde{\omega} C^T C^{ab} \dot{P}_1
  + \ddot{C}^T C^{ab} P_1 
  + 2 \dot{C}^T C^{ab} \dot{P}_1 
  + C^T C^{ab} \ddot{P}_1 
\right) 
- \dot{F_1^*} 
= 0 \\
% eqn 2.
- \dot{C}^T C^{ab} M_1
- C^T C^{ab} \dot{M}_1
- \dot{\bar{m}_1^{-}}
+ \frac{\Delta L}{2} 
\left(
  \dot{\tilde{\omega}} C^T C^{ab} H_1
  + \tilde{\omega} \dot{C}^T C^{ab} H_1
  + \tilde{\omega} C^T C^{ab} \dot{H}_1
  + \ddot{C}^T C^{ab} H_1
  + 2\dot{C}^T C^{ab} \dot{H}_1
  + C^T C^{ab} \ddot{H}_1
  + \dot{C}^T C^{ab} 
  ( 
    \tilde{V}_1 P_1 
    - (\tilde{e}_1 + \tilde{\gamma}_1) F_1
  )
  + C^T C^{ab} 
  (
    \dot{\tilde{V}}_1 P_1 
    + \tilde{V}_1 \dot{P}_1 
    - \dot{\tilde{\gamma}}_1 F_1
    - (\tilde{e}_1 + \tilde{\gamma}_1) \dot{F}_1
  )
\right) 
- \dot{M_1^*} 
= 0 \\
% eqn 3.
\dot{u}_1 
- \frac{\Delta L}{2} 
\left(
  \dot{C}^T C^{ab} (e_1 + \gamma_1)
  + C^T C^{ab} \dot{\gamma}_1
\right) 
- \dot{\hat{u}_1} 
= 0 \\
% eqn 4.
\dot{\theta}_1
- \frac{\Delta L}{2} 
\left(
  - Q_a^{-1} \dot{Q_a} Q_a^{-1} C^{ab} \kappa_1
  + Q_a^{-1} C^{ab} \dot{\kappa}_1
\right)
- \dot{\hat{\theta}_1} 
= 0 \\
% eqn 5.
\dot{C}^T C^{ab} F_1
+ C^T C^{ab} \dot{F}_1
- \dot{\bar{f}_1^{+}} 
+ \frac{\Delta L}{2} 
\left( 
  \dot{\tilde{\omega}} C^T C^{ab} P_1
  + \tilde{\omega} \dot{C}^T C^{ab} P_1
  + \tilde{\omega} C^T C^{ab} \dot{P}_1
  + \ddot{C}^T C^{ab} P_1 
  + 2 \dot{C}^T C^{ab} \dot{P}_1 
  + C^T C^{ab} \ddot{P}_1 
\right) 
- \dot{F_1^*} 
= 0 \\
% eqn 6.
\dot{C}^T C^{ab} M_1
C^T C^{ab} \dot{M}_1
- \dot{\bar{m}_1^{+}}
+ \frac{\Delta L}{2} 
\left(
  \dot{\tilde{\omega}} C^T C^{ab} H_1
  + \tilde{\omega} \dot{C}^T C^{ab} H_1
  + \tilde{\omega} C^T C^{ab} \dot{H}_1
  + \ddot{C}^T C^{ab} H_1
  + 2\dot{C}^T C^{ab} \dot{H}_1
  + C^T C^{ab} \ddot{H}_1
  + \dot{C}^T C^{ab} 
  ( 
    \tilde{V}_1 P_1 
    - (\tilde{e}_1 + \tilde{\gamma}_1) F_1
  )
  + C^T C^{ab} 
  (
    \dot{\tilde{V}}_1 P_1 
    + \tilde{V}_1 \dot{P}_1 
    - \dot{\tilde{\gamma}}_1 F_1
    - (\tilde{e}_1 + \tilde{\gamma}_1) \dot{F}_1
  )
\right) 
- \dot{M_1^*} 
= 0 \\
% eqn 7. 
- \dot{u}_1 
- \frac{\Delta L}{2} 
\left(
  \dot{C}^T C^{ab} (e_1 + \gamma_1)
  + C^T C^{ab} \dot{\gamma}_1
\right) 
- \dot{\hat{u}_1} 
= 0 \\
% eqn 8.
- \dot{\theta}_1
- \frac{\Delta L}{2} 
\left(
  - Q_a^{-1} \dot{Q_a} Q_a^{-1} C^{ab} \kappa_1
  + Q_a^{-1} C^{ab} \dot{\kappa}_1
\right)
- \dot{\hat{\theta}_1} 
= 0 \\
% eqn 9.
\dot{u}_1 = C^T C^{ab} V_1 - v_1 - \tilde{\omega_a} u_1 \\
% eqn 10.
\dot{\theta} = Q_a^{-1} C^{ab} \left(\Omega_1 - C^{ba} C \omega_a \right)  \\
% eqn 11.
\dot{V} = \dot{V} \\
% eqn 12.
\dot{\Omega} = \dot{\Omega}
$$

We would like to be able to construct an explicit expression for the state variable rates.  To do so, we isolate terms with the state rates $\dot{u}$, $\dot{\theta}$, $\dot{F}$, $\dot{M}$, $\dot{V}$, $\dot{\Omega}$, $\dot{F_1}$, $\dot{M_1}$, $\dot{\hat{u}}_2$, or $\dot{\hat{\theta}}_2$ on the left hand side.  Note that since we have explicit expressions for $\dot{\theta}$ and $\ddot{\theta}$, these terms can be treated as lower order terms which can be solved explicitly. 

$$
% eqn 1.
- C^T C^{ab} \dot{F}_1
+ \frac{\Delta L}{2} C^T C^{ab} M^{11}_1 \ddot{V}_1
+ \frac{\Delta L}{2} C^T C^{ab} M^{12}_1 \ddot{\Omega}_1
- \dot{F_1^*} 
= 
\dot{C}^T C^{ab} F_1
+ \dot{\bar{f}_1^{-}} 
- \frac{\Delta L}{2} 
\left( 
  \dot{\tilde{\omega}} C^T C^{ab} P_1
  + \tilde{\omega} \dot{C}^T C^{ab} P_1
  + \tilde{\omega} C^T C^{ab} \dot{P}_1
  + \ddot{C}^T C^{ab} P_1 
  + 2 \dot{C}^T C^{ab} \dot{P}_1
\right) 
\\
% eqn 2.
+ \frac{\Delta L}{2} C^T C^{ab} \left(\tilde{F_1} S^{11}_1 - (\tilde{e}_1 + \tilde{\gamma}_1) \right) \dot{F}_1 
+ C^T C^{ab} \left(\frac{\Delta L}{2} \tilde{F_1} S^{12}_1
- I \right) \dot{M}_1
+ \frac{\Delta L}{2} C^T C^{ab} M^{21}_1 \ddot{V}_1
+ \frac{\Delta L}{2} C^T C^{ab} M^{22}_1 \ddot{\Omega}_1
- \dot{M_1^*} 
= 
\dot{C}^T C^{ab} M_1
+ \dot{\bar{m}_1^{-}}
- \frac{\Delta L}{2} 
\left(
  \dot{\tilde{\omega}} C^T C^{ab} H_1
  + \tilde{\omega} \dot{C}^T C^{ab} H_1
  + \tilde{\omega} C^T C^{ab} \dot{H}_1
  + \ddot{C}^T C^{ab} H_1
  + 2\dot{C}^T C^{ab} \dot{H}_1
  + \dot{C}^T C^{ab} 
  ( 
    \tilde{V}_1 P_1 
    - (\tilde{e}_1 + \tilde{\gamma}_1) F_1
  )
  + C^T C^{ab} 
  (
    \dot{\tilde{V}}_1 P_1 
    + \tilde{V}_1 \dot{P}_1 
  )
\right) 
\\
% eqn 3.
\dot{u}_1 
- \frac{\Delta L}{2} C^T C^{ab} S^{11}_1 \dot{F}_1
- \frac{\Delta L}{2} C^T C^{ab} S^{12}_1 \dot{M}_1
= 
\frac{\Delta L}{2} \dot{C}^T C^{ab} (e_1 + \gamma_1)
 \\
% eqn 4.
\dot{\theta}_1
- \frac{\Delta L}{2} Q_a^{-1} C^{ab} S^{21}_1 \dot{F}_1
- \frac{\Delta L}{2} Q_a^{-1} C^{ab} S^{22}_1 \dot{M}_1
- \dot{\hat{\theta}_1} 
= 
- \frac{\Delta L}{2} Q_a^{-1} \dot{Q_a} Q_a^{-1} C^{ab} \kappa_1 
- \dot{\hat{u}_1} 
\\
% eqn 5.
C^T C^{ab} \dot{F}_1
+ \frac{\Delta L}{2} C^T C^{ab} M^{11}_1 \ddot{V}_1 
+ \frac{\Delta L}{2} C^T C^{ab} M^{12}_1 \ddot{\Omega}_1 
= 
- \dot{C}^T C^{ab} F_1
+ \dot{\bar{f}_1^{+}} 
- \frac{\Delta L}{2} 
\left( 
    \dot{\tilde{\omega}} C^T C^{ab} P_1
  + \tilde{\omega} \dot{C}^T C^{ab} P_1
  + \tilde{\omega} C^T C^{ab} \dot{P}_1
  + \ddot{C}^T C^{ab} P_1 
  + 2 \dot{C}^T C^{ab} \dot{P}_1 
\right) 
+ \dot{F_1^*} 
\\
% eqn 6.
\frac{\Delta L}{2} C^T C^{ab} \tilde{F_1} S^{11}_1 \dot{F}_1
- \frac{\Delta L}{2} C^T C^{ab} (\tilde{e}_1 + \tilde{\gamma}_1) \dot{F}_1
+ \frac{\Delta L}{2} C^T C^{ab} \tilde{F_1} S^{12}_1 \dot{M}_1
+ C^T C^{ab} \dot{M}_1
+ \frac{\Delta L}{2} C^T C^{ab} S^{21}_1 \ddot{V}_1
+ \frac{\Delta L}{2} C^T C^{ab} S^{22}_1 \ddot{\Omega}_1
= 
- \dot{C}^T C^{ab} M_1
+ \dot{\bar{m}_1^{+}}
- \frac{\Delta L}{2} 
\left(
  \dot{\tilde{\omega}} C^T C^{ab} H_1
  + \tilde{\omega} \dot{C}^T C^{ab} H_1
  + \tilde{\omega} C^T C^{ab} \dot{H}_1
  + \ddot{C}^T C^{ab} H_1
  + 2\dot{C}^T C^{ab} \dot{H}_1
  + \dot{C}^T C^{ab} 
  ( 
    \tilde{V}_1 P_1 
    - (\tilde{e}_1 + \tilde{\gamma}_1) F_1
  )
  + C^T C^{ab} 
  (
    \dot{\tilde{V}}_1 P_1 
    + \tilde{V}_1 \dot{P}_1 
  )
\right) 
+ \dot{M_1^*} 
\\
% eqn 7. 
- \dot{u}_1 
- \frac{\Delta L}{2} C^T C^{ab} S^{11}_1 \dot{F}_1
- \frac{\Delta L}{2} C^T C^{ab} S^{12}_1 \dot{M}_1
- \dot{\hat{u}_1} 
= 
\frac{\Delta L}{2} \dot{C}^T C^{ab} (e_1 + \gamma_1)
\\
% eqn 8.
- \frac{\Delta L}{2} Q_a^{-1} C^{ab} S^{21}_1 \dot{F}_1
- \frac{\Delta L}{2} Q_a^{-1} C^{ab} S^{22}_1 \dot{M}_1
- \dot{\hat{\theta}_1} 
= 
\dot{\theta}_1
- \frac{\Delta L}{2} Q_a^{-1} \dot{Q_a} Q_a^{-1} C^{ab} \kappa_1
\\
$$

We can then form a matrix problem which can be solved to obtain explicit expressions for the state rates in terms of the state variables.

$$
\begin{bmatrix}
0 & 0 & - C^T C^{ab} & 0 & 0 & 0 & \frac{\Delta L}{2} C^T C^{ab} M^{11}_1 & \frac{\Delta L}{2} C^T C^{ab} M^{12}_1  & -I & 0 & 0 & 0 \\
0 & 0 & \frac{\Delta L}{2} C^T C^{ab} \left(\tilde{F_1} S^{11}_1 - (\tilde{e}_1 + \tilde{\gamma}_1)\right) & -C^T C^{ab} \left(I - \frac{\Delta L}{2} \tilde{F_1} S^{12}_1\right) & 0 & 0 & \frac{\Delta L}{2} C^T C^{ab} M^{21}_1 & \frac{\Delta L}{2} C^T C^{ab} M^{22}_1  & 0  & -I & 0 & 0 \\
0 & 0 & \frac{\Delta L}{2} C^T C^{ab} S^{11}_1 & \frac{\Delta L}{2} C^T C^{ab} S^{12}_1 & 0 & 0 & 0 & 0 & 0  & 0 & 0 & 0 \\
0 & 0 & - \frac{\Delta L}{2} Q_a^{-1} C^{ab} S^{21}_1 & - \frac{\Delta L}{2} Q_a^{-1} C^{ab} S^{22}_1 & 0 & 0  & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & C^T C^{ab} & 0 & 0 & 0 & \frac{\Delta L}{2} C^T C^{ab} M^{11}_1 & \frac{\Delta L}{2} C^T C^{ab} M^{12}_1 & 0  & 0 & 0 & 0 \\
0 & 0 & \frac{\Delta L}{2} C^T C^{ab} \left( \tilde{F_1} S^{11}_1 - (\tilde{e}_1 + \tilde{\gamma}_1) \right) & C^T C^{ab} \left(I+\frac{\Delta L}{2} \tilde{F_1} S^{12}_1\right)
 & 0 & 0 & \frac{\Delta L}{2} C^T C^{ab} M^{21}_1 & \frac{\Delta L}{2} C^T C^{ab} M^{22}_1 & 0 & 0 & 0  & 0 \\
0 & 0 & - \frac{\Delta L}{2} C^T C^{ab} S^{11}_1 & - \frac{\Delta L}{2} C^T C^{ab} S^{12}_1 & 0 & 0 & 0 & 0 & 0 & 0  & -I & 0 \\
0 & 0 & -\frac{\Delta L}{2} Q_a^{-1} C^{ab} S^{21}_1 & -\frac{\Delta L}{2} Q_a^{-1} C^{ab} S^{22}_1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & -I \\
I & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0  & 0 & 0 & 0 \\
0 & I & 0 & 0 & 0 & 0 & 0 & 0 & 0  & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & I & 0 & 0 & 0 & 0  & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & I & 0 & 0 & 0  & 0 & 0 & 0 \\
\end{bmatrix}
\begin{bmatrix}
\dot{u} \\ \dot{\theta} \\ \dot{F} \\ \dot{M} \\ \dot{V} \\ \dot{\Omega} \\ \ddot{V} \\ \ddot{\Omega} \\ \dot{F^*_1} \\ \dot{M^*_1} \\ \dot{\hat{u}}_2 \\ \dot{\hat{\theta}}_2
\end{bmatrix} =
\begin{bmatrix}
\dot{C}^T C^{ab} F_1
+ \dot{\bar{f}_1^{-}} 
- \frac{\Delta L}{2} 
\left( 
  \dot{\tilde{\omega}} C^T C^{ab} P_1
  + \tilde{\omega} \dot{C}^T C^{ab} P_1
  + \tilde{\omega} C^T C^{ab} \dot{P}_1
  + \ddot{C}^T C^{ab} P_1 
  + 2 \dot{C}^T C^{ab} \dot{P}_1
\right) 
\\ 
\dot{C}^T C^{ab} M_1
+ \dot{\bar{m}_1^{-}}
- \frac{\Delta L}{2} 
\left(
  \dot{\tilde{\omega}} C^T C^{ab} H_1
  + \tilde{\omega} \dot{C}^T C^{ab} H_1
  + \tilde{\omega} C^T C^{ab} \dot{H}_1
  + \ddot{C}^T C^{ab} H_1
  + 2\dot{C}^T C^{ab} \dot{H}_1
  + \dot{C}^T C^{ab} 
  ( 
    \tilde{V}_1 P_1 
    - (\tilde{e}_1 + \tilde{\gamma}_1) F_1
  )
  + C^T C^{ab} 
  (
    \dot{\tilde{V}}_1 P_1 
    + \tilde{V}_1 \dot{P}_1 
  )
\right)  
\\ 
\frac{\Delta L}{2} \dot{C}^T C^{ab} (e_1 + \gamma_1) \\ 
- \frac{\Delta L}{2} Q_a^{-1} \dot{Q_a} Q_a^{-1} C^{ab} \kappa_1 
- \dot{\hat{u}_1} 
\\ 
- \dot{C}^T C^{ab} F_1
+ \dot{\bar{f}_1^{+}} 
- \frac{\Delta L}{2} 
\left( 
    \dot{\tilde{\omega}} C^T C^{ab} P_1
  + \tilde{\omega} \dot{C}^T C^{ab} P_1
  + \tilde{\omega} C^T C^{ab} \dot{P}_1
  + \ddot{C}^T C^{ab} P_1 
  + 2 \dot{C}^T C^{ab} \dot{P}_1 
\right) 
+ \dot{F_1^*} 
\\ 
- \dot{C}^T C^{ab} M_1
+ \dot{\bar{m}_1^{+}}
- \frac{\Delta L}{2} 
\left(
  \dot{\tilde{\omega}} C^T C^{ab} H_1
  + \tilde{\omega} \dot{C}^T C^{ab} H_1
  + \tilde{\omega} C^T C^{ab} \dot{H}_1
  + \ddot{C}^T C^{ab} H_1
  + 2\dot{C}^T C^{ab} \dot{H}_1
  + \dot{C}^T C^{ab} 
  ( 
    \tilde{V}_1 P_1 
    - (\tilde{e}_1 + \tilde{\gamma}_1) F_1
  )
  + C^T C^{ab} 
  (
    \dot{\tilde{V}}_1 P_1 
    + \tilde{V}_1 \dot{P}_1 
  )
\right) 
+ \dot{M_1^*} 
\\ 
\frac{\Delta L}{2} \dot{C}^T C^{ab} (e_1 + \gamma_1)
\\
- \frac{\Delta L}{2} Q_a^{-1} \dot{Q_a} Q_a^{-1} C^{ab} \kappa_1
\\ 
C^T C^{ab} V_i - v_i - \tilde{\omega_a} u_i \\
Q_a^{-1} C^{ab} \left(\Omega_i - C^{ba} C \omega_a \right)  \\
\dot{V} \\ 
\dot{\Omega}
\end{bmatrix}
 
$$