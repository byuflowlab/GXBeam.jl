# Examples

These examples show how to use the various analysis capabilities of GEBT.jl.  Some of the examples also provide a verification/validation for the implementation of geometrically exact beam theory in GEBT.jl.


```@contents
Pages = ["examples.md"]
Depth = 3
```

## Linear Analysis Examples

## Nonlinear Static Analysis Examples



### Cantilever Subjected to a Constant Moment

This problem is a common benchmark problem for the geometrically nonlinear analysis of beams and has an analytical solution.

```@example

using GEBT
using Plots
pyplot()

L = 12 # inches
h = w = 1 # inches
E = 30e6 # lb/in^4 Young's Modulus

A = h*w
Iyy = w*h^3/12
Izz = w^3*h/12

# bending moment (applied at end)
# note that solutions for λ > 1.8 do not converge
λ = [0.0, 0.4, 0.8, 1.2, 1.6, 1.8]
m = pi*E*Iyy/L
M = λ*m

# create points
nelem = 16
x = range(0, L, length=nelem+1)
y = zero(x)
z = zero(x)
points = [[x[i],y[i],z[i]] for i = 1:length(x)]

# index of first and last endpoints of each beam element
start = 1:nelem
stop = 2:nelem+1

# compliance matrix for each beam element
compliance = fill(Diagonal([1/(E*A), 0, 0, 0, 1/(E*Iyy), 1/(E*Izz)]), nelem)

# create assembly of interconnected nonlinear beams
assembly = Assembly(points, start, stop, compliance=compliance)

# pre-initialize system storage
static = true
keep_points = [1, nelem+1] # points that we request are included in the system of equations
n_tf = 0 # number of time functions
system = System(assembly, preserved_points, static, n_tf)

# run an analysis for each prescribed bending moment

states = Vector{AssemblyState{Float64}}(undef, length(M))

for i = 1:length(M)

    # create dictionary of prescribed conditions
    prescribed_conditions = Dict(
        # fixed left side
        1 => PrescribedConditions(
            force_dof = [false, false, false, false, false, false]
            ),
        # moment on right side
        nelem+1 => PrescribedConditions(
            follower = [0, 0, 0, 0, 0, M[i]]
        )
    )

    static_analysis!(system, assembly, prescribed_conditions=prescribed_conditions)

    states[i] = AssemblyState(system, assembly, prescribed_conditions=prescribed_conditions)

end

# analytical solution (ρ = E*I/M)
analytical(x, ρ) = ifelse(ρ == Inf, zeros(3), [ρ*sin(x/ρ)-x, ρ*(1-cos(x/ρ)), 0])

# Plot Results

# set up the plot
plot(
    xlim = (-0.25, 1.1),
    xticks = -0.25:0.25:1.0,
    xlabel = "x/L",
    ylim = (-0.05, 0.8),
    yticks = 0.0:0.2:0.8,
    ylabel = "y/L",
    aspect_ratio = 1.0,
    grid = false,
    overwrite_figure=false
    )

# create dummy legend entries for GEBT and Analytical
scatter!([NaN, NaN], [NaN, NaN], color=:black, label="GEBT")
plot!([NaN, NaN], [NaN, NaN], color=:black, label="Analytical")

# plot the data
for i = 1:length(M)
    # GEBT
    x = [assembly.points[ipoint][1] + states[i].points[ipoint].u[1] for ipoint = 1:length(assembly.points)]
    y = [assembly.points[ipoint][2] + states[i].points[ipoint].u[2] for ipoint = 1:length(assembly.points)]
    scatter!(x/L, y/L, label="", color = i)

    # Analytical
    x0 = range(0, L, length=100)
    deflection = analytical.(x0, E*Iyy/M[i])
    x = (x0 + getindex.(deflection, 1))
    y = getindex.(deflection, 2)
    plot!(x/L, y/L, label="λ=$(λ[i])", color=i)
end

# show the plot
plot!(show=true)

savefig("cantilever.png"); nothing # hide
```

![](cantilever.png)

### Bending of a Curved Beam in 3D Space

This problem is also a common benchmark problem for the geometrically exact bending of nonlinear beams, but does not have an analytical solution.

```

# problem constants
R = 100
L = R*pi/4 # inches
h = w = 1 # inches
E = 1e7 # psi Young's Modulus
ν = 0.0
G = E/(2*(1+ν))

# beam starting point and curvature
r = [0, 0, 0]
frame = [0 -1 0; 1 0 0; 0 0 1]
k = [0, 0, -1/R]

# cross section properties
A = h*w
Ay = A
Az = A
Iyy = w*h^3/12
Izz = w^3*h/12
J = Iyy + Izz

# discretize the beam
nelem = 16
ΔL, xp, xm, Cab = discretize_beam(L, r, nelem, Cab=frame, k = k)

# force
P = 600 # lbs

# index of left and right endpoints of each beam element
pt1 = 1:nelem
pt2 = 2:nelem+1

# compliance matrix for each beam element
compliance = fill(Diagonal([1/(E*A), 1/(G*Ay), 1/(G*Az), 1/(G*J), 1/(E*Iyy), 1/(E*Izz)]), nelem)

# create assembly of interconnected nonlinear beams
assembly = Assembly(xp, pt1, pt2, compliance=compliance, frames=Cab,
    lengths=ΔL, midpoints=xm)

# create dictionary of prescribed conditions
prescribed_conditions = Dict(
    # fixed left endpoint
    1 => PrescribedConditions(
        force_dof = [false, false, false, false, false, false]
        ),
    # force on right endpoint
    nelem+1 => PrescribedConditions(
        value = [0, 0, P, 0, 0, 0]
    )
)

system, converged = static_analysis(assembly, prescribed_conditions=prescribed_conditions)

state = AssemblyState(system, assembly, prescribed_conditions=prescribed_conditions)

println("Tip Displacement: ", state.points[end].u)
println("Tip Displacement (Bathe and Bolourch): [-13.4, -23.5, 53.4]")

# write a file that can be visualized in ParaView
write_vtk("curved", assembly, state)

```

![](curved.png)

The calculated tip displacements match those found by Bathe and Bolourch closely, thus verifying our GEBT implementation.


## Nonlinear Steady State Examples

## Nonlinear Stability Analysis Examples

## Nonlinear Time Marching Analysis Examples
