# Examples

These examples show how to use the various analysis capabilities of GEBT.jl.  Some of the examples also provide a verification/validation for the implementation of geometrically exact beam theory in GEBT.jl.


```@contents
Pages = ["examples.md"]
Depth = 3
```

## Linear Analysis of a Cantilever Subjected to a Uniform Distributed Load

```@example

using GEBT, LinearAlgebra

nelem = 12

# create points
n1 = n3 = div(nelem, 3)
n2 = nelem - n1 - n3
x1 = range(0, 0.3, length=n1+1)
x2 = range(0.3, 0.7, length=n2+1)
x3 = range(0.7, 1.0, length=n3+1)
x = vcat(x1, x2[2:end], x3[2:end])
y = zero(x)
z = zero(x)
points = [[x[i],y[i],z[i]] for i = 1:length(x)]

# index of endpoints for each beam element
start = 1:nelem
stop = 2:nelem+1

# create compliance matrix for each beam element
compliance = fill(Diagonal([2.93944738387698E-10, 0, 0, 4.69246721094557E-08, 6.79584e-8, 1.37068861370898E-09]), nelem)

# create assembly
assembly = Assembly(points, start, stop, compliance=compliance)

# set prescribed conditions (fixed right endpoint)
prescribed_conditions = Dict(
    nelem+1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0)
)

# create distributed load
distributed_loads = Dict()
for ielem in n1+1:n1+n2
    distributed_loads[ielem] = DistributedLoads(assembly, ielem; fz = (s) -> 1000)
end

system, converged = static_analysis(assembly, prescribed_conditions=prescribed_conditions,
    distributed_loads=distributed_loads, linear=true)

state = AssemblyState(system, assembly, prescribed_conditions=prescribed_conditions)

# Plot Results

using Plots
pyplot()

# set up the plot
plot(
    xlim = (0.0, 1.0),
    xticks = 0.0:0.2:1.0,
    xlabel = "x (m)",
    ylim = (0, 3e-6),
    yticks = 0.0:5e-7:3e-6,
    ylabel = "Deflection (m)",
    grid = false,
    overwrite_figure=false
    )

# GEBT
x = [assembly.points[ipoint][1] + state.points[ipoint].u[1] for ipoint = 1:length(assembly.points)]
y = [state.points[ipoint].u[3] for ipoint = 1:length(assembly.points)]
plot!(x, y, markershape=:circle, label="")

savefig("cantilever1.png"); nothing # hide

```

![](cantilever1.png)

## Linear Analysis of a Beam Under a Linear Distributed Load

```@example

using GEBT, LinearAlgebra

nelem = 16

# create points
x = range(0, 1, length=nelem+1)
y = zero(x)
z = zero(x)
points = [[x[i],y[i],z[i]] for i = 1:length(x)]

# index of endpoints for each beam element
start = 1:nelem
stop = 2:nelem+1

# create compliance matrix for each beam element
compliance = fill(Diagonal([2.93944738387698E-10, 0, 0, 4.69246721094557E-08, 6.79584e-8, 1.37068861370898E-09]), nelem)

# create assembly
assembly = Assembly(points, start, stop, compliance=compliance)

# set prescribed conditions
prescribed_conditions = Dict(
    # simply supported left endpoint
    1 => PrescribedConditions(uz=0),
    # clamped right endpoint
    nelem+1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0)
)

# create distributed load
distributed_loads = Dict()
for i = 1:nelem
    distributed_loads[i] = DistributedLoads(assembly, i; s1=x[i],
        s2=x[i+1], fz = (s) -> 1000*s)
end

system, converged = static_analysis(assembly, prescribed_conditions=prescribed_conditions,
    distributed_loads=distributed_loads, linear=true)

state = AssemblyState(system, assembly, prescribed_conditions=prescribed_conditions)

# Plot Results

using Plots
pyplot()

# set up the plot
plot(
    xlim = (0.0, 1.0),
    xticks = 0.0:0.2:1.0,
    xlabel = "x (m)",
    ylim = (0, 1.7e-7),
    yticks = 0.0:5e-8:1.5e-7,
    ylabel = "Deflection (m)",
    grid = false,
    overwrite_figure=false
    )

# GEBT
x = [assembly.points[ipoint][1] + state.points[ipoint].u[1] for ipoint = 1:length(assembly.points)]
y = [state.points[ipoint].u[3] for ipoint = 1:length(assembly.points)]
plot!(x, y, markershape=:circle, label="")

savefig("cantilever2.png"); nothing # hide

```

![](cantilever2.png)


## Nonlinear Analysis of a Cantilever Subjected to a Constant Moment

This problem is a common benchmark problem for the geometrically nonlinear analysis of beams and has an analytical solution.

```@example

using GEBT, LinearAlgebra

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

# index of endpoints of each beam element
start = 1:nelem
stop = 2:nelem+1

# compliance matrix for each beam element
compliance = fill(Diagonal([1/(E*A), 0, 0, 0, 1/(E*Iyy), 1/(E*Izz)]), nelem)

# create assembly of interconnected nonlinear beams
assembly = Assembly(points, start, stop, compliance=compliance)

# pre-initialize system storage
static = true
keep_points = [1, nelem+1] # points that we request are included in the system of equations
system = System(assembly, preserved_points, static)

# run an analysis for each prescribed bending moment

states = Vector{AssemblyState{Float64}}(undef, length(M))

for i = 1:length(M)

    # create dictionary of prescribed conditions
    prescribed_conditions = Dict(
        # fixed left side
        1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
        # moment on right side
        nelem+1 => PrescribedConditions(Mz = M[i])
    )

    static_analysis!(system, assembly, prescribed_conditions=prescribed_conditions)

    states[i] = AssemblyState(system, assembly, prescribed_conditions=prescribed_conditions)

end

# Plot Results

using Plots
pyplot()

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

savefig("cantilever3.png"); nothing # hide
```

![](cantilever3.png)

## Nonlinear Analysis of the Bending of a Curved Beam in 3D Space

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
    1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
    # force on right endpoint
    nelem+1 => PrescribedConditions(Fz=P)
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
