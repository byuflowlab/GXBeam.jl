# # [Flexible Pendulum Simulation](@id pendulum-simulation)
#
# In this example we simulate the motion of a flexible cantilever beam and pendulum.
# 
# ![](../assets/pendulum-drawing.svg)
#
#-
#md # !!! tip
#md #     This example is also available as a Jupyter notebook:
#md #     [`pendulum-simulation.ipynb`](@__NBVIEWER_ROOT_URL__/examples/pendulum-simulation.ipynb).
#-

using Plots, GXBeam, DifferentialEquations, LinearAlgebra

L = 1 # m

## create points
nelem = 10
x = range(0, L, length=nelem + 1)
y = zeros(nelem + 1)
z = zeros(nelem + 1)
points = [[x[i],y[i],z[i]] for i = 1:length(x)]

## index of endpoints of each beam element
start = 1:nelem
stop = 2:nelem+1

## chosen damping coefficient
μ = 0.01 # s

## compliance matrices
compliance = fill(Diagonal([1e-6, 1e-6, 1e-6, 1e-6, 1/0.15, 1e-6]), nelem)

## mass matrices
mass = fill(Diagonal([0.15, 0.15, 0.15, 1e-6, 1e-6, 1e-6]), nelem)

## damping coefficients
damping = fill([μ, μ, μ, μ, μ, μ], nelem)

## create assembly of interconnected nonlinear beams
assembly = Assembly(points, start, stop; compliance=compliance, mass=mass, damping=damping)

## prescribed conditions
prescribed_conditions = Dict(
    ## fixed left side, with rigid body motion about the y-axis
    1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, My=0, theta_z=0),
)

## define simulation time
tspan = (0.0, 0.68)

## define named tuple with parameters
p = (; prescribed_conditions=prescribed_conditions, gravity=[0.0, 0.0, -9.8])

## run initial condition analysis to get a consistent set of initial conditions
system, converged = initial_condition_analysis(assembly, tspan[1]; 
    prescribed_conditions = prescribed_conditions,
    structural_damping = false)

## construct ODEProblem
prob = ODEProblem(system, assembly, tspan, p; 
    constant_mass_matrix = false,
    structural_damping = false)

## solve ODE
sol = solve(prob, Rodas4())

## times at which to save the simulation outputs
t = 0.0:0.01:0.68

## post-process solution
history = [AssemblyState(system, assembly, sol(t[i]); prescribed_conditions) for i in eachindex(t)]

## visualize solution
write_vtk("pendulum-simulation", assembly, history, t)
