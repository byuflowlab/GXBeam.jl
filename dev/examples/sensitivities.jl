using GXBeam, LinearAlgebra
import ForwardDiff # for forward-mode automatic differentiation
import ReverseDiff # for reverse-mode automatic differentiation
using BenchmarkTools # for benchmarking function performance

L = 12 # inches
h = w = 1 # inches
E = 30e6 # lb/in^4 Young's Modulus

A = h*w
Iyy = w*h^3/12
Izz = w^3*h/12

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

# construct objective function
objfun1 = (p) -> begin

    # non-dimensional tip moment
    λ = p[1]

    # dimensionalized tip moment
    m = pi*E*Iyy/L
    M = λ*m

    # prescribed conditions
    prescribed_conditions = Dict(
        # fixed left side
        1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
        # moment on right side
        nelem+1 => PrescribedConditions(Mz = M)
    )

    # initialize internal storage with the correct type
    system = StaticSystem(eltype(p), assembly)

    # perform static analysis
    system, state, converged = static_analysis!(system, assembly; prescribed_conditions)

    # return the desired outputs
    return [state.points[end].u[1], state.points[end].u[2]]
end

# compute sensitivities using ForwardDiff with λ = 1.0
@btime ForwardDiff.jacobian(objfun1, [1.0])

# compute sensitivities using ReverseDiff with λ = 1.0
@btime ReverseDiff.jacobian(objfun1, [1.0])

# construct pfunc to overwrite prescribed conditions
pfunc = (p, t) -> begin

    # non-dimensional tip moment
    λ = p[1]

    # dimensionalized tip moment
    m = pi*E*Iyy/L
    M = λ*m

    # create dictionary of prescribed conditions
    prescribed_conditions = Dict(
        # fixed left side
        1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
        # moment on right side
        nelem+1 => PrescribedConditions(Mz = M)
    )

    # return named tuple with new arguments
    return (; prescribed_conditions=prescribed_conditions)
end

# construct objective function
objfun2 = (p) -> begin

    # perform static analysis
    system, state, converged = static_analysis(assembly; pfunc, p)

    # return the desired outputs
    return [state.points[end].u[1], state.points[end].u[2]]
end

# compute sensitivities using ForwardDiff with λ = 1.0
@btime ForwardDiff.jacobian(objfun2, [1.0])

# compute sensitivities using ReverseDiff with λ = 1.0
@btime ReverseDiff.jacobian(objfun2, [1.0])

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

