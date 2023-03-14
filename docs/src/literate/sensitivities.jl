# # [Sensitivity Analysis](@id sensitivities)
#
# In order to facilitate gradient-based optimization, all linear solves in GXBeam have been
# overloaded to support automatic differentiation using the ImplicitAD package.  Automatic
# differentiation using ForwardDiff and ReverseDiff should therefore work with without any
# major headaches, as long as you initialize the internal storage with the correct type.
# (e.g. `system = DynamicSystem(TF, assembly)` where `TF` is an appropriate floating point
# type).
#
# For example, consider the [Cantilever with a Tip Moment](@id tipmoment) example.  Suppose
# we were interested in the sensitivity of tip x and y-displacement with respect to
# the nondimensional tip moment ``\lambda`` when ``\lambda=1``.  These sensitivites may
# be computed as follows:

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

## create points
nelem = 16
x = range(0, L, length=nelem+1)
y = zero(x)
z = zero(x)
points = [[x[i],y[i],z[i]] for i = 1:length(x)]

## index of endpoints of each beam element
start = 1:nelem
stop = 2:nelem+1

## compliance matrix for each beam element
compliance = fill(Diagonal([1/(E*A), 0, 0, 0, 1/(E*Iyy), 1/(E*Izz)]), nelem)

## create assembly of interconnected nonlinear beams
assembly = Assembly(points, start, stop, compliance=compliance)

## construct objective function
objfun1 = (p) -> begin

    ## non-dimensional tip moment
    λ = p[1]

    ## dimensionalized tip moment
    m = pi*E*Iyy/L
    M = λ*m

    ## prescribed conditions
    prescribed_conditions = Dict(
        ## fixed left side
        1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
        ## moment on right side
        nelem+1 => PrescribedConditions(Mz = M)
    )

    ## initialize internal storage with the correct type
    system = StaticSystem(eltype(p), assembly)

    ## perform static analysis
    system, state, converged = static_analysis!(system, assembly; prescribed_conditions)

    ## return the desired outputs
    return [state.points[end].u[1], state.points[end].u[2]]
end

## compute sensitivities using ForwardDiff with λ = 1.0
@btime ForwardDiff.jacobian(objfun1, [1.0])

#-

## compute sensitivities using ReverseDiff with λ = 1.0
@btime ReverseDiff.jacobian(objfun1, [1.0])

# Advanced users, however, may wish to use overloaded versions of each nonlinear solve
# in order to further decrease the total computational costs associated with obtaining
# design sensitivites.  Overloading the nonlinear solver also significantly reduces the
# memory requirements associated with using ReverseDiff. Using these overloads, however,
# requires that the user provide the parameter function `parameters = pfunc(p, t)` and
# associated parameter vector `p`.  As described in the documentation for each analysis
# type, the `pfunc` function returns a named tuple which contains updated arguments for
# the analysis, based on the contents of the parameter vector `p` and the current time `t`.

## construct pfunc to overwrite prescribed conditions
pfunc = (p, t) -> begin

    ## non-dimensional tip moment
    λ = p[1]

    ## dimensionalized tip moment
    m = pi*E*Iyy/L
    M = λ*m

    ## create dictionary of prescribed conditions
    prescribed_conditions = Dict(
        ## fixed left side
        1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
        ## moment on right side
        nelem+1 => PrescribedConditions(Mz = M)
    )

    ## return named tuple with new arguments
    return (; prescribed_conditions=prescribed_conditions)
end

## construct objective function
objfun2 = (p) -> begin

    ## perform static analysis
    system, state, converged = static_analysis(assembly; pfunc, p)

    ## return the desired outputs
    return [state.points[end].u[1], state.points[end].u[2]]
end

## compute sensitivities using ForwardDiff with λ = 1.0
@btime ForwardDiff.jacobian(objfun2, [1.0])

#-

## compute sensitivities using ReverseDiff with λ = 1.0
@btime ReverseDiff.jacobian(objfun2, [1.0])
