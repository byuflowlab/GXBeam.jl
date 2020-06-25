module GEBT

using LinearAlgebra
using StaticArrays
using FLOWMath

abstract type AbstractAnalysisType end
struct Steady <: AbstractSystemType end
struct Unsteady <: AbstractSystemType end
struct InitialStep <: AbstractSystemType end
struct TimeMarching{TF} <: AbstractSystemType
	dt::TF
end

"""
    tilde(x)

Constructs the cross product operator
"""
tilde(x) = SMatrix{3,3}(0, -x[3], x[2], x[3], 0, -x[1], -x[2], x[1], 0)

"""
    wiener_milenkovic

Constructs the wiener milenkovic rotation matrix
"""
function wiener_milenkovic(c)
	c0 = 2 - c'*c/8
	return ((c0^2 - c'*c)*I - 2*c0*tilde(c) + 2*c*c')/(4-c0)^2
end

const I3 = @SMatrix [1 0 0; 0 1 0; 0 0 1]

function static_analysis(assembly, prescribed_conditions, distributed_loads)

end

function steady_state_analysis(x0, ω0, v0, niter, nstep)

end

function time_domain_analysis(x0, ω0, v0, V, Ωt, Vt, niter, nstep)

end

function eigenvalue_analysis(x0, ω0, v0, niter, nstep, nev)

end

function apply_load(; follower=false)

end

function element_jacobian()

end


function element_mass()

end

"""
    Evaluates nonlinear equations for a point
"""
function point_rhs()

end

function point_jacobian()

end

function member_rhs()

end

function member_jacobian()

end

function load_jacobian()

end

function update_follower()

end

function transform_follower()

end

function follower_jacobian()

end

function load_integration(ds, se)

end

function direction_cosine_transpose()

end

function direction_cosine_derivative()

end

function system_jacobian()

end

function system_rhs()

end



end
