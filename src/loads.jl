"""
	TimeFunction{TF}

Defines a vector that defines a function for various time steps
"""
struct TimeFunction{TF} <: AbstractArray{TF, 1}
	indices::Vector{TF}
	values::Vector{TF}
	function TimeFunction{TF}(indices, values) where TF
		@assert indices[1] == 1
		@assert issorted(indices)
		@assert length(indices) == length(values)
		return new(indices, values)
	end
end

TimeFunction(indices, values) = TimeFunction{eltype(values)}(indices, values)

Base.size(tf::TimeFunction) = size(tf.indices)
Base.getindex(tf::TimeFunction{TF}, i::Int) where TF = tf.values[searchsortedlast(tf.indices, i)]
Base.eltype(::TimeFunction{TF}) where TF = TF
Base.eltype(::Type{TimeFunction{TF}}) where TF = TF

"""
    PrescribedConditions{T}

Describes the forces, moments, displacements, and/or rotations prescribed at a
point.
"""
struct PrescribedConditions{T}
    force_dof::SVector{6, Bool}
    value::SVector{6, T}
    follower::SVector{6, T}
	value_tf::SVector{6, Int}
	follower_tf::SVector{6, Int}
end
Base.eltype(::PrescribedConditions{T}) where T = T

"""
    PrescribedConditions(; kwargs...)

Construct an object of type PrescribedConditions which stores the prescribed
conditions for a point.

# Keyword Arguments
 - force_dof: Collection of six booleans indicating whether forces/moments and/or
 	displacements/rotations are prescribed. True/false flags in each of the six slots correspond to the
    following degrees of freedom - (F[1]/u[1], F[2]/u[2], F[3]/u[3], M[1]/θ[1],
    M[2]/θ[2], M[3]/θ[3]).  The degree of freedom that are not prescribed is
    solved for as part of the system of equations.
 - value: Prescribed non-follower force/moment and/or displacement/rotations for
 	each of the six degrees of freedom.
 - follower: Follower force and moment applied when displacement/rotation variables
 	are not prescribed.  Note that both follower and non-follower forces can be
	applied simultaneously.
-  value_tf: Time function for prescribed non-follower force/moment and/or
    displacement/rotations.
 - follower_tf: Time function for follower force and moment
"""
function PrescribedConditions(;
	force_dof = (@SVector ones(Bool, 6)),
	value = (@SVector zeros(6)),
	follower = (@SVector zeros(6)),
	value_tf = (@SVector zeros(Int, 6)),
	follower_tf = (@SVector zeros(Int, 6)))

	T = promote_type(eltype(value), eltype(follower))

	return PrescribedConditions{T}(SVector{6, Bool}(force_dof),
		SVector{6, T}(value), SVector{6, T}(follower),
		SVector{6, Int}(value_tf), SVector{6, Int}(follower_tf))
end

"""
    DistributedLoads{T}

Describes the distributed forces and moments applied on a beam element.  These
forces/moments are stored as values evaluated at the 6 Chebyshev-Gauss
quadrature points described by x_i = -cos((2*i-1)/(2*n)*pi).

# Fields (all with shape (3,6))
 - f: Non-follower distributed force pre-evaluated at gaussian quadrature points
 - m: Non-follower distributed moment pre-evaluated at gaussian quadrature points
 - f_follower: Follower distributed force pre-evaluated at gaussian quadrature points
 - m_follower: Follower distributed force pre-evaluated at gaussian quadrature points
 - f_tf: Time function for each value which defines the non-follower distributed force
 - m_tf: Time function for each value which defines the non-follower distributed moment
 - f_follower_tf: Time function for each value which defines the follower distributed force
 - m_follower_tf: Time function for each value which defines the follower distributed moments

# Quadrature Points: (defined on interval from -1 to 1, with -1 as the left
side and 1 as the right side of the beam element
of the beam element
 - x_1: -0.9659258262890682
 - x_2: -0.7071067811865475
 - x_3: -0.25881904510252085
 - x_4: 0.25881904510252074
 - x_5: 0.7071067811865476
 - x_6: 0.9659258262890683
"""
struct DistributedLoads{T}
    f::SMatrix{3, N_GAUSS, T, 18}
	m::SMatrix{3, N_GAUSS, T, 18}
	f_follower::SMatrix{3, N_GAUSS, T, 18}
	m_follower::SMatrix{3, N_GAUSS, T, 18}
	f_tf::SMatrix{3, N_GAUSS, Int, 18}
	m_tf::SMatrix{3, N_GAUSS, Int, 18}
	f_follower_tf::SMatrix{3, N_GAUSS, Int, 18}
	m_follower_tf::SMatrix{3, N_GAUSS, Int, 18}
end
Base.eltype(::DistributedLoads{T}) where T = T

function DistributedLoads(;
	f = (@SMatrix zeros(3, N_GAUSS)),
	m = (@SMatrix zeros(3, N_GAUSS)),
	f_follower = (@SMatrix zeros(3, N_GAUSS)),
	m_follower = (@SMatrix zeros(3, N_GAUSS)),
	f_tf = (@SMatrix zeros(Int, 3, N_GAUSS)),
	m_tf = (@SMatrix zeros(Int, 3, N_GAUSS)),
	f_follower_tf = (@SMatrix zeros(Int, 3, N_GAUSS)),
	m_follower_tf = (@SMatrix zeros(Int, 3, N_GAUSS)))

	T = promote_type(eltype(f), eltype(m), eltype(f_follower), eltype(m_follower))

	return DistributedLoads{T}(
		SMatrix{3, N_GAUSS}(f), SMatrix{3, N_GAUSS}(m),
		SMatrix{3, N_GAUSS}(f_follower), SMatrix{3, N_GAUSS}(m_follower),
		SMatrix{3, N_GAUSS, Int}(f_tf), SMatrix{3, N_GAUSS, Int}(m_tf),
		SMatrix{3, N_GAUSS, Int}(f_follower_tf), SMatrix{3, N_GAUSS, Int}(m_follower_tf),
		)
end

"""
    quadrature_points(x1, x2)

Return Chebyshev-Gauss quadrature points for a beam element that starts at
x1 and stops at x2.
"""
quadrature_points(x1, x2) = (1 .+ X_CHEBYSHEV)*(x2-x1)/2 .+ x1

"""
	integrate_element_loads(ΔL, Ct, distributed_load)

Return the integrated loads (`f1, m1, f2, m2`) on each element given the
element length (`ΔL`), rotation matrix (`Ct`), and distributed load.
"""
function integrate_element_loads(ΔL, Ct, distributed_load, time_function_values)

	# apply time functions to forces/moments
	f = distributed_load.f .* time_function_values[distributed_load.f_tf]
	m = distributed_load.m .* time_function_values[distributed_load.m_tf]
	f_follower = distributed_load.f_follower .* time_function_values[distributed_load.f_follower_tf]
	m_follower = distributed_load.m_follower .* time_function_values[distributed_load.m_follower_tf]

    # integrate non-follower loads using Gauss-quadrature
    f1 = ΔL*f*W1
    m1 = ΔL*m*W1
    f2 = ΔL*f*W2
    m2 = ΔL*m*W2

	# integrate follower loads using Gauss-quadrature
    f1_follower = ΔL*f_follower*W1
    m1_follower = ΔL*m_follower*W1
    f2_follower = ΔL*f_follower*W2
    m2_follower = ΔL*m_follower*W2

	# add follower and non-follower loads together
	f1 += Ct*f1_follower
	m1 += Ct*m1_follower
	f2 += Ct*f2_follower
	m2 += Ct*m2_follower

    return f1, f2, m1, m2
end

"""
	follower_load_jacobians(ΔL, Ct_θ1, Ct_θ2, Ct_θ3, distributed_load, time_function_values)

Return the jacobian of the follower loads with respect to θ.
"""
function follower_load_jacobians(ΔL, Ct_θ1, Ct_θ2, Ct_θ3, distributed_load, time_function_values)

	# apply time functions to forces/moments
	f_follower = distributed_load.f_follower .* time_function_values[distributed_load.f_follower_tf]
	m_follower = distributed_load.m_follower .* time_function_values[distributed_load.m_follower_tf]

	# integrate follower loads using Gauss-quadrature
    f1_follower = ΔL*f_follower*W1
    m1_follower = ΔL*m_follower*W1
    f2_follower = ΔL*f_follower*W2
    m2_follower = ΔL*m_follower*W2

	# get jacobian of follower loads w.r.t. θ
	f1_θ = matrix_jacobian_product(Ct_θ1, Ct_θ2, Ct_θ3, f1_follower)
	m1_θ = matrix_jacobian_product(Ct_θ1, Ct_θ2, Ct_θ3, f2_follower)
	f2_θ = matrix_jacobian_product(Ct_θ1, Ct_θ2, Ct_θ3, m1_follower)
	m2_θ = matrix_jacobian_product(Ct_θ1, Ct_θ2, Ct_θ3, m2_follower)

    return f1_θ, f2_θ, m1_θ, m2_θ
end
