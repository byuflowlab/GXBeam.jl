"""
    PrescribedCondition{T}

Describes the forces, moments, displacements, and/or rotations prescribed at a
point.

# Fields
 - force: Six booleans indicating whether forces/moments or displacements/rotations
    are prescribed.  True/false flags in each of the six slots correspond to the
    following degrees of freedom - (F[1]/u[1], F[2]/u[2], F[3]/u[3], M[1]/θ[1],
    M[2]/θ[2], M[3]/θ[3]).  The degree of freedom that are not prescribed is
    solved for as part of the system of equations.
 - value: Prescribed non-follower force/moment or displacement for each of the
    six prescribed conditions.
 - follower: Follower force and moment. Only applied when displacement/rotation
    variables are not prescribed
"""
struct PrescribedCondition{T}
    force::NTuple{6, Bool}
    value::NTuple{6, T}
    follower::NTuple{6, T}
end

const N_GAUSS = 6
const X_GAUSS = @SVector [-cos((2*i-1)/(2*N_GAUSS)*pi) for i = 1:N_GAUSS]
const W_GAUSS = @SVector [pi/N_GAUSS for i = 1:N_GAUSS]
const W1 = (1 .- X_GAUSS) ./ 4 .* W_GAUSS
const W2 = (1 .+ X_GAUSS) ./ 4 .* W_GAUSS

"""
    DistributedLoad{T}

Describes the distributed forces and moments applied on a beam element.  These
forces/moments are stored as values evaluated at the 6 Chebyshev-Gauss
quadrature points described by x_i = -cos((2*i-1)/(2*n)*pi).

# Fields (all with shape (3,6))
 - f: Non-follower distributed force pre-evaluated at gaussian quadrature points
 - m: Non-follower distributed moment pre-evaluated at gaussian quadrature points
 - f_follower: Follower distributed force pre-evaluated at gaussian quadrature points
 - m_follower: Follower distributed force pre-evaluated at gaussian quadrature points

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
struct DistributedLoad{T}
    f::SMatrix{3, N_GAUSS, T, 18}
	m::SMatrix{3, N_GAUSS, T, 18}
	f_follower::SMatrix{3, N_GAUSS, T, 18}
	m_follower::SMatrix{3, N_GAUSS, T, 18}
end

"""
    quadrature_points(x1, x2)

Return Chebyshev-Gauss quadrature points for a beam element that starts at
x1 and stops at x2.
"""
quadrature_points(x1, x2) = (1 .+ X_CHEBYSHEV)*(x2-x1)/2 .+ x1

"""
	element_load(ΔL, Ct, distributed_load)

Return the integrated loads (`f1, m1, f2, m2`) on each element given the
element length (`ΔL`), rotation matrix (`Ct`), and distributed load.
"""
function element_load(ΔL, Ct, distributed_load)

    # integrate non-follower loads using Gauss-quadrature
    f1 = ΔL*distributed_load.f*W1
    m1 = ΔL*distributed_load.m*W1
    f2 = ΔL*distributed_load.f*W2
    m2 = ΔL*distributed_load.m*W2

	# integrate follower loads using Gauss-quadrature
    f1_follower = ΔL*distributed_load.f_follower*W1
    m1_follower = ΔL*distributed_load.m_follower*W1
    f2_follower = ΔL*distributed_load.f_follower*W2
    m2_follower = ΔL*distributed_load.m_follower*W2

	# add follower and non-follower loads together
	f1 += Ct*f1_follower
	m1 += Ct*m1_follower
	f2 += Ct*f2_follower
	m2 += Ct*m2_follower

    return f1, f2, m1, m2
end
