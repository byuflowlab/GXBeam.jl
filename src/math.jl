# unit vectors and identity matrix
const e1 = SVector(1, 0, 0)
const e2 = SVector(0, 1, 0)
const e3 = SVector(0, 0, 1)
const I3 = @SMatrix [1 0 0; 0 1 0; 0 0 1]

"""
    tilde(x)

Construct the cross product operator matrix
"""
tilde(x) = @SMatrix [0 -x[3] x[2]; x[3] 0 -x[1]; -x[2] x[1] 0]

"""
    wiener_milenkovic

Construct the Wiener-Milenkovic rotation matrix
"""
function wiener_milenkovic(c)
	c0 = 2 - c'*c/8

	return 1/(4-c0)^2*(@SMatrix [
		c0^2 + c[1]^2 - c[2]^2 - c[3]^2  2*(c[1]*c[2] + c0*c[3])          2*(c[1]*c[3]-c0*c[2]);
		2*(c[1]*c[2] - c0*c[3])          c0^2 - c[1]^2 + c[2]^2 - c[3]^2  2*(c[2]*c[3] + c0*c[1]);
		2*(c[1]*c[3] + c0*c[2])          2*(c[2]*c[3] - c0*c[1])          c0^2 - c[1]^2 - c[2]^2 + c[3]^2
		]
	)
end

"""
    wiener_milenkovic_jacobian(c)
	wiener_milenkovic_jacobian(C, c)

Return the gradient of the Wiener-Milenkovic rotation matrix w.r.t. `c`
"""
wiener_milenkovic_jacobian

wiener_milenkovic_jacobian(c) = wiener_milenkovic_jacobian(wiener_milenkovic(c), c)

function wiener_milenkovic_jacobian(C, c)
	c0 = 2 - c'*c/8
	c0_c = -c/4

	tmp = 1/(4-c0)^2
	tmp_c = 2/(4-c0)^3*c0_c

	C_c1 = tmp_c[1]*C/tmp + tmp*(@SMatrix [
	    2*c0*c0_c[1] + 2*c[1]    2*(c[2] + c0_c[1]*c[3]) 2*(c[3] - c0_c[1]*c[2]);
	 	2*(c[2] - c0_c[1]*c[3])  2*c0*c0_c[1] - 2*c[1]   2*(c0_c[1]*c[1] + c0);
 		2*(c[3] + c0_c[1]*c[2]) -2*(c0_c[1]*c[1] + c0)   2*c0*c0_c[1] - 2*c[1]
		]
	)

	C_c2 = tmp_c[2]*C/tmp + tmp*(@SMatrix [
	    2*c0*c0_c[2] - 2*c[2]   2*(c[1] + c0_c[2]*c[3]) -2*(c0_c[2]*c[2] + c0);
		2*(c[1] - c0_c[2]*c[3]) 2*c0*c0_c[2] + 2*c[2]    2*(c[3] + c0_c[2]*c[1]);
 		2*(c0_c[2]*c[2] + c0)   2*(c[3] - c0_c[2]*c[1])  2*c0*c0_c[2] - 2*c[2]
		]
	)

	C_c3 = tmp_c[3]*C/tmp + tmp*(@SMatrix [
		 2*c0*c0_c[3] - 2*c[3]   2*(c0_c[3]*c[3] + c0)    2*(c[1] - c0_c[3]*c[2]);
		-2*(c0_c[3]*c[3] + c0)   2*c0*c0_c[3] - 2*c[3]    2*(c[2] + c0_c[3]*c[1]);
		 2*(c[1] + c0_c[3]*c[2]) 2*(c[2] - c0_c[3]*c[1])  2*c0*c0_c[3] + 2*c[3]
		]
	)

	return C_c1, C_c2, C_c3
end

"""
    get_Q(c)

Return the matrix Q as defined in the paper "Geometrically nonlinear analysis
of composite beams using Wiener-Milenković parameters" by Qi Wang and Wenbin Yu.
"""
function get_Q(c)
	c0 = 2 - c'*c/8 # TODO: check whether this is a mistake: 2 - c'*c
	return 1/(4-c0)^2*((4 - 1/4*c'*c)*I - 2*tilde(c) + 1/2*c*c')
end

"""
    get_Q_jacobian(c)

Return the jacobian of the matrix Q as defined in the paper "Geometrically
nonlinear analysis of composite beams using Wiener-Milenković parameters" by
Qi Wang and Wenbin Yu.
"""
get_Q_jacobian(c) = get_Q_jacobian(get_Q(c), c)

function get_Q_jacobian(Q, c)

	c0 = 2 - c'*c/8
	c0_c = -c/4

	tmp = 1/(4-c0)^2
	tmp_c = 2/(4-c0)^3*c0_c

	Q_c1 = tmp_c[1]*Q/tmp + tmp*(-c[1]/2*I - 2*tilde(e1) + (e1*c' + c*e1')/2)

	Q_c2 = tmp_c[2]*Q/tmp + tmp*(-c[2]/2*I - 2*tilde(e2) + (e2*c' + c*e2')/2)

	Q_c3 = tmp_c[3]*Q/tmp + tmp*(-c[3]/2*I - 2*tilde(e3) + (e3*c' + c*e3')/2)

	return Q_c1, Q_c2, Q_c3
end

"""
    get_Qinv(c)

Return the inverse of the matrix Q as defined in the paper "Geometrically
nonlinear analysis of composite beams using Wiener-Milenković parameters" by
Qi Wang and Wenbin Yu.
"""
get_Qinv(c) =  (1 - 1/16*c'*c)*I3 + 1/2*tilde(c) + 1/8*c*c'

"""
    get_Qinv_jacobian(c)

Return the inverse of the matrix Q as defined in the paper "Geometrically
nonlinear analysis of composite beams using Wiener-Milenković parameters" by
Qi Wang and Wenbin Yu.
"""
function get_Qinv_jacobian(c)

	Qinv_c1 = -c[1]/8*I3 + 1/2*tilde(e1) + 1/8*(e1*c' + c*e1')
	Qinv_c2 = -c[2]/8*I3 + 1/2*tilde(e2) + 1/8*(e2*c' + c*e2')
	Qinv_c3 = -c[3]/8*I3 + 1/2*tilde(e3) + 1/8*(e3*c' + c*e3')

	return Qinv_c1, Qinv_c2, Qinv_c3
end

"""
	matrix_jacobian_product(A_1, A_2, A_3, b)

Return a matrix which is the product of the jacobian of a rotation matrix `A`
and an arbitrary vector b.  The jacobian of the rotation matrix is stored in
A_1, A_2, and A_3 which correspond to the partials with respect to θ1, θ2, and
θ3 respectively.
"""
matrix_jacobian_product(A_1, A_2, A_3, b) = hcat(A_1*b, A_2*b, A_3*b)

"""
	curve_length(r1, r2, k)

Calculate the length of a curve given its endpoints and its curvature vector
"""
function curve_length(r1, r2, k)
	kn = sqrt(k'*k)
	r = r2 - r1
	rn = sqrt(r'*r)
	if kn == 0 || (k[2] == 0 && k[3] == 0)
		ΔL = rn
	else
		if k[1] == 0
			ΔL = 2*asin(kn*rn/2)/kn
		else
			lower_bound = rn
			upper_bound = rn*kn/k
			residual = (L) -> (2*(kn^2-k[1]^2)*(1-cos(kn*L)) + k[1]^2*(kn*L)^2)/kn^4 - rnorm
			ΔL = fzero(residual, lower_bound, upper_bound)
		end
	end
	return ΔL
end

"""
	discretize_beam(L, r, ndiv;
		Cab = SMatrix{3,3}(Diagonal(ones(3))),
		k = (@SVector zeros(3)))

Divide a beam into `ndiv` uniformly spaced beam elements.  `L` is the beam length,
`r` is the starting point, `Cab` is the rotation matrix at the starting
point, and `k` is the curvature vector of the beam.

Return the lengths, endpoints, midpoints, and rotation matrices of the beam elements.
"""
function discretize_beam(L, r, ndiv; Cab = I3, k = (@SVector zeros(3)))

	# discretize length into ndiv beam elements
	ΔL = fill(L/ndiv, ndiv)
	sr = cumsum(ΔL)
	sm = sr .- 0.5*L/ndiv
	sp = vcat(0, sr)

	# precompute some curvature quantities
	kkt = k*k'
	ktilde = tilde(k)
	kn = sqrt(k'*k)

	if k == zero(k)
		triads = fill(Cab, ndiv)
		xp = [r + L*i/(ndiv)*Cab*e1 for i = 0:ndiv]
		xm = [r + (L*i/(ndiv) - L/(2*ndiv))*Cab*e1 for i = 1:ndiv]
	else
		triads = curve_triad.(Ref(Cab), Ref(kkt), Ref(ktilde), Ref(kn), sm)
		xp = curve_coordinates.(Ref(r), Ref(Cab), Ref(kkt), Ref(ktilde), Ref(kn), sp)
		xm = curve_coordinates.(Ref(r), Ref(Cab), Ref(kkt), Ref(ktilde), Ref(kn), sm)
	end

	return ΔL, xp, xm, triads
end

curve_triad(Cab, k, s) = curve_triad(Cab, k*k', tilde(k), sqrt(k'*k), s)
curve_triad(Cab, kkt, ktilde, kn, s) = SMatrix{3,3}(Cab*((I - kkt/kn^2)*cos(kn*s) +
	ktilde/kn*sin(kn*s) + kkt/kn^2))

curve_coordinates(r, Cab, k, s) = curve_coordinates(r, Cab, k*k', tilde(k), sqrt(k'*k), s)
curve_coordinates(r, Cab, kkt, ktilde, kn, s) = r + SVector{3}(Cab*((I/kn - kkt/kn^2)*sin(kn*s) + ktilde/kn^2*(1-cos(kn*s)) + kkt/kn^2*s)*e1)
