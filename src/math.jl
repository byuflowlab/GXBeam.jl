# unit vectors and identity matrix
const e1 = SVector(1, 0, 0)
const e2 = SVector(0, 1, 0)
const e3 = SVector(0, 0, 1)
const I3 = @SMatrix [1 0 0; 0 1 0; 0 0 1]

"""
    tilde(x)

Construct the cross product operator matrix
"""
@inline tilde(x) = @SMatrix [0 -x[3] x[2]; x[3] 0 -x[1]; -x[2] x[1] 0]

"""
    wiener_milenkovic(c)

Construct a Wiener-Milenkovic rotation matrix
"""
@inline function wiener_milenkovic(c)

	c0 = 2 - c'*c/8

	return 1/(4-c0)^2*(@SMatrix [
		c0^2 + c[1]^2 - c[2]^2 - c[3]^2  2*(c[1]*c[2] + c0*c[3])          2*(c[1]*c[3]-c0*c[2]);
		2*(c[1]*c[2] - c0*c[3])          c0^2 - c[1]^2 + c[2]^2 - c[3]^2  2*(c[2]*c[3] + c0*c[1]);
		2*(c[1]*c[3] + c0*c[2])          2*(c[2]*c[3] - c0*c[1])          c0^2 - c[1]^2 - c[2]^2 + c[3]^2
		]
	)
end

"""
    get_C(c)

Alias for `wiener_milenkovic(c)`
"""
@inline get_C(c) = wiener_milenkovic(c)

"""
	get_C_t([C, ] c, c_t)

Calculate the derivative of a Wiener-Milenkovic rotation matrix `C` with respect
to the scalar parameter `t`. `c_t` is the derivative of `c` with respect to `t`.
"""
@inline get_C_t(c, c_t) = get_C_t(rotation(c, model), c, c_t)

@inline function get_C_t(C, c, c_t)

	c0 = 2 - c'*c/8
	c0dot = -c'*c_t/4
	tmp = 1/(4-c0)^2
	tmpdot = 2/(4-c0)^3*c0dot

	Cdot = tmpdot*C/tmp + tmp*(@SMatrix [
		2*c0*c0dot + 2*c[1]*c_t[1] - 2*c[2]*c_t[2] - 2*c[3]*c_t[3] (
		2*(c_t[1]*c[2] + c[1]*c_t[2] + c0dot*c[3] + c0*c_t[3])) (
		2*(c_t[1]*c[3] + c[1]*c_t[3] - c0dot*c[2] - c0*c_t[2]));
		2*(c_t[1]*c[2] + c[1]*c_t[2] - c0dot*c[3] - c0*c_t[3]) (
		2*c0*c0dot - 2*c[1]*c_t[1] + 2*c[2]*c_t[2] - 2*c[3]*c_t[3]) (
		2*(c_t[2]*c[3] + c[2]*c_t[3] + c0dot*c[1] + c0*c_t[1]));
		2*(c_t[1]*c[3] + c[1]*c_t[3] + c0dot*c[2] + c0*c_t[2]) (
		2*(c_t[2]*c[3] + c[2]*c_t[3] - c0dot*c[1] - c0*c_t[1])) (
		2*c0*c0dot - 2*c[1]*c_t[1] - 2*c[2]*c_t[2] + 2*c[3]*c_t[3])]
	)

	return Cdot
end

"""
	get_C_c([C, ] c)

Calculate the derivative of the Wiener-Milenkovic rotation matrix `C` with respect
to each of the rotation parameters in `c`.
"""
@inline get_C_c(c) = get_C_c(rotation(c, model), c)

@inline function get_C_c(C, c)

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
    get_C_cdot([C, ] c)

Calculate the derivative of the time derivative of the Wiener-Milenkovic rotation
matrix `C` with respect to each of the time derivatives of `c`. Used for
constructing the "mass" matrix for eigenvalue computations.
"""
get_C_cdot

@inline get_C_cdot(c) = get_C_cdot(rotation(c), c)

@inline function get_C_cdot(C, c)

	c0 = 2 - c'*c/8
	c0dot_cdot = -c/4
	tmp = 1/(4-c0)^2
	tmpdot_cdot = 2/(4-c0)^3*c0dot_cdot

	C_cdot1 = tmpdot_cdot[1]*C/tmp + tmp*(@SMatrix [
		2*c0*c0dot_cdot[1] + 2*c[1]   2*(c[2] + c[3]*c0dot_cdot[1]) 2*(c[3] - c[2]*c0dot_cdot[1]);
		2*(c[2] - c[3]*c0dot_cdot[1]) 2*c0*c0dot_cdot[1] - 2*c[1]   2*(c[1]*c0dot_cdot[1] + c0);
		2*(c[3] + c[2]*c0dot_cdot[1]) 2*(-c[1]*c0dot_cdot[1] - c0)  2*c0*c0dot_cdot[1] - 2*c[1]]
	)

	C_cdot2 = tmpdot_cdot[2]*C/tmp + tmp*(@SMatrix [
		2*c0*c0dot_cdot[2] - 2*c[2]   2*(c[1] + c[3]*c0dot_cdot[2]) 2*(-c[2]*c0dot_cdot[2] - c0);
		2*(c[1] - c[3]*c0dot_cdot[2]) 2*c0*c0dot_cdot[2] + 2*c[2]   2*(c[3] + c[1]*c0dot_cdot[2]);
		2*(c[2]*c0dot_cdot[2] + c0)   2*(c[3] - c[1]*c0dot_cdot[2]) 2*c0*c0dot_cdot[2]- 2*c[2]]
	)

	C_cdot3 = tmpdot_cdot[3]*C/tmp + tmp*(@SMatrix [
		2*c0*c0dot_cdot[3] - 2*c[3] 2*(c[3]*c0dot_cdot[3] + c0) 2*(c[1] - c[2]*c0dot_cdot[3]);
		2*(-c[3]*c0dot_cdot[3] - c0) 2*c0*c0dot_cdot[3] - 2*c[3] 2*(c[2] + c[1]*c0dot_cdot[3]);
		2*(c[1] + c[2]*c0dot_cdot[3]) 2*(c[2] - c[1]*c0dot_cdot[3]) 2*c0*c0dot_cdot[3] + 2*c[3]]
	)

	return C_cdot1, C_cdot2, C_cdot3
end

"""
    get_Q(c)

Calculate the matrix Q as defined in the paper "Geometrically nonlinear analysis
of composite beams using Wiener-Milenković parameters" by Qi Wang and Wenbin Yu.
"""
@inline function get_Q(c)
	c0 = 2 - c'*c/8
	return 1/(4-c0)^2*((4 - 1/4*c'*c)*I - 2*tilde(c) + 1/2*c*c')
end

"""
    get_Q_c(c)
	get_Q_c(Q, c)

Calculate the derivative of the matrix `Q` with respect to each of the rotation
parameters in `c`.
"""
@inline get_Q_c(c) = get_Q_c(get_Q(c), c)

@inline function get_Q_c(Q, c)

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

Calculate the matrix inverse `Qinv` as defined in the paper "Geometrically
nonlinear analysis of composite beams using Wiener-Milenković parameters" by
Qi Wang and Wenbin Yu.
"""
@inline get_Qinv(c) =  (1 - 1/16*c'*c)*I3 + 1/2*tilde(c) + 1/8*c*c'

"""
    get_Qinv_c(c)

Calculate the derivative of the matrix inverse `Qinv` with respect to each of the rotation
parameters in `c`.
"""
@inline function get_Qinv_c(c)

	Qinv_c1 = -c[1]/8*I3 + 1/2*tilde(e1) + 1/8*(e1*c' + c*e1')
	Qinv_c2 = -c[2]/8*I3 + 1/2*tilde(e2) + 1/8*(e2*c' + c*e2')
	Qinv_c3 = -c[3]/8*I3 + 1/2*tilde(e3) + 1/8*(e3*c' + c*e3')

	return Qinv_c1, Qinv_c2, Qinv_c3
end

"""
	mul3(A_1, A_2, A_3, b)

Returns the product of a 3x3x3 tensor represented by `A_1`, `A_2`, and `A_3` with
the vector `b`.
"""
@inline mul3(A_1, A_2, A_3, b) = hcat(A_1*b, A_2*b, A_3*b)
