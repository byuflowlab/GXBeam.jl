# unit vectors and identity matrix
const e1 = SVector(1, 0, 0)
const e2 = SVector(0, 1, 0)
const e3 = SVector(0, 0, 1)
const I3 = @SMatrix [1 0 0; 0 1 0; 0 0 1]

"""
    tilde(x)

Constructs the cross product operator
"""
tilde(x) = @SMatrix [0 -x[3] x[2]; x[3] 0 -x[1]; -x[2] x[1] 0]

"""
    wiener_milenkovic

Constructs the Wiener-Milenkovic rotation matrix
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

Returns the gradient of the Wiener-Milenkovic rotation matrix w.r.t. `c`
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

function Q(c)
	c0 = 2 - c'*c/8
	return 1/(4-c0)^2*((4 - 1/4*c'*c)*I - 2*tilde(c) + 1/2*c*c')
end

Q_jacobian(c) = Q_jacobian(Q(c), c)

function Q_jacobian(Q, c)

	c0 = 2 - c'*c/8
	c0_c = -c/4

	tmp = 1/(4-c0)^2
	tmp_c = 2/(4-c0)^3*c0_c

	Q_c1 = tmp_c[1]*Q/tmp + tmp*(-c[1]/2*I - 2*tilde(e1) + (e1*c' + c*e1')/2)

	Q_c2 = tmp_c[2]*Q/tmp + tmp*(-c[2]/2*I - 2*tilde(e2) + (e2*c' + c*e2')/2)

	Q_c3 = tmp_c[3]*Q/tmp + tmp*(-c[3]/2*I - 2*tilde(e3) + (e3*c' + c*e3')/2)

	return Q_c1, Q_c2, Q_c3
end

Qinv(c) =  (1 - 1/16*c'*c)*I3 + 1/2*tilde(c) + 1/8*c*c'

function Qinv_jacobian(c)

	Qinv_c1 = -c[1]/8*I3 + 1/2*tilde(e1) + 1/8*(e1*c' + c*e1')
	Qinv_c2 = -c[2]/8*I3 + 1/2*tilde(e2) + 1/8*(e2*c' + c*e2')
	Qinv_c3 = -c[3]/8*I3 + 1/2*tilde(e3) + 1/8*(e3*c' + c*e3')

	return Qinv_c1, Qinv_c2, Qinv_c3
end

# function curve_length(r1, r2, k)
# 	kn = sqrt(k'*k)
# 	r = r2 - r1
# 	rn = sqrt(r'*r)
# 	if kn == 0 || (k[2] == 0 && k[3] == 0)
# 		ΔL = rn
# 	else
# 		if k[1] == 0
# 			ΔL = 2*asin(kn*rn/2)/kn
# 		else
# 			lower_bound = rn
# 			upper_bound = rn*kn/k
# 			residual = (L) -> (2*(kn^2-k[1]^2)*(1-cos(kn*L)) + k[1]^2*(kn*L)^2)/kn^4 - rnorm
# 			ΔL = fzero(residual, lower_bound, upper_bound)
# 		end
# 	return ΔL
# end
