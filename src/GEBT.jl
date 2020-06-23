module GEBT

"""
    element_equations()

Evaluates nonlinear equations for an element using shape functions of the lowest
possible order in order to avoid integration.

See "GEBT: A general-purpose nonlinear analysis tool for composite beams" by Wenbin Yu

The Rodriguez parameters have been swapped out for Wiener-Milenkovic parameters
"""
function element_equations()

    f_u1 = -C'*Cab*F
    f_u2 =  C'*Cab*F
    f_ψ1 = -C'*Cab*M - ΔL/2*cross(e1 + γ, F)
    f_ψ2 =  C'*Cab*M - ΔL/2*cross(e1 + γ, F)
    f_F1 =  u - ΔL/2*(C'*Cab*(e1 + γ) - Cab*e1)
    f_F2 = -u - ΔL/2*(C'*Cab*(e1 + γ) - Cab*e1)
    f_M1 =  θ - ΔL/2*(θ'*θ*I/16 + θtilde/2 + θ*θ'/8)*Cab*κ
    f_M2 = -θ - ΔL/2*(θ'*θ*I/16 + θtilde/2 + θ*θ'/8)*Cab*κ

    if loads
        f1, f2, m1, m2 = element_loads()
        f_u1 -= f1
        f_u2 -= f2
        f_ψ1 -= m1
        f_ψ2 -= m2
    end

    if !steady
        f_u1 += cross(ωa, ΔL/2*C'*Cab*P)
        f_u2 += cross(ωa, ΔL/2*C'*Cab*P)
        if time_marching
            f_u1 += 2/dt*ΔL/2*C'*Cab*P
            f_u2 += 2/dt*ΔL/2*C'*Cab*P
        end

        f_ψ1 = cross(ωa, ΔL/2*C'*Cab*H) + ΔL/2*C'*Cab*cross(V, P)
        f_ψ2 = cross(ωa, ΔL/2*C'*Cab*H) + ΔL/2*C'*Cab*cross(V, P)
        if time_marching
            m_ψ1 += 2/dt*ΔL/2*C'*Cab*H
            m_ψ2 += 2/dt*ΔL/2*C'*Cab*H
        end

        f_P  = C'*Cab*V - v - cross(ωa, u)
		if time_marching
			f_P -= 2/dt*u
		end

		f_H = Ω - Cba*C*ωa

    end

    if initial_step || time_marching
        f_u1 += ΔL/2*CtCabPdot
        f_u2 += ΔL/2*CtCabPdot
        f_ψ1 += ΔL/2*CtCabHdot
        f_ψ2 += ΔL/2*CtCabHdot

        f_P -= udot
        f_H -= Cab'*(((4 - θ'*θ/4 - 2*θtilde + θ*θ'/2)*θdot)/(4-(2-θ'*θ))^2)
	end

	return f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2, f_P, f_H
end

"""
    boundary_equations()

Evaluates the equations at the start or end of a beam.

See "GEBT: A general-purpose nonlinear analysis tool for composite beams" by Wenbin Yu
"""
function boundary_equations(side, xpt, internal_displacement, internal_force, displacement)

	displacement_residual = internal_displacement + side*xpt
	force_residual = internal_force - side*force

	return displacement_residual, force_residual
end

"""
    connection_equations()

Evaluates the equations at the start or end of a beam.

See "GEBT: A general-purpose nonlinear analysis tool for composite beams" by Wenbin Yu
"""
function connection_equations(side, xpt, internal_displacement, internal_force, force)

	displacement_residual = internal_displacement - side*xpt
	force_residual = internal_force + force

	return displacement_residual, force_residual
end

function system_residual!(b, irow_pnt, iconn, xpt)

	ndof = ifelse(steady, 12, 18)
	initialized = zeros(Bool, npnt)

	for ibeam = 1:nbeam
		# solve equations for the first element
		f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2, f_P, f_H = element_equations()

		# set internal forces/displacements at right side of the initial node
		fp = (f_u1[1], f_u1[2], f_u1[3], f_ψ1[1], f_ψ1[2], f_ψ1[3],
			f_F1[1], f_F1[2], f_F1[3], f_M1[1], f_M1[2], f_M1[3])

		# check if node has been initialized
		if initialized[ifrst[ibeam]]
			# add to equilibrium equations
			irow = irow_pnt[beam.ifrst] - 1
			b[irow+1:irow+12] .+= fp
			# add compatability equations
			irow = irow_pnt[beam.ifrst] - 1
			b[irow+1:irow+6] = fp - xpt
		else
			# add equilibrium equation
			irow = irow_pnt[beam.ifrst] - 1
			b[irow+1:irow+12] = fp + prescribed_condition[ipnt]
		end

		# construct residual equations for this element
		if !steady
			irow =
			b[irow+1:irow+3] .= f_P
			b[irow+4:irow+6] .= f_H
		end

		# move to next node
		fm = (f_u2[1], f_u2[2], f_u2[3], f_ψ2[1], f_ψ2[2], f_ψ2[3],
			f_F2[1], f_F2[2], f_F2[3], f_M2[1], f_M2[2], f_M2[3])

		for i = 2:nelem
			# solve equations for next element
			f_u1, f_u2, f_ψ1, f_ψ2, f_F1, f_F2, f_M1, f_M2, f_P, f_H = element_equations()

			# set internal forces/displacements at right side of the node
			fp = (f_u1[1], f_u1[2], f_u1[3], f_ψ1[1], f_ψ1[2], f_ψ1[3],
				f_F1[1], f_F1[2], f_F1[3], f_M1[1], f_M1[2], f_M1[3])

			# construct residual equations for this node
			b[irow+1:irow+12] .= fm + fp
			irow += 12

			# construct residual equations for this element
			if !steady
				b[irow+1:irow+3] .= f_P
				b[irow+4:irow+6] .= f_H
				irow += 6
			end

			# move to next node/element
			fm = (f_u2[1], f_u2[2], f_u2[3], f_ψ2[1], f_ψ2[2], f_ψ2[3],
				f_F2[1], f_F2[2], f_F2[3], f_M2[1], f_M2[2], f_M2[3])
		end

		# check if node has been initialized
		if initialized[ilast[ibeam]]
			# add to equilibrium equations
			irow =
			b[irow+1:irow+12] .+= fm
			# add compatability equations
			irow =
			b[irow+1:irow+6] = fm
		else
			# add equilibrium equation
			irow =
			b[irow+1:irow+12] = fm + prescribed_condition[ipnt]
		end
	end
end

function element_jacobian()

end

function boundary_jacobian()

end

function connection_jacobian()

end

function system_jacobian()

end




function static_analysis(, niter, nstep)

end

function steady_state_analysis(Ω, V, niter, nstep)

end

function time_domain_analysis(Ω, V, Ωt, Vt, niter, nstep)

end

function eigenvalue_analysis(Ω, V, niter, nstep, nev)

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
