using LinearAlgebra: Symmetric, det
using FLOWMath: trapz, linear

MaterialPlane(E1, E2, G12, nu12, rho) = MaterialPlane(E1, E2, G12, nu12, rho, 1.0, 1.0, 1.0, 1.0, 1.0)
MaterialPlane(E1, E2, G12, nu12, rho, S1t, S1c, S2t, S2c, S12) = Material(E1, E2, 0.0, G12, 0.0, 0.0, nu12, 0.0, 0.0, rho,
                                                                            S1t, S1c, S2t, S2c, 0.0, 0.0, S12, 0.0, 0.0)

# """
#     Lamina(t, theta, materials)

# A ply

# **Arguments**
# - `t::float`: thickness of ply
# - `theta::float`: corresponding orientation (rad)
# - `material::Material`: corresponding material
# """
# struct Lamina{TF1, TF2, TM}
#     t::TF1
#     theta::TF2
#     material::TM
# end


"""
    Qbar(lamina)

Computes the lamina stiffness matrix at some arbitrary orientation theta
Transforms a specially orthotropic lamina from principal axis to
an arbitrary axis defined by the ply orientation.
"""
function Qbar(lamina)

    mat = lamina.material
    E1 = mat.E1
    E2 = mat.E2
    nu12 = mat.nu12
    nu21 = nu12*E2/E1
    G12 = mat.G12
    delta = 1.0/(1 - nu12*nu21)

    c = cos(lamina.theta)
    s = sin(lamina.theta)
    c2 = c*c
    s2 = s*s
    cs = c*s

    Q = Symmetric([E1*delta  nu12*E2*delta  0.0;
                  nu12*E2*delta  E2*delta  0.0;
                  0.0  0.0  G12])
    Ts = [c2 s2 -2*cs;
          s2 c2 2*cs;
          cs -cs c2-s2]

    Te = [c2 s2 -cs;
          s2 c2 cs;
          2*cs -2*cs c2-s2]

    Qbar = Symmetric(Ts * Q * Ts')

    return Qbar, Ts, Te
end

function zspacing(laminate)

    n = length(laminate)

    # compute z vector
    z = zeros(n+1)
    for i = 1:n
        z[i+1] = z[i] + laminate[i].t
    end
    h = z[end] - z[1]
    mid = h/2.0  # compute midpoint
    z .-= mid  # recenter at midpoint

    return z, h
end

function stiffnessmatrix(laminate, z)
    A = Symmetric(zeros(3, 3))
    B = Symmetric(zeros(3, 3))
    D = Symmetric(zeros(3, 3))

    n = length(laminate)
    for k = 1:n
        Q, _, _ = Qbar(laminate[k])
        A += Q*(z[k+1] - z[k])
        B += Q*(z[k+1]^2 - z[k]^2)/2
        D += Q*(z[k+1]^3 - z[k]^3)/3
    end

    return A, B, D
end

function compliancematrix(A, B, D)
    Ainv = inv(A)
    Hinv = Symmetric(inv(D - B*Ainv*B))
    alpha = Symmetric(Ainv + Ainv*B*Hinv*B*Ainv)
    beta = -Ainv*B*Hinv
    delta = Hinv

    return alpha, beta, delta
end

function laminatecompliance(laminate)
    z, h = zspacing(laminate)
    A, B, D = stiffnessmatrix(laminate, z)
    alpha, beta, delta = compliancematrix(A, B, D)
    return alpha, beta, delta
end

function strains(alpha, beta, delta, z, forces)
    C = [alpha beta; beta' delta]

    alleps = C*forces
    epsilonbar = alleps[1:3]
    kappa = alleps[4:6]

    # setup new z vector at top and bottom of each ply
    nz = 2*(length(z)-1)
    zvec = zeros(nz)
    zvec[1] = z[1]
    zvec[end] = z[end]
    j = 2
    for i = 2:length(z)-1
        zvec[j] = z[i]
        zvec[j+1] = z[i]
        j += 2
    end

    epsilonp = zeros(3, nz)
    for i = 1:3
        epsilonp[i, :] = epsilonbar[i] .+ kappa[i]*zvec
    end

    return epsilonbar, kappa, zvec, epsilonp
end

function stresses(laminate, epsilonp)

    n = length(laminate)
    sigmap = zeros(3, 2*n)
    sigma = zeros(3, 2*n)
    epsilon = zeros(3, 2*n)

    i = 1
    for k = 1:n
        Q, Ts, Te = Qbar(laminate[k])
        # compute for top and bottom
        sigmap[:, i] = Q*epsilonp[:, i]
        sigmap[:, i+1] = Q*epsilonp[:, i+1]
        sigma[:, i] = Te' * sigmap[:, i]
        sigma[:, i+1] = Te' * sigmap[:, i+1]
        epsilon[:, i] = Ts' * epsilonp[:, i]
        epsilon[:, i+1] = Ts' * epsilonp[:, i+1]
        i += 2
    end

    return sigmap, sigma, epsilon
end



function clt(laminate, forces)

    z, h = zspacing(laminate)
    A, B, D = stiffnessmatrix(laminate, z)
    alpha, beta, delta = compliancematrix(A, B, D)
    epsilonbar, kappa, zvec, epsilonp = strains(alpha, beta, delta, z, forces)
    sigmap, sigma, epsilon = stresses(laminate, epsilonp)
    failure = tsai_wu_plane(sigma, laminate)

    return sigma, epsilon, failure
end



function tsai_hill(sigma, strength)
    S1t = strength.S1t
    S1c = strength.S1c
    S2t = strength.S2t
    S2c = strength.S2c
    S12 = strength.S12

    _, n = size(sigma)
    failure = zeros(n)  # fails if > 1
    for i = 1:n
        if sigma[1, i] >= 0.0
            S1 = S1t
        else
            S1 = S1c
        end
        if sigma[2, i] >= 0.0
            S2 = S2t
        else
            S2 = S2c
        end
        failure[i] = sigma[1, i]^2/S1^2 + sigma[2, i]^2/S2^2 + sigma[3, i]^2/S12^2 - sigma[1, i]*sigma[2, i]/S1^2
    end

    return failure
end

function tsai_wu_plane(sigma, laminate)
    # S1t = strength.S1t
    # S1c = strength.S1c
    # S2t = strength.S2t
    # S2c = strength.S2c
    # S12 = strength.S12

    # _, n = size(sigma)
    n = length(laminate)
    failure = zeros(2*n)  # failure if > 1
    i = 1
    for k = 1:n
        (; S1t, S1c, S2t, S2c, S12) = laminate[k].material
        for j = i:i+1
            failure[j] = sigma[1, j]^2/(S1t*S1c) + sigma[2, j]^2/(S2t*S2c) - sigma[1, j]*sigma[2, j]/sqrt(S1t*S1c*S2t*S2c) +
                sigma[1, j]*(1/S1t - 1/S1c) + sigma[2, j]*(1/S2t - 1/S2c) + sigma[3, j]^2/S12^2
        end
        i += 2
    end

    return failure
end

struct BeamSection{VL, VF}
    laminate::VL  #Vector{Lamina}
    y::VF  # Vector{Float}
    z::VF  # Vector{Float}
end

function centroid(beamsections)

    m = length(beamsections)

    ynum = 0.0
    znum = 0.0
    den = 0.0
    for i = 1:m
        sec = beamsections[i]
        yp = sec.y
        zp = sec.z
        n = length(yp)

        alpha, beta, delta = laminatecompliance(sec.laminate)
        D = alpha[1, 1]*delta[1, 1] - beta[1, 1]^2

        s = 0.0
        ybar = 0.0
        zbar = 0.0
        ca = 0.0
        sa = 0.0
        for j = 2:n
            ds = sqrt((yp[j] - yp[j-1])^2 + (zp[j] - zp[j-1])^2)
            ym = (yp[j] + yp[j-1])/2
            zm = (zp[j] + zp[j-1])/2
            s += ds
            ca += (yp[j] - yp[j-1])  # / ds * ds
            sa += (zp[j] - zp[j-1])  # / ds * ds
            ybar += ym * ds
            zbar += zm * ds
        end

        ynum += delta[1, 1]/D*ybar + beta[1, 1]/D*sa
        znum += delta[1, 1]/D*zbar - beta[1, 1]/D*ca
        den += delta[1, 1]/D*s
    end

    yc = ynum/den
    zc = znum/den


    return yc, zc
end

function beamstiffnessold(beamsections)

    yc, zc = centroid(beamsections)

    m = length(beamsections)

    EA = 0.0
    EIyy = 0.0
    EIzz = 0.0
    EIyz = 0.0
    A = 0.0
    GJden = 0.0
    for i = 1:m
        sec = beamsections[i]
        yp = sec.y
        zp = sec.z
        n = length(yp)

        alpha, beta, delta = laminatecompliance(sec.laminate)  # different sections should generally have different laminates for a closed section (otherwise combine into one section), so no loss in efficiency here
        a11hat = alpha[1, 1] - beta[1, 2]^2/delta[2, 2]
        b11hat = beta[1, 1] - beta[1, 2]*delta[1, 2]/delta[2, 2]
        d11hat = delta[1, 1] - delta[1, 2]^2/delta[2, 2]
        Dhat = a11hat*d11hat - b11hat^2
        D = alpha[1, 1]*delta[1, 1] - beta[1, 1]^2
        alpha_nu = alpha[3, 3] - beta[3, 3]^2/delta[3, 3]

        s = 0.0
        z2 = 0.0
        y2 = 0.0
        yz = 0.0
        zca = 0.0
        zsa = 0.0
        yca = 0.0
        ysa = 0.0
        c2 = 0.0
        s2 = 0.0
        cs = 0.0
        Asub = 0.0
        for j = 2:n
            ds = sqrt((yp[j] - yp[j-1])^2 + (zp[j] - zp[j-1])^2)
            ca = (yp[j] - yp[j-1])/ds
            sa = (zp[j] - zp[j-1])/ds
            ym = (yp[j] + yp[j-1])/2
            zm = (zp[j] + zp[j-1])/2
            y = ym - yc
            z = zm - zc
            s += ds
            z2 += z^2 * ds + ds^3/12*sa^2
            y2 += y^2 * ds + ds^3/12*ca^2
            yz += y*z * ds + ds^3/12*sa*ca
            zca += z*ca * ds
            zsa += z*sa * ds
            yca += y*ca * ds
            ysa += y*sa * ds
            c2 += ca*ca * ds
            s2 += sa*sa * ds
            cs += ca*sa * ds
            Asub += (zp[j-1] + zp[j]) * (yp[j-1] - yp[j])  # if closed section
        end

        EA += d11hat/Dhat * s  # for a closed section
        EIyy += delta[1, 1]/D*z2 - 2*beta[1, 1]/D*zca + alpha[1, 1]/D*c2
        EIzz += delta[1, 1]/D*y2 + 2*beta[1, 1]/D*ysa + alpha[1, 1]/D*s2
        EIyz += delta[1, 1]/D*yz + beta[1, 1]/D*(zsa - yca) - alpha[1, 1]/D*cs
        GJden += s*alpha_nu
        A += abs(Asub)/2.0
    end

    GJ = 4*A^2/GJden

    return EA, EIyy, EIzz, EIyz, GJ
end

# this one is for an arbitrary beam. not necessarily orthotropic
function beamstiffness(sections; closedsection=true)

    m = length(sections)

    Pbar = zeros(4, 4)
    I = zeros(2, 4)
    F = zeros(2, 2)
    A = 0.0

    for i = 1:m
        sec = sections[i]
        yp = sec.y
        zp = sec.z
        n = length(yp)

        alpha, beta, delta = laminatecompliance(sec.laminate)
        a = Symmetric([alpha[1, 1] beta[1, 1] beta[1, 3];
                 beta[1, 1] delta[1, 1] delta[1, 3];
                 beta[1, 3] delta[1, 3] delta[3, 3]])
        atinv = (a[2, 2]*a[3, 3] - a[2, 3]^2) / det(a)
        # Atilde = inv(atilde)  # TODO: we only need (1, 1) component so don't need to invert everything

        Asub = 0.0
        for k = 2:n
            ybar = (yp[k-1] + yp[k])/2
            zbar = (zp[k-1] + zp[k])/2
            b = sqrt((yp[k] - yp[k-1])^2 + (zp[k] - zp[k-1])^2)
            ca = (yp[k] - yp[k-1])/b
            sa = (zp[k] - zp[k-1])/b
            Asub += (zp[k-1] + zp[k]) * (yp[k-1] - yp[k])  # if closed section

            Rk = [1.0 zbar ybar 0.0;
                0.0 ca -sa 0.0;
                0.0 sa ca 0.0;
                0 0 0 1]
            omega = 1.0/b*Symmetric([
                    alpha[1, 1] beta[1, 1] 0.0 -beta[1, 3]/2.0;
                    beta[1, 1] delta[1, 1] 0.0 -delta[1, 3]/2.0;
                    0 0 12.0/(atinv*b^2) 0;
                    -beta[1, 3]/2.0 -delta[1, 3]/2.0 0 delta[3, 3]/4.0])
            # omegainv = inv(omega)
            Pbar += Rk'*(omega\Rk)

            I1 = [alpha[1, 3] beta[3, 1] 0.0 -beta[3, 3]/2.0;
                  beta[1, 2] delta[1, 2] 0.0 -delta[2, 3]/2.0]
            I += I1*(omega\Rk)  # repeated, could cache

            F1 = [alpha[3, 3] beta[3, 2];
                 beta[3, 2] delta[2, 2]]
            F += b*F1 - I1*(omega\I1')

        end

        if i <= 10   #TODO: temporary hack
            A += abs(Asub)/2.0
        end
    end
    L = -I
    L[1, 4] += 2*A
    if closedsection
        Pbar += L'*(F\L)
    end
    Wbar = inv(Symmetric(Pbar))
    cent = -Symmetric([Wbar[2, 2] Wbar[2, 3]; Wbar[2, 3] Wbar[3, 3]]) \ [Wbar[1, 2]; Wbar[1, 3]]
    zc = cent[1]; yc = cent[2]

    Rb = [1.0 0 0 0;
          zc 1 0 0;
          yc 0 1 0;
          0 0 0 1]
    W = Symmetric(Rb'*Wbar*Rb)
    P = inv(W)

    # s, S, ysc, zsc = shearflow(sections, Wbar, F, L, yc, zc)
    return W, P, yc, zc
end

function qopen(eta, y1, y2, z1, z2, m1, m2)
    ys = y1 + eta*(y2 - y1)
    zs = z1 + eta*(z2 - z1)
    ybarsub = (y1 + ys)/2
    zbarsub = (z1 + zs)/2
    bsub = sqrt((ys - y1)^2 + (zs - z1)^2)
    qo = m1*zbarsub*bsub + m2*ybarsub*bsub
    return qo
end

function shearflow(sections, P, yc, zc; closedsection=true, npts=20)

    EIyy = P[2, 2]
    EIzz = P[3, 3]
    EIyz = P[2, 3]

    eta = range(0.0, 1.0, length=npts)

    function factors(Vy, Vz)
        factor1 = (-EIzz*Vz + EIyz*Vy)/(EIyy*EIzz - EIyz^2)
        factor2 = (EIyz*Vz - EIyy*Vy)/(EIyy*EIzz - EIyz^2)
        return factor1, factor2
    end
    factor1y, factor2y = factors(1.0, 0.0)
    factor1z, factor2z = factors(0.0, 1.0)


    m = length(sections)


    qcy = 0.0
    qcz = 0.0
    if closedsection
        # ------ compute qc -------------
        qoprevy = 0.0
        qoprevz = 0.0
        qcy_num = 0.0
        qcy_den = 0.0
        qcz_num = 0.0
        qcz_den = 0.0
        for i = 1:m
            sec = sections[i]
            yp = sec.y .- yc  # relative to centroid
            zp = sec.z .- zc
            n = length(yp)

            alpha, beta, delta = laminatecompliance(sec.laminate)
            a11 = alpha[1, 1]
            a_nu = alpha[3, 3] - beta[3, 3]^2/delta[3, 3]

            m1y = factor1y/a11
            m2y = factor2y/a11
            m1z = factor1z/a11
            m2z = factor2z/a11

            for k = 2:n
                b = sqrt((yp[k] - yp[k-1])^2 + (zp[k] - zp[k-1])^2)

                qo(s, qoprev, m1, m2) = qoprev + qopen(s, yp[k-1], yp[k], zp[k-1], zp[k], m1, m2)

                # integraly, _ = quadgk(s -> qo(s, qoprevy, m1y, m2y), 0.0, b)
                # integralz, _ = quadgk(s -> qo(s, qoprevz, m1z, m2z), 0.0, b)
                integraly = b*trapz(eta, qo.(eta, qoprevy, m1y, m2y))
                integralz = b*trapz(eta, qo.(eta, qoprevz, m1z, m2z))
                qcy_num += a_nu*integraly
                qcy_den += a_nu*b
                qcz_num += a_nu*integralz
                qcz_den += a_nu*b

                qoprevy = qo(1.0, qoprevy, m1y, m2y)
                qoprevz = qo(1.0, qoprevz, m1z, m2z)
            end
        end

        qcy = -qcy_num/qcy_den
        qcz = -qcz_num/qcz_den
    end


    # ------ compute ysc, zsc, s, S -------------
    qoprevy = 0.0
    qoprevz = 0.0
    ysc = 0.0
    zsc = 0.0
    syy = 0.0
    szz = 0.0
    syz = 0.0

    for i = 1:m
        sec = sections[i]
        yp = sec.y .- yc  # relative to centroid
        zp = sec.z .- zc
        n = length(yp)

        alpha, beta, delta = laminatecompliance(sec.laminate)
        a11 = alpha[1, 1]
        a_nu = alpha[3, 3] - beta[3, 3]^2/delta[3, 3]

        m1y = factor1y/a11
        m2y = factor2y/a11
        m1z = factor1z/a11
        m2z = factor2z/a11

        for k = 2:n
            b = sqrt((yp[k] - yp[k-1])^2 + (zp[k] - zp[k-1])^2)

            qo(s, qoprev, m1, m2) = qoprev + qopen(s, yp[k-1], yp[k], zp[k-1], zp[k], m1, m2)
            qy(s) = qcy + qo(s, qoprevy, m1y, m2y)
            qz(s) = qcz + qo(s, qoprevz, m1z, m2z)
            function p(s)
                ys = yp[k-1] + s*(yp[k] - yp[k-1])
                zs = zp[k-1] + s*(zp[k] - zp[k-1])
                return sqrt(ys^2 + zs^2)
            end

            # integraly, _ = quadgk(s -> qy(s)*p(s), 0.0, b)
            # integralz, _ = quadgk(s -> qz(s)*p(s), 0.0, b)
            integraly = b*trapz(eta, @. qy(eta)*p(eta))
            integralz = b*trapz(eta, @. qz(eta)*p(eta))

            zsc -= integraly
            ysc += integralz


            # integral, _ = quadgk(s -> qy(s)^2, 0.0, b)
            integral = b*trapz(eta, @. qy(eta)^2)
            syy += a_nu*integral
            # integral, _ = quadgk(s -> qz(s)^2, 0.0, b)
            integral = b*trapz(eta, @. qz(eta)^2)
            szz += a_nu*integral
            # integral, _ = quadgk(s -> qy(s)*qz(s), 0.0, b)
            integral = b*trapz(eta, @. qy(eta)*qz(eta))
            syz += a_nu*integral

            qoprevy = qo(1.0, qoprevy, m1y, m2y)
            qoprevz = qo(1.0, qoprevz, m1z, m2z)
        end
    end

    s = Symmetric([syy syz;
                   syz szz])
    S = inv(s)

    return s, S, ysc, zsc
end



function shearflow_general_attempt(sections, Wbar, F, L, yc, zc)

    m = length(sections)

    qoy = 0.0
    qoz = 0.0
    qcy = 0.0
    qcz = 0.0
    eta = 0.0
    AMy = zeros(2, 2)
    bvy = zeros(2)
    AMz = zeros(2, 2)
    bvz = zeros(2)
    for i = 1:m
        sec = sections[i]
        yp = sec.y
        zp = sec.z
        n = length(yp)

        alpha, beta, delta = laminatecompliance(sec.laminate)
        muk = Symmetric([alpha[1, 1] beta[1, 1] beta[1, 3];
                            beta[1, 1] delta[1, 1] delta[1, 3]
                            beta[1, 3] delta[1, 3] delta[3, 3]])
        nuk = [alpha[1, 3] beta[1, 2]
            beta[3, 1] delta[1, 2]
            beta[3, 3] delta[2, 3]]

        for k = 2:n
            ybar = (yp[k-1] + yp[k])/2
            zbar = (zp[k-1] + zp[k])/2
            b = sqrt((yp[k] - yp[k-1])^2 + (zp[k] - zp[k-1])^2)
            ca = (yp[k] - yp[k-1])/b
            sa = (zp[k] - zp[k-1])/b

            # open section shear flow
            eta += b
            Rk = [1.0 zbar ybar 0.0;
                0.0 ca -sa 0.0;
                0.0 sa ca 0.0;
                0 0 0 1]
            Reta = [1.0 0 eta 0;
                    0 1 0 0;
                    0 0 0 -2]

            M = (muk\(Reta*Rk - nuk*(F\L)))*Wbar

            # --------
            Vy = 1.0; Vz = 0.0
            qoy = shearflowsub!(M, Vy, Vz, b, alpha, beta, delta, qoy, AMy, bvy)

            Vy = 0.0; Vz = 1.0
            qoz = shearflowsub!(M, Vy, Vz, b, alpha, beta, delta, qoz, AMz, bvz)
        end
    end

    qcy = 0.0  #(AMy\bvy)[1]
    qcz = 0.0  #(AMz\bvz)[1]

    eta = 0.0
    qoy = 0.0
    qoz = 0.0
    ysc = 0.0
    zsc = 0.0
    syy = 0.0
    szz = 0.0
    syz = 0.0
    for i = 1:m
        sec = sections[i]
        yp = sec.y
        zp = sec.z
        n = length(yp)

        alpha, beta, delta = laminatecompliance(sec.laminate)
        muk = Symmetric([alpha[1, 1] beta[1, 1] beta[1, 3];
                            beta[1, 1] delta[1, 1] delta[1, 3]
                            beta[1, 3] delta[1, 3] delta[3, 3]])
        nuk = [alpha[1, 3] beta[1, 2]
            beta[3, 1] delta[1, 2]
            beta[3, 3] delta[2, 3]]

        alpha_nu = alpha[3, 3] - beta[3, 3]^2/delta[3, 3]

        for k = 2:n
            ybar = (yp[k-1] + yp[k])/2
            zbar = (zp[k-1] + zp[k])/2
            b = sqrt((yp[k] - yp[k-1])^2 + (zp[k] - zp[k-1])^2)
            ca = (yp[k] - yp[k-1])/b
            sa = (zp[k] - zp[k-1])/b

            # open section shear flow
            eta += b
            Rk = [1.0 zbar ybar 0.0;
                0.0 ca -sa 0.0;
                0.0 sa ca 0.0;
                0 0 0 1]
            Reta = [1.0 0 eta 0;
                    0 1 0 0;
                    0 0 0 -2]

            M = (muk\(Reta*Rk - nuk*(F\L)))*Wbar

            Vy = 1.0; Vz = 0.0
            a = M*[0.0; Vz; Vy; 0]
            dN = a[1]
            qoy -= dN*b
            qy = qoy + qcy

            Vy = 0.0; Vz = 1.0
            a = M*[0.0; Vz; Vy; 0]
            dN = a[1]
            qoz -= dN*b
            qz = qoz + qcz

            # shear center
            p = sqrt((ybar - yc)^2 + (zbar - zc)^2)
            ysc += qz*p*b
            zsc -= qy*p*b

            # stiffnesses
            syy += alpha_nu*qy^2*b
            szz += alpha_nu*qz^2*b
            syz += alpha_nu*qy*qz*b
        end
    end

    s = Symmetric([syy syz;
                   syz szz])
    S = inv(s)

    return s, S, ysc, zsc
end

function shearflowsub!(M, Vy, Vz, b, alpha, beta, delta, qo, AM, bv)

    a = M*[0.0; Vz; Vy; 0]
    dN = a[1]
    qo -= dN*b

    # closed section due to Vy/Vz
    AM .-= b*Symmetric([alpha[3, 3] beta[2, 3]
                      beta[2, 3] delta[2, 2]])
    bv .+= qo*b*[alpha[3, 3]; beta[2, 3]]

    return qo
end

function fullstiffnessmatrix(P, S)
    K = zeros(6, 6)
    K[2:3, 2:3] .= S
    idx = [1, 5, 6, 4]
    j = 1
    for i in idx
        K[i, idx] .= P[j, :]
        j += 1
    end

    return K
end

function rotatestiffnessmatrix(K, r, theta)
    R = [1.0 0 0;
        0 cos(theta) -sin(theta);
        0 sin(theta) cos(theta)]
    p = [0.0 -r[3] r[2];
        r[3] 0 -r[1];
        -r[2] r[1] 0]
    HinvT = [R p*R; zeros(3, 3) R]
    Hinv = transpose(HinvT)

    Kp = Hinv*K*HinvT
    return Kp
end


# # ---------- tests ------------
# using Test

# E1 = 129000.0e6  # note must be typo on this book
# E2 = 11000.0e6
# G12 = 6600.0e6
# nu12 = 0.28
# rho = 1.0
# mat = OrthotropicMaterial(E1, E2, nu12, G12, rho, "carbonfiber")

# t = 8e-3
# theta = [45 -45 0 90 90 0 -45 45]*pi/180
# laminate = Lamina.(t/8, theta, Ref(mat))

# Q1, _, _ = Qbar(laminate[1])
# @test isapprox(Q1[1, 1]/1e6, 43385.924, atol=1e-3)
# @test isapprox(Q1[1, 2]/1e6, 30185.924, atol=1e-3)
# @test isapprox(Q1[1, 3]/1e6, 29698.543, atol=1e-3)
# @test isapprox(Q1[2, 2]/1e6, 43385.924, atol=1e-3)
# @test isapprox(Q1[2, 3]/1e6, 29698.543, atol=1e-3)
# @test isapprox(Q1[3, 3]/1e6, 33685.195, atol=1e-3)
# @test Q1[1, 2] == Q1[2, 1]
# @test Q1[1, 3] == Q1[3, 1]
# @test Q1[2, 3] == Q1[3, 2]

# Q2, _, _ = Qbar(laminate[2])
# @test isapprox(Q2[1, 1]/1e6, 43385.924, atol=1e-3)
# @test isapprox(Q2[1, 2]/1e6, 30185.924, atol=1e-3)
# @test isapprox(Q2[1, 3]/1e6, -29698.543, atol=1e-3)
# @test isapprox(Q2[2, 2]/1e6, 43385.924, atol=1e-3)
# @test isapprox(Q2[2, 3]/1e6, -29698.543, atol=1e-3)
# @test isapprox(Q2[3, 3]/1e6, 33685.195, atol=1e-3)
# @test Q2[1, 2] == Q2[2, 1]
# @test Q2[1, 3] == Q2[3, 1]
# @test Q2[2, 3] == Q2[3, 2]

# Q3, _, _ = Qbar(laminate[3])
# @test isapprox(Q3[1, 1]/1e6, 129868.204, atol=1e-3)
# @test isapprox(Q3[1, 2]/1e6, 3100.729, atol=1e-3)
# @test isapprox(Q3[1, 3]/1e6, 0.0, atol=1e-3)
# @test isapprox(Q3[2, 2]/1e6, 11074.033, atol=1e-3)
# @test isapprox(Q3[2, 3]/1e6, 0.0, atol=1e-3)
# @test isapprox(Q3[3, 3]/1e6, 6600.0, atol=1e-3)

# Q4, _, _ = Qbar(laminate[4])
# @test isapprox(Q4[1, 1]/1e6, 11074.033, atol=1e-3)
# @test isapprox(Q4[1, 2]/1e6, 3100.729, atol=1e-3)
# @test isapprox(Q4[1, 3]/1e6, 0.0, atol=1e-3)
# @test isapprox(Q4[2, 2]/1e6, 129868.204, atol=1e-3)
# @test isapprox(Q4[2, 3]/1e6, 0.0, atol=1e-3)
# @test isapprox(Q4[3, 3]/1e6, 6600.0, atol=1e-3)

# z, h = zspacing(laminate)
# A, B, D = stiffnessmatrix(laminate, z)

# @test isapprox(A[1, 1]/1e3, 455428.170, atol=1e-3)
# @test isapprox(A[1, 2]/1e3, 133146.612, atol=1e-3)
# @test isapprox(A[1, 3]/1e3, 0.0, atol=1e-3)
# @test isapprox(A[2, 2]/1e3, 455428.170, atol=1e-3)
# @test isapprox(A[2, 3]/1e3, 0.0, atol=1e-3)
# @test isapprox(A[3, 3]/1e3, 161140.779, atol=1e-3)

# @test isapprox(B[1, 1], 0.0, atol=1e-10)
# @test isapprox(B[1, 2], 0.0, atol=1e-10)
# @test isapprox(B[1, 3], 0.0, atol=1e-10)
# @test isapprox(B[2, 2], 0.0, atol=1e-10)
# @test isapprox(B[2, 3], 0.0, atol=1e-10)
# @test isapprox(B[3, 3], 0.0, atol=1e-10)

# @test isapprox(D[1, 1]*1e3, 2233175.466, atol=1e-3)
# @test isapprox(D[1, 2]*1e3, 1143478.381, atol=1e-3)
# @test isapprox(D[1, 3]*1e3, 356382.514, atol=1e-3)
# @test isapprox(D[2, 2]*1e3, 1757998.781, atol=1e-3)
# @test isapprox(D[2, 3]*1e3, 356382.514, atol=1e-3)
# @test isapprox(D[3, 3]*1e3, 1292780.601, atol=1e-3)

# alpha, beta, delta = compliancematrix(A, B, D)

# @test isapprox(alpha[1, 1]/1e-10, 24.009482, atol=1e-6)
# @test isapprox(alpha[1, 2]/1e-10, -7.019287, atol=1e-6)
# @test isapprox(alpha[1, 3]/1e-10, 0.0, atol=1e-6)
# @test isapprox(alpha[2, 2]/1e-10, 24.009482, atol=1e-6)
# @test isapprox(alpha[2, 3]/1e-10, 0.0, atol=1e-6)
# @test isapprox(alpha[3, 3]/1e-10, 62.057538, atol=1e-6)

# @test isapprox(beta[1, 1], 0.0, atol=1e-10)
# @test isapprox(beta[1, 2], 0.0, atol=1e-10)
# @test isapprox(beta[1, 3], 0.0, atol=1e-10)
# @test isapprox(beta[2, 2], 0.0, atol=1e-10)
# @test isapprox(beta[2, 3], 0.0, atol=1e-10)
# @test isapprox(beta[3, 3], 0.0, atol=1e-10)

# @test isapprox(delta[1, 1]/1e-4, 6.771890, atol=1e-6)
# @test isapprox(delta[1, 2]/1e-4, -4.264613, atol=1e-6)
# @test isapprox(delta[1, 3]/1e-4, -0.691184, atol=1e-6)
# @test isapprox(delta[2, 2]/1e-4, 8.710637, atol=1e-6)
# @test isapprox(delta[2, 3]/1e-4, -1.225641, atol=1e-6)
# @test isapprox(delta[3, 3]/1e-4, 8.263678, atol=1e-6)

# Nx = 1000e3
# forces = [Nx; 0.0; 0.0; 0.0; 0.0; 0.0]

# epsilonbar, kappa, zvec, epsilonp = strains(alpha, beta, delta, z, forces)

# @test isapprox(epsilonbar[1]*1e3, 2.400948, atol=1e-6)
# @test isapprox(epsilonbar[2]*1e3, -0.701929, atol=1e-6)
# @test isapprox(epsilonbar[3]*1e3, 0.0, atol=1e-10)
# @test isapprox(kappa[1], 0.0, atol=1e-10)
# @test isapprox(kappa[2], 0.0, atol=1e-10)
# @test isapprox(kappa[3], 0.0, atol=1e-10)


# sigmap, sigma, epsilon = stresses(laminate, epsilonp)

# using PyPlot
# close("all"); pygui(true)

# figure()
# plot(sigmap[1, :]/1e6, zvec)
# figure()
# plot(epsilonp[1, :]/1e6, zvec)

# @test isapprox(sigma[1, 1]/1e6, 112.958, atol=1e-3)
# @test isapprox(sigma[1, 2]/1e6, 112.958, atol=1e-3)
# @test isapprox(sigma[1, 3]/1e6, 112.958, atol=1e-3)
# @test isapprox(sigma[1, 4]/1e6, 112.958, atol=1e-3)
# @test isapprox(sigma[1, 5]/1e6, 309.630, atol=1e-3)
# @test isapprox(sigma[1, 6]/1e6, 309.630, atol=1e-3)
# @test isapprox(sigma[1, 7]/1e6, -83.714, atol=1e-3)
# @test isapprox(sigma[1, 8]/1e6, -83.714, atol=1e-3)
# @test isapprox(sigma[1, 9]/1e6, -83.714, atol=1e-3)
# @test isapprox(sigma[1, 10]/1e6, -83.714, atol=1e-3)
# @test isapprox(sigma[1, 11]/1e6, 309.630, atol=1e-3)
# @test isapprox(sigma[1, 12]/1e6, 309.630, atol=1e-3)
# @test isapprox(sigma[1, 13]/1e6, 112.958, atol=1e-3)
# @test isapprox(sigma[1, 14]/1e6, 112.958, atol=1e-3)
# @test isapprox(sigma[1, 15]/1e6, 112.958, atol=1e-3)
# @test isapprox(sigma[1, 16]/1e6, 112.958, atol=1e-3)

# @test isapprox(sigma[2, 1]/1e6, 12.042, atol=1e-3)
# @test isapprox(sigma[2, 2]/1e6, 12.042, atol=1e-3)
# @test isapprox(sigma[2, 3]/1e6, 12.042, atol=1e-3)
# @test isapprox(sigma[2, 4]/1e6, 12.042, atol=1e-3)
# @test isapprox(sigma[2, 5]/1e6, -0.328, atol=1e-3)
# @test isapprox(sigma[2, 6]/1e6, -0.328, atol=1e-3)
# @test isapprox(sigma[2, 7]/1e6, 24.412, atol=1e-3)
# @test isapprox(sigma[2, 8]/1e6, 24.412, atol=1e-3)
# @test isapprox(sigma[2, 9]/1e6, 24.412, atol=1e-3)
# @test isapprox(sigma[2, 10]/1e6, 24.412, atol=1e-3)
# @test isapprox(sigma[2, 11]/1e6, -0.328, atol=1e-3)
# @test isapprox(sigma[2, 12]/1e6, -0.328, atol=1e-3)
# @test isapprox(sigma[2, 13]/1e6, 12.042, atol=1e-3)
# @test isapprox(sigma[2, 14]/1e6, 12.042, atol=1e-3)
# @test isapprox(sigma[2, 15]/1e6, 12.042, atol=1e-3)
# @test isapprox(sigma[2, 16]/1e6, 12.042, atol=1e-3)

# @test isapprox(sigma[3, 1]/1e6, -20.479, atol=1e-3)
# @test isapprox(sigma[3, 2]/1e6, -20.479, atol=1e-3)
# @test isapprox(sigma[3, 3]/1e6, 20.479, atol=1e-3)
# @test isapprox(sigma[3, 4]/1e6, 20.479, atol=1e-3)
# @test isapprox(sigma[3, 5]/1e6, 0.0, atol=1e-3)
# @test isapprox(sigma[3, 6]/1e6, 0.0, atol=1e-3)
# @test isapprox(sigma[3, 7]/1e6, 0.0, atol=1e-3)
# @test isapprox(sigma[3, 8]/1e6, 0.0, atol=1e-3)
# @test isapprox(sigma[3, 9]/1e6, 0.0, atol=1e-3)
# @test isapprox(sigma[3, 10]/1e6, 0.0, atol=1e-3)
# @test isapprox(sigma[3, 11]/1e6, 0.0, atol=1e-3)
# @test isapprox(sigma[3, 12]/1e6, 0.0, atol=1e-3)
# @test isapprox(sigma[3, 13]/1e6, 20.479, atol=1e-3)
# @test isapprox(sigma[3, 14]/1e6, 20.479, atol=1e-3)
# @test isapprox(sigma[3, 15]/1e6, -20.479, atol=1e-3)
# @test isapprox(sigma[3, 16]/1e6, -20.479, atol=1e-3)

# @test isapprox(epsilon[1, 1]*1e3, 0.849510, atol=1e-6)
# @test isapprox(epsilon[1, 2]*1e3, 0.849510, atol=1e-6)
# @test isapprox(epsilon[1, 3]*1e3, 0.849510, atol=1e-6)
# @test isapprox(epsilon[1, 4]*1e3, 0.849510, atol=1e-6)
# @test isapprox(epsilon[1, 5]*1e3, 2.400948, atol=1e-6)
# @test isapprox(epsilon[1, 6]*1e3, 2.400948, atol=1e-6)
# @test isapprox(epsilon[1, 7]*1e3, -0.701929, atol=1e-6)
# @test isapprox(epsilon[1, 8]*1e3, -0.701929, atol=1e-6)
# @test isapprox(epsilon[1, 9]*1e3, -0.701929, atol=1e-6)
# @test isapprox(epsilon[1, 10]*1e3, -0.701929, atol=1e-6)
# @test isapprox(epsilon[1, 11]*1e3, 2.400948, atol=1e-6)
# @test isapprox(epsilon[1, 12]*1e3, 2.400948, atol=1e-6)
# @test isapprox(epsilon[1, 13]*1e3, 0.849510, atol=1e-6)
# @test isapprox(epsilon[1, 14]*1e3, 0.849510, atol=1e-6)
# @test isapprox(epsilon[1, 15]*1e3, 0.849510, atol=1e-6)
# @test isapprox(epsilon[1, 16]*1e3, 0.849510, atol=1e-6)

# @test isapprox(epsilon[2, 1]*1e3, 0.849510, atol=1e-6)
# @test isapprox(epsilon[2, 2]*1e3, 0.849510, atol=1e-6)
# @test isapprox(epsilon[2, 3]*1e3, 0.849510, atol=1e-6)
# @test isapprox(epsilon[2, 4]*1e3, 0.849510, atol=1e-6)
# @test isapprox(epsilon[2, 5]*1e3, -0.701929, atol=1e-6)
# @test isapprox(epsilon[2, 6]*1e3, -0.701929, atol=1e-6)
# @test isapprox(epsilon[2, 7]*1e3, 2.400948, atol=1e-6)
# @test isapprox(epsilon[2, 8]*1e3, 2.400948, atol=1e-6)
# @test isapprox(epsilon[2, 9]*1e3, 2.400948, atol=1e-6)
# @test isapprox(epsilon[2, 10]*1e3, 2.400948, atol=1e-6)
# @test isapprox(epsilon[2, 11]*1e3, -0.701929, atol=1e-6)
# @test isapprox(epsilon[2, 12]*1e3, -0.701929, atol=1e-6)
# @test isapprox(epsilon[2, 13]*1e3, 0.849510, atol=1e-6)
# @test isapprox(epsilon[2, 14]*1e3, 0.849510, atol=1e-6)
# @test isapprox(epsilon[2, 15]*1e3, 0.849510, atol=1e-6)
# @test isapprox(epsilon[2, 16]*1e3, 0.849510, atol=1e-6)

# @test isapprox(epsilon[3, 1]*1e3, -3.102877, atol=1e-6)
# @test isapprox(epsilon[3, 2]*1e3, -3.102877, atol=1e-6)
# @test isapprox(epsilon[3, 3]*1e3, 3.102877, atol=1e-6)
# @test isapprox(epsilon[3, 4]*1e3, 3.102877, atol=1e-6)
# @test isapprox(epsilon[3, 5]*1e3, 0.0, atol=1e-6)
# @test isapprox(epsilon[3, 6]*1e3, 0.0, atol=1e-6)
# @test isapprox(epsilon[3, 7]*1e3, 0.0, atol=1e-6)
# @test isapprox(epsilon[3, 8]*1e3, 0.0, atol=1e-6)
# @test isapprox(epsilon[3, 9]*1e3, 0.0, atol=1e-6)
# @test isapprox(epsilon[3, 10]*1e3, 0.0, atol=1e-6)
# @test isapprox(epsilon[3, 11]*1e3, 0.0, atol=1e-6)
# @test isapprox(epsilon[3, 12]*1e3, 0.0, atol=1e-6)
# @test isapprox(epsilon[3, 13]*1e3, 3.102877, atol=1e-6)
# @test isapprox(epsilon[3, 14]*1e3, 3.102877, atol=1e-6)
# @test isapprox(epsilon[3, 15]*1e3, -3.102877, atol=1e-6)
# @test isapprox(epsilon[3, 16]*1e3, -3.102877, atol=1e-6)

# forces = [0.0; 0.0; 0.0; 1000; 0.0; 0.0]

# epsilonbar, kappa, zvec, epsilonp = strains(alpha, beta, delta, z, forces)

# @test isapprox(epsilonbar[1], 0.0, atol=1e-6)
# @test isapprox(epsilonbar[2], 0.0, atol=1e-6)
# @test isapprox(epsilonbar[3], 0.0, atol=1e-10)
# @test isapprox(kappa[1], 0.677189, atol=1e-6)
# @test isapprox(kappa[2], -0.426461, atol=1e-6)
# @test isapprox(kappa[3], -0.069118, atol=1e-6)

# sigmap, sigma, epsilon = stresses(laminate, epsilonp)


# figure()
# plot(sigmap[1, :]/1e6, zvec*1e3)
# figure()
# plot(epsilonp[1, :]*1e3, zvec*1e3)

# @test isapprox(sigma[1, 1]/1e6, -49.154, atol=1e-3)
# @test isapprox(sigma[1, 2]/1e6, -36.866, atol=1e-3)
# @test isapprox(sigma[1, 3]/1e6, -63.151, atol=1e-3)
# @test isapprox(sigma[1, 4]/1e6, -42.101, atol=1e-3)
# @test isapprox(sigma[1, 5]/1e6, -173.246, atol=1e-3)
# @test isapprox(sigma[1, 6]/1e6, -86.623, atol=1e-3)
# @test isapprox(sigma[1, 7]/1e6, 53.284, atol=1e-3)
# @test isapprox(sigma[1, 8]/1e6, 0.0, atol=1e-3)
# @test isapprox(sigma[1, 9]/1e6, 0.0, atol=1e-3)
# @test isapprox(sigma[1, 10]/1e6, -53.284, atol=1e-3)
# @test isapprox(sigma[1, 11]/1e6, 86.623, atol=1e-3)
# @test isapprox(sigma[1, 12]/1e6, 173.246, atol=1e-3)
# @test isapprox(sigma[1, 13]/1e6, 42.101, atol=1e-3)
# @test isapprox(sigma[1, 14]/1e6, 63.151, atol=1e-3)
# @test isapprox(sigma[1, 15]/1e6, 36.866, atol=1e-3)
# @test isapprox(sigma[1, 16]/1e6, 49.154, atol=1e-3)

# @test isapprox(sigma[2, 1]/1e6, -8.210, atol=1e-3)
# @test isapprox(sigma[2, 2]/1e6, -6.158, atol=1e-3)
# @test isapprox(sigma[2, 3]/1e6, -4.504, atol=1e-3)
# @test isapprox(sigma[2, 4]/1e6, -3.003, atol=1e-3)
# @test isapprox(sigma[2, 5]/1e6, 5.246, atol=1e-3)
# @test isapprox(sigma[2, 6]/1e6, 2.623, atol=1e-3)
# @test isapprox(sigma[2, 7]/1e6, -6.177, atol=1e-3)
# @test isapprox(sigma[2, 8]/1e6, 0.0, atol=1e-3)
# @test isapprox(sigma[2, 9]/1e6, 0.0, atol=1e-3)
# @test isapprox(sigma[2, 10]/1e6, 6.177, atol=1e-3)
# @test isapprox(sigma[2, 11]/1e6, -2.623, atol=1e-3)
# @test isapprox(sigma[2, 12]/1e6, -5.246, atol=1e-3)
# @test isapprox(sigma[2, 13]/1e6, 3.003, atol=1e-3)
# @test isapprox(sigma[2, 14]/1e6, 4.504, atol=1e-3)
# @test isapprox(sigma[2, 15]/1e6, 6.158, atol=1e-3)
# @test isapprox(sigma[2, 16]/1e6, 8.210, atol=1e-3)

# @test isapprox(sigma[3, 1]/1e6, 29.136, atol=1e-3)
# @test isapprox(sigma[3, 2]/1e6, 21.852, atol=1e-3)
# @test isapprox(sigma[3, 3]/1e6, -21.852, atol=1e-3)
# @test isapprox(sigma[3, 4]/1e6, -14.568, atol=1e-3)
# @test isapprox(sigma[3, 5]/1e6, 0.912, atol=1e-3)
# @test isapprox(sigma[3, 6]/1e6, 0.456, atol=1e-3)
# @test isapprox(sigma[3, 7]/1e6, -0.456, atol=1e-3)
# @test isapprox(sigma[3, 8]/1e6, 0.0, atol=1e-3)
# @test isapprox(sigma[3, 9]/1e6, 0.0, atol=1e-3)
# @test isapprox(sigma[3, 10]/1e6, 0.456, atol=1e-3)
# @test isapprox(sigma[3, 11]/1e6, -0.456, atol=1e-3)
# @test isapprox(sigma[3, 12]/1e6, -0.912, atol=1e-3)
# @test isapprox(sigma[3, 13]/1e6, 14.568, atol=1e-3)
# @test isapprox(sigma[3, 14]/1e6, 21.852, atol=1e-3)
# @test isapprox(sigma[3, 15]/1e6, -21.852, atol=1e-3)
# @test isapprox(sigma[3, 16]/1e6, -29.136, atol=1e-3)


# # ----- 4.3 --------
# theta = [45 -45 0 90 0 90 45 -45]*pi/180
# laminate = Lamina.(t/8, theta, Ref(mat))

# forces = [1000*1e3; 0.0; 0.0; 0.0; 0.0; 0.0]

# z, h = zspacing(laminate)
# A, B, D = stiffnessmatrix(laminate, z)
# alpha, beta, delta = compliancematrix(A, B, D)
# epsilonbar, kappa, zvec, epsilonp = strains(alpha, beta, delta, z, forces)
# sigmap, sigma, epsilon = stresses(laminate, epsilonp)

# @test isapprox(A[1, 1]/1e3, 455428.170, atol=1e-3)
# @test isapprox(A[1, 2]/1e3, 133146.612, atol=1e-3)
# @test isapprox(A[1, 3]/1e3, 0.0, atol=1e-3)
# @test isapprox(A[2, 2]/1e3, 455428.170, atol=1e-3)
# @test isapprox(A[2, 3]/1e3, 0.0, atol=1e-3)
# @test isapprox(A[3, 3]/1e3, 161140.779, atol=1e-3)

# @test isapprox(B[1, 1], -118794.171, atol=1e-3)
# @test isapprox(B[1, 2], 1.455e-11, atol=1e-3)
# @test isapprox(B[1, 3], -59397.086, atol=1e-3)
# @test isapprox(B[2, 2], 118794.171, atol=1e-3)
# @test isapprox(B[2, 3], -59397.086, atol=1e-3)
# @test isapprox(B[3, 3], 0.0, atol=1e-10)

# @test isapprox(D[1, 1]*1e3, 1995587.124, atol=1e-3)
# @test isapprox(D[1, 2]*1e3, 1143478.381, atol=1e-3)
# @test isapprox(D[1, 3]*1e3, 0.0, atol=1e-3)
# @test isapprox(D[2, 2]*1e3, 1995587.124, atol=1e-3)
# @test isapprox(D[2, 3]*1e3, 0.0, atol=1e-3)
# @test isapprox(D[3, 3]*1e3, 1292780.601, atol=1e-3)

# @test isapprox(alpha[1, 1]/1e-10, 24.562273, atol=1e-6)
# @test isapprox(alpha[1, 2]/1e-10, -6.911749, atol=1e-6)
# @test isapprox(alpha[1, 3]/1e-10, 0.445254, atol=1e-6)
# @test isapprox(alpha[2, 2]/1e-10, 24.562273, atol=1e-6)
# @test isapprox(alpha[2, 3]/1e-10, -0.445254, atol=1e-6)
# @test isapprox(alpha[3, 3]/1e-10, 62.948045, atol=1e-6)

# @test isapprox(beta[1, 1]/1e-7, 1.834321, atol=1e-6)
# @test isapprox(beta[1, 2]/1e-7, -0.626374, atol=1e-6)
# @test isapprox(beta[1, 3]/1e-7, 0.810957, atol=1e-6)
# @test isapprox(beta[2, 2]/1e-7, -1.834321, atol=1e-6)
# @test isapprox(beta[2, 3]/1e-7, 0.810957, atol=1e-6)
# @test isapprox(beta[3, 3]/1e-7, 0.0, atol=1e-6)

# @test isapprox(delta[1, 1]/1e-4, 7.677865, atol=1e-6)
# @test isapprox(delta[1, 2]/1e-4, -4.400777, atol=1e-6)
# @test isapprox(delta[1, 3]/1e-4, 0.113057, atol=1e-6)
# @test isapprox(delta[2, 2]/1e-4, 7.677865, atol=1e-6)
# @test isapprox(delta[2, 3]/1e-4, -0.113057, atol=1e-6)
# @test isapprox(delta[3, 3]/1e-4, 7.809784, atol=1e-6)

# @test isapprox(epsilonbar[1]*1e3, 2.456227, atol=1e-6)
# @test isapprox(epsilonbar[2]*1e3, -0.691175, atol=1e-6)
# @test isapprox(epsilonbar[3]*1e3, 0.044525, atol=1e-6)
# @test isapprox(kappa[1], 0.183432, atol=1e-6)
# @test isapprox(kappa[2], -0.062637, atol=1e-6)
# @test isapprox(kappa[3], 0.081096, atol=1e-6)

# @test isapprox(sigma[1, 1]/1e6, 67.486, atol=1e-3)
# @test isapprox(sigma[1, 2]/1e6, 80.657, atol=1e-3)
# @test isapprox(sigma[1, 3]/1e6, 105.854, atol=1e-3)
# @test isapprox(sigma[1, 4]/1e6, 108.745, atol=1e-3)
# @test isapprox(sigma[1, 5]/1e6, 269.587, atol=1e-3)
# @test isapprox(sigma[1, 6]/1e6, 293.215, atol=1e-3)
# @test isapprox(sigma[1, 7]/1e6, -74.580, atol=1e-3)
# @test isapprox(sigma[1, 8]/1e6, -82.146, atol=1e-3)
# @test isapprox(sigma[1, 9]/1e6, 316.843, atol=1e-3)
# @test isapprox(sigma[1, 10]/1e6, 340.470, atol=1e-3)
# @test isapprox(sigma[1, 11]/1e6, -89.711, atol=1e-3)
# @test isapprox(sigma[1, 12]/1e6, -97.277, atol=1e-3)
# @test isapprox(sigma[1, 13]/1e6, 146.513, atol=1e-3)
# @test isapprox(sigma[1, 14]/1e6, 159.684, atol=1e-3)
# @test isapprox(sigma[1, 15]/1e6, 123.199, atol=1e-3)
# @test isapprox(sigma[1, 16]/1e6, 126.090, atol=1e-3)

# @test isapprox(sigma[2, 1]/1e6, 10.201, atol=1e-3)
# @test isapprox(sigma[2, 2]/1e6, 10.734, atol=1e-3)
# @test isapprox(sigma[2, 3]/1e6, 9.149, atol=1e-3)
# @test isapprox(sigma[2, 4]/1e6, 10.328, atol=1e-3)
# @test isapprox(sigma[2, 5]/1e6, 0.212, atol=1e-3)
# @test isapprox(sigma[2, 6]/1e6, 0.087, atol=1e-3)
# @test isapprox(sigma[2, 7]/1e6, 23.220, atol=1e-3)
# @test isapprox(sigma[2, 8]/1e6, 25.057, atol=1e-3)
# @test isapprox(sigma[2, 9]/1e6, -0.038, atol=1e-3)
# @test isapprox(sigma[2, 10]/1e6, -0.163, atol=1e-3)
# @test isapprox(sigma[2, 11]/1e6, 26.894, atol=1e-3)
# @test isapprox(sigma[2, 12]/1e6, 28.731, atol=1e-3)
# @test isapprox(sigma[2, 13]/1e6, 13.398, atol=1e-3)
# @test isapprox(sigma[2, 14]/1e6, 13.931, atol=1e-3)
# @test isapprox(sigma[2, 15]/1e6, 16.225, atol=1e-3)
# @test isapprox(sigma[2, 16]/1e6, 17.405, atol=1e-3)



# # ------- 4.4 --------

# theta = [45 0 30 -45]*pi/180
# laminate = Lamina.(t/4, theta, Ref(mat))

# forces = [1000*1e3; 500e3; 0.0; 0.0; 0.0; 0.0]

# S1t = 1950.0e6
# S1c = 1480.0e6
# S2t = 48.0e6
# S2c = 200.0e6
# S12 = 79.0e6
# strength = CompositeStrength(S1t, S1c, S2t, S2c, S12)

# sigma, epsilon, wu = clt(laminate, forces, strength)

# @test isapprox(sigma[1, 1]/1e6, 18.776, atol=1e-3)
# # @test isapprox(sigma[1, 2]/1e6, 18.776, atol=1e-3)  # typo in text, repeated value
# @test isapprox(sigma[1, 3]/1e6, 238.813, atol=1e-3)
# @test isapprox(sigma[1, 4]/1e6, 229.609, atol=1e-3)
# @test isapprox(sigma[1, 5]/1e6, 149.134, atol=1e-3)
# @test isapprox(sigma[1, 6]/1e6, 229.851, atol=1e-3)
# @test isapprox(sigma[1, 7]/1e6, 214.113, atol=1e-3)
# @test isapprox(sigma[1, 8]/1e6, 13.857, atol=1e-3)

# @test isapprox(epsilon[3, 1]*1e3, 2.663312, atol=1e-6)
# @test isapprox(epsilon[3, 2]*1e3, 1.790924, atol=1e-6)
# @test isapprox(epsilon[3, 3]*1e3, -4.138210, atol=1e-6)
# @test isapprox(epsilon[3, 4]*1e3, -1.996388, atol=1e-6)
# @test isapprox(epsilon[3, 5]*1e3, -0.202719, atol=1e-6)
# @test isapprox(epsilon[3, 6]*1e3, 0.112682, atol=1e-6)
# @test isapprox(epsilon[3, 7]*1e3, -0.046148, atol=1e-6)
# @test isapprox(epsilon[3, 8]*1e3, 0.826241, atol=1e-6)


# failure = tsai_hill(sigma, strength)
# R = sqrt.(1.0 ./ failure)

# @test isapprox(R[1], 0.683, atol=1e-3)
# @test isapprox(R[2], 0.880, atol=1e-3)
# @test isapprox(R[3], 1.001, atol=1e-3)
# @test isapprox(R[4], 1.346, atol=1e-3)
# @test isapprox(R[5], 1.214, atol=1e-3)
# @test isapprox(R[6], 1.999, atol=1e-3)
# @test isapprox(R[7], 1.929, atol=1e-3)
# @test isapprox(R[8], 1.828, atol=1e-3)

# R = [0.681; 0.901; 1.041; 1.467; 1.296; 2.325; 2.216; 1.836]
# sigmaR = copy(sigma)
# sigmaR[1, :] .*= R
# sigmaR[2, :] .*= R
# sigmaR[3, :] .*= R
# failure = tsai_wu(sigmaR, strength)

# @test isapprox(failure[1], 1.0, atol=1e-2)
# @test isapprox(failure[2], 1.0, atol=1e-2)
# @test isapprox(failure[3], 1.0, atol=1e-2)
# @test isapprox(failure[4], 1.0, atol=1e-2)
# @test isapprox(failure[5], 1.0, atol=1e-2)
# @test isapprox(failure[6], 1.0, atol=1e-2)
# @test isapprox(failure[7], 1.0, atol=1e-2)
# @test isapprox(failure[8], 1.0, atol=1e-2)

# forces = [0.0; 0.0; 0.0; 1000; 500; 0.0]

# sigma, epsilon, wu = clt(laminate, forces, strength)

# failure = tsai_hill(sigma, strength)
# R = sqrt.(1.0 ./ failure)

# @test isapprox(R[1], 5.324, atol=1e-2)
# @test isapprox(R[2], 9.990, atol=3e-3)
# # @test isapprox(R[3], 5.973, atol=1e-3)  # typo?  all stresses and strains match
# @test isapprox(R[4], 10.557, atol=1e-3)
# @test isapprox(R[5], 20.539, atol=1e-3)
# @test isapprox(R[6], 9.977, atol=1e-3)
# @test isapprox(R[7], 4.698, atol=1e-3)
# @test isapprox(R[8], 2.523, atol=1e-3)

# R = [6.727; 10.909; 9.786; 14.021; 19.030; 10.114; 4.463; 2.602]
# sigmaR = copy(sigma)
# sigmaR[1, :] .*= R
# sigmaR[2, :] .*= R
# sigmaR[3, :] .*= R
# failure = tsai_wu(sigmaR, strength)

# @test isapprox(failure[1], 1.0, atol=1e-2)
# @test isapprox(failure[2], 1.0, atol=1e-2)
# # @test isapprox(failure[3], 1.0, atol=1e-2)  # again there appears an error in the text for this station.  stresses checkout.
# @test isapprox(failure[4], 1.0, atol=1e-2)
# @test isapprox(failure[5], 1.0, atol=1e-2)
# @test isapprox(failure[6], 1.0, atol=1e-2)
# @test isapprox(failure[7], 1.0, atol=1e-2)
# @test isapprox(failure[8], 1.0, atol=1e-2)



# # ------- Kollar example 6.2 ------

# E1 = 148e9
# E2 = 9.65e9
# G12 = 4.55e9
# nu12 = 0.3
# rho = 1.0
# m0 = OrthotropicMaterial(E1, E2, nu12, G12, rho, "m0")
# t0 = 0.1e-3

# E1 = 16.39e9
# E2 = 16.39e9
# G12 = 38.19e9
# nu12 = 0.801
# rho = 1.0
# m45 = OrthotropicMaterial(E1, E2, nu12, G12, rho, "m45")
# t45 = 0.2e-3


# t = [t45*ones(2); t0*ones(12); t45*ones(2)]
# theta = [0; 0; zeros(12); 0; 0]*pi/180
# mat = [fill(m45, 2); fill(m0, 12); fill(m45, 2)]
# laminate = Lamina.(t, theta, mat)

# # alpha, beta, delta = laminatecompliance(laminate)

# # b1 = BeamSection(laminate, [0.0; 50e-3], [61e-3; 61e-3])
# # b2 = BeamSection(laminate, [25e-3; 25e-3], [0; 60e-3])
# b1 = BeamSection(laminate, [0.0; 50e-3], [0.0; 0])
# b2 = BeamSection(laminate, [25e-3; 25e-3], [-1e-3; -61e-3])
# sections = [b1; b2]

# yc, zc = centroid(sections)

# @test isapprox(yc, 25e-3, atol=1e-4)
# @test isapprox(61e-3+zc, 0.0441, atol=1e-4)

# EA, EIyy, EIzz, EIyz, GJ = beamstiffnessold(sections)

# @test isapprox(EA/1e6, 21.22, atol=1e-2)  # same for open/closed if B = 0
# @test isapprox(EIyy/1e3, 8.530, atol=1e-3)


# # ----------- example 6.3 ------------

# b1 = BeamSection(laminate, [0.0; 52e-3], [0.0; 0.0])
# b2 = BeamSection(laminate, [51e-3; 51e-3], [1e-3; 69e-3])
# b3 = BeamSection(laminate, [52e-3; 0.0], [70e-3; 70e-3])
# b4 = BeamSection(laminate, [1e-3; 1e-3], [69e-3; 1e-3])
# sections = [b1; b2; b3; b4]
# EA, EIyy, EIzz, EIyz, GJ = beamstiffnessold(sections)

# @test isapprox(EIyy/1e3, 34.692, atol=1e-3)
# @test isapprox(EIzz/1e3, 20.924, atol=1e-3)

# # the above gives exact with book for EI because we defined geometry consistent with the formulas, but is approximate with GJ.  to match GJ exactly need closed path below.
# b = BeamSection(laminate, [0.0; 50e-3; 50e-3; 0.0; 0.0], [0.0; 0.0; 70e-3; 70e-3; 0.0])
# EA, EIyy, EIzz, EIyz, GJ = beamstiffnessold([b])
# @test isapprox(EIyy/1e3, 34.692, atol=1e-1)  # looser tolerances for these - exact above.
# @test isapprox(EIzz/1e3, 20.924, atol=1e-1)  # looser tolerances for these
# @test isapprox(GJ/1e3, 7.352, atol=1e-3)

# W, P, _ = beamstiffness(sections)
# @test isapprox(P[1, 1]/1e6, 46.303, atol=1e-3)
# @test isapprox(P[2, 2]/1e3, 34.692, atol=1e-3)
# @test isapprox(P[3, 3]/1e3, 20.924, atol=1e-3)
# @test isapprox(P[1, 2]/1e3, 0.0, atol=1e-6)
# @test isapprox(P[1, 3]/1e3, 0.0, atol=1e-6)
# @test isapprox(P[1, 4]/1e3, 0.0, atol=1e-6)
# @test isapprox(P[2, 1]/1e3, 0.0, atol=1e-6)
# @test isapprox(P[2, 3]/1e3, 0.0, atol=1e-6)
# @test isapprox(P[2, 4]/1e3, 0.0, atol=1e-6)
# @test isapprox(P[3, 4]/1e3, 0.0, atol=1e-6)
# W, P, _ = beamstiffness([b])

# @test isapprox(P[4, 4]/1e3, 7.3, atol=1e-1)  # looser tolerance - Area is a bit ambiguous


# # ---------- example 6.5 ------------------

# t = [t0*ones(10); t0*ones(10)]
# theta = [zeros(10); 45*ones(10)]*pi/180
# mat = [fill(m0, 10); fill(m0, 10)]
# laminate = Lamina.(t, theta, mat)

# b = BeamSection(laminate, [0.0; 50e-3; 50e-3; 0.0; 0.0], [0.0; 0.0; 70e-3; 70e-3; 0.0])
# # b = BeamSection(laminate, [-25; 25; 25; -25; -25]*1e-3, [-35; -35; 35; 35; -35]*1e-3)

# W, P = beamstiffness([b])
# @test isapprox(W[1, 1]/1e-6, 0.02576, atol=1e-5)
# @test isapprox(W[1, 4]/1e-6, -0.4237, atol=1e-4)
# @test isapprox(W[2, 2]/1e-6, 33.91, atol=1e-2)
# @test isapprox(W[3, 3]/1e-6, 55.63, atol=1e-2)
# @test isapprox(W[4, 4]/1e-6, 250.56, atol=1e-2)

# # TODO: add checks for yc/zc - although these are explicitly checked by EA


# # # (just an internal check on shear flow)
# # # -------- example 6.2 -------
# # t = [t45*ones(2); t0*ones(12); t45*ones(2)]
# # theta = [0; 0; zeros(12); 0; 0]*pi/180
# # mat = [fill(m45, 2); fill(m0, 12); fill(m45, 2)]
# # laminate = Lamina.(t, theta, mat)

# # b1 = BeamSection(laminate, [25e-3; 25e-3], [-60e-3; 0.0])
# # b2 = BeamSection(laminate, [0.0; 50e-3], [0.0; 0.0])
# # sections = [b1; b2]
# # W, P, yc, zc = beamstiffness(sections)
# # s, S, ysc, zsc = shearflow(sections, P, yc, zc)



# # ----- example 6.6 / Appendix A.8 ------
# t = [t45*ones(2); t0*ones(12); t45*ones(2)]
# theta = [0; 0; zeros(12); 0; 0]*pi/180
# mat = [fill(m45, 2); fill(m0, 12); fill(m45, 2)]
# laminate = Lamina.(t, theta, mat)

# df = 49e-3
# d = 62e-3
# b = BeamSection(laminate, [0.0; df; df; 0.0], [0.0; 0.0; d; d])
# sections = [b]
# W, P, yc, zc = beamstiffness(sections, closedsection=false)
# s, S, ysc, zsc = shearflow(sections, P, yc, zc, closedsection=false)

# a11 = 5.18e-9
# bf = 50e-3
# e = 3*bf^2/a11 / ( (6*bf + d)/a11 )
# ysca = e + df - yc  # TODO: shear center does not seem to agree.  book may be using the approximate formula though...

# alpha, beta, delta = laminatecompliance(laminate)
# alpha_nu = alpha[3, 3] - beta[3, 3]^2/delta[3, 3]
# deltay = (df - yc)/yc
# rhoy = 3/5*(8 - 9*deltay + 3*deltay^2)/(2 - deltay)^2
# qyu = 1/(2*df)*(3 - 3*deltay)/(2 - deltay)
# gamma = 1 + 1/6*d/df
# syy = rhoy/2*alpha_nu/df + d/3*qyu^2*alpha_nu
# szz = alpha_nu/d + 2/3*alpha_nu*df/(d*gamma)^2
# syz = 0.0

# @test isapprox(syy/1e-7, s[1, 1]/1e-7, atol=0.01)
# @test isapprox(szz/1e-7, s[2, 2]/1e-7, atol=0.1)
# @test isapprox(syz/1e-7, s[1, 2]/1e-7, atol=1e-6)

# # ---- appendix A.7
# df = 2.0
# d = 4.0
# b = BeamSection(laminate, [0.0; df; df; 0.0; 0.0], [0.0; 0.0; d; d; 0.0])
# W, P, yc, zc = beamstiffness([b])
# s, S, ysc, zsc = shearflow([b], P, yc, zc)

# gammaz = 1.0 + 1/3*d/df
# gammay = 1.0 + 1/3*df/d
# szz = alpha_nu/(2*d) + alpha_nu*df/(6*d^2*gammaz^2)
# syy = alpha_nu/(2*df) + alpha_nu*d/(6*df^2*gammay^2)
# syz = 0.0

# @test isapprox(syy/1e-8, s[1, 1]/1e-8, atol=0.01)
# @test isapprox(szz/1e-8, s[2, 2]/1e-8, atol=0.05)
# @test isapprox(syz/1e-8, s[1, 2]/1e-8, atol=1e-6)


# # -------- Precomp test ------------
# using PreComp

# chord = 1.0
# twist = 0.0
# twist_prime = 0.0
# le_loc = 0.0
# xnode = [0.0, 0.00116, 0.0083, 0.0206, 0.0377, 0.0592, 0.0847, 0.11409, 0.14685, 0.18266, 0.22111, 0.26177, 0.30418, 0.34829, 0.39439, 0.44237, 0.49169, 0.54177, 0.59199, 0.64174, 0.69037, 0.78169, 0.82312, 0.86095, 0.8946, 0.9238, 0.94879, 0.96963, 0.98582, 0.99632, 1.0, 0.99606, 0.98364, 0.96197, 0.93154, 0.89348, 0.84887, 0.79901, 0.74529, 0.68912, 0.63189, 0.57485, 0.51909, 0.46552, 0.41483, 0.36753, 0.32394, 0.28389, 0.24712, 0.21322, 0.18122, 0.15072, 0.12193, 0.0952, 0.0708, 0.0492, 0.031, 0.0164, 0.00607, 0.00048]
# ynode = [0.0, 0.00674, 0.0185, 0.0307, 0.043, 0.0551, 0.0667, 0.0777, 0.0877, 0.0965, 0.103972, 0.109823, 0.1136479, 0.1149801, 0.113861, 0.1108122, 0.1062507, 0.1004443, 0.0937, 0.0861, 0.0697, 0.0611, 0.0526, 0.0441, 0.0357, 0.0271, 0.0186, 0.0108, 0.00473, 0.00113, 0.0, 0.00101, 0.00378, 0.00717, 0.00961, 0.0102, 0.00824, 0.00345, -0.00428, -0.0148, -0.0276, -0.0423, -0.0579, -0.0736, -0.0886, -0.1017592, -0.1121818, -0.1189823, -0.1217499, -0.1202246, -0.114725, -0.1063324, -0.0958, -0.0837, -0.0706, -0.0569, -0.043, -0.0293, -0.0162, -0.00435]
# E1 = [3.7e10, 1.03e10, 10.0, 1.03e10, 1.0e7]
# E2 = [9.0e9, 1.03e10, 10.0, 1.03e10, 1.0e7]
# G12 = [4.0e9, 8.0e9, 1.0, 8.0e9, 200000.0]
# nu12 = [0.28, 0.3, 0.3, 0.3, 0.3]
# rho = [1860.0, 1830.0, 1830.0, 1664.0, 128.0]
# xsec_nodeU = [0.0, 1.0]
# n_laminaU = [7]
# n_pliesU = [1, 1, 17, 38, 0, 37, 16]
# t_lamU = [0.000381, 0.00051, 0.00053, 0.00053, 0.003125, 0.00053, 0.00053]
# tht_lamU = [0.0, 0.0, 20.0, 30.0, 0.0, 30.0, 20.0]
# mat_lamU = [3, 4, 2, 1, 5, 1, 2]
# xsec_nodeL =  [0.0, 1.0]
# n_laminaL = n_laminaU
# n_pliesL = n_pliesU
# t_lamL = t_lamU
# tht_lamL = tht_lamU
# mat_lamL = mat_lamU
# loc_web= [0.0]
# n_laminaW = [0]
# n_pliesW = [0]
# t_lamW = [0.0]
# tht_lamW = [0.0]
# mat_lamW = [0]

# input = PreComp.Input(chord, twist, twist_prime, le_loc, xnode, ynode, E1, E2, G12, nu12, rho,
#     xsec_nodeU, n_laminaU, n_pliesU, t_lamU, tht_lamU, mat_lamU,
#     xsec_nodeL, n_laminaL, n_pliesL, t_lamL, tht_lamL, mat_lamL,
#     loc_web, n_laminaW, n_pliesW, t_lamW, tht_lamW, mat_lamW)

# output = PreComp.properties(input)

# # println("ei_flap = ", output.ei_flap)
# # println("ei_lag = ", output.ei_lag)
# # println("gj = ", output.gj)
# # println("ea = ", output.ea)
# # println("s_fl = ", output.s_fl)
# # println("s_af = ", output.s_af)
# # println("s_al = ", output.s_al)
# # println("s_ft = ", output.s_ft)
# # println("s_lt = ", output.s_lt)
# # println("s_at = ", output.s_at)
# # println("x_sc = ", output.x_sc)
# # println("y_sc = ", output.y_sc)
# # println("x_tc = ", output.x_tc)
# # println("y_tc = ", output.y_tc)
# # println("mass = ", output.mass)
# # println("flap_iner = ", output.flap_iner)
# # println("lag_iner = ", output.lag_iner)
# # println("tw_iner_d = ", output.tw_iner_d)
# # println("x_cm = ", output.x_cm)
# # println("y_cm = ", output.y_cm)

# xu = [0.0, 0.00116, 0.0083, 0.0206, 0.0377, 0.0592, 0.0847, 0.11409, 0.14685, 0.18266, 0.22111, 0.26177, 0.30418, 0.34829, 0.39439, 0.44237, 0.49169, 0.54177, 0.59199, 0.64174, 0.69037, 0.78169, 0.82312, 0.86095, 0.8946, 0.9238, 0.94879, 0.96963, 0.98582, 0.99632, 1.0]
# xl = [0.0, 0.00048, 0.00607, 0.0164, 0.031, 0.0492, 0.0708, 0.0952, 0.12193, 0.15072, 0.18122, 0.21322, 0.24712, 0.28389, 0.32394, 0.36753, 0.41483, 0.46552, 0.51909, 0.57485, 0.63189, 0.68912, 0.74529, 0.79901, 0.84887, 0.89348, 0.93154, 0.96197, 0.98364, 0.99606, 1.0]
# yu = [0.0, 0.00674, 0.0185, 0.0307, 0.043, 0.0551, 0.0667, 0.0777, 0.0877, 0.0965, 0.103972, 0.109823, 0.1136479, 0.1149801, 0.113861, 0.1108122, 0.1062507, 0.1004443, 0.0937, 0.0861, 0.0697, 0.0611, 0.0526, 0.0441, 0.0357, 0.0271, 0.0186, 0.0108, 0.00473, 0.00113, 0.0]
# yl = [0.0, -0.00435, -0.0162, -0.0293, -0.043, -0.0569, -0.0706, -0.0837, -0.0958, -0.1063324, -0.114725, -0.1202246, -0.1217499, -0.1189823, -0.1121818, -0.1017592, -0.0886, -0.0736, -0.0579, -0.0423, -0.0276, -0.0148, -0.00428, 0.00345, 0.00824, 0.0102, 0.00961, 0.00717, 0.00378, 0.00101, 0.0]
# xaf = [xl; xu[end-1:-1:2]]
# yaf = [yl; yu[end-1:-1:2]]

# allmats = OrthotropicMaterial.(E1, E2, nu12, G12, rho, "name")
# mats = allmats[mat_lamU]
# laminate = Lamina.(n_pliesU.*t_lamU, tht_lamU*pi/180, mats)
# b = BeamSection(laminate, xaf, yaf)
# W, P, yc, zc = beamstiffness([b])
# EA2, EIyy2, EIzz2, EIyz2, GJ2 = beamstiffnessold([b])
# s, S, ysc, zsc = shearflow([b], P, yc, zc, npts=2)

# println(P[1, 1], " ", output.ea, " ", abs(P[1, 1] - output.ea) / output.ea * 100, "% error")
# println(P[2, 2], " ", output.ei_flap, " ", abs(P[2, 2] - output.ei_flap) / output.ei_flap * 100, "% error")
# println(P[3, 3], " ", output.ei_lag, " ", abs(P[3, 3] - output.ei_lag) / output.ei_lag * 100, "% error")
# println(P[4, 4], " ", output.gj, " ", abs(P[4, 4] - output.gj) / output.gj * 100, "% error")
# println(P[2, 3], " ", output.s_fl, " ", abs(P[2, 3] - output.s_fl) / output.s_fl * 100, "% error")
# println(P[1, 4], " ", output.s_at, " ", abs(P[1, 4] - output.s_at) / output.s_at * 100, "% error")
# println(yc, " ", output.y_tc, " ", (yc - output.y_tc) / output.y_tc * 100, "% error")
# println(zc, " ", output.x_tc, " ", (zc - output.x_tc) / output.x_tc * 100, "% error")

# println(P[1, 2], " ", output.s_af, " ", P[1, 2] / output.s_af)
# println(P[1, 3], " ", output.s_al, " ", P[1, 3] / output.s_al)
# println(P[2, 4], " ", output.s_ft, " ", P[2, 4] / output.s_ft)
# println(P[3, 4], " ", output.s_lt, " ", P[3, 4] / output.s_lt)
# println(yc+ysc, " ", output.y_sc, " ", (yc+ysc) / output.y_sc)
# println(zc+zsc, " ", output.x_sc, " ", (zc+zsc) / output.x_sc)


# # ----------- analyswift paper -----------------
# # -------- circular ring case ---------

# R = 0.3
# E = 73e9
# nu = 0.33
# rho = 2800.0
# tratio = 1.0/15

# t = tratio*2*R
# G = E/(2*(1 + nu))

# n = 400
# theta = range(0, 2*pi, length=n)
# y = (R-t/2)*cos.(theta)
# z = (R-t/2)*sin.(theta)

# mat = OrthotropicMaterial(E, E, nu, G, rho)
# laminate = Lamina(t, 0.0, mat)
# b = BeamSection([laminate], y, z)
# W, P, yc, zc = beamstiffness([b])

# @test isapprox(P[1, 1], 5.137e9, rtol=.001)
# @test isapprox(P[2, 2], 2.024e8, rtol=.005)
# @test isapprox(P[3, 3], 2.024e8, rtol=.005)
# @test isapprox(P[4, 4], 1.553e8, rtol=.02)


# theta = reverse(range(0, pi, length=201))
# yu = 0.5*cos.(theta)
# zu = 0.5*sin.(theta)
# theta = reverse(range(-pi, 0, length=201))
# yl = 0.5*cos.(theta)
# zl = 0.5*sin.(theta)
# y = [yu; yl[2:end-1]]
# z = [zu; zl[2:end-1]]

# y .-= y[1]

# chord = 2*(R-t/2)
# twist = 0.0
# twist_prime = 0.0
# le_loc = 0.0
# xnode = y
# ynode = z
# E1 = [E]
# E2 = [E]
# G12 = [G]
# nu12 = [nu]
# rhov = [rho]
# xsec_nodeU = [0.0, 1.0]
# n_laminaU = [1]
# n_pliesU = [1]
# t_lamU = [t]
# tht_lamU = [0.0]
# mat_lamU = [1]
# xsec_nodeL =  xsec_nodeU
# n_laminaL = n_laminaU
# n_pliesL = n_pliesU
# t_lamL = t_lamU
# tht_lamL = tht_lamU
# mat_lamL = mat_lamU
# loc_web= [0.0]
# n_laminaW = [0]
# n_pliesW = [0]
# t_lamW = [0.0]
# tht_lamW = [0.0]
# mat_lamW = [0]

# input = PreComp.Input(chord, twist, twist_prime, le_loc, xnode, ynode, E1, E2, G12, nu12, rhov,
#     xsec_nodeU, n_laminaU, n_pliesU, t_lamU, tht_lamU, mat_lamU,
#     xsec_nodeL, n_laminaL, n_pliesL, t_lamL, tht_lamL, mat_lamL,
#     loc_web, n_laminaW, n_pliesW, t_lamW, tht_lamW, mat_lamW)

# output = PreComp.properties(input)

# # output.ei_flap + (output.x_sc-output.x_tc)^2*output.ea

# println("---- errors -------")
# EAexact = 5.137e9
# EAvabs = 5.128e9
# println("EA (mine): ", round((P[1, 1]/EAexact-1)*100, digits=3), "%, EA (PreComp): ", round((output.ea/EAexact-1)*100, digits=3), "%, EA (VABS): ", round((EAvabs/EAexact-1)*100, digits=3), "%")
# EIexact = 2.024e8
# EIvabs = 2.022e8
# println("EI (mine): ", round((P[2, 2]/EIexact-1)*100, digits=2), "%, EI (PreComp): ", round((output.ei_flap/EIexact-1)*100, digits=2), "%, EI (VABS): ", round((EIvabs/EIexact-1)*100, digits=2), "%")
# GJexact = 1.553e8
# GJvabs = 1.544e8
# println("GJ (mine): ", round((P[4, 4]/GJexact-1)*100, digits=2), "%, GJ (PreComp): ", round((output.gj/GJexact-1)*100, digits=2), "%, GJ (VABS): ", round((GJvabs/GJexact-1)*100, digits=2), "%")


# tratio = 1.0/3

# t = tratio*2*R

# n = 400
# theta = range(0, 2*pi, length=n)
# y = (R-t/2)*cos.(theta)
# z = (R-t/2)*sin.(theta)

# mat = OrthotropicMaterial(E, E, nu, G, rho)
# laminate = Lamina(t, 0.0, mat)
# b = BeamSection([laminate], y, z)
# W, P, yc, zc = beamstiffness([b])
# s, S, ysc, zsc = shearflow([b], P, yc, zc, npts=2)

# @test isapprox(P[1, 1], 1.835e10, rtol=.001)
# @test isapprox(P[2, 2], 4.587e8, rtol=.14)
# @test isapprox(P[3, 3], 4.587e8, rtol=.14)
# @test isapprox(P[4, 4], 3.519e8, rtol=.05)

# @test isapprox(S[1, 1], 4.682e9, rtol=.14)
# @test isapprox(S[2, 2], 4.682e9, rtol=.14)

# chord = 2*(R-t/2)
# t_lamU = [t]
# t_lamL = t_lamU

# input = PreComp.Input(chord, twist, twist_prime, le_loc, xnode, ynode, E1, E2, G12, nu12, rhov,
#     xsec_nodeU, n_laminaU, n_pliesU, t_lamU, tht_lamU, mat_lamU,
#     xsec_nodeL, n_laminaL, n_pliesL, t_lamL, tht_lamL, mat_lamL,
#     loc_web, n_laminaW, n_pliesW, t_lamW, tht_lamW, mat_lamW)

# output = PreComp.properties(input)


# println("---- errors -------")
# EAexact = 1.835e10
# EAvabs = 1.834e10
# println("EA (mine): ", round((P[1, 1]/EAexact-1)*100, digits=3), "%, EA (PreComp): ", round((output.ea/EAexact-1)*100, digits=3), "%, EA (VABS): ", round((EAvabs/EAexact-1)*100, digits=3), "%")
# EIexact = 4.587e08
# EIvabs = 4.586e08
# println("EI (mine): ", round((P[2, 2]/EIexact-1)*100, digits=2), "%, EI (PreComp): ", round((output.ei_flap/EIexact-1)*100, digits=2), "%, EI (VABS): ", round((EIvabs/EIexact-1)*100, digits=2), "%")
# GJexact = 3.519e08
# GJvabs = 3.515e08
# println("GJ (mine): ", round((P[4, 4]/GJexact-1)*100, digits=2), "%, GJ (PreComp): ", round((output.gj/GJexact-1)*100, digits=2), "%, GJ (VABS): ", round((GJvabs/GJexact-1)*100, digits=2), "%")

# tratios = [1/15, 1/7.5, 1/5, 2/7.5, 1/3]
# EAvec = zeros(5)
# EIvec = zeros(5)
# GJvec = zeros(5)

# for i = 1:5
#     t = tratios[i]*2*R

#     y = (R-t/2)*cos.(theta)
#     z = (R-t/2)*sin.(theta)

#     laminate = Lamina(t, 0.0, mat)
#     b = BeamSection([laminate], y, z)
#     W, P, yc, zc = beamstiffness([b])
#     EAvec[i] = P[1, 1]
#     EIvec[i] = P[2, 2]
#     GJvec[i] = P[4, 4]
# end

# EAexact = [5.137e9, 9.540e9, 1.321e10, 1.615e10, 1.835e10]
# EIexact = [2.024e8, 3.301e8, 4.042e8, 4.424e8, 4.587e8]
# GJexact = [1.553e8, 2.532e8, 3.101e8, 3.394e8, 3.519e8]
# EAerror = @. abs(EAvec/EAexact - 1)*100
# EIerror = @. abs(EIvec/EIexact - 1)*100
# GJerror = @. abs(GJvec/GJexact - 1)*100

# figure()
# plot(tratios, EAerror)
# plot(tratios, EIerror)
# plot(tratios, GJerror)
# plot([0, 0.35], [0, 0], "k--")
# xlim([0, 0.35])
# ylim([-10, 60])
# # savefig("/Users/aning/Downloads/stiffnesserror.pdf")

# # -------- highly heterogenous section ---------

# E1 = 206.843e9
# nu1 = 0.49
# rho1 = 1068.69
# G1 = E1/(2*(1 + nu1))
# mat1 = OrthotropicMaterial(E1, E1, nu1, G1, rho1)

# E2 = E1*1e-12
# nu2 = nu1*1e-12
# rho2 = rho1*1e-12
# G2 = E2/(2*(1 + nu2))
# mat2 = OrthotropicMaterial(E2, E2, nu2, G2, rho2)

# t = 1.524e-3
# laminate1 = Lamina(t, 0.0, mat1)
# laminate2 = Lamina(t, 0.0, mat2)

# b1 = BeamSection([laminate1], [12.7, 0, 0, 25.4]*1e-3, [12.7, 12.7, -12.7, -12.7]*1e-3)

# W, P, yc, zc = beamstiffness([b1], closedsection=false)
# s, S, ysc, zsc = shearflow([b1], P, yc, zc, closedsection=false)

# K = fullstiffnessmatrix(P, S)
# theta = 0.0
# r = [0.0, ysc, zsc]
# Kp = rotatestiffnessmatrix(K, r, theta)

# EIyexact = 2.463e3
# EIzexact = 3.542e3
# GJexact = 4.918
# EAexact = 1.9056e7
# S34exact = -6.263e2
# S13exact = 1.053e5
# S14exact = -2.191e5
# println("EA error: ", round(abs(P[1, 1]/EAexact-1)*100, digits=2), "%")
# println("EIy error: ", round(abs(P[2, 2]/EIyexact-1)*100, digits=2), "%")
# println("EIy error: ", round(abs((P[2, 2] + zsc^2*P[1, 1])/EIyexact-1)*100, digits=2), "%")
# println("EIz error: ", round(abs(P[3, 3]/EIzexact-1)*100, digits=2), "%")
# println("EIz error: ", round(abs((P[3, 3] + ysc^2*P[1, 1])/EIzexact-1)*100, digits=2), "%")
# println("GJ error: ", round(abs(P[4, 4]/GJexact-1)*100, digits=2), "%")

# println("S34 error: ", round(abs(P[2, 3]/S34exact-1)*100, digits=2), "%")
# println("S34 error: ", round(abs((P[2, 3] + ysc*zsc*P[1, 1])/S34exact-1)*100, digits=2), "%")
# println("S13 error: ", round(abs(P[1, 2]/S13exact-1)*100, digits=2), "%")
# println("S13 error: ", round(abs((P[1, 2] - zsc*P[1, 1])/S13exact-1)*100, digits=2), "%")
# println("S14 error: ", round(abs(P[1, 3]/S14exact-1)*100, digits=2), "%")
# println("S14 error: ", round(abs((P[1, 3] + ysc*P[1, 1])/S14exact-1)*100, digits=2), "%")

# # ------ composite pipe -----------

# # E1 = 141.963e9
# # E2 = 9.79056e9
# # nu12 = 0.42
# # G12 = 59.9844e9
# E1 = 20.59e6
# E2 = 1.42e6
# nu12 = 0.42
# G12 = 0.87e6
# rho = 1.0
# mat1 = OrthotropicMaterial(E1, E2, nu12, G12, rho)

# # t = [2.54e-3, 2.54e-3]
# t = [0.1, 0.1]
# theta = [0.0, 90.0]*pi/180
# laminate1 = Lamina.(t, theta, Ref(mat1))

# theta = [-45, 45]*pi/180
# laminate2 = Lamina.(t, theta, Ref(mat1))

# # s = 50.8e-3
# # r = 10.16e-3
# s = 2.0
# r = 0.4
# b1 = BeamSection(laminate1, [-s/2, s/2], [-r, -r])
# # figure(); plot([-s/2, s/2], [-r, -r])
# theta = range(-pi/2, pi/2, length=20)
# y = r*cos.(theta)
# z = r*sin.(theta)
# y .+= s/2
# # plot(y, z)
# b2 = BeamSection(laminate2, y, z)
# b3 = BeamSection(laminate1, [s/2, -s/2], [r, r])
# # plot([s/2, -s/2], [r, r])
# theta = range(pi/2, 3*pi/2, length=20)
# y = r*cos.(theta)
# z = r*sin.(theta)
# y .-= s/2
# b4 = BeamSection(laminate2, y, z)
# # plot(y, z)

# beams = [b1, b2, b3, b4]

# W, P, yc, zc = beamstiffness(beams)
# s, S, ysc, zsc = shearflow(beams, P, yc, zc)

# K = fullstiffnessmatrix(P, S)
# theta = 0.0
# r = [0.0, ysc, zsc]
# Kp = rotatestiffnessmatrix(K, r, theta)

# println("here ------")
# abs(Kp[1, 1]/4.621e7 - 1)*100
# abs(Kp[2, 2]/3.489e6 - 1)*100
# abs(Kp[3, 3]/1.463e6 - 1)*100
# abs(Kp[4, 4]/1.971e3 - 1)*100
# abs(Kp[5, 5]/5.402e3 - 1)*100
# abs(Kp[6, 6]/1.547e4 - 1)*100
# abs(Kp[1, 4]/1.111e4 - 1)*100
# abs(Kp[2, 5]/-9.251e2 - 1)*100
# abs(Kp[3, 6]/-5.859e3 - 1)*100

# println("K11 = ", round((Kp[1, 1]/1.03892e7 - 1)*100, digits=2), "%")
# println("K22 = ", round((Kp[2, 2]/7.85310e5 - 1)*100, digits=2), "%")
# println("K33 = ", round((Kp[3, 3]/3.29279e5 - 1)*100, digits=2), "%")
# println("K14 = ", round((Kp[1, 4]/9.84575e4 - 1)*100, digits=2), "%")
# println("K25 = ", round((Kp[2, 5]/-8.21805e3 - 1)*100, digits=2), "%")
# println("K36 = ", round((Kp[3, 6]/-5.20981e4 - 1)*100, digits=2), "%")
# println("K44 = ", round((Kp[4, 4]/6.87275e5 - 1)*100, digits=2), "%")
# println("K55 = ", round((Kp[5, 5]/1.88238e6 - 1)*100, digits=2), "%")
# println("K66 = ", round((Kp[6, 6]/5.38987e6 - 1)*100, digits=2), "%")


# # -------- simple airfoil -----

# E1 = 206.843e9
# nu1 = 0.49
# rho1 = 1068.69
# G1 = E1/(2*(1 + nu1))
# mat1 = OrthotropicMaterial(E1, E1, nu1, G1, rho1)

# t = 1.524e-3
# theta = 0.0
# laminate = Lamina(t, theta, mat1)
# laminate1 = [laminate]

# r = 12.7e-3
# r -= t/2
# l = 38.1e-3 - t/2
# alpha = asin(r/l)
# d = l*cos(alpha)

# theta = range(pi/2 - alpha, 3*pi/2 + alpha, length=60)
# y = r*cos.(theta)
# z = r*sin.(theta)
# b1 = BeamSection(laminate1, y, z)
# b2 = BeamSection(laminate1, [r*sin(alpha), l], [-r*cos(alpha), 0.0])
# b3 = BeamSection(laminate1, [l, r*sin(alpha)], [0.0, r*cos(alpha)])
# profile = [b1, b2, b3]
# # figure(); plot(y, z); plot([r*sin(alpha), l], [-r*cos(alpha), 0.0]); plot([l, r*sin(alpha)], [0.0, r*cos(alpha)])

# W, P, yc, zc = beamstiffness(profile)
# s, S, ysc, zsc = shearflow(profile, P, yc, zc)
# K = fullstiffnessmatrix(P, S)

# theta = 0.0
# r = [0.0, ysc, zsc]
# Kp = rotatestiffnessmatrix(K, r, theta)

# abs(Kp[1, 1]/3.566e7 - 1)*100
# abs(Kp[2, 2]/8.252e6 - 1)*100
# abs(Kp[3, 3]/2.444e6 - 1)*100
# abs(Kp[4, 4]/1.762e3 - 1)*100
# abs(Kp[5, 5]/2.101e3 - 1)*100
# abs(Kp[6, 6]/1.113e4 - 1)*100
# abs(Kp[3, 4]/2.382e3 - 1)*100
# abs(Kp[1, 6]/-3.394e5 - 1)*100

# # -------- wind turbine blade -------

# chord = 1.9
# web = [0.15, 0.5]
# nodes = [0.0, 0.0016, 0.0041, 0.1147, 0.5366, 1.0]
# mat = Vector{OrthotropicMaterial}(undef, 5)
# mat[1] = OrthotropicMaterial(3.7e10, 9.0e9, 0.28, 4.0e9, 1.0)
# mat[2] = OrthotropicMaterial(1.03e10, 1.03e10, 0.30, 8.0e9, 1.0)
# mat[3] = OrthotropicMaterial(10.0, 10.0, 0.30, 1.0, 1.0)
# mat[4] = OrthotropicMaterial(1.03e10, 1.03e10, 0.30, 8.0e9, 1.0)
# mat[5] = OrthotropicMaterial(1.0e7, 1.0e7, 0.30, 2.0e5, 1.0)

# t = [0.000381, 0.00051, 18*0.00053]
# theta = [0, 0, 20]*pi/180
# segment1 = Lamina.(t, theta, [mat[3], mat[4], mat[2]])
# t = [0.000381, 0.00051, 33*0.00053]
# theta = [0, 0, 20]*pi/180
# segment3 = Lamina.(t, theta, [mat[3], mat[4], mat[2]])
# t = [0.000381, 0.00051, 17*0.00053, 38*0.00053, 1*0.003125, 37*0.00053, 16*0.00053]
# theta = [0, 0, 20, 30, 0, 30, 20]*pi/180
# idx = [3, 4, 2, 1, 5, 1, 2]
# segment4 = Lamina.(t, theta, mat[idx])
# t = [0.000381, 0.00051, 17*0.00053, 0.003125, 16*0.00053]
# theta = [0, 0, 20, 0, 0]*pi/180
# idx = [3, 4, 2, 5, 2]
# segment5 = Lamina.(t, theta, mat[idx])
# t = [38*0.00053, 0.003125, 38*0.00053]
# theta = [0, 0, 0]*pi/180
# idx = [1, 5, 1]
# webs = Lamina.(t, theta, mat[idx])

# yn = [1.00000000, 0.99619582, 0.98515158, 0.96764209, 0.94421447, 0.91510964, 0.88074158, 0.84177999, 0.79894110, 0.75297076, 0.70461763, 0.65461515, 0.60366461, 0.55242353, 0.50149950, 0.45144530, 0.40276150, 0.35589801, 0.31131449, 0.26917194, 0.22927064, 0.19167283, 0.15672257, 0.12469599, 0.09585870, 0.07046974, 0.04874337, 0.03081405, 0.01681379, 0.00687971, 0.00143518, 0.00053606, 0.00006572, 0.00001249, 0.00023032, 0.00079945, 0.00170287, 0.00354717, 0.00592084, 0.01810144, 0.03471169, 0.05589286, 0.08132751, 0.11073805, 0.14391397, 0.18067874, 0.22089879, 0.26433734, 0.31062190, 0.35933893, 0.40999990, 0.46204424, 0.51483073, 0.56767889, 0.61998250, 0.67114514, 0.72054815, 0.76758733, 0.81168064, 0.85227225, 0.88883823, 0.92088961, 0.94797259, 0.96977487, 0.98607009, 0.99640466, 1.00000000]
# zn = [0.00000000, 0.00017047, 0.00100213, 0.00285474, 0.00556001, 0.00906779, 0.01357364, 0.01916802, 0.02580144, 0.03334313, 0.04158593, 0.05026338, 0.05906756, 0.06766426, 0.07571157, 0.08287416, 0.08882939, 0.09329359, 0.09592864, 0.09626763, 0.09424396, 0.09023579, 0.08451656, 0.07727756, 0.06875796, 0.05918984, 0.04880096, 0.03786904, 0.02676332, 0.01592385, 0.00647946, 0.00370956, 0.00112514, -0.00046881, -0.00191488, -0.00329201, -0.00470585, -0.00688469, -0.00912202, -0.01720842, -0.02488211, -0.03226730, -0.03908459, -0.04503763, -0.04986836, -0.05338180, -0.05551392, -0.05636585, -0.05605816, -0.05472399, -0.05254383, -0.04969990, -0.04637175, -0.04264894, -0.03859653, -0.03433153, -0.02996944, -0.02560890, -0.02134397, -0.01726049, -0.01343567, -0.00993849, -0.00679919, -0.00402321, -0.00180118, -0.00044469, 0.00000000]


# i0 = argmin(yn)
# yu = yn[1:i0]
# zu = zn[1:i0]

# ynew = yu
# znew = zu
# for i = 2:5
#     ix = findlast(ynew .> nodes[i])

#     ynew = [ynew[1:ix]; nodes[i]; ynew[ix+1:end]]
#     znew = [znew[1:ix]; linear(reverse(yu), reverse(zu), nodes[i]); znew[ix+1:end]]
# end

# idx = zeros(Int64, 6)
# idx[1] = length(ynew)
# for i = 2:5
#     idx[i] = findfirst(ynew .== nodes[i])
# end
# idx[6] = 1

# ynew *= chord
# znew *= chord

# # # add half thickness
# # segments = [segment1, segment1, segment3, segment4, segment5]
# # nx = length(ynew)
# # zt = zeros(nx)
# # for i = 1:nx
# #     segidx = searchsortedfirst(nodes, ynew[i]) - 1
# #     if  segidx == 6
# #         segidx = 5
# #     end
# #     seg = segments[segidx]
# #     zt[i] = znew[i]
# #     T = 0.0
# #     for j = 1:length(seg)
# #         T += seg[j].t
# #     end
# #     zt[i] -= T/2.0
# # end



# b1 = BeamSection(segment5, ynew[idx[6]:idx[5]], znew[idx[6]:idx[5]])
# b2 = BeamSection(segment4, ynew[idx[5]:idx[4]], znew[idx[5]:idx[4]])
# b3 = BeamSection(segment3, ynew[idx[4]:idx[3]], znew[idx[4]:idx[3]])
# b4 = BeamSection(segment1, ynew[idx[3]:idx[2]], znew[idx[3]:idx[2]])
# b5 = BeamSection(segment1, ynew[idx[2]:idx[1]], znew[idx[2]:idx[1]])

# # figure()
# # plot(ynew[idx[6]:idx[5]], znew[idx[6]:idx[5]])
# # plot(ynew[idx[5]:idx[4]], znew[idx[5]:idx[4]])
# # plot(ynew[idx[4]:idx[3]], znew[idx[4]:idx[3]])
# # plot(ynew[idx[3]:idx[2]], znew[idx[3]:idx[2]])
# # plot(ynew[idx[2]:idx[1]], znew[idx[2]:idx[1]])

# yl = yn[i0:end]
# zl = zn[i0:end]

# ynew = yl
# znew = zl
# for i = 2:5
#     ix = findlast(ynew .< nodes[i])

#     ynew = [ynew[1:ix]; nodes[i]; ynew[ix+1:end]]
#     znew = [znew[1:ix]; linear(yl, zl, nodes[i]); znew[ix+1:end]]
# end

# idx[1] = 1
# for i = 2:5
#     idx[i] = findfirst(ynew .== nodes[i])
# end
# idx[6] = length(ynew)

# # figure()
# # plot(ynew[idx[1]:idx[2]], znew[idx[1]:idx[2]])
# # plot(ynew[idx[2]:idx[3]], znew[idx[2]:idx[3]])
# # plot(ynew[idx[3]:idx[4]], znew[idx[3]:idx[4]])
# # plot(ynew[idx[4]:idx[5]], znew[idx[4]:idx[5]])
# # plot(ynew[idx[5]:idx[6]], znew[idx[5]:idx[6]])

# ynew *= chord
# znew *= chord


# # # add half thickness
# # zt = zeros(nx)
# # for i = 1:nx
# #     segidx = searchsortedfirst(nodes, ynew[i]) - 1
# #     if segidx == 6
# #         segidx = 5
# #     end
# #     seg = segments[segidx]
# #     zt[i] = znew[i]
# #     T = 0.0
# #     for j = 1:length(seg)
# #         T += seg[j].t
# #     end
# #     zt[i] += T/2.0
# # end

# b6 = BeamSection(segment1, ynew[idx[6]:idx[5]], znew[idx[6]:idx[5]])
# b7 = BeamSection(segment1, ynew[idx[5]:idx[4]], znew[idx[5]:idx[4]])
# b8 = BeamSection(segment3, ynew[idx[4]:idx[3]], znew[idx[4]:idx[3]])
# b9 = BeamSection(segment4, ynew[idx[3]:idx[2]], znew[idx[3]:idx[2]])
# b10 = BeamSection(segment5, ynew[idx[2]:idx[1]], znew[idx[2]:idx[1]])

# webseg = Vector{BeamSection}(undef, 2)
# for i = 1:2
#     wzu = linear(reverse(yu), reverse(zu), web[i])
#     wzl = linear(yl, zl, web[i])
#     webseg[i] = BeamSection(webs, chord*[web[i], web[i]], chord*[wzl, wzu])
#     webseg[i] = BeamSection(webs, chord*[web[i], web[i]], chord*[wzl, wzu])
# end

# profile = [b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, webseg[1], webseg[2]]

# W, P, yc, zc = beamstiffness(profile)

# K = fullstiffnessmatrix(P, zeros(2, 2))

# theta = 0.0
# r = [0.0, 0.031, 0.040]
# Kp = rotatestiffnessmatrix(K, r, theta)

# abs(Kp[1, 1]/2.389e9 - 1)*100
# # abs(Kp[2, 2]/8.252e6 - 1)*100
# # abs(Kp[3, 3]/2.444e6 - 1)*100
# abs(Kp[4, 4]/2.167e7 - 1)*100
# abs(Kp[5, 5]/1.97e7 - 1)*100
# abs(Kp[6, 6]/4.406e8 - 1)*100
# abs(Kp[1, 4]/-3.382e7 - 1)*100
# abs(Kp[1, 5]/-2.627e7 - 1)*100
# abs(Kp[1, 6]/-4.736e8 - 1)*100
# abs(Kp[4, 5]/-6.279e4 - 1)*100
# abs(Kp[4, 6]/1.430e6 - 1)*100
# abs(Kp[5, 6]/1.209e7 - 1)*100