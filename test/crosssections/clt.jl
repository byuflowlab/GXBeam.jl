using GXBeam.GXBeamCS
using Test

#------- Foundations of Classical Laminate Theory, Andreas Öchsner -----
# -------- 4.2. Problem 1 ------------

E1 = 129000.0e6  # note must be typo on this book
E2 = 11000.0e6
G12 = 6600.0e6
nu12 = 0.28
rho = 1.0
S1t = 1950.0e6
S1c = 1480e6
S2t = 48e6
S2c = 200e6
S12 = 79e6
mat = MaterialPlane(E1, E2, G12, nu12, rho, S1t, S1c, S2t, S2c, S12)

t = 8e-3
theta = [45 -45 0 90 90 0 -45 45]*pi/180
laminate = Layer.(Ref(mat), t/8, theta)

Q1, _, _ = GXBeamCS.Qbar(laminate[1])
@test isapprox(Q1[1, 1]/1e6, 43385.924, atol=1e-3)
@test isapprox(Q1[1, 2]/1e6, 30185.924, atol=1e-3)
@test isapprox(Q1[1, 3]/1e6, 29698.543, atol=1e-3)
@test isapprox(Q1[2, 2]/1e6, 43385.924, atol=1e-3)
@test isapprox(Q1[2, 3]/1e6, 29698.543, atol=1e-3)
@test isapprox(Q1[3, 3]/1e6, 33685.195, atol=1e-3)
@test Q1[1, 2] == Q1[2, 1]
@test Q1[1, 3] == Q1[3, 1]
@test Q1[2, 3] == Q1[3, 2]

Q2, _, _ = GXBeamCS.Qbar(laminate[2])
@test isapprox(Q2[1, 1]/1e6, 43385.924, atol=1e-3)
@test isapprox(Q2[1, 2]/1e6, 30185.924, atol=1e-3)
@test isapprox(Q2[1, 3]/1e6, -29698.543, atol=1e-3)
@test isapprox(Q2[2, 2]/1e6, 43385.924, atol=1e-3)
@test isapprox(Q2[2, 3]/1e6, -29698.543, atol=1e-3)
@test isapprox(Q2[3, 3]/1e6, 33685.195, atol=1e-3)
@test Q2[1, 2] == Q2[2, 1]
@test Q2[1, 3] == Q2[3, 1]
@test Q2[2, 3] == Q2[3, 2]

Q3, _, _ = GXBeamCS.Qbar(laminate[3])
@test isapprox(Q3[1, 1]/1e6, 129868.204, atol=1e-3)
@test isapprox(Q3[1, 2]/1e6, 3100.729, atol=1e-3)
@test isapprox(Q3[1, 3]/1e6, 0.0, atol=1e-3)
@test isapprox(Q3[2, 2]/1e6, 11074.033, atol=1e-3)
@test isapprox(Q3[2, 3]/1e6, 0.0, atol=1e-3)
@test isapprox(Q3[3, 3]/1e6, 6600.0, atol=1e-3)

Q4, _, _ = GXBeamCS.Qbar(laminate[4])
@test isapprox(Q4[1, 1]/1e6, 11074.033, atol=1e-3)
@test isapprox(Q4[1, 2]/1e6, 3100.729, atol=1e-3)
@test isapprox(Q4[1, 3]/1e6, 0.0, atol=1e-3)
@test isapprox(Q4[2, 2]/1e6, 129868.204, atol=1e-3)
@test isapprox(Q4[2, 3]/1e6, 0.0, atol=1e-3)
@test isapprox(Q4[3, 3]/1e6, 6600.0, atol=1e-3)

z, h = GXBeamCS.zspacing(laminate)
A, B, D = GXBeamCS.stiffnessmatrix(laminate, z)

@test isapprox(A[1, 1]/1e3, 455428.170, atol=1e-3)
@test isapprox(A[1, 2]/1e3, 133146.612, atol=1e-3)
@test isapprox(A[1, 3]/1e3, 0.0, atol=1e-3)
@test isapprox(A[2, 2]/1e3, 455428.170, atol=1e-3)
@test isapprox(A[2, 3]/1e3, 0.0, atol=1e-3)
@test isapprox(A[3, 3]/1e3, 161140.779, atol=1e-3)

@test isapprox(B[1, 1], 0.0, atol=1e-10)
@test isapprox(B[1, 2], 0.0, atol=1e-10)
@test isapprox(B[1, 3], 0.0, atol=1e-10)
@test isapprox(B[2, 2], 0.0, atol=1e-10)
@test isapprox(B[2, 3], 0.0, atol=1e-10)
@test isapprox(B[3, 3], 0.0, atol=1e-10)

@test isapprox(D[1, 1]*1e3, 2233175.466, atol=1e-3)
@test isapprox(D[1, 2]*1e3, 1143478.381, atol=1e-3)
@test isapprox(D[1, 3]*1e3, 356382.514, atol=1e-3)
@test isapprox(D[2, 2]*1e3, 1757998.781, atol=1e-3)
@test isapprox(D[2, 3]*1e3, 356382.514, atol=1e-3)
@test isapprox(D[3, 3]*1e3, 1292780.601, atol=1e-3)

alpha, beta, delta = GXBeamCS.compliancematrix(A, B, D)

@test isapprox(alpha[1, 1]/1e-10, 24.009482, atol=1e-6)
@test isapprox(alpha[1, 2]/1e-10, -7.019287, atol=1e-6)
@test isapprox(alpha[1, 3]/1e-10, 0.0, atol=1e-6)
@test isapprox(alpha[2, 2]/1e-10, 24.009482, atol=1e-6)
@test isapprox(alpha[2, 3]/1e-10, 0.0, atol=1e-6)
@test isapprox(alpha[3, 3]/1e-10, 62.057538, atol=1e-6)

@test isapprox(beta[1, 1], 0.0, atol=1e-10)
@test isapprox(beta[1, 2], 0.0, atol=1e-10)
@test isapprox(beta[1, 3], 0.0, atol=1e-10)
@test isapprox(beta[2, 2], 0.0, atol=1e-10)
@test isapprox(beta[2, 3], 0.0, atol=1e-10)
@test isapprox(beta[3, 3], 0.0, atol=1e-10)

@test isapprox(delta[1, 1]/1e-4, 6.771890, atol=1e-6)
@test isapprox(delta[1, 2]/1e-4, -4.264613, atol=1e-6)
@test isapprox(delta[1, 3]/1e-4, -0.691184, atol=1e-6)
@test isapprox(delta[2, 2]/1e-4, 8.710637, atol=1e-6)
@test isapprox(delta[2, 3]/1e-4, -1.225641, atol=1e-6)
@test isapprox(delta[3, 3]/1e-4, 8.263678, atol=1e-6)

Nx = 1000e3
forces = [Nx; 0.0; 0.0; 0.0; 0.0; 0.0]

epsilonbar, kappa, zvec, epsilonp = GXBeamCS.strains(alpha, beta, delta, z, forces)

@test isapprox(epsilonbar[1]*1e3, 2.400948, atol=1e-6)
@test isapprox(epsilonbar[2]*1e3, -0.701929, atol=1e-6)
@test isapprox(epsilonbar[3]*1e3, 0.0, atol=1e-10)
@test isapprox(kappa[1], 0.0, atol=1e-10)
@test isapprox(kappa[2], 0.0, atol=1e-10)
@test isapprox(kappa[3], 0.0, atol=1e-10)


sigmap, sigma, epsilon = GXBeamCS.stresses(laminate, epsilonp)

# using PyPlot
# close("all"); pygui(true)

# figure()
# plot(sigmap[1, :]/1e6, zvec)
# figure()
# plot(epsilonp[1, :]/1e6, zvec)

@test isapprox(sigma[1, 1]/1e6, 112.958, atol=1e-3)
@test isapprox(sigma[1, 2]/1e6, 112.958, atol=1e-3)
@test isapprox(sigma[1, 3]/1e6, 112.958, atol=1e-3)
@test isapprox(sigma[1, 4]/1e6, 112.958, atol=1e-3)
@test isapprox(sigma[1, 5]/1e6, 309.630, atol=1e-3)
@test isapprox(sigma[1, 6]/1e6, 309.630, atol=1e-3)
@test isapprox(sigma[1, 7]/1e6, -83.714, atol=1e-3)
@test isapprox(sigma[1, 8]/1e6, -83.714, atol=1e-3)
@test isapprox(sigma[1, 9]/1e6, -83.714, atol=1e-3)
@test isapprox(sigma[1, 10]/1e6, -83.714, atol=1e-3)
@test isapprox(sigma[1, 11]/1e6, 309.630, atol=1e-3)
@test isapprox(sigma[1, 12]/1e6, 309.630, atol=1e-3)
@test isapprox(sigma[1, 13]/1e6, 112.958, atol=1e-3)
@test isapprox(sigma[1, 14]/1e6, 112.958, atol=1e-3)
@test isapprox(sigma[1, 15]/1e6, 112.958, atol=1e-3)
@test isapprox(sigma[1, 16]/1e6, 112.958, atol=1e-3)

@test isapprox(sigma[2, 1]/1e6, 12.042, atol=1e-3)
@test isapprox(sigma[2, 2]/1e6, 12.042, atol=1e-3)
@test isapprox(sigma[2, 3]/1e6, 12.042, atol=1e-3)
@test isapprox(sigma[2, 4]/1e6, 12.042, atol=1e-3)
@test isapprox(sigma[2, 5]/1e6, -0.328, atol=1e-3)
@test isapprox(sigma[2, 6]/1e6, -0.328, atol=1e-3)
@test isapprox(sigma[2, 7]/1e6, 24.412, atol=1e-3)
@test isapprox(sigma[2, 8]/1e6, 24.412, atol=1e-3)
@test isapprox(sigma[2, 9]/1e6, 24.412, atol=1e-3)
@test isapprox(sigma[2, 10]/1e6, 24.412, atol=1e-3)
@test isapprox(sigma[2, 11]/1e6, -0.328, atol=1e-3)
@test isapprox(sigma[2, 12]/1e6, -0.328, atol=1e-3)
@test isapprox(sigma[2, 13]/1e6, 12.042, atol=1e-3)
@test isapprox(sigma[2, 14]/1e6, 12.042, atol=1e-3)
@test isapprox(sigma[2, 15]/1e6, 12.042, atol=1e-3)
@test isapprox(sigma[2, 16]/1e6, 12.042, atol=1e-3)

@test isapprox(sigma[3, 1]/1e6, -20.479, atol=1e-3)
@test isapprox(sigma[3, 2]/1e6, -20.479, atol=1e-3)
@test isapprox(sigma[3, 3]/1e6, 20.479, atol=1e-3)
@test isapprox(sigma[3, 4]/1e6, 20.479, atol=1e-3)
@test isapprox(sigma[3, 5]/1e6, 0.0, atol=1e-3)
@test isapprox(sigma[3, 6]/1e6, 0.0, atol=1e-3)
@test isapprox(sigma[3, 7]/1e6, 0.0, atol=1e-3)
@test isapprox(sigma[3, 8]/1e6, 0.0, atol=1e-3)
@test isapprox(sigma[3, 9]/1e6, 0.0, atol=1e-3)
@test isapprox(sigma[3, 10]/1e6, 0.0, atol=1e-3)
@test isapprox(sigma[3, 11]/1e6, 0.0, atol=1e-3)
@test isapprox(sigma[3, 12]/1e6, 0.0, atol=1e-3)
@test isapprox(sigma[3, 13]/1e6, 20.479, atol=1e-3)
@test isapprox(sigma[3, 14]/1e6, 20.479, atol=1e-3)
@test isapprox(sigma[3, 15]/1e6, -20.479, atol=1e-3)
@test isapprox(sigma[3, 16]/1e6, -20.479, atol=1e-3)

@test isapprox(epsilon[1, 1]*1e3, 0.849510, atol=1e-6)
@test isapprox(epsilon[1, 2]*1e3, 0.849510, atol=1e-6)
@test isapprox(epsilon[1, 3]*1e3, 0.849510, atol=1e-6)
@test isapprox(epsilon[1, 4]*1e3, 0.849510, atol=1e-6)
@test isapprox(epsilon[1, 5]*1e3, 2.400948, atol=1e-6)
@test isapprox(epsilon[1, 6]*1e3, 2.400948, atol=1e-6)
@test isapprox(epsilon[1, 7]*1e3, -0.701929, atol=1e-6)
@test isapprox(epsilon[1, 8]*1e3, -0.701929, atol=1e-6)
@test isapprox(epsilon[1, 9]*1e3, -0.701929, atol=1e-6)
@test isapprox(epsilon[1, 10]*1e3, -0.701929, atol=1e-6)
@test isapprox(epsilon[1, 11]*1e3, 2.400948, atol=1e-6)
@test isapprox(epsilon[1, 12]*1e3, 2.400948, atol=1e-6)
@test isapprox(epsilon[1, 13]*1e3, 0.849510, atol=1e-6)
@test isapprox(epsilon[1, 14]*1e3, 0.849510, atol=1e-6)
@test isapprox(epsilon[1, 15]*1e3, 0.849510, atol=1e-6)
@test isapprox(epsilon[1, 16]*1e3, 0.849510, atol=1e-6)

@test isapprox(epsilon[2, 1]*1e3, 0.849510, atol=1e-6)
@test isapprox(epsilon[2, 2]*1e3, 0.849510, atol=1e-6)
@test isapprox(epsilon[2, 3]*1e3, 0.849510, atol=1e-6)
@test isapprox(epsilon[2, 4]*1e3, 0.849510, atol=1e-6)
@test isapprox(epsilon[2, 5]*1e3, -0.701929, atol=1e-6)
@test isapprox(epsilon[2, 6]*1e3, -0.701929, atol=1e-6)
@test isapprox(epsilon[2, 7]*1e3, 2.400948, atol=1e-6)
@test isapprox(epsilon[2, 8]*1e3, 2.400948, atol=1e-6)
@test isapprox(epsilon[2, 9]*1e3, 2.400948, atol=1e-6)
@test isapprox(epsilon[2, 10]*1e3, 2.400948, atol=1e-6)
@test isapprox(epsilon[2, 11]*1e3, -0.701929, atol=1e-6)
@test isapprox(epsilon[2, 12]*1e3, -0.701929, atol=1e-6)
@test isapprox(epsilon[2, 13]*1e3, 0.849510, atol=1e-6)
@test isapprox(epsilon[2, 14]*1e3, 0.849510, atol=1e-6)
@test isapprox(epsilon[2, 15]*1e3, 0.849510, atol=1e-6)
@test isapprox(epsilon[2, 16]*1e3, 0.849510, atol=1e-6)

@test isapprox(epsilon[3, 1]*1e3, -3.102877, atol=1e-6)
@test isapprox(epsilon[3, 2]*1e3, -3.102877, atol=1e-6)
@test isapprox(epsilon[3, 3]*1e3, 3.102877, atol=1e-6)
@test isapprox(epsilon[3, 4]*1e3, 3.102877, atol=1e-6)
@test isapprox(epsilon[3, 5]*1e3, 0.0, atol=1e-6)
@test isapprox(epsilon[3, 6]*1e3, 0.0, atol=1e-6)
@test isapprox(epsilon[3, 7]*1e3, 0.0, atol=1e-6)
@test isapprox(epsilon[3, 8]*1e3, 0.0, atol=1e-6)
@test isapprox(epsilon[3, 9]*1e3, 0.0, atol=1e-6)
@test isapprox(epsilon[3, 10]*1e3, 0.0, atol=1e-6)
@test isapprox(epsilon[3, 11]*1e3, 0.0, atol=1e-6)
@test isapprox(epsilon[3, 12]*1e3, 0.0, atol=1e-6)
@test isapprox(epsilon[3, 13]*1e3, 3.102877, atol=1e-6)
@test isapprox(epsilon[3, 14]*1e3, 3.102877, atol=1e-6)
@test isapprox(epsilon[3, 15]*1e3, -3.102877, atol=1e-6)
@test isapprox(epsilon[3, 16]*1e3, -3.102877, atol=1e-6)

forces = [0.0; 0.0; 0.0; 1000; 0.0; 0.0]

epsilonbar, kappa, zvec, epsilonp = GXBeamCS.strains(alpha, beta, delta, z, forces)

@test isapprox(epsilonbar[1], 0.0, atol=1e-6)
@test isapprox(epsilonbar[2], 0.0, atol=1e-6)
@test isapprox(epsilonbar[3], 0.0, atol=1e-10)
@test isapprox(kappa[1], 0.677189, atol=1e-6)
@test isapprox(kappa[2], -0.426461, atol=1e-6)
@test isapprox(kappa[3], -0.069118, atol=1e-6)

sigmap, sigma, epsilon = GXBeamCS.stresses(laminate, epsilonp)


# figure()
# plot(sigmap[1, :]/1e6, zvec*1e3)
# figure()
# plot(epsilonp[1, :]*1e3, zvec*1e3)

@test isapprox(sigma[1, 1]/1e6, -49.154, atol=1e-3)
@test isapprox(sigma[1, 2]/1e6, -36.866, atol=1e-3)
@test isapprox(sigma[1, 3]/1e6, -63.151, atol=1e-3)
@test isapprox(sigma[1, 4]/1e6, -42.101, atol=1e-3)
@test isapprox(sigma[1, 5]/1e6, -173.246, atol=1e-3)
@test isapprox(sigma[1, 6]/1e6, -86.623, atol=1e-3)
@test isapprox(sigma[1, 7]/1e6, 53.284, atol=1e-3)
@test isapprox(sigma[1, 8]/1e6, 0.0, atol=1e-3)
@test isapprox(sigma[1, 9]/1e6, 0.0, atol=1e-3)
@test isapprox(sigma[1, 10]/1e6, -53.284, atol=1e-3)
@test isapprox(sigma[1, 11]/1e6, 86.623, atol=1e-3)
@test isapprox(sigma[1, 12]/1e6, 173.246, atol=1e-3)
@test isapprox(sigma[1, 13]/1e6, 42.101, atol=1e-3)
@test isapprox(sigma[1, 14]/1e6, 63.151, atol=1e-3)
@test isapprox(sigma[1, 15]/1e6, 36.866, atol=1e-3)
@test isapprox(sigma[1, 16]/1e6, 49.154, atol=1e-3)

@test isapprox(sigma[2, 1]/1e6, -8.210, atol=1e-3)
@test isapprox(sigma[2, 2]/1e6, -6.158, atol=1e-3)
@test isapprox(sigma[2, 3]/1e6, -4.504, atol=1e-3)
@test isapprox(sigma[2, 4]/1e6, -3.003, atol=1e-3)
@test isapprox(sigma[2, 5]/1e6, 5.246, atol=1e-3)
@test isapprox(sigma[2, 6]/1e6, 2.623, atol=1e-3)
@test isapprox(sigma[2, 7]/1e6, -6.177, atol=1e-3)
@test isapprox(sigma[2, 8]/1e6, 0.0, atol=1e-3)
@test isapprox(sigma[2, 9]/1e6, 0.0, atol=1e-3)
@test isapprox(sigma[2, 10]/1e6, 6.177, atol=1e-3)
@test isapprox(sigma[2, 11]/1e6, -2.623, atol=1e-3)
@test isapprox(sigma[2, 12]/1e6, -5.246, atol=1e-3)
@test isapprox(sigma[2, 13]/1e6, 3.003, atol=1e-3)
@test isapprox(sigma[2, 14]/1e6, 4.504, atol=1e-3)
@test isapprox(sigma[2, 15]/1e6, 6.158, atol=1e-3)
@test isapprox(sigma[2, 16]/1e6, 8.210, atol=1e-3)

@test isapprox(sigma[3, 1]/1e6, 29.136, atol=1e-3)
@test isapprox(sigma[3, 2]/1e6, 21.852, atol=1e-3)
@test isapprox(sigma[3, 3]/1e6, -21.852, atol=1e-3)
@test isapprox(sigma[3, 4]/1e6, -14.568, atol=1e-3)
@test isapprox(sigma[3, 5]/1e6, 0.912, atol=1e-3)
@test isapprox(sigma[3, 6]/1e6, 0.456, atol=1e-3)
@test isapprox(sigma[3, 7]/1e6, -0.456, atol=1e-3)
@test isapprox(sigma[3, 8]/1e6, 0.0, atol=1e-3)
@test isapprox(sigma[3, 9]/1e6, 0.0, atol=1e-3)
@test isapprox(sigma[3, 10]/1e6, 0.456, atol=1e-3)
@test isapprox(sigma[3, 11]/1e6, -0.456, atol=1e-3)
@test isapprox(sigma[3, 12]/1e6, -0.912, atol=1e-3)
@test isapprox(sigma[3, 13]/1e6, 14.568, atol=1e-3)
@test isapprox(sigma[3, 14]/1e6, 21.852, atol=1e-3)
@test isapprox(sigma[3, 15]/1e6, -21.852, atol=1e-3)
@test isapprox(sigma[3, 16]/1e6, -29.136, atol=1e-3)


# ----- 4.3. Problem 2 --------
theta = [45 -45 0 90 0 90 45 -45]*pi/180
laminate = Layer.(Ref(mat), t/8, theta)

forces = [1000*1e3; 0.0; 0.0; 0.0; 0.0; 0.0]

z, h = GXBeamCS.zspacing(laminate)
A, B, D = GXBeamCS.stiffnessmatrix(laminate, z)
alpha, beta, delta = GXBeamCS.compliancematrix(A, B, D)
epsilonbar, kappa, zvec, epsilonp = GXBeamCS.strains(alpha, beta, delta, z, forces)
sigmap, sigma, epsilon = GXBeamCS.stresses(laminate, epsilonp)

@test isapprox(A[1, 1]/1e3, 455428.170, atol=1e-3)
@test isapprox(A[1, 2]/1e3, 133146.612, atol=1e-3)
@test isapprox(A[1, 3]/1e3, 0.0, atol=1e-3)
@test isapprox(A[2, 2]/1e3, 455428.170, atol=1e-3)
@test isapprox(A[2, 3]/1e3, 0.0, atol=1e-3)
@test isapprox(A[3, 3]/1e3, 161140.779, atol=1e-3)

@test isapprox(B[1, 1], -118794.171, atol=1e-3)
@test isapprox(B[1, 2], 1.455e-11, atol=1e-3)
@test isapprox(B[1, 3], -59397.086, atol=1e-3)
@test isapprox(B[2, 2], 118794.171, atol=1e-3)
@test isapprox(B[2, 3], -59397.086, atol=1e-3)
@test isapprox(B[3, 3], 0.0, atol=1e-10)

@test isapprox(D[1, 1]*1e3, 1995587.124, atol=1e-3)
@test isapprox(D[1, 2]*1e3, 1143478.381, atol=1e-3)
@test isapprox(D[1, 3]*1e3, 0.0, atol=1e-3)
@test isapprox(D[2, 2]*1e3, 1995587.124, atol=1e-3)
@test isapprox(D[2, 3]*1e3, 0.0, atol=1e-3)
@test isapprox(D[3, 3]*1e3, 1292780.601, atol=1e-3)

@test isapprox(alpha[1, 1]/1e-10, 24.562273, atol=1e-6)
@test isapprox(alpha[1, 2]/1e-10, -6.911749, atol=1e-6)
@test isapprox(alpha[1, 3]/1e-10, 0.445254, atol=1e-6)
@test isapprox(alpha[2, 2]/1e-10, 24.562273, atol=1e-6)
@test isapprox(alpha[2, 3]/1e-10, -0.445254, atol=1e-6)
@test isapprox(alpha[3, 3]/1e-10, 62.948045, atol=1e-6)

@test isapprox(beta[1, 1]/1e-7, 1.834321, atol=1e-6)
@test isapprox(beta[1, 2]/1e-7, -0.626374, atol=1e-6)
@test isapprox(beta[1, 3]/1e-7, 0.810957, atol=1e-6)
@test isapprox(beta[2, 2]/1e-7, -1.834321, atol=1e-6)
@test isapprox(beta[2, 3]/1e-7, 0.810957, atol=1e-6)
@test isapprox(beta[3, 3]/1e-7, 0.0, atol=1e-6)

@test isapprox(delta[1, 1]/1e-4, 7.677865, atol=1e-6)
@test isapprox(delta[1, 2]/1e-4, -4.400777, atol=1e-6)
@test isapprox(delta[1, 3]/1e-4, 0.113057, atol=1e-6)
@test isapprox(delta[2, 2]/1e-4, 7.677865, atol=1e-6)
@test isapprox(delta[2, 3]/1e-4, -0.113057, atol=1e-6)
@test isapprox(delta[3, 3]/1e-4, 7.809784, atol=1e-6)

@test isapprox(epsilonbar[1]*1e3, 2.456227, atol=1e-6)
@test isapprox(epsilonbar[2]*1e3, -0.691175, atol=1e-6)
@test isapprox(epsilonbar[3]*1e3, 0.044525, atol=1e-6)
@test isapprox(kappa[1], 0.183432, atol=1e-6)
@test isapprox(kappa[2], -0.062637, atol=1e-6)
@test isapprox(kappa[3], 0.081096, atol=1e-6)

@test isapprox(sigma[1, 1]/1e6, 67.486, atol=1e-3)
@test isapprox(sigma[1, 2]/1e6, 80.657, atol=1e-3)
@test isapprox(sigma[1, 3]/1e6, 105.854, atol=1e-3)
@test isapprox(sigma[1, 4]/1e6, 108.745, atol=1e-3)
@test isapprox(sigma[1, 5]/1e6, 269.587, atol=1e-3)
@test isapprox(sigma[1, 6]/1e6, 293.215, atol=1e-3)
@test isapprox(sigma[1, 7]/1e6, -74.580, atol=1e-3)
@test isapprox(sigma[1, 8]/1e6, -82.146, atol=1e-3)
@test isapprox(sigma[1, 9]/1e6, 316.843, atol=1e-3)
@test isapprox(sigma[1, 10]/1e6, 340.470, atol=1e-3)
@test isapprox(sigma[1, 11]/1e6, -89.711, atol=1e-3)
@test isapprox(sigma[1, 12]/1e6, -97.277, atol=1e-3)
@test isapprox(sigma[1, 13]/1e6, 146.513, atol=1e-3)
@test isapprox(sigma[1, 14]/1e6, 159.684, atol=1e-3)
@test isapprox(sigma[1, 15]/1e6, 123.199, atol=1e-3)
@test isapprox(sigma[1, 16]/1e6, 126.090, atol=1e-3)

@test isapprox(sigma[2, 1]/1e6, 10.201, atol=1e-3)
@test isapprox(sigma[2, 2]/1e6, 10.734, atol=1e-3)
@test isapprox(sigma[2, 3]/1e6, 9.149, atol=1e-3)
@test isapprox(sigma[2, 4]/1e6, 10.328, atol=1e-3)
@test isapprox(sigma[2, 5]/1e6, 0.212, atol=1e-3)
@test isapprox(sigma[2, 6]/1e6, 0.087, atol=1e-3)
@test isapprox(sigma[2, 7]/1e6, 23.220, atol=1e-3)
@test isapprox(sigma[2, 8]/1e6, 25.057, atol=1e-3)
@test isapprox(sigma[2, 9]/1e6, -0.038, atol=1e-3)
@test isapprox(sigma[2, 10]/1e6, -0.163, atol=1e-3)
@test isapprox(sigma[2, 11]/1e6, 26.894, atol=1e-3)
@test isapprox(sigma[2, 12]/1e6, 28.731, atol=1e-3)
@test isapprox(sigma[2, 13]/1e6, 13.398, atol=1e-3)
@test isapprox(sigma[2, 14]/1e6, 13.931, atol=1e-3)
@test isapprox(sigma[2, 15]/1e6, 16.225, atol=1e-3)
@test isapprox(sigma[2, 16]/1e6, 17.405, atol=1e-3)



# ------- 4.4. Problem 3 --------

theta = [45 0 30 -45]*pi/180
laminate = Layer.(Ref(mat), t/4, theta)

forces = [1000*1e3; 500e3; 0.0; 0.0; 0.0; 0.0]

# S1t = 1950.0e6
# S1c = 1480.0e6
# S2t = 48.0e6
# S2c = 200.0e6
# S12 = 79.0e6
# strength = CompositeStrength(S1t, S1c, S2t, S2c, S12)

sigma, epsilon, wu = GXBeamCS.clt(laminate, forces)

@test isapprox(sigma[1, 1]/1e6, 18.776, atol=1e-3)
# @test isapprox(sigma[1, 2]/1e6, 18.776, atol=1e-3)  # typo in text, repeated value
@test isapprox(sigma[1, 3]/1e6, 238.813, atol=1e-3)
@test isapprox(sigma[1, 4]/1e6, 229.609, atol=1e-3)
@test isapprox(sigma[1, 5]/1e6, 149.134, atol=1e-3)
@test isapprox(sigma[1, 6]/1e6, 229.851, atol=1e-3)
@test isapprox(sigma[1, 7]/1e6, 214.113, atol=1e-3)
@test isapprox(sigma[1, 8]/1e6, 13.857, atol=1e-3)

@test isapprox(epsilon[3, 1]*1e3, 2.663312, atol=1e-6)
@test isapprox(epsilon[3, 2]*1e3, 1.790924, atol=1e-6)
@test isapprox(epsilon[3, 3]*1e3, -4.138210, atol=1e-6)
@test isapprox(epsilon[3, 4]*1e3, -1.996388, atol=1e-6)
@test isapprox(epsilon[3, 5]*1e3, -0.202719, atol=1e-6)
@test isapprox(epsilon[3, 6]*1e3, 0.112682, atol=1e-6)
@test isapprox(epsilon[3, 7]*1e3, -0.046148, atol=1e-6)
@test isapprox(epsilon[3, 8]*1e3, 0.826241, atol=1e-6)


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

R = [0.681; 0.901; 1.041; 1.467; 1.296; 2.325; 2.216; 1.836]
sigmaR = copy(sigma)
sigmaR[1, :] .*= R
sigmaR[2, :] .*= R
sigmaR[3, :] .*= R
failure = GXBeamCS.tsai_wu_plane(sigmaR, laminate)

@test isapprox(failure[1], 1.0, atol=1e-2)
@test isapprox(failure[2], 1.0, atol=1e-2)
@test isapprox(failure[3], 1.0, atol=1e-2)
@test isapprox(failure[4], 1.0, atol=1e-2)
@test isapprox(failure[5], 1.0, atol=1e-2)
@test isapprox(failure[6], 1.0, atol=1e-2)
@test isapprox(failure[7], 1.0, atol=1e-2)
@test isapprox(failure[8], 1.0, atol=1e-2)

forces = [0.0; 0.0; 0.0; 1000; 500; 0.0]

sigma, epsilon, wu = GXBeamCS.clt(laminate, forces)

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

R = [6.727; 10.909; 9.786; 14.021; 19.030; 10.114; 4.463; 2.602]
sigmaR = copy(sigma)
sigmaR[1, :] .*= R
sigmaR[2, :] .*= R
sigmaR[3, :] .*= R
failure = GXBeamCS.tsai_wu_plane(sigmaR, laminate)

@test isapprox(failure[1], 1.0, atol=1e-2)
@test isapprox(failure[2], 1.0, atol=1e-2)
# @test isapprox(failure[3], 1.0, atol=1e-2)  # again there appears an error in the text for this station.  stresses checkout.
@test isapprox(failure[4], 1.0, atol=1e-2)
@test isapprox(failure[5], 1.0, atol=1e-2)
@test isapprox(failure[6], 1.0, atol=1e-2)
@test isapprox(failure[7], 1.0, atol=1e-2)
@test isapprox(failure[8], 1.0, atol=1e-2)
