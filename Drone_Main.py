import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import pi, g
from Drone_RigidFW2 import RigidFW2
from Drone_BEMT_Hover import BEMT_Hover

# Initial parameters
N = 2
R = 0.06
A = pi * R**2
r_segn = np.linspace(0.1, 1, 30)  # Array 1D
theta0 = np.radians(20)  # [rad]
theta_tw = np.radians(-9)
theta = theta0 + theta_tw * (r_segn - r_segn[0])

c = 0.015  # [m]
solidity = N * c / (pi * R)
N_rotors = 4
W = 0.246 * g / N_rotors   # [N]
mu = 0
Clalpha = 2 * pi           # [rad^-1]
Cd0 = 0.02
height = 2000              # [m]
rho_inf = 1.0065           # [Kg/m^3]

Tc, _, _, _, _ = BEMT_Hover(N, r_segn, solidity, theta, Clalpha, Cd0, mu)
print("The value of Tc is:", Tc)
Omega_h = np.sqrt(W / (rho_inf * R**2 * A * Tc))  # [rad/s]
print("The value of Omega_h is:", Omega_h)

Pc_i = Tc**(3/2) / np.sqrt(2)
Pc0 = solidity * Cd0 / 8
RPM_h = Omega_h / (2 * pi) * 60  # [RPM]
lambda_i = np.sqrt(Tc / 2)
w_h = lambda_i * Omega_h * R

P_ih = Pc_i * rho_inf * Omega_h**3 * R**3 * A
Pc_0h = Pc0
Vinf_vec = np.linspace(0, 16, 30)  #  [m/s]
f_A = 0.07
f = f_A * (4 * A)

linear_inflow = 0
results = RigidFW2(Vinf_vec, N, r_segn, c, theta, W, w_h, height, R, A,
                   P_ih, f, Pc_0h, Clalpha, Cd0, Omega_h, N_rotors, linear_inflow)

(P_i, P_0, P_fus, _, P_new, P_new_i, P_new_0, P_new_fus, alpha_vec,
 Omega_new, alpha_e_vec, Mach_vec, dTc, U_P, dH_dpsi_0, dQ_dpsi_0,
 Tc_vec, Y_vec, Qc_vec) = results

Vinf_vec = Vinf_vec * 1.94384  # [kts]

# Graphics BEMT vs IT
#  P_tot 
plt.figure(1)
plt.plot(Vinf_vec, 4 * P_new, '-^', markersize=6, color='k', label='BEMT')
plt.plot(Vinf_vec, 4 * (P_i + P_0 + P_fus / N_rotors), '-s', markersize=6, color='k', label='IT')
plt.legend()
plt.xlabel(r'$V_\infty \, [kts]$', fontsize=12)
plt.ylabel(r'$P_{TOT} \, [W]$', fontsize=12)
plt.grid(which='minor')

# P_fus
plt.figure(2)
plt.plot(Vinf_vec, P_new_fus, '-^', markersize=6, color='k', label='BEMT')
plt.plot(Vinf_vec, P_fus, '-s', markersize=6, color='k', label='IT')
plt.legend()
plt.xlabel(r'$V_\infty \, [kts]$', fontsize=12)
plt.ylabel(r'$P_{fus} \, [W]$', fontsize=12)
plt.grid(which='minor')

# P_profile
plt.figure(3)
plt.plot(Vinf_vec, P_new_0, '-^', markersize=6, color='k', label='BEMT')
plt.plot(Vinf_vec, P_0, '-s', markersize=6, color='k', label='IT')
plt.legend()
plt.xlabel(r'$V_\infty \, [kts]$', fontsize=12)
plt.ylabel(r'$P_{0} \, [W]$', fontsize=12)
plt.grid(which='minor')

# P_induced
plt.figure(4)
plt.plot(Vinf_vec, P_new_i, '-^', markersize=6, color='k', label='BEMT')
plt.plot(Vinf_vec, P_i, '-s', markersize=6, color='k', label='IT')
plt.legend()
plt.xlabel(r'$V_\infty \, [kts]$', fontsize=12)
plt.ylabel(r'$P_{i} \, [W]$', fontsize=12)
plt.grid(which='minor')

# Tc
plt.figure(5)
plt.plot(Vinf_vec, Tc_vec, 'k')
plt.xlabel(r'$V_\infty \, [kts]$', fontsize=12)
plt.ylabel(r'$T_{c}$', fontsize=12)
plt.grid(which='minor')

# Qc
plt.figure(6)
plt.plot(Vinf_vec, Qc_vec, 'k')
plt.xlabel(r'$V_\infty \, [kts]$', fontsize=12)
plt.ylabel(r'$Q_{c}$', fontsize=12)
plt.grid(which='minor')

plt.show()
