# DRONE IN FORWARD FLIGHT F2ROTOR
#
# Description: This is the test case function in order to determine power curve contributions, thrust and torque
#              coefficients with Blade Element Momentum Theory (BEMT) and simple impulse theory (IT) for the DJI-Mini4k
#              Drone.
#
# References:  Renato Tognaccini - Lezioni per il corso di Aerodinamica dell'Ala Rotante - Eliche rotori ed aeromotori
#              con un'introduzione all'aerodinamica instazionaria" - a.a. 2023/2024 - vsn 2.04.
#
# Authors:     Pasquale De Riso, Samuel Filosa.
# Rotary Wing Aerodynamics Course, Prof. Renato Tognaccini.
# University of Naples Federico II.
# Academic Year 2023-2024.
#
# Date: 11/28/2024.
#
# Version: 1.0.0.



import numpy as np
import matplotlib.pyplot as plt
from ambiance import Atmosphere
from scipy.constants import pi, g
from Drone_functions import RigidFW2
from Drone_functions import BEMT_Hover

# Initial parameters.
N = 2                                                                                                # Number of blades.
R = 0.06                                                                                        # Disk rotor radius [m].
A = pi * R**2                                                                                   # Disk rotor area [m^2].
r_segn = np.linspace(0.1, 1, 30)                                      # Dimensionless disk rotor radius.
theta0 = np.radians(20)                                                                        # Collective pitch [rad].
theta_tw = np.radians(-9)                                                                           # Twist pitch [rad].
theta = theta0 + theta_tw * (r_segn - r_segn[0])                                                      # Pitch law [rad].

c = 0.015                                                                                     # Blade element chord [m].
solidity = N * c / (pi * R)                                                                            # Rotor solidity.
N_rotors = 4                                                                                         # Number of rotors.
W = 0.246 * g / N_rotors                                                                      # Single rotor weight [N].
mu = 0                                                                                                  # Advance ratio.
Clalpha = 2 * pi                                                                            # Lift slope curve [rad^-1].
Cd0 = 0.02                                                                                 # Parasitic drag coefficient.
height = 2000                                                                                            # Altitude [m].
atmosphere = Atmosphere(height)
rho_inf = atmosphere.density[0]                                                                  # Air density [kg/m^3].

Tc, _, _, _, _ = BEMT_Hover(N, r_segn, solidity, theta, Clalpha, Cd0, mu)                              # Tc calculation.
Omega_h = np.sqrt(W / (rho_inf * R**2 * A * Tc))                                       # Hover angular velocity [rad/s].

Pc_i = Tc**(3/2) / np.sqrt(2)                                                   # Power induced coefficient, eq. (6.17).
Pc0 = solidity * Cd0 / 8                                                      # Power parasitic coefficient, eq. (6.17).
lambda_i = np.sqrt(Tc / 2)                                                       # Axial interference factor, eq. (6.2).
w_h = lambda_i * Omega_h * R                                                          # Induced velocity in hover [m/s].

P_ih = Pc_i * rho_inf * Omega_h**3 * R**3 * A                                # Induced power in hovering, eq. (2.5) [W].
Pc_0h = Pc0
Vinf_vec = np.linspace(0, 16, 30)                                          # Asymptotic air speed [m/s].
f_A = 0.07                                                                       # Dimensionless equivalent wetted area.
f = f_A * (4 * A)                                                                        # Equivalent wetted area [m^2].

linear_inflow = 0                                                                                    # Unknown function.
results = RigidFW2(Vinf_vec, N, r_segn, c, theta, W, w_h, height, R, A,                      # RigidFW2 function output.
                   P_ih, f, Pc_0h, Clalpha, Cd0, Omega_h, N_rotors, linear_inflow)

(P_i, P_0, P_fus, _, P_new, P_new_i, P_new_0, P_new_fus, alpha_vec,
 Omega_new, alpha_e_vec, Mach_vec, dTc, U_P, dH_dpsi_0, dQ_dpsi_0,
 Tc_vec, Y_vec, Qc_vec) = results                                                     # RigidFW2 function output saving.

Vinf_vec = Vinf_vec * 1.94384                                                                     # From [m/s] to [kts].

# Graphics BEMT vs IT.
#  P_tot.
plt.figure(1)
plt.plot(Vinf_vec, 4 * P_new, '-^', markersize=6, color='k', label='BEMT')
plt.plot(Vinf_vec, 4 * (P_i + P_0 + P_fus / N_rotors), '-s', markersize=6, color='k', label='IT')
plt.legend()
plt.xlabel(r'$V_\infty \, [kts]$', fontsize=12)
plt.ylabel(r'$P_{TOT} \, [W]$', fontsize=12)
plt.grid(which='minor')

# P_fus.
plt.figure(2)
plt.plot(Vinf_vec, P_new_fus, '-^', markersize=6, color='k', label='BEMT')
plt.plot(Vinf_vec, P_fus, '-s', markersize=6, color='k', label='IT')
plt.legend()
plt.xlabel(r'$V_\infty \, [kts]$', fontsize=12)
plt.ylabel(r'$P_{fus} \, [W]$', fontsize=12)
plt.grid(which='minor')

# P_profile.
plt.figure(3)
plt.plot(Vinf_vec, P_new_0, '-^', markersize=6, color='k', label='BEMT')
plt.plot(Vinf_vec, P_0, '-s', markersize=6, color='k', label='IT')
plt.legend()
plt.xlabel(r'$V_\infty \, [kts]$', fontsize=12)
plt.ylabel(r'$P_{0} \, [W]$', fontsize=12)
plt.grid(which='minor')

# P_induced.
plt.figure(4)
plt.plot(Vinf_vec, P_new_i, '-^', markersize=6, color='k', label='BEMT')
plt.plot(Vinf_vec, P_i, '-s', markersize=6, color='k', label='IT')
plt.legend()
plt.xlabel(r'$V_\infty \, [kts]$', fontsize=12)
plt.ylabel(r'$P_{i} \, [W]$', fontsize=12)
plt.grid(which='minor')

# Tc.
plt.figure(5)
plt.plot(Vinf_vec, Tc_vec, 'k')
plt.xlabel(r'$V_\infty \, [kts]$', fontsize=12)
plt.ylabel(r'$T_{c}$', fontsize=12)
plt.grid(which='minor')

# Qc.
plt.figure(6)
plt.plot(Vinf_vec, Qc_vec, 'k')
plt.xlabel(r'$V_\infty \, [kts]$', fontsize=12)
plt.ylabel(r'$Q_{c}$', fontsize=12)
plt.grid(which='minor')

plt.show()
