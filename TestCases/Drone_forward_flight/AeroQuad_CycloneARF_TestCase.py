# DRONE IN FORWARD FLIGHT F2ROTOR
#
# Description: This is the test case function in order to determine power curve contributions, thrust and torque
#              coefficients with Blade Element Momentum Theory (BEMT) and simple impulse theory (IT) for the AeroQuad
#              Cyclone ARF Drone.
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
from Drone_functions import Hover_state

N = 2                                                                                                # Number of blades.
R = 0.1524                                                                                           # Rotor radius [m].
theta0 = np.deg2rad(21.5)                                                                      # Root pitch angle [rad].
theta_end = np.deg2rad(11.1)                                                                    # Tip pitch angle [rad].
c_r = 0.031                                                                                            # Root chord [m].
c_t = 0.012                                                                                             # Tip chord [m].
W = 2*9.81                                                                                                 # Weight [N].

A = pi * R**2                                                                                   # Disk rotor area [m^2].
r_segn = np.linspace(0.1, 1, 20)                                      # Dimensionless disk rotor radius.
theta_tw = (theta_end-theta0)/(r_segn[-1]-r_segn[0])                                          # Twist pitch angle [rad].
c_tw = (c_t-c_r)/(r_segn[-1]-r_segn[0])                                                               # Twist chord [m].
theta = theta0 + theta_tw*(r_segn-r_segn[0])                                                          # Pitch law [rad].
c = c_r + c_tw*(r_segn-r_segn[0])                                                                       # Chord law [m].

solidity = N * c / (pi * R)                                                                            # Rotor solidity.
N_rotors = 4                                                                                         # Number of rotors.

W = W/N_rotors                                                                                # Single rotor weight [N].
mu = 0                                                                                                  # Advance ratio.
Clalpha = 2 * pi                                                                            # Lift slope curve [rad^-1].
Cd0 = 0.02                                                                                 # Parasitic drag coefficient.
height = 0                                                                                               # Altitude [m].
atmosphere = Atmosphere(height)                                        # Introducing atmosphere model, according to S.I.
rho_inf = atmosphere.density[0]                                                                  # Air density [kg/m^3].


Tc, Qc, Pc_i_hover, Pc_0_hover, _, dTcdr_segn_PR, dPcdr_segn, phi_hover, alpha_hover, lambda_i_hover = \
    Hover_state(N, r_segn, solidity, theta, Clalpha, Cd0, mu)           # Tc, Qc, Pc_i_hover and Pc_0_hover calculation.
Omega_h = np.sqrt(W / (rho_inf * R**2 * A * Tc))                                       # Hover angular velocity [rad/s].
dT_drsegn = dTcdr_segn_PR*rho_inf*Omega_h**2*pi*R**4/R/N                                                 # From Tc to T.
dQ_drsegn = dPcdr_segn*rho_inf*Omega_h**2*pi*R**5/R/N                                                    # From Qc to Q.
Q_hover = Qc*rho_inf*Omega_h**2*pi*R**5                                                       # Torque in hovering [Nm].

Pc_i = Pc_i_hover                                                               # Power induced coefficient, eq. (6.17).
Pc0  = Pc_0_hover                                                             # Power parasitic coefficient, eq. (6.17).
lambda_i = np.sqrt(Tc / 2)                                                       # Axial interference factor, eq. (6.2).
w_h = lambda_i * Omega_h * R                                                          # Induced velocity in hover [m/s].

P_ih = Pc_i * rho_inf * Omega_h**3 * R**3 * A                                # Induced power in hovering, eq. (2.5) [W].
Pc_0h = Pc0                                                                   # Parasitic power coefficient in hovering.
Vinf_vec = np.linspace(0, 14, 30)                                          # Asymptotic air speed [m/s].
f_A = 0.16                                                                       # Dimensionless equivalent wetted area.
f = f_A * (4*A)                                                                          # Equivalent wetted area [m^2].

results = RigidFW2(Vinf_vec, N, r_segn, c, theta, W, w_h, height, R, A,                      # RigidFW2 function output.
                   P_ih, f, Pc_0h, Clalpha, Cd0, Omega_h, N_rotors)

(P_i, P_0, P_fus, P, P_new, P_new_i, P_new_0, P_new_fus, alpha_vec,
 Omega_new, alpha_e_vec, dTc, U_P, dH_dpsi_0, dQ_dpsi_0,
 T_vec, Y_vec, Q_vec, dT_dr, dQ_dr) = results                                         # RigidFW2 function output saving.

Omega_new_RPM = Omega_new*60/(2*pi)                                                             # From [rad/s] to [RPM].

# Graphics.
#  P_tot.
plt.figure(1)
plt.plot(Vinf_vec, 4*P_new, '-^', markersize=6, color='k', label='BEMT')
plt.plot(Vinf_vec, 4*P, '-s', markersize=6, color='k', label='IT')
plt.legend()
plt.xlabel(r'$V_\infty \, [m/s]$', fontsize=12)
plt.ylabel(r'$P_{TOT} \, [W]$', fontsize=12)
plt.minorticks_on()
plt.grid(which='major', linestyle='-', linewidth=0.75)                                                     # Major grid.
plt.grid(which='minor', linestyle='--', linewidth=0.5)                                                     # Minor grid.

# P_fus.
plt.figure(2)
plt.plot(Vinf_vec, P_new_fus, '-^', markersize=6, color='k', label='BEMT')
plt.plot(Vinf_vec, P_fus, '-s', markersize=6, color='k', label='IT')
plt.legend()
plt.xlabel(r'$V_\infty \, [m/s]$', fontsize=12)
plt.ylabel(r'$P_{fus} \, [W]$', fontsize=12)
plt.minorticks_on()
plt.grid(which='major', linestyle='-', linewidth=0.75)
plt.grid(which='minor', linestyle='--', linewidth=0.5)

# P_profile.
plt.figure(3)
plt.plot(Vinf_vec, P_new_0, '-^', markersize=6, color='k', label='BEMT')
plt.plot(Vinf_vec, P_0, '-s', markersize=6, color='k', label='IT')
plt.legend()
plt.xlabel(r'$V_\infty \, [m/s]$', fontsize=12)
plt.ylabel(r'$P_{0} \, [W]$', fontsize=12)
plt.minorticks_on()
plt.grid(which='major', linestyle='-', linewidth=0.75)
plt.grid(which='minor', linestyle='--', linewidth=0.5)

# P_induced.
plt.figure(4)
plt.plot(Vinf_vec, P_new_i, '-^', markersize=6, color='k', label='BEMT')
plt.plot(Vinf_vec, P_i, '-s', markersize=6, color='k', label='IT')
plt.legend()
plt.xlabel(r'$V_\infty \, [m/s]$', fontsize=12)
plt.ylabel(r'$P_{i} \, [W]$', fontsize=12)
plt.minorticks_on()
plt.grid(which='major', linestyle='-', linewidth=0.75)
plt.grid(which='minor', linestyle='--', linewidth=0.5)

# Thrust.
plt.figure(5)
plt.plot(Vinf_vec, T_vec, 'k')
plt.xlabel(r'$V_\infty \, [m/s]$', fontsize=12)
plt.ylabel(r'$T \, [N]$', fontsize=12)
plt.minorticks_on()
plt.grid(which='major', linestyle='-', linewidth=0.75)
plt.grid(which='minor', linestyle='--', linewidth=0.5)
# Calculate y limits
y_min = min(T_vec) - 0.3
y_max = max(T_vec) + 0.3
plt.ylim(y_min, y_max)

# Torque.
plt.figure(6)
plt.plot(Vinf_vec, Q_vec, 'k')
plt.xlabel(r'$V_\infty \, [m/s]$', fontsize=12)
plt.ylabel(r'$Q \, [Nm]$', fontsize=12)
plt.minorticks_on()
plt.grid(which='major', linestyle='-', linewidth=0.75)
plt.grid(which='minor', linestyle='--', linewidth=0.5)

# Omega_new.
plt.figure(7)
plt.plot(Vinf_vec, Omega_new_RPM, 'k')
plt.xlabel(r'$V_\infty \, [m/s]$', fontsize=12)
plt.ylabel(r'$\Omega \, [RPM]$', fontsize=12)
plt.minorticks_on()
plt.grid(which='major', linestyle='-', linewidth=0.75)
plt.grid(which='minor', linestyle='--', linewidth=0.5)

# dT_dr.
plt.figure(8)
plt.plot(r_segn, dT_dr[0,:,0], '-^', markersize=6, color='k', label='BEMT')
plt.plot(r_segn, dT_drsegn, '-s', markersize=6, color='k', label='HOVER_STATE')
plt.legend()
plt.xlabel(r'$\bar{r}$', fontsize=12)
plt.ylabel(r'$dT/{d\bar{r}}$', fontsize=12)
plt.minorticks_on()
plt.grid(which='major', linestyle='-', linewidth=0.75)
plt.grid(which='minor', linestyle='--', linewidth=0.5)

# dQ_dr.
#dQ_dr[0,:] = dQ_dr[0,:]*((Omega_h/Omega_new[0])**2)
plt.figure(9)
plt.plot(r_segn, dQ_dr[0,:,0], '-^', markersize=6, color='k', label='BEMT')
plt.plot(r_segn, dQ_drsegn, '-s', markersize=6, color='k', label='HOVER_STATE')
plt.legend()
plt.xlabel(r'$\bar{r}$', fontsize=12)
plt.ylabel(r'$dQ/{d\bar{r}}$', fontsize=12)
plt.minorticks_on()
plt.grid(which='major', linestyle='-', linewidth=0.75)
plt.grid(which='minor', linestyle='--', linewidth=0.5)

# Effective alpha vs r_bar.
plt.figure(10)
plt.plot(r_segn,alpha_e_vec[0,:,0],'-^', markersize=6, color='k', label='BEMT')
plt.plot(r_segn,alpha_hover, '-s', markersize=6, color='k', label='HOVER_STATE')
plt.legend()
plt.xlabel(r'$\bar{r}$', fontsize=12)
plt.ylabel(r'$\alpha_e \, [rad]$', fontsize=12)
plt.minorticks_on()
plt.grid(which='major', linestyle='-', linewidth=0.75)
plt.grid(which='minor', linestyle='--', linewidth=0.5)

#Effective alpha vs V_inf.
plt.figure(11)
plt.plot(Vinf_vec,alpha_e_vec[:,5,0],'k')
plt.xlabel(r'$V_\infty \, [m/s]$', fontsize=12)
plt.ylabel(r'$\alpha_e \, [rad]$', fontsize=12)
plt.minorticks_on()
plt.grid(which='major', linestyle='-', linewidth=0.75)
plt.grid(which='minor', linestyle='--', linewidth=0.5



plt.show()
