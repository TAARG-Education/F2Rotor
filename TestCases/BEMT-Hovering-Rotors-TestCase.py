# BEMT COUNTER-ROTATING ROTORS IN HOVERING
#
# Description:
# This function is part of the set of functions that computes the performances of 
# a counter-rotating configuration of two rotors in hovering.
# Specifically, this function calls the functions 'TeoriaElemPalaRotoreUpper' and 'TeoriaElemPalaRotoreLower'
# and stores the output values; The trends of both the upper and lower rotors are exploited. 
#
# The Hypothesis assumed are:
# - Inviscid Flow.
# - Incompressible Flow.
# - It is assumed to be working in the linear trend of the lift curve, whose slope is assumed to be 2*pi.
# - a' = 0. 
# - phi << 1.
# - The lower rotor does not influence by any mean the upper rotor. 
# - The Influence of the upper rotor is consistent. An interference factor is assumed: r_segn_interf = 0.81. 
#   The interference factor accounts the elements of the lower rotor's blades that are affected by the wake of the upper rotor.
#   This values may be changed by the user if needed. 
#
# Authors of the Function: Giovanni Perez, Catello Meglio.
# Rotary Wing Aerodynamics Course, Prof. Renato Tognaccini.
# University of Naples Federico II.
# Academic year 2023/2024.
#

import numpy as np
from scipy.integrate import trapezoid
import matplotlib.pyplot as plt
from TeoriaElemPalaRotoreUpper import TeoriaElemPalaRotoreUpper
from TeoriaElemPalaRotoreLower import TeoriaElemPalaRotoreLower

# Example variables
r_segn = np.linspace(0.1, 1, 100) # Adimensionalized element blase distance
N      = 3 # Number of blades
solidity = np.linspace(0.05, 0.1, 100) # Solidity vector
theta = np.deg2rad(np.linspace(0, 10, 100)) # Twist angle vector.

Clalpha = 2*np.pi 
Cd0 = 0.01 
r_segn_interf = 0.81 # Distance of the element of the blades that are affected by the wake of the upper rotor.

Results_Upper = TeoriaElemPalaRotoreUpper(r_segn, N, solidity, theta, Clalpha, Cd0) 
Results_Lower = TeoriaElemPalaRotoreLower(r_segn, r_segn_interf, N, solidity, theta, Clalpha, Cd0, Results_Upper[-1]) 

# Plots
fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
fig1.suptitle('Prandtl Correction Factor', fontsize=16)

ax1.plot(r_segn, Results_Upper[13], color='b')
ax1.set_title('Prandtl Correction Factor for the Upper rotor')
ax1.set_xlabel('r')
ax1.set_ylabel('F')
ax1.grid(True, linestyle='--', linewidth=0.5) 

ax2.plot(r_segn, Results_Lower[13], color='b')
ax2.set_title('Prandtl Correction Factor for the Lower rotor')
ax2.set_xlabel('r')
ax2.set_ylabel('F')
ax2.grid(True, linestyle='--', linewidth=0.5) 

plt.tight_layout(rect=[0, 0, 1, 0.95])  

fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
fig2.suptitle('Elemental Thrust Coefficient per unit length', fontsize=16)

ax1.plot(r_segn, Results_Upper[0], color='b')
ax1.plot(r_segn, Results_Upper[1], color='r')
ax1.set_title('Thrust Coefficient per unit length for the Upper rotor')
ax1.set_xlabel('r')
ax1.set_ylabel('dTcdr')
#ax1.legend('dTcdr','dTcdr_Pr')
ax1.grid(True, linestyle='--', linewidth=0.5) 

ax2.plot(r_segn, Results_Lower[0], color='b')
ax2.plot(r_segn, Results_Lower[1], color='r')
ax2.set_title('Thrust Coefficient per unit length for the Lower rotor')
ax2.set_xlabel('r')
ax2.set_ylabel('dTcdr')
#ax2.legend('dTcdr','dTcdr_Pr')
ax2.grid(True, linestyle='--', linewidth=0.5) 

plt.tight_layout(rect=[0, 0, 1, 0.95])  

fig3, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
fig3.suptitle('Elemental Torque Coefficient per unit length', fontsize=16)

ax1.plot(r_segn, Results_Upper[2], color='b')
ax1.plot(r_segn, Results_Upper[3], color='r')
ax1.set_title('Torque Coefficient per unit length for the Upper rotor')
ax1.set_xlabel('r')
ax1.set_ylabel('dQcdr')
#ax1.legend('dTcdr','dTcdr_Pr')
ax1.grid(True, linestyle='--', linewidth=0.5) 

ax2.plot(r_segn, Results_Lower[2], color='b')
ax2.plot(r_segn, Results_Lower[3], color='r')
ax2.set_title('Torque Coefficient per unit length for the Lower rotor')
ax2.set_xlabel('r')
ax2.set_ylabel('dQcdr')
#ax2.legend('dTcdr','dTcdr_Pr')
ax2.grid(True, linestyle='--', linewidth=0.5) 

plt.tight_layout(rect=[0, 0, 1, 0.95])  

fig3, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
fig3.suptitle('Elemental inducted Power', fontsize=16)

ax1.plot(r_segn, Results_Upper[4], color='b')
ax1.set_title('Induced power for the upper rotor')
ax1.set_xlabel('r')
ax1.set_ylabel('dP_i')
#ax1.legend('dTcdr','dTcdr_Pr')
ax1.grid(True, linestyle='--', linewidth=0.5) 

ax2.plot(r_segn, Results_Lower[4], color='b')
ax2.set_title('Induced power for the Lower rotor')
ax2.set_xlabel('r')
ax2.set_ylabel('dP_i')
#ax2.legend('dTcdr','dTcdr_Pr')
ax2.grid(True, linestyle='--', linewidth=0.5) 

plt.tight_layout(rect=[0, 0, 1, 0.95])  

fig4, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
fig4.suptitle('Elemental parasite Power', fontsize=16)

ax1.plot(r_segn, Results_Upper[5], color='b')
ax1.set_title('parasite power for the upper rotor')
ax1.set_xlabel('r')
ax1.set_ylabel('dP_0')
#ax1.legend('dTcdr','dTcdr_Pr')
ax1.grid(True, linestyle='--', linewidth=0.5) 

ax2.plot(r_segn, Results_Lower[5], color='b')
ax2.set_title('parasite power for the Lower rotor')
ax2.set_xlabel('r')
ax2.set_ylabel('dP_0')
#ax2.legend('dTcdr','dTcdr_Pr')
ax2.grid(True, linestyle='--', linewidth=0.5) 

plt.tight_layout(rect=[0, 0, 1, 0.95])  

plt.show()

# Print the values of the coefficients 
print('Upper rotor variables values')
variables_upper = {
    'Tc_PR_u': Results_Upper[6],
    'Tc_u': Results_Upper[7],
    'Qc_u': Results_Upper[8],
    'Qc_PR_u': Results_Upper[9],
    'Pc_i_u': Results_Upper[10],
    'Pc_0_u': Results_Upper[11],
    'FM_u': Results_Upper[12],
}

for nome, valore in variables_upper.items():
    print(f"{nome} = {valore}")

print('Lower rotor variables values')
variables_lower = {
'Tc_PR_l': Results_Lower[6],
'Tc_l': Results_Lower[7],
'Qc_l': Results_Lower[8],
'Qc_PR_l': Results_Lower[9],
'Pc_i_l': Results_Lower[10],
'Pc_0_l': Results_Lower[11],
'FM_l': Results_Lower[12],
}

for nome, valore in variables_lower.items():
    print(f"{nome} = {valore}")
