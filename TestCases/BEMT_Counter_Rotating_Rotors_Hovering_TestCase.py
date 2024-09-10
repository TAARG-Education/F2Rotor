# BEMT COUNTER-ROTATING ROTORS IN HOVERING
#
# Description:
# This function calls the functions 'BEMT_Counter_Rotating_Rotors_Hovering'. 
# Given the inputs parameters requested by the function, it stores the resulting values and plots the outputs. 
#
# Authors of the Function: Giovanni Perez, Catello Meglio.
# Rotary Wing Aerodynamics Course, Prof. Renato Tognaccini.
# University of Naples Federico II.
# Academic year 2023/2024.

import numpy as np
import matplotlib.pyplot as plt
from BEMT_Counter_Rotating_Rotors_Hovering import BEMT_Counter_Rotating_Rotors_Hovering

# Number of points with which discretize the domain
N_points = 100

# Example variables for the upper rotor
r_segn_u = np.linspace(0.1, 1, N_points) # Adimensionalized element blase distance
N_u      = 2 # Number of blades
solidity_u = 0.027*np.ones(N_points) # Solidity vector
theta_u = np.deg2rad(np.linspace(0, 10, N_points)) # Twist angle vector.
Clalpha_u = 2*np.pi 
Cd0_u = 0.011 

# Example variables for the lower rotor
r_segn_l = np.linspace(0.1, 1, N_points) # Adimensionalized element blase distance
N_l      = 2 # Number of blades
solidity_l = 0.027*np.ones(N_points) # Solidity vector
theta_l = np.deg2rad(np.linspace(0, 10, N_points)) # Twist angle vector.
Clalpha_l = 2*np.pi 
Cd0_l = 0.011 

r_segn_interf = 0.81 # Distance of the element of the blades that are affected by the wake of the upper rotor.

# Results are considered as lists so the notation works as follow:
# Results[k][i]: 
# 'k' selects the output of the function (0-14) ('-1' selects the last element of the array and is equivalent to 14)
# 'i' selects the upper and lower rotors (0-1)
Results = BEMT_Counter_Rotating_Rotors_Hovering(r_segn_u, r_segn_l, N_u, N_l, solidity_u, solidity_l, theta_u, theta_l, Clalpha_u, Clalpha_l, Cd0_u, Cd0_l, r_segn_interf, N_points)

# Plots
fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
fig1.suptitle('Prandtl Correction Factor', fontsize=16)

ax1.plot(Results[0][0], Results[-1][0], color='b')
ax1.set_title('Prandtl Correction Factor for the Upper rotor')
ax1.set_xlabel('r')
ax1.set_ylabel('F')
ax1.grid(True, linestyle='--', linewidth=0.5) 

ax2.plot(Results[0][1], Results[-1][1], color='b')
ax2.set_title('Prandtl Correction Factor for the Lower rotor')
ax2.set_xlabel('r')
ax2.set_ylabel('F')
ax2.grid(True, linestyle='--', linewidth=0.5) 

fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
fig2.suptitle('Elemental Thrust Coefficient per unit length', fontsize=16)

ax1.plot(Results[0][0], Results[1][0], color='b')
ax1.plot(Results[0][0], Results[2][0], color='r')
ax1.set_title('Thrust Coefficient per unit length for the Upper rotor')
ax1.set_xlabel('r')
ax1.set_ylabel('dTcdr')
ax1.grid(True, linestyle='--', linewidth=0.5) 

ax2.plot(Results[0][1], Results[1][1], color='b')
ax2.plot(Results[0][1], Results[2][1], color='r')
ax2.set_title('Thrust Coefficient per unit length for the Lower rotor')
ax2.set_xlabel('r')
ax2.set_ylabel('dTcdr')
ax2.grid(True, linestyle='--', linewidth=0.5) 

plt.tight_layout(rect=[0, 0, 1, 0.95])  

fig3, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
fig3.suptitle('Elemental Torque Coefficient per unit length', fontsize=16)

ax1.plot(Results[0][0], Results[3][0], color='b')
ax1.plot(Results[0][0], Results[4][0], color='r')
ax1.set_title('Torque Coefficient per unit length for the Upper rotor')
ax1.set_xlabel('r')
ax1.set_ylabel('dQcdr')
ax1.grid(True, linestyle='--', linewidth=0.5) 

ax2.plot(Results[0][1], Results[3][1], color='b')
ax2.plot(Results[0][1], Results[4][1], color='r')
ax2.set_title('Torque Coefficient per unit length for the Lower rotor')
ax2.set_xlabel('r')
ax2.set_ylabel('dQcdr')
ax2.grid(True, linestyle='--', linewidth=0.5) 

plt.tight_layout(rect=[0, 0, 1, 0.95]) 

fig3, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
fig3.suptitle('Elemental induced Power', fontsize=16)

ax1.plot(Results[0][0], Results[5][0], color='b')
ax1.set_title('Induced power for the upper rotor')
ax1.set_xlabel('r')
ax1.set_ylabel('dP_i')
ax1.grid(True, linestyle='--', linewidth=0.5) 

ax2.plot(Results[0][1], Results[5][1], color='b')
ax2.set_title('Induced power for the Lower rotor')
ax2.set_xlabel('r')
ax2.set_ylabel('dP_i')
ax2.grid(True, linestyle='--', linewidth=0.5) 

plt.tight_layout(rect=[0, 0, 1, 0.95])  
 
fig4, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
fig4.suptitle('Elemental parasite Power', fontsize=16)

ax1.plot(Results[0][0], Results[6][0], color='b')
ax1.set_title('parasite power for the upper rotor')
ax1.set_xlabel('r')
ax1.set_ylabel('dP_0')
ax1.grid(True, linestyle='--', linewidth=0.5) 

ax2.plot(Results[0][1], Results[6][1], color='b')
ax2.set_title('parasite power for the Lower rotor')
ax2.set_xlabel('r')
ax2.set_ylabel('dP_0')
ax2.grid(True, linestyle='--', linewidth=0.5) 

plt.tight_layout(rect=[0, 0, 1, 0.95])  

plt.show()

# Print the values of the coefficients 
print('Upper rotor variables values')
variables_upper = {
    'Tc_PR_u': Results[7][0],
    'Tc_u': Results[8][0],
    'Qc_u': Results[9][0],
    'Qc_PR_u': Results[10][0],
    'Pc_i_u': Results[11][0],
    'Pc_0_u': Results[12][0],
    'FM_u': Results[13][0],
}

for nome, valore in variables_upper.items():
    print(f"{nome} = {valore}")

print('Lower rotor variables values')
variables_lower = {
    'Tc_PR_l': Results[7][1],
    'Tc_l': Results[8][1],
    'Qc_l': Results[9][1],
    'Qc_PR_l': Results[10][1],
    'Pc_i_l': Results[11][1],
    'Pc_0_l': Results[12][1],
    'FM_l': Results[13][1],
}

for nome, valore in variables_lower.items():
    print(f"{nome} = {valore}")
