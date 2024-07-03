from Helicopter_Class import Helicopter
from TailRotor import *
import numpy as np
import matplotlib.pyplot as plt

'''
This test case aims to validate the function "TailRotor.py," which calculates 
the power coefficients of a helicopter's tail rotor. Specifically, the induced 
power coefficient (Cp_tr_i), the profile power coefficient (Cp_tr_0), and their 
sum (Cp_tr = Cp_tr_0 + Cp_tr_i). To achieve this, once the coefficients values have
been obtained, they have been dimensionalized to power. This has been done for a 
speed variation from 0 m/s to 88 m/s (170 kts) to match well-known power curves. 
To obtain the torque from the main rotor, a dedicated function has been created 
("TorqueCoefficient"), which calculates it at an altitude of 0 m 
(air density rho = 1.225 kg/m}^3).

Input data used are from the Agusta A109 Helicopter.

Authors: Fabio Beltratti, Giuseppe Russo

'''

# Main Rotor Input Data
Omega_r_mr = 221.59         # Main rotor tip speed [m/s]
R_mr = 5.4986               # Main rotor radius [m]
c_mr = 0.3353               # Main rotor chord [m]
N_mr = 4                    # Main rotor blades number [-]

# Tail Rotor Input Data
Omega_r_tr = 221.59         # Tail rotor tip speed [m/s]
R_tr = 1.015                # Tail rotor radius [m]
c_tr = 0.2012               # Tail rotor chord [m]
N_tr = 2                    # Tail rotor blades number [-]
l_tr = 7.1015               # Tail rotor arm [m]

# Helicopter Class Definition (0 values are related to not used input in this test case)
helicopter = Helicopter(Omega_r_mr, R_mr, 0, c_mr, N_mr, 0, 0, 0, 0, Omega_r_tr, R_tr, 0, c_tr, N_tr, l_tr)

# Other Input Data
rho = 1.225                 # S/L air density [kg/m^3]
mass = 2601.42              # helicopter TO mass [kg]         

# Arrays for speed and powers
V_infty_values = range(89)  # V_infty from 0 to 88 m/s (170 kts)
P_tr_i_values = []          # Induced Power Array
P_tr_0_values = []          # Profile Power Array
P_tr_values = []            # Total Power Array

# Torque Coefficient Function Definition
def TorqueCoefficient(Helicopter, mass, V_infty, rho):

    Omega_r_mr = Helicopter.Omega_r_mr            # Main rotor tip speed [m/s]
    R_mr = Helicopter.R_mr                        # Main rotor radius [m]
    A_mr = np.pi*R_mr**2                          # Main rotor area [m^2]
    N_mr = Helicopter.N_mr                        # Main rotor blades number [-]                 
    c_mr = Helicopter.c_mr                        # Main rotor chord [m]                                

    CD_mr = 0.0121                                # Main rotor drag coefficient [-]   
    fA = 0.009                                    # Main rotor friction coefficient * area[-]  
    f = fA*A_mr                                   # Main rotor friction coefficient [-]
    mu_mr = V_infty / Omega_r_mr                  # Main rotor advance ratio (hp: AoA << 1)

    # Induced input ratio: induced inflow ratio
    lambda_it = np.sqrt(- V_infty ** 2 / 2 + 0.5 * np.sqrt(V_infty ** 4 
                                                     + 4*(mass * 9.81 / (2 * rho * A_mr)) ** 2)) / Omega_r_mr  
    k_mr = 1.20                                   # induced power losses factor
    
    # Powers calculation
    P_i_mr = k_mr * lambda_it * Omega_r_mr * mass * 9.81    # main rotor induced power 
    P_prf_mr = 1 / 8 * rho * N_mr * c_mr * R_mr * CD_mr * Omega_r_mr** 3 * (1 + 4.7 * mu_mr ** 2) # main rotor profile power
    P_0_mr = 0.5* rho * f * V_infty ** 3 # main rotor parasite power

    P_tot = P_i_mr + P_0_mr + P_prf_mr  # main rotor required power
    
    # Main rotor torque and torque coefficient calculation
    Q = P_tot * R_mr/Omega_r_mr
    Qc = Q / (rho * Omega_r_mr ** 2 * R_mr * A_mr)

    return Qc

# Calculate power for each speed value
for V_infty in V_infty_values:
    Q_c = TorqueCoefficient(helicopter, mass, V_infty, rho)
    Cp_tr_i, Cp_tr_0, Cp_tr = TailRotor(helicopter, V_infty, Q_c)
    
    P_tr_i_values.append(Cp_tr_i * rho * Omega_r_tr**3 * np.pi * R_tr**2)  # tail rotor induced power
    P_tr_0_values.append(Cp_tr_0 * rho * Omega_r_tr**3 * np.pi * R_tr**2)  # tail rotor parasite power
    P_tr_values.append(Cp_tr * rho * Omega_r_tr**3 * np.pi * R_tr**2)      # tail rotor total power

print(f"Qc: {Q_c}")
print(f"Cp_tr_i: {Cp_tr_i}, Cp_tr_0: {Cp_tr_0}, Cp_tr: {Cp_tr}")

# Plotting all three powers in one graph with LaTeX formatting
plt.figure(figsize=(10, 6))

# Plot of Induced Power vs Speed
plt.plot(V_infty_values, P_tr_i_values, label=r'$P_{tr_i}$', color='b')

# Plot of Profile Power vs Speed
plt.plot(V_infty_values, P_tr_0_values, label=r'$P_{tr_0}$', color='r')

# Plot of Total Power vs Speed
plt.plot(V_infty_values, P_tr_values, label=r'$P_{tr}$', color='g')

# Formatting the plot
plt.xlabel(r'$V_{\infty}$ (m/s)')
plt.ylabel('Power (W)')
plt.title('Tail Rotor Power')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

