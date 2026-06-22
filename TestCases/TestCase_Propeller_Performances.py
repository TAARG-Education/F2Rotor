# TEST CASE FOR EVALUATING PROPELLER COEFFICIENTS
# This test case demonstrates how to evaluate a propeller's coefficients using the provided 'aero' and 'geometry' classes,
# as well as the functions from 'bemt.py'.
#
# Authors: Simone Orazzo, Leonardo Parisi, A. Dario Marotta

import Geometry 
import matplotlib.pyplot as plt
import numpy as np
import bemt
import aero

# INPUT DATA
R = 1.60                                                                # Blade radius (m)
R_hub = 0.152                                                           # Hub radius (m)  
N = 2                                                                   # Number of blades
RPM = 1000                                                              # revolutions per minute

# stations along the blade (r/R)
r_R_known = [ 0.0940,  0.1587,    0.1905,    0.2857,    0.3810,    0.4762,  0.5714,    0.6667,    0.7619,    0.8571,    0.9524,    1.0000]
# chord distribution   (c/R)
c_R_known = [ 0.0615, 0.1106,    0.1168,    0.1283 ,   0.1335,    0.1338,    0.1297,  0.1195,    0.1037,    0.0827,    0.0594,       0.02]
# pitch distribution measured from the chord (deg)
beta_known = [42.4,   42.4,      39.1,      35.1,      30.2,      26.7,      23.9,   21.75,     20.01,     18.8,      17.3,      16] 
sweep_known = [0 for b in beta_known]
beta75 = 15.5    #pitch angle imposed at 75% along the blade (deg)

# airfoils along the blade: NACA 4 and 5-digits or 'custom'
airfoil_known = ['0012','0012', '0012',  '0012', '0012','0012','0012', '0012','0012','0012','0012','0012']

pitch = 0.0                                                             # measured from nominal pitch (deg)
v_J=np.linspace(0.152, 0.856, 80)                                       # range of J=V/nD
interp_kind = 'cubic'
dx = 0.0254                                                             # new spacing between the stations used in numerical integration (m)
z = 0                                                                   # altitude (km)

CT_1 = []
CP_1 = []
ETA_1 = []

CT_2 = []
CP_2 = []
ETA_2 = []

CT_3 = []
CP_3 = []
temp = []
ETA_3 = []

xrot_params = {
    'alpha_0_lift': np.deg2rad(-3.1),
    'Cl_alpha': 6.28,
    'Cl_max': 1.5,
    'Cl_min':-0.5,
    'Cd': 0.01,
    'Ka': 0.87
}

for J in v_J:

    # Class iniziatization
    g = Geometry.Geometry( R, R_hub, N, RPM, r_R_known, c_R_known, beta_known,sweep_known, beta75, airfoil_known, pitch,interp_kind)

    # Defining aerodynamic method: method 1 
    # NOTE: as in aero version 1.1, sweep is implemented only in the aerodynamic method 1
    Aero = aero.Aerodynamics(1, True, True, xrot_params)

    #Computing C_T, C_P and eta and storing them in a vector
    ct_1, cp_1=bemt.BEMT_timp(z, J, dx, g, Aero)

    CT_1.append(ct_1)
    CP_1.append(cp_1)
    if ct_1 >= 0:
        ETA_1.append(J*ct_1/cp_1)
    else:
        ETA_1.append(np.nan)

    ct_2, cp_2=bemt.BEMT_tvorpd(z, J, dx, g, Aero)

    CT_2.append(ct_2)
    CP_2.append(cp_2)
    if ct_2 >= 0:
        ETA_2.append(J*ct_2/cp_2)
    else:
        ETA_2.append(np.nan)

    ct_3, cp_3 =bemt.BEMT_tvor(z, J, dx, g, Aero)

    CT_3.append(ct_3)
    CP_3.append(cp_3)
    if ct_3 >= 0:
        ETA_3.append(J*ct_3/cp_3)
    else:
        ETA_3.append(np.nan)


# Array conversion
Ct_1=np.array(CT_1)
Cp_1=np.array(CP_1)
Eta_1=np.array(ETA_1)

Ct_2=np.array(CT_2)
Cp_2=np.array(CP_2)
Eta_2=np.array(ETA_2)

Ct_3=np.array(CT_3)
Cp_3=np.array(CP_3)
Eta_3=np.array(ETA_3)


#Plot
plt.figure(figsize=(5, 6))
plt.plot(v_J, Ct_1, '--', label='Ct: timp', color='k')
plt.plot(v_J, Ct_2, '-.', label='Ct: tvorpd', color='k')
plt.plot(v_J, Ct_3, '-', label='Ct: tvor', color='k')
plt.plot(v_J, Cp_1, '--', label='Cp: timp', color='r')
plt.plot(v_J, Cp_2, '-.', label='Cp: tvorpd', color='r')
plt.plot(v_J, Cp_3, '-', label='Cp: tvor', color='r')
plt.xlabel('J=V/nD')
plt.legend()
plt.grid(True)

plt.figure(figsize=(5,6))
plt.plot(v_J, Eta_1, '--', label='$\eta$: timp', color='k')
plt.plot(v_J, Eta_2, '-.', label='$\eta$: tvorpd', color='k')
plt.plot(v_J, Eta_3, '-', label='$\eta$: tvor', color='k')

plt.xlabel('J=V/nD')
plt.legend()
plt.grid(True)
plt.show()