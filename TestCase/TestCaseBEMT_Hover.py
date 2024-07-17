from BEMT_Hover import *
from Helicopter_Class import Helicopter
import numpy as np


# This test case aims to validate the function "BEMT_Hover.py", which uses the Blade Element
# Momentum Theory for an helicopter in hover to carachterize the hover (or slow climbing)
# condition for an helicopter. Specifically, given the structural carachteristics of the
# helicopter and the aerodynamic conditions, it calculate the Figure of Merit, the Prandtl
# correction function and the Thrust and Power distribution along the rotor's blade and the 
# relative coefficients.
# For this Test Case an helicopter class has been defined using the structural properties
# of the helicopter Augusta A109, whose image can be seen in the documentation.

# Authors: Francesco Gervasio

# Input data

V_infty = 0            # Ascensional flight speed (m/s):
Omega_R_mr = 224.2067  # Main rotor tip speed (m/s)
R_mr = 5.502           # Main rotor radius (m)
R_hub = 0.05502        # Rotor hub radius (m)
c_mr = 0.3355          # Rotor's chord (m)
N_mr = 4               # Number of blades of the helicopter
theta0 = 70            # Pitch angle at the root (°)
theta_tw = -6          # Blade twist (tip minus root incidence) (°)
Cla = 6.28             # Lift curve slope of the airfoil (1/rad)
Cd0 = 0.0150           # Profile drag coefficient of the airfoil


# Helicopter Definition
Helicopter = Helicopter(Omega_R_mr, R_mr, R_hub, c_mr, N_mr, theta0, theta_tw, Cla, Cd0)

# Calculation of the coefficients
[Tc_PR, Pc_i, Pc_0, Qc, dTcdr_bar_PR, dPcdr_bar, FM, F_Prandtl] = BEMT_Hover(Helicopter, V_infty)

# Printing the results
print("The Thrust coefficient (considering the Prandtl correction function) Tc is:" + str(Tc_PR))
print("The Torque (Power) coefficient Qc is: " + str(Qc))
print("The induced Power (Torque) coefficient Pc_i is: " + str(Pc_i))
print("The parasite Power (Torque) coefficient Pc_0 is: " + str(Pc_0))
print("The Figure of Merit FM is: " + str(FM))

