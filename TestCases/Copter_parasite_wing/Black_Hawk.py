
#This test case script calculates aerodynamic coefficients for the horizontal stabilizer of a simulated Sirosky UH-60 Black Hawk using Python.
#The script utilizes parameters such as wing span, wing chord, wing thickness ratio, and aspect ratio specific to the copter selected, all defined in a configuration dictionary (config).
#Additionally, flight parameters such as cruise speed, altitude, speed of sound, kinematic viscosity, and density are specified to compute the Mach number and Reynolds number necessary for aerodynamic calculations.

#Author: Eduardo Duraccio
#Date: 17/07/2024
#Version: 1.02

import matplotlib.pyplot as plt
import math
import numpy as np
from Copter_parasite_wing import parasite_wing
from Helicopter_HS_Class import Helicopter_HS

b_w = 4.38      # Wing span [m] 
c_w = 0.95      # wing chord [m]  
tau = 0.14      # Wing thickness ratio 
AR = 4.61       # wing aspect ratio 
R_mr = 8.18     # Radius of main rotor [m] 
R_LS =1.08      # lift-induced skin friction factor
Hf =1.1         #hinge factor
r = 0.89        #compressibility correction factor
gamma = 1.4     # ratio of specific heats for air
K =0.0009       #admissible roughness
R_w_b = 1.07    # interference factor wing-fusolage
ala_fissa = False # fixed or rotary wing
profilo = "4-series" # type of airfoil
turbolentflow = True # laminar or turbolent flow


# Additional parameters needed for calculation of flight condition
v = 49.5833  # 50% of maximum speed  [m/s]
height = 2000     # Altitude [m] 
a = 332.529       # Speed of sound [m/s]
nu = 1.7148e-5    # Kinematic viscosity [m/s]
rho = 1.0065      # Density [kg/m³]


helicopter_HS = Helicopter_HS(tau=tau, b_w=b_w, AR=AR, R_LS=R_LS, Hf=Hf, r=r, gamma=gamma, k=K, R_mr=R_mr, 
                              c_w=c_w, R_w_b=R_w_b, ala_fissa= ala_fissa, profilo= profilo, turbolentflow= turbolentflow)
helicopter_HS.set_M_inf(v, a)
helicopter_HS.set_Re_inf(v, rho, nu)






# Call the wing_comp function with the simulated aircraft instance
f_Awing, f_wing = parasite_wing(helicopter_HS)

# Print the calculated coefficients for the wing
print("M_inf:", helicopter_HS.M_inf)
print("R_inf:", helicopter_HS.Re_inf)
print("Drag coefficient for the wing:", f_Awing)
print("Dimensional drag coefficient for the wing:", f_wing)
