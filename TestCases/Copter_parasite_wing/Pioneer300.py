#This test case script calculates aerodynamic coefficients for the wing of a simulated Pioneer 300 aircraft using Python.
#The script utilizes parameters such as wing span, wing chord, wing thickness ratio, and aspect ratio specific to the Pioneer 300, all defined in a configuration dictionary (config).
#Additionally, flight parameters such as cruise speed, altitude, speed of sound, kinematic viscosity, and density are specified to compute the Mach number and Reynolds number necessary for aerodynamic calculations.

#Author: Eduardo Duraccio
#Date: 17/07/2024
#Version: 1.02

import matplotlib.pyplot as plt
import math
import numpy as np

import matplotlib.pyplot as plt
import math
import numpy as np
from Copter_parasite_wing import parasite_wing
from Helicopter_Class import Helicopter_HS


b_w = 8.10    # Wing span [m] 
c_w = 1.2     # Rectangular wing chord [m]  
tau = 0.12    # Wing thickness ratio 
AR = 5.96      # Rectangular wing aspect ratio 
R_mr = 0    # Radius of main rotor [m] 
R_LS =1.08    # lift-induced skin friction factor
Hf =1.1       #hinge factor
r = 0.89       #compressibility correction factor
gamma = 1.4     # ratio of specific heats for air
K =0.0009      #admissible roughness
R_w_b = 1.07    #interference factor wing-fusolage
ala_fissa = True # fixed or rotary wing
profilo = "4-series" # type of airfoil
turbolentflow = True    #laminar or turbolent flow

# Additional parameters needed for calculation of flight condition 
v = 37.5  # Cruise speed at 75% power in m/s
height = 2000    # Altitude in meters
a = 332.5        # Speed of sound in m/s
nu = 1.7e-5      # Kinematic viscosity in m²/s
rho = 1.0065     # Density in kg/m³


helicopter_HS = Helicopter_HS(tau=tau, b_w=b_w, AR=AR, R_LS=R_LS, Hf=Hf, r=r, gamma=gamma, k=K, 
                              R_mr=R_mr, c_w=c_w, R_w_b=R_w_b, ala_fissa= ala_fissa, profilo= profilo, turbolentflow= turbolentflow)
helicopter_HS.set_M_inf(v, a)
helicopter_HS.set_Re_inf(v, rho, nu)

# Call the parasite_wing function with helicopter_HS istance
f_Awing, f_wing = parasite_wing(helicopter_HS)

# Print the calculated coefficients for the wing
print("M_inf:", helicopter_HS.M_inf)
print("R_inf:", helicopter_HS.Re_inf)
print("Drag coefficient for the wing:", f_Awing)
print("Dimensional drag coefficient for the wing:", f_wing)

