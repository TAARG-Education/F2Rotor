"""
This test case script calculates aerodynamic coefficients for the horizontal stabilizer of a simulated Sirosky UH-60 Black Hawk using Python.
The script utilizes parameters such as wing span, wing chord, wing thickness ratio, and aspect ratio specific to the copter selected, all defined in a configuration dictionary (config).
Additionally, flight parameters such as cruise speed, altitude, speed of sound, kinematic viscosity, and density are specified to compute the Mach number and Reynolds number necessary for aerodynamic calculations.

Author: Eduardo Duraccio
    Date: 17/07/2024
    Version: 1.02
"""
import matplotlib.pyplot as plt
import math
import numpy as np


def wing_comp(self):
    b = self.config["Lift-Compound"][0]["b_w"]  # wing span
    c = self.config["Lift-Compound"][0]["c_w"]  # rectangular wing chord
    tau = self.config["Lift-Compound"][0]["tau"]  # wing thickness
    AR_w = self.config["Lift-Compound"][0]["AR"]  # aspect ratio of rectangular wing
    r_mr = self.config["Lift-Compound"][0]["r_mr"] # radius of main rotor 
    R_LS = self.config["Lift-Compound"][0]["R_LS"] # lift-induced skin friction factor 
    Hf = self.config["Lift-Compound"][0]["Hf"] # hinge factor
    r = self.config["Lift-Compound"][0]["r"] # compressibility correction factor
    gamma = self.config["Lift-Compound"][0]["gamma"] # ratio of specific heats for air
    K = self.config["Lift-Compound"][0]["K"] # admissible roughness
    M_oo = self.M_inf  # freestream Mach number
    Re_oo = self.Re_inf  # freestream Reynolds number

    S = b**2/AR_w # wing area
    S_rotor= 3.14*r_mr**2 # reference area (rotor area)
    Swet = 2 * S  # wetted area of the wing (both sides)
    FF = 1 + 1.68 * (tau / c) + 3 * (tau / c)**2  # form factor for 4-series airfoil
    IF = R_LS * Hf # interference factor for horizontal stabilizer 

    
    t = (1 + r * (gamma - 1) / 2 * M_oo**2)**(-1)  # thermodynamic temperature ratio factor
    F = 1 + 0.03916 * M_oo**2 * t  # compressibility correction factor

    Re_L1 = Re_oo  # initial Reynolds number condition
    K1 = 37.587 + 4.615 * M_oo + 2.949 * M_oo**2 + 4.132 * M_oo**3  # # admissible roughness
    Re_L2 = K1 * (c / K)**(1.0489)  # second Reynolds number condition
    Re_L = min([Re_L1, Re_L2])  # select the smaller Reynolds number condition
    Cf = t * F**2 * 0.430 / (math.log10(Re_L * t**1.67 * F))**2.56  # skin friction coefficient for turbolent flow


    f_Awing = (Cf * Swet / S_rotor) * FF * IF  # drag coefficient for the wing
    f_wing = f_Awing * S_rotor  # dimensional drag coefficient for the wing
    f_Awing = round(f_Awing, 5)  # round drag coefficient for the wing
    f_wing = round(f_wing, 4)  # round dimensional drag coefficient for the wing
    return f_Awing, f_wing  # return dimensional and non dimensional drag coefficient



# In this case study it's used Sirosky UH-60 Black Hawk
config = {
    "Lift-Compound": [
        {
            "b_w": 4.38,    # Wing span [m] 
            "c_w": 0.95,     # Rectangular wing chord [m]  
            "tau": 0.14,    # Wing thickness ratio 
            "AR": 4.61,      # Rectangular wing aspect ratio 
            "r_mr": 8.18,    # Radius of main rotor [m] 
            "R_LS":1.08,    # lift-induced skin friction factor
            "Hf":1.1,       #hinge factor
            "r":0.89,       #compressibility correction factor
            "gamma":1.4,     # ratio of specific heats for air
            "K":0.0009      #admissible roughness
        }
    ]
}

# Additional parameters needed for calculation of flight condition considering 50% of maximum spped and selected height
v = 49.5833  # 50% of maximum speed  [m/s]
height = 2000    # Altitude [m] 
a = 332.529       # Speed of sound [m/s]
nu = 1.7148e-5      # Kinematic viscosity [m/s]
rho = 1.0065     # Density [kg/m³]

# Calculate the Mach number
M_inf = v / a 

# Calculate the Reynolds number
Re_inf = rho * v * config["Lift-Compound"][0]["c_w"] / nu  # Reynolds number

# Define a simulated aircraft class to hold configuration and flight parameters
class SimulatedAircraft:
    def __init__(self, config):
        self.config = config
        self.M_inf = M_inf
        self.Re_inf = Re_inf

# Create an instance of SimulatedAircraft with the configuration
aircraft = SimulatedAircraft(config)

# Call the wing_comp function with the simulated aircraft instance
f_Awing, f_wing = wing_comp(aircraft)

# Print the calculated coefficients for the wing
print("M_inf:", M_inf)
print("R_inf:", Re_inf)
print("Drag coefficient for the wing:", f_Awing)
print("Dimensional drag coefficient for the wing:", f_wing)
