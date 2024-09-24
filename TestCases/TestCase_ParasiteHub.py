#This is a test case for the Parasite Hub function which calculates the parasite drag area
#of the hub, assuming it is a circular section cylinder.

#Author: Arianna Palumbo
#Rotary Wing Aerodynamics course, prof. Renato Tognaccini
#University of Naples Federico II
#Academic Year 2023-2024

import math as ma
from parasite_hub import hub


config = {
    "ParasiteArea": [
        {}, {}, {},
        {
            "R": 0.4,  # Radius of the hub (m)
            "Lcs": 0.7,  # Height of the hub (m)
            "alpha": ma.pi/2  # Angle of attack (radians)
        }
    ]
}
S = 35  # Reference Area (m^2)

# Calling the function for the calculation of the parasite drag area of the hub
f_Acs, f_cs = hub(config, S)

# Results printing
print(f"f_Acs: {f_Acs}")
print(f"f_cs: {f_cs}")
