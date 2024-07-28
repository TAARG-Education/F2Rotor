#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#x  Test cases for validating the function 'parasite_fus.py' to calculate the equivalent wetted area of conventional and compound helicopters.                             x
#x   - For the validation related to the conventional helicopter case, Aérospatiale SA 341 Gazelle is considered:                                                          x
#x     https://www.thisdayinaviation.com/tag/aerospatiale-sa-341-gazelle/                                                                                                  x                                                                                                    
#x   - For the validation related to the compound helicopter case, the Lockheed AH-56A Cheyenne is considered:                                                             x
#x      dimensions: https://www.aviastar.org/helicopters_eng/lok_cheyenne.php                                                                                              x
#x      max speed (SL): https://www.armyaviationmagazine.com/images/archive/backissues/1980/80_06_07.pdf                                                                   x
#x                                                                                                                                                                         x
#x To use the 'parasite_fus.py' function, it is necessary to instantiate the DataImport class by assigning a value to the following attributes:                            x                                                                                                      
#x   1) Re_inf: asymptotic Reynolds number based on fuselage length                                                                                                        x                                                           
#x   2) M_inf: asymptotic Mach number                                                                                                                                      x
#x   3) config: geometric dimensions of the fuselage (length L and cross-sectional dimensions h and d)                                                                     x
#x   4) S: reference area (rotor disk area)                                                                                                                                x                                     
#x   5) helicopter_type: helicopter type (conventional or compound). We use the instance “conventional” to indicate that the studied helicopter is conventional;           x
#x      “compound” to indicate that it is a compound helicopter.                                                                                                           x                                               
#x                                                                                                                                                                         x
#x  The max speed at sea level is considered for calculating the asymptotic Mach number (clearly different choices are possible).                                          x                                                                                                      
#x                                                                                                                                                                         x
#x Author: Iole Paolillo                                                                                                                                                   x
#x Latest update: 21/07/2024                                                                                                                                               x
#x Version: 1.0                                                                                                                                                            x                                                                                                                                                                                                                                                                                                                                                                   
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

import math as ma
from ambiance import Atmosphere
from parasite_fus import *       # Importing all objects defined in the module

# Input parameters to be used to study the Aérospatiale SA 341 Gazelle helicopter (input parameters to be used to study the Lockheed AH-56A Cheyenne compound helicopter are shown in parentheses)
helicopter_type = "conventional" # (for Lockheed AH-56A Cheyenne use: "compound")       
R_mr = 10.500 / 2                 # [m] main rotor radius (for Lockheed AH-56A Cheyenne use: 15.62 / 2 m)                                                        
S = ma.pi * R_mr**2               # [m^2] rotor disk area
config = {
    "ParasiteArea": [{"h":2.16, "d": 1.4, "L": 9.533}]   # fuselage dimensions Bell 505 [m] (for Lockheed AH-56A Cheyenne use: "h":3.39, "d": 2.15, "L": 16.63)
}
L_fuselage = config["ParasiteArea"][0]["L"]                
V_max = 86.43                                            # max speed SL [m/s] (= 168 kts) (for Lockheed AH-56A Cheyenne use: 113.178 m/s = 220 knots)

# Standard atmospheric parameters at sea level
atmosphere    = Atmosphere(0)                              # set sea level
rho_sea_level = atmosphere.density[0]                      # air density at sea level in kg/m^3
mu_sea_level  = atmosphere.dynamic_viscosity[0]            # air dynamic viscosity at sea level in kg/(m·s)
a_sea_level   = atmosphere.speed_of_sound[0]               # speed of sound at sea level in m/s
# Calculation of Re_inf and M_inf
Re_inf = (rho_sea_level * V_max * L_fuselage) / mu_sea_level
M_inf  = V_max / a_sea_level

# Instantiate the DataImport class with the calculated values
mydata = DataImport(Re_inf, M_inf, config, S, helicopter_type)
# Print results
result_inter_fact = inter_fact(mydata)
if result_inter_fact is not None:
    f_Af, f_f = fuselage(mydata)
    print(f"The correlation factor is: {result_inter_fact}")
    print(f"f_Af: {f_Af}, f_f: {f_f}")
else:
    print("Failed to calculate the correlation factor due to input errors.")
