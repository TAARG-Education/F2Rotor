from Helicopter_Class import Helicopter
import numpy as np
from GroundEffect import *

'''
This test case aims to validate the function "GroundEffect.py"

Authors: Andrea Iuliano, Michele Follo
'''
# Input data
R_mr = 10  # Main rotor radius [m]
Height = 4 # Ground distance [m]

# Helicopter Class Definition (0 values are related to not used input in this test case)                                                          
helicopter = Helicopter(0, R_mr, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

# Calculation of the correction coefficients
T_ratio_IGE, P_ratio_IGE = ground_effect(helicopter, Height)
print("T_ratio_IGE: {}, P_ratio_IGE: {}".format(T_ratio_IGE, P_ratio_IGE))