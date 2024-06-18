import numpy as np
def evaluate_range(V_inf, W_fuel, SFC, Pn):

"""
Function that evaluates the autonomy with a fixed amount of fuel at different speeds of flight. 
The function requires as input the array of the speed of flight, the fuel weight, the specific fuel consumption and the required power. 
The input data values ​​must be in m/s, N, kg/W*h and W, respectively.
The function returns a tuple containing the range, maximum range speed and maximum range. Indicated respectively as Range, V_max_range, max_range.
The output data is provided in km/h for the speed and in km for the range.

    Args:
    V_inf (array): Velocity of flight (m/s).
    W_fuel (float): Weight of fuel (N).
    SFC (float): Specific Fuel Consuption (Kg/Wh).
    Pn (array): Power required (W).

    Returns:
    tuple:  Range array (Range), the velocity of max range (V_max_range),
           and the max range (max_range).


The theory behind this code can be found in: Giovanni Di Giorgio - "Lezioni integrative dell'insegnamento di Aerodinamica dell'Ala Rotante" - a.a. 2023/2024 - pages 87.

    Authors: Fabio Pisano, Beniamino Ferrara
    Date: 13/06/2024
    Version: 1.00
""" 


    Range = (3.6 * V_inf) * W_fuel / (SFC * Pn) #Calculation of the range in km
    Range = np.round(Range, 3)                  #Rounding the range to three decimal places
    max_range_index = np.argmax(Range)          #Localization of the index corresponding to the range
    max_range = np.max(Range)                   #Extracting the maximum value of the range array
    
    
    #Inverse formula for speed corresponding to maximum range
    V_max_range = max_range * (SFC * Pn[max_range_index]) / W_fuel
    V_max_range = np.round(V_max_range, 3)     #Rounding the range to three decimal places
    
    return Range, V_max_range, max_range

