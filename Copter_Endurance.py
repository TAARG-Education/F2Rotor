    #Function that evaluates the endurance with a fixed amount of fuel at different speeds of flight. 
    #The function requires as input the array of the Speed of Flight, the Fuel Weight, the Specific Fuel Consumption and the Required Power. 
    #The input data values ​​must be in m/s, Kg, kg/KW*h and KW, respectively.
    #The function returns the maximum Time of Flight for each Speed.
    #The output data is provided in min.
    
    #Authors: Davide Sergio, Carmine Marra
    #Date: 23/09/2024
    #Version: 1.00
    
    #The theory behind this code can be found in: Giovanni Di Giorgio - "Lezioni integrative dell'insegnamento di Aerodinamica dell'Ala Rotante" - a.a. 2023/2024 - pages 87.

import numpy as np

def Endurance(Pn, Wfuel, Voo):
    
    """
    Input:
    Pn    = Required Power[KWatt]
    Wfuel = Fuel Weight [Kg]
    Voo   = Cruise Speed Array [m/s]

    Output:
    Endurance: Time of flight Array for any given speed.
    
    """

    P_s = Pn            # Power [KW]
    Wfuel_kg = Wfuel    # Fuel Weight [Kg]
    SFC = 0.458         # Specific Fuel Consumption (Constant)

    time = np.zeros(len(Voo))  # Array initialization for flight times

    for i, v in enumerate(Voo):
        t_i = Wfuel_kg / (SFC * P_s/60)  # Endurance [min]
        time = t_i
        i +=1
    return time

