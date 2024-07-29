#Copter climb F2ROTOR
#The climb performance of a helicopter is determined by its ability to convert available power into vertical ascent.
#The primary factors influencing this performance are the available power, the power required by the rotors, the maximum takeoff weight (MTOW), and the horizontal flight speed.
#- Available Power (Pd): This is the total power output provided by the helicopter's engines.
#- Rotor Power (Pn): This includes the power consumed by both the main rotor and the tail rotor.
#- Maximum Takeoff Weight (MTOW): This is the maximum weight at which the helicopter is certified to take off.
#- Horizontal Flight Speed (Voo): This is the speed at which the helicopter is moving horizontally.

#The theory behind this code can be found in: Giovanni Di Giorgio - "Lezioni integrative dell'insegnamento di Aerodinamica dell'Ala Rotante" - a.a. 2023/2024 - page 89.

#Authors: Alessandra Di Rauso, Sara Vitale
    #Date: 16/06/2024
    #Version: 1.00

import numpy as np
def Helicopter_properties(Pd, Pn, MTOW, Voo):
    g = 9.81  
    '''
    Acceleration due to gravity
    '''

    '''
     Function to check and replace zero values for Voo
     '''
    def check_velocity(value, default=1):
        return value if value != 0 else default

    '''
    Check and replace Voo if it is zero
    '''
    Voo = check_velocity(Voo)

    properties = {
        "Pd": Pd,      
        "Pn": Pn,     
        "MTOW": MTOW,  
        "Voo": Voo,    
        "g": g        

        '''
        Available power Pd, Main and tail rotor power Pn, Maximum takeoff weight MTOW, Horizontal speed (fixed) Voo, Acceleration due to gravity g
        '''
    }

    def climb():
        '''
         Calculate the rate of climb (ROC) in meters per second (m/s)
         '''
        ROC = (Pd - Pn) / (MTOW * g)

        '''
         Calculate the climb angle gamma in radians
         '''
        gamma = np.arctan(ROC / Voo)

        '''
         Convert the ROC from meters per second (m/s) to feet per minute (ft/min)
         '''
        ROC_ft_min = ROC * 3.28 * 60

        '''
         Convert the climb angle gamma from radians to degrees
         '''
        gamma_deg = np.rad2deg(gamma)

        '''
         Return the rate of climb in feet per minute and the climb angle in degrees
         '''
        return ROC_ft_min, gamma_deg
    
    '''
     Return the properties and the climb function
     '''
    return properties, climb

'''
 Function to run tests and display results
 '''
def test_Helicopter_climb(test_cases):
    for Pd, Pn, MTOW, Voo in test_cases:
        properties, climb = Helicopter_properties(Pd, Pn, MTOW, Voo)
        ROC, gamma = climb()
        print(f"Pd: {Pd}, Pn: {Pn}, MTOW: {MTOW}, Voo: {Voo}")
        print("Rate of Climb (ft/min):", ROC)
        print("Climb Angle (degrees):", gamma)
        print("---------")


