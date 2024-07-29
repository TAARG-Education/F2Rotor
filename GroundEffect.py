# This function computes the ground effect correction coefficients 
# for the performances of the main rotor

# Steps and parameters involved:
    # 1. Accept in input:
    #     - Helicopter class with main and tail rotor parameters 
    #     - Height: ground distance 
    #     - Radius: main rotor radius
    # 2. Calculate the ground effect correction coefficient 
    #     for the main rotor thrust.
    # 3. Calculate the ground effect correction coefficient 
    #     for the main rotor needed power.
    # 4. Return the ground effect correction coefficients.

# INPUT: Radius, Height
# OUTPUT: T_ratio_IGE, P_ratio_IGE

# References: 
# - "Lezioni di aerodinamica dell'ala rotante - Eliche rotori ed aeromotori 
#     con un'introduzione all'aerodinamica instazionaria",
#     Renato Tognaccini. From page 88 to 89. 
# - "The Effect of the Ground on a Helicopter Rotor in Forward Flight" - 1957, 
#     I. C. Cheeseman and W. E. Bennett. 

# Authors: Andrea Iuliano, Michele Follo
# Latest update: 15/07/2024
# Version: 1.0



import numpy as np
 
# Definizione della funzione ground_effect
def ground_effect(Radius, Height): 

    '''
    This function computes the ground effect correction coefficients 
    for the performances of the main rotor

    INPUT: Radius, Height
    OUTPUT: T_ratio_IGE, P_ratio_IGE

    Authors: Andrea Iuliano, Michele Follo
    Latest update: 15/07/2024
    Version: 1.0
    '''

    # Thrust Correction Coefficient                       
    T_ratio_IGE = 1 / (1 - (Radius / (4 * Height))**2) 
    K_G = 1 / T_ratio_IGE
    # Power Correction Coefficient
    P_ratio_IGE = K_G                                
    return T_ratio_IGE, P_ratio_IGE