import numpy as np

'''
    This function computes the ground effect correction coefficients for the performances of the main rotor

    Steps and parameters involved:
    1. Accept in input:
        - Helicopter class with main and tail rotor parameters 
        - Height: ground distance 
    2. Extract main rotor radius R_mr from the Helicopter class
    3. Calculate the ground effect correction coefficient for the main rotor thrust.
    4. Calculate the ground effect correction coefficient for the main rotor needed power.
    5. Return the ground effect correction coefficients.

    INPUT: Helicopter class, Height
    OUTPUT: T_ratio_IGE, P_ratio_IGE

    References: 
    - "Lezioni di aerodinamica dell'ala rotante - Eliche rotori ed aeromotori con un'introduzione all'aerodinamica instazionaria",
       Renato Tognaccini. From page 88 to 89. 
    - "The Effect of the Ground on a Helicopter Rotor in Forward Flight" - 1957, I. C. Cheeseman and W. E. Bennett. 

    Authors: Andrea Iuliano, Michele Follo
    Latest update: 13/07/2024
    Version: 1.0
'''
 
# Definizione della funzione ground_effect
def ground_effect(Helicopter, Height):
    R_mr = Helicopter.R_mr                           # Main Rotor Radius
    T_ratio_IGE = 1 / (1 - (R_mr / (4 * Height))**2) # Thrust Correction Coefficient
    K_G = 1 / T_ratio_IGE
    P_ratio_IGE = K_G                                # Power Correction Coefficient
    return T_ratio_IGE, P_ratio_IGE