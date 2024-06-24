import numpy as np # type: ignore

def TailRotor(Helicopter, V_infty, Qc, k_tr = 1.20):
      
    '''
    This function computes the induced power, profile power, and total power
    coefficients needed for the tail rotor operation. The tail rotor counteracts the torque
    produced by the main rotor, ensuring directional stability and control.

    Steps and parameters involved:
    1. Accept in input:
        - Helicopter class with main and tail rotor parameters 
        - V_infty 
        - the main rotor torque coefficient (Qc) --> one can retrive it from Helicopter_power_advance function
        - k_tr which is the factor related to the induced power losses (default value: 1.20)
    2. Extract necessary tail rotor configuration parameters from the helicopter's config dictionary:
        - R_tr: Tail rotor radius
        - c_tr: Tail rotor chord
        - N_tr: Number of tail rotor blades
        - Omega_tr: Tail rotor tip speed
        - b_tr: Arm relative to the helicopter's center of gravity for the calculation of counter-torque
    3. Calculate the tail rotor disk area (A_tr).
    4. Define constants for induced power correction (k_tr) and profile drag coefficient (CD_tr).
    5. Calculate the advance ratio (mu_tr) based on forward velocity (V_infty) and tip speed.
    6. Determine the dimensionless ratio of tail rotor thrust to rho.
    7. Compute the inflow ratio (lambda_it) using a quadratic equation.
    8. Calculate the induced power (Cp_tr_i) and profile power (Cp_tr_0) coefficients.
    9. Sum the induced and profile power coefficients to get the total tail rotor power coefficient (Cp_tr).
    10. Return the induced power, profile power, and total power coefficients.

    INPUT: Qc, Helicopter class, V_infty, k_tr
    OUTPUT: Cp_tr_i, Cp_tr_0, Cp_tr

    References: 
    - "Lezioni di aerodinamica dell'ala rotante - Eliche rotori ed aeromotori con un'introduzione all'aerodinamica instazionaria",
       Renato Tognaccini. From page 97 to 100. 
    - "Lezioni integrative dell'insegnamento di Aerodinamica dell'Ala Rotante" - a.a. 2023/2024, Giovanni Di Giorgio. 

    Authors: Fabio Beltratti, Giuseppe Russo
    Latest update: 21/06/2024
    Version: 1.0
    '''
    
    # Main rotor parameters 
    Omega_R_mr = Helicopter.Omega_r_mr                            # Main rotor tip speed
    R_mr = Helicopter.R_mr                                        # Main rotor radius 
    A_mr = np.pi * R_mr ** 2                                      # Main rotor area

    # Tail rotor parameters
    R_tr = Helicopter.R_tr                                        # Tail rotor radius
    c_tr = Helicopter.c_tr                                        # Tail rotor chord
    N_tr = Helicopter.N_tr                                        # Tail rotor blades number
    Omega_r_tr = Helicopter.Omega_r_tr                            # Tail rotor tip speed                 
    l_tr = Helicopter.l_tr                                        # Tail rotor torque arm
    A_tr = np.pi * R_tr ** 2                                      # Tail rotor area 

    CD_tr = 0.01                                                  # Tail rotor profile drag coefficient                    
    mu_tr = V_infty / Omega_r_tr                                  # Tail rotor advance ratio (hp: AoA << 1)

    T_tr_over_rho = Omega_R_mr ** 2 * A_mr * R_mr * Qc / l_tr     # Tail rotor thrust over density

    # Induced input ratio
    lambda_it = np.sqrt(- V_infty ** 2 / 2 + 0.5 * np.sqrt(V_infty ** 4 
                     + 4 * (T_tr_over_rho / (2 * A_tr)) ** 2)) / Omega_r_tr   
    
    # Power Coefficients calculation
    Cp_tr_i = k_tr * lambda_it * R_tr * Qc / l_tr
    Cp_tr_0 = 1 / 8 * N_tr * c_tr * R_tr * CD_tr * (1 + 4.7 * mu_tr ** 2) / A_tr
    Cp_tr = Cp_tr_0 + Cp_tr_i

    return Cp_tr_i, Cp_tr_0, Cp_tr
